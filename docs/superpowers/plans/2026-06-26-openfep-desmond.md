# openfep-desmond Implementation Plan
**Date:** 2026-06-26  
**Spec:** `docs/superpowers/specs/2026-06-25-openfep-desmond-design.md`

---

## Overview

Build `openfep` — an open-source FEP toolkit that drives Schrodinger's Desmond MD engine without requiring Schrodinger's Python stack at runtime. Supports RBFE and ABFE via MSJ generation, job submission, and MBAR-based analysis.

**Forcefield:** OPLS4 (OPLSVersion.F17, string `"OPLS4"` in MSJ)  
**Input:** Single complex `.mae` with `s_fep_struc_tag` CT properties  
**Production stage name:** `lambda_hopping` (not `simulate`)

---

## Phases and Tasks

### Phase 1 — Package Scaffold + Core Types

#### Task 1 — Package skeleton and constants

Create the package layout and all shared constants before any other code is written.

**Files to create:**
- `openfep/__init__.py` — version + public re-exports
- `openfep/constants.py` — all numeric and string constants

**Constants to define in `openfep/constants.py`:**

```python
# Forcefield
FORCEFIELD = "OPLS4"           # OPLSVersion.F17 string value in MSJ

# Buffer widths (Å)
COMPLEX_BUFFER_WIDTH = 5.0
NET_CHARGE_COMPLEX_BUFFER_WIDTH = 8.0
SOLVENT_BUFFER_WIDTH = 10.0
VACUUM_BUFFER_WIDTH = 100.0

# HMR
HMR_TIMESTEPS = [0.004, 0.004, 0.008]   # ps
HMR_MIGRATION_INTERVAL = 0.024           # ps

# Salt
MIN_CHARGED_SALT_CONC = 0.15            # M

# GCMC
GCMC_SOLVENT_VDW_SCALE_FACTOR = 0.75

# ENERGY_GROUP output
ENEGRP_NAME_TEMPLATE = "$JOBNAME$[_replica$REPLICA$]_enegrp.dat"
ENEGRP_FIRST = 0.0   # ps
ENEGRP_INTERVAL = 1.2  # ps

# s_fep_struc_tag values
class StrucTag:
    RECEPTOR = "receptor"
    LIGAND   = "ligand"
    COMPLEX  = "complex"
    SOLVENT  = "solvent"
    MEMBRANE = "membrane"

# ABFE lambda counts (from fep_schedule.py:282-301)
# Default (neutral)
ABFE_COMPLEX_LAMBDAS = 68
ABFE_SOLVENT_LAMBDAS = 60
ABFE_RESTRAINED_COMPLEX_LAMBDAS = 80
ABFE_RESTRAINED_SOLVENT_LAMBDAS = 68
# Charged
ABFE_CHARGED_COMPLEX_LAMBDAS = 108
ABFE_CHARGED_SOLVENT_LAMBDAS = 60
ABFE_CHARGED_RESTRAINED_COMPLEX_LAMBDAS = 128
ABFE_CHARGED_RESTRAINED_SOLVENT_LAMBDAS = 68

# RBFE lambda count (default; adaptive override allowed)
RBFE_DEFAULT_LAMBDAS = 12

# Backbone ASL (receptor restraints)
BACKBONE_ASL = (
    "((protein and backbone) or "
    "(nucleic_acids and nucleic_backbone)) and not a.ele H"
)
RECEPTOR_ASL = "protein or nucleic_acids"
```

**Keep-struc-tag maps** (mirroring `desmond/msj_generator/keep_struc_tags.py`):

```python
# FepLeg tag selection per fep_type / leg_type
KEEP_STRUC_TAGS = {
    "rbfe": {
        "complex": [StrucTag.RECEPTOR, StrucTag.MEMBRANE,
                    StrucTag.SOLVENT, StrucTag.LIGAND],
        "solvent": [StrucTag.LIGAND],
        "vacuum":  [StrucTag.LIGAND],
    },
    "abfe": {
        "md":      [StrucTag.RECEPTOR, StrucTag.SOLVENT,
                    StrucTag.MEMBRANE, StrucTag.LIGAND],
        "complex": [StrucTag.RECEPTOR, StrucTag.MEMBRANE,
                    StrucTag.SOLVENT, StrucTag.COMPLEX],
        "solvent": [StrucTag.LIGAND],
    },
}
```

---

#### Task 2 — MAE CT parser

Parse a single `.mae` file with multiple CT blocks and extract per-CT metadata.

**File:** `openfep/mae_parser.py`

**API:**
```python
@dataclass
class CTBlock:
    struc_tag: str       # value of s_fep_struc_tag property
    atom_count: int
    raw_block: str       # the full CT text

def parse_mae(path: str | Path) -> list[CTBlock]:
    """Return all CT blocks in order, each annotated with struc_tag."""
```

**Implementation notes:**
- Split on `f_m_ct {` boundaries (MAE format is plain text)
- Read `s_fep_struc_tag` from the CT-level property block
- CTs without `s_fep_struc_tag` are silently skipped
- No Schrodinger Python stack required — plain string parsing

**Tests:** `tests/test_mae_parser.py` with a synthetic two-CT `.mae` fixture (receptor + ligand).

---

### Phase 2 — FEP Graph

#### Task 3 — Atom mappers (Kartograf + Lomap)

Two concrete mapper implementations behind a common abstract base.

**Files:**
- `openfep/mapping/base.py` — `AtomMapper` ABC with `map(mol_a, mol_b) -> AtomMapping`
- `openfep/mapping/kartograf_mapper.py` — wraps `kartograf.atom_mapper.KartografAtomMapper`
- `openfep/mapping/lomap_mapper.py` — wraps `lomap.MCS`

**`AtomMapping` dataclass:**
```python
@dataclass
class AtomMapping:
    mol_a: Chem.Mol
    mol_b: Chem.Mol
    a_to_b: dict[int, int]    # atom-index mapping
    score: float              # Lomap-compatible 0–1 score
```

**Kartograf mapper:**
- Uses 3D geometry: `KartografAtomMapper(atom_map_hydrogens=False)`
- Score = Lomap MCS Tanimoto of the mapped core

**Lomap mapper:**
- `lomap.MCS` with `max3d=1.0` threshold; score = MCS Lomap score
- Fallback when Kartograf unavailable

---

#### Task 4 — FEP perturbation network

Build the transformation graph from a list of molecules, score edges, identify optimal perturbation network.

**File:** `openfep/mapping/network.py`

**API:**
```python
class FEPNetwork:
    def __init__(self, molecules: list[Chem.Mol], mapper: AtomMapper): ...
    def build(self) -> nx.Graph:
        """All-pairs mapping; nodes=molecules, edges weighted by score."""
    def optimal_network(self) -> nx.Graph:
        """Lomap-style MST-based star/hub network; prune low-score edges."""
    def rbfe_pairs(self) -> list[tuple[Chem.Mol, Chem.Mol, AtomMapping]]:
        """Return ordered list of ligand pairs for RBFE."""
```

**Scoring:** Edge weight = 1 - score (so MST maximises mapping quality).  
**Pruning:** Drop edges with score < 0.2 (Lomap default threshold).

---

### Phase 3 — System Preparer

#### Task 5 — RBFE system preparer

Prepare the `.fmp` graph file and configure the `fep_mapper` MSJ stage for each RBFE edge.

`fep_mapper` is a Desmond MSJ stage (`desmond/stage/fep_mapper.py:490`), not a subprocess utility. It internally calls `run -FROM scisol fep_mapper.py`. The dual-topology merge happens inside `multisim` at runtime. openfep does not produce a pre-merged dual-topology `.mae` — instead it emits the `.fmp` graph file and a `main.msj` with the `fep_mapper` stage pre-configured.

**File:** `openfep/prep/rbfe_prep.py`

**Workflow:**
1. Parse input `.mae` → `CTBlock` list via `mae_parser.parse_mae()`
2. Write Kartograf atom mapping to `{edge}_atom_mapping.txt` (format: `atom_a_idx atom_b_idx` per line)
3. Build `{jobname}.fmp` FEP graph file: edge list + structure references (required by `fep_mapper` stage)
4. Detect net charge of ligand pair → select protocol:
   - Neutral: standard complex + solvent legs
   - Charged: cubic box, `neutralize_system=True`, `min_salt_conc=0.15`, `make_alchemical_water=True`
5. Compute `receptor` count for `fep_mapper` stage: `1 + bool(has_solvent_ct) + bool(has_membrane_ct)`
6. Return `RBFESystem(fmp_path, atom_mapping_path, receptor_count, is_charged, buffer_width)`

The `fep_mapper` stage block is emitted in `main.msj` by the MSJ writer (Task 10), not here.

**Buffer width selection:**
```python
buffer = NET_CHARGE_COMPLEX_BUFFER_WIDTH if is_charged else COMPLEX_BUFFER_WIDTH
```

---

#### Task 6 — ABFE system preparer

Prepare the 3-phase (md + complex + solvent) ABFE input.

**File:** `openfep/prep/abfe_prep.py`

**Workflow:**
1. Parse input `.mae` → filter CTs per `KEEP_STRUC_TAGS["abfe"]`
2. Detect net charge → set `add_alchemical_ions=True` if charged (ABFE uses ions, not alchemical water)
3. Detect membrane presence (any CT with `struc_tag == "membrane"`) → set GCMC flag
4. Compute Boresch restraint geometry:
   - Select 3 receptor atoms + 3 ligand atoms from 3D coordinates
   - Return as `BoreschRestraint` dataclass (6 distances + angles + dihedrals)
5. Return `ABFESystem(mae_path, boresch, is_charged, has_membrane, lambda_counts)`

**Lambda count selection (from constants):**
```python
if has_boresch_restraints:
    complex_lambdas = ABFE_RESTRAINED_COMPLEX_LAMBDAS   # 80
    solvent_lambdas = ABFE_RESTRAINED_SOLVENT_LAMBDAS   # 68
else:
    complex_lambdas = ABFE_COMPLEX_LAMBDAS              # 68
    solvent_lambdas = ABFE_SOLVENT_LAMBDAS              # 60
```

---

### Phase 4 — Lambda Scheduler

#### Task 7 — Lambda presets

**File:** `openfep/lambda_scheduler/presets.py`

**RBFE presets** — fixed window counts, linear λ spacing:
```python
def rbfe_lambdas(n: int = RBFE_DEFAULT_LAMBDAS) -> list[float]:
    """Uniform spacing 0..1 with softcore endpoints."""
```

**ABFE Gibbs schedules** — NOT simple lambda lists. Schrodinger uses full Gibbs weight vectors stored in SEA-format MSJ files at `$MMSHARE/data/desmond/abfep/`:
- `{leg}_schedule.msj` (default, neutral)
- `{leg}_schedule_restrained.msj`
- `{leg}_schedule_chg.msj` (charged)
- `{leg}_schedule_restrained_chg.msj`

```python
def load_abfe_schedule(mmshare_dir: Path, leg: str,
                       restrained: bool = False,
                       charged: bool = False) -> str:
    """
    Read the SEA-format schedule MSJ text from $MMSHARE/data/desmond/abfep/.
    Returns the raw .schedule block text to be injected verbatim into the
    subjob MSJ via _set_custom_lambda_schedule() equivalent.
    leg: "complex" | "solvent"
    """
```

Parse the returned SEA text for `.schedule.weights.N.val` to get the window count (do not hardcode it — read it from the file). The file also contains `.flexible[0].A.val` Gibbs weight vectors needed by Desmond's lambda-hopping engine.

**Solvent restrained-window count** (for REST scaling end_win):
```python
def abfe_solvent_restrained_windows(mmshare_dir: Path,
                                    restrained: bool,
                                    charged: bool) -> int:
    """
    Returns count of non-zero flexible[0].A.val entries from solvent schedule.
    Passed as end_win to _set_solute_tempering_temperature_ladder().
    Returns 0 when restraint_enabled=False.
    """
```

---

#### Task 8 — Adaptive lambda scheduler (a3fe)

Pilot-run–based redistribution using √Var[dG/dλ] criterion.

**File:** `openfep/lambda_scheduler/adaptive.py`

**Algorithm:**
1. Run short pilot simulation (1 ns) at n_init=5 windows
2. Parse `_enegrp.dat` from each window → compute dG/dλ gradient
3. Compute variance profile √Var[dG/dλ] across λ
4. Redistribute N_final windows to equidistribute variance
5. Re-run full simulation at redistributed windows

**API:**
```python
class A3FEScheduler:
    def __init__(self, n_final: int, pilot_ns: float = 1.0): ...
    def fit(self, pilot_results: list[EnergySeries]) -> list[float]:
        """Return redistributed lambda list."""
```

This task is **optional for MVP**; presets from Task 7 unblock all downstream tasks.

---

### Phase 5 — MSJ Generator

#### Task 9 — MSJ stage builder (shared primitives)

Low-level building blocks for assembling MSJ stage blocks.

**File:** `openfep/msj/stages.py`

**Stage constructors (each returns a dict that serialises to MSJ SEA format):**

```python
def build_assign_forcefield(forcefield: str = FORCEFIELD,
                            custom_charge_mode: str = "assign") -> dict: ...

def build_solvate_pocket(buffer: float, is_charged: bool = False) -> dict: ...

def build_minimize(restraints_on_heavy: bool = True) -> dict: ...

def build_simulate(time_ns: float, ensemble: str,
                   polarization_restraint: str | None = None,
                   is_last_equil: bool = False) -> dict:
    """
    ensemble: "NPT" | "NVT" | "muVT"
    polarization_restraint: only set when OPLS4 (F17).
      - "full"  for all equilibration stages except the last
      - "decay" for the last equilibration stage
    """

def build_lambda_hopping(time_ns: float, n_lambda: int,
                         lambda_values: list[float],
                         energy_groups: bool = True) -> dict:
    """
    Production stage — always named 'lambda_hopping' (not 'simulate').
    energy_groups=True always (deviation from Schrodinger default of complex-only).
    """

def build_load_restraints_from_structure() -> dict: ...

def build_trim() -> dict:
    """Erases *-in.cms.gz, *.cms.gz, *-out.tgz after production."""
```

**Polarization restraint gating (CRITICAL — OPLS4 only):**
```python
# Only inject polarization_restraints when forcefield == "OPLS4" (F17)
# All equilibration stages except last: polarization_restraints = full
# Last equilibration stage:            polarization_restraints = decay
# Production (lambda_hopping):         no polarization_restraints
```

**HMR injection** (applied to all simulate + lambda_hopping stages):
```python
HMR_BLOCK = {
    "hydrogen_mass_repartition": True,
    "timestep": HMR_TIMESTEPS,      # [0.004, 0.004, 0.008] ps
    "migration_interval": HMR_MIGRATION_INTERVAL,  # 0.024 ps
}
```

**Energy group block** (always-on — deliberate deviation):
```python
ENERGY_GROUP_BLOCK = {
    "name": ENEGRP_NAME_TEMPLATE,
    "first": ENEGRP_FIRST,      # 0.0 ps
    "interval": ENEGRP_INTERVAL,  # 1.2 ps
}
```

**GCMC stage** (only when membrane present):
```python
def build_gcmc(scale_solvent_vdw: float = GCMC_SOLVENT_VDW_SCALE_FACTOR) -> dict: ...
```

---

#### Task 10 — RBFE MSJ writer

Assemble the 2-file RBFE MSJ pair (complex.msj + solvent.msj or vacuum.msj).

**File:** `openfep/msj/rbfe_writer.py`

**Stage sequence for each leg:**

```
assign_forcefield      (custom_charge_mode="assign" for complex; "keep" for solvent)
solvate_pocket         (buffer from RBFESystem)
[GCMC]                 (complex only, if membrane)
minimize               (x2: heavy-atom restraints on, then off)
simulate NPT           (pre-equil 1: polarization_restraints="full")
simulate NPT           (pre-equil 2: polarization_restraints="full")
simulate NPT           (pre-equil 3: polarization_restraints="full")
simulate NPT           (last equil:  polarization_restraints="decay")
lambda_hopping         (production; energy_groups=True)
trim
```

**Number of equilibration stages:** Mirror the actual Schrodinger RBFE template — the spec sketch showed 2 stages but the source template has more. Implementer should count equilibration stages in `desmond/msj_generator/fep/small_molecule.py` and match exactly.

**`assign_custom_charge.mode`:**
- Complex leg: `"assign"` (charges assigned from OPLS4)
- Solvent leg: `"keep"` (reuses charges from complex leg output)

**`fep_mapper` stage in `main.msj`** (first stage, before leg subjobs):
```
fep_mapper {
  graph_file    = "{jobname}.fmp"
  atom_mapping  = "{edge}_atom_mapping.txt"
  receptor      = {receptor_count}   # 1 + bool(solvent_ct) + bool(membrane_ct)
}
```
The `main.msj` then dispatches the leg subjob MSJs (`complex.msj`, `solvent.msj`, optionally `vacuum.msj`) in sequence after `fep_mapper` completes.

**API:**
```python
def write_rbfe_msj(system: RBFESystem, output_dir: Path,
                   lambdas: list[float]) -> dict[str, Path]:
    """Write main.msj + complex.msj + solvent.msj; return path dict."""
```

---

#### Task 11 — ABFE MSJ writer

Assemble the 4-file ABFE MSJ set: `main.msj`, `md.msj`, `complex.msj`, `solvent.msj`.

**File:** `openfep/msj/abfe_writer.py`

**`md.msj` stage sequence (MD pre-equil + FEP primer):**
```
assign_forcefield      (custom_charge_mode="keep")
solvate_pocket
minimize               (x2)
simulate NPT           (polarization_restraints="full")
simulate NPT           (polarization_restraints="full")
simulate NPT           (polarization_restraints="full")
simulate NPT           (last equil: polarization_restraints="decay")
FepAbsoluteBindingFepPrimer    ← FINAL STAGE of md.msj
```

> **Critical:** `FepAbsoluteBindingFepPrimer` is the FINAL stage of `md.msj`, NOT of `complex.msj`.  
> Source: `desmond/msj_generator/fep/absolute_binding.py:231` (`generate_md_msj` calls `_set_fep_absolute_binding_fep_primer()`).

**`complex.msj` stage sequence:**
```
load_restraints_from_structure   ← NOT FepPrimer
assign_forcefield      (custom_charge_mode="keep")
solvate_pocket
minimize
simulate NPT           (polarization_restraints="full")
simulate NPT           (last equil: polarization_restraints="decay")
lambda_hopping         (n_lambda from ABFESystem; energy_groups=True)
trim
```

**`assign_custom_charge.mode` for ABFE legs (CRITICAL fix from review):**
- `md` leg: `"keep"` (charges from input structure)
- `complex` leg: `"keep"` (NOT `"assign"` — source: `absolute_binding.py:328`, `CUSTOM_CHARGE_MODE.KEEP`)
- `solvent` leg: `"keep"`

**`solvent.msj` stage sequence:**
```
assign_forcefield      (custom_charge_mode="keep")
solvate_pocket         (buffer=SOLVENT_BUFFER_WIDTH)
minimize
simulate NPT           (polarization_restraints="full")
simulate NPT           (last equil: polarization_restraints="decay")
lambda_hopping         (energy_groups=True)
trim
```

**REST scaling — solvent leg (restrained protocol only):**
Call `_set_solute_tempering_temperature_ladder(msj, start_win=0, end_win=end_win)` where `end_win = lambda_windows - abfe_solvent_restrained_windows(mmshare_dir, restrained, charged)`. This omits REST temperature scaling on the restrained λ-windows (where the ligand is still held in its bound pose). When `restrained=False`, `abfe_solvent_restrained_windows()` returns 0 so REST is not applied at all. Source: `absolute_binding.py:298-302`, `subjob_msj.py:194-200`.

**ABFE charged protocol:**
- Load schedule from `{leg}_schedule_chg.msj` or `{leg}_schedule_restrained_chg.msj`
- Window counts: 108/60 (unrestrained) or 128/68 (restrained) — read from schedule file, do NOT hardcode
- `add_alchemical_ions=True` (not `make_alchemical_water`)
- `neutralize_system=True`, cubic box, salt ≥ 0.15M

**`main.msj`:** Orchestrates md.msj → complex.msj + solvent.msj (parallel legs).

**API:**
```python
def write_abfe_msj(system: ABFESystem, output_dir: Path) -> dict[str, Path]:
    """Return {'main': ..., 'md': ..., 'complex': ..., 'solvent': ...}"""
```

---

### Phase 6 — Job Runner

#### Task 12 — Local job runner

Execute `multisim` for a single FEP job on the local machine.

**File:** `openfep/runner/local.py`

**API:**
```python
class LocalRunner:
    def __init__(self, schrodinger_dir: Path): ...

    def run_rbfe(self, pair_mae: Path, msj_files: dict[str, Path],
                 output_dir: Path, wait: bool = True) -> subprocess.Popen: ...

    def run_abfe(self, complex_mae: Path, msj_files: dict[str, Path],
                 output_dir: Path, wait: bool = True) -> subprocess.Popen: ...
```

**`multisim` invocation:**
```bash
$SCHRODINGER/utilities/multisim \
    -JOBNAME <jobname> \
    -m <main.msj> \
    -i <input.mae> \
    -o <output.cms> \
    -HOST localhost \
    -maxjob 4 \
    -WAIT
```

---

#### Task 13 — HPC job runner (SLURM/PBS)

Submit FEP jobs to a cluster; poll status; collect outputs.

**File:** `openfep/runner/hpc.py`

**API:**
```python
class HPCRunner:
    def __init__(self, schrodinger_dir: Path,
                 scheduler: str = "slurm",   # or "pbs"
                 queue: str = "gpu",
                 gpu_count: int = 1): ...

    def submit(self, msj: Path, input_mae: Path, output_dir: Path) -> str:
        """Submit job; return job ID."""

    def status(self, job_id: str) -> str:
        """Return 'running' | 'done' | 'failed'."""

    def wait(self, job_ids: list[str], poll_interval: int = 60): ...
```

**SLURM wrapper template:**
```bash
#!/bin/bash
#SBATCH --gres=gpu:{gpu_count}
#SBATCH --partition={queue}
$SCHRODINGER/utilities/multisim -JOBNAME {jobname} -m {msj} -i {input} -o {output} -HOST localhost -WAIT
```

---

### Phase 7 — Analysis

#### Task 14 — `_enegrp.dat` parser

Parse Desmond's energy group output files for MBAR input.

**File:** `openfep/analysis/enegrp_parser.py`

**Two line formats from `desmond/ene_utils.py`:**

1. **Inline key=value** (frame-level energies):
   ```
   time=0.000000 en=1.87e+03 pair_elec=-2.34e+02 pair_vdw=1.23e-01
   ```
   Parse with: `re.findall(r'(\w+)=([\d.eE+\-]+)', line)`

2. **Multi-value** (named per-quantity lines):
   ```
   Dispersion_Correction (0.000000) -5.67e+00
   ```
   Parse with: `re.match(r'(\S+)\s+\(([^)]+)\)\s+([\d.eE+\-]+)', line)`

**FEP-relevant fields:**
- Inline: `pair_elec`, `pair_vdw`
- Multi-value: `Dispersion_Correction`, `Self_Energy_Correction`, `Net_Charge_Correction`

**API:**
```python
@dataclass
class EnergySample:
    time: float
    pair_elec: float
    pair_vdw: float
    dispersion_correction: float
    self_energy_correction: float
    net_charge_correction: float

def parse_enegrp(path: Path) -> list[EnergySample]:
    """Parse one _enegrp.dat file; return all samples in time order."""

def load_lambda_series(enegrp_dir: Path,
                       pattern: str = "*_enegrp.dat") -> dict[float, list[EnergySample]]:
    """
    Load all _enegrp.dat files in a directory.
    Keys are lambda values parsed from filename.
    """
```

**Tests:** `tests/test_enegrp_parser.py` with synthetic inline + multi-value fixture lines.

---

#### Task 15 — TI / MBAR estimator

Compute ΔG and uncertainty from energy group data.

**File:** `openfep/analysis/estimator.py`

**Library:** `alchemlyb` (wraps `pymbar`)

**API:**
```python
class MBAREstimator:
    def __init__(self, temperature: float = 298.15): ...

    def fit(self, lambda_series: dict[float, list[EnergySample]]) -> None:
        """Build u_kn matrix; run MBAR via alchemlyb."""

    @property
    def delta_g(self) -> float: ...       # kcal/mol

    @property
    def uncertainty(self) -> float: ...   # kcal/mol (1σ)

    def overlap_matrix(self) -> np.ndarray:
        """MBAR overlap matrix for convergence diagnostics."""
```

**Convergence checks:**
- Overlap matrix diagonal ≥ 0.05 (adjacent windows must overlap)
- Warn if any λ has fewer than 100 samples after decorrelation

---

#### Task 16 — RBFE and ABFE result wrappers

High-level result objects with ΔΔG, cycle closure, and per-leg decomposition.

**Files:**
- `openfep/analysis/rbfe_result.py`
- `openfep/analysis/abfe_result.py`

**RBFE:**
```python
@dataclass
class RBFEResult:
    ligand_a: str
    ligand_b: str
    ddg_complex: float    # ΔG_complex (kcal/mol)
    ddg_solvent: float    # ΔG_solvent (kcal/mol)
    ddg_binding: float    # ddg_complex - ddg_solvent
    uncertainty: float
    overlap_ok: bool

def cycle_closure(results: list[RBFEResult]) -> dict[str, float]:
    """Compute cycle closure errors for each closed triangle in network."""
```

**ABFE:**
```python
@dataclass
class ABFEResult:
    ligand: str
    dg_complex: float     # absolute ΔG, complex leg
    dg_solvent: float     # absolute ΔG, solvent leg
    dg_binding: float     # dg_complex - dg_solvent
    boresch_correction: float   # analytical Boresch free energy correction
    dg_corrected: float   # dg_binding + boresch_correction
    uncertainty: float
```

**Boresch analytical correction:**  
Implement `boresch_correction(restraint: BoreschRestraint, temperature: float) -> float` per Boresch et al. 2003.

**Convergence for ABFE** (not cycle closure — ABFE uses forward/backward decomposition):
```python
def abfe_convergence(result: ABFEResult,
                     forward_fraction: float = 0.5) -> dict:
    """Compare ΔG from first/last half of trajectory."""

def abfe_energy_decomposition(lambda_series: dict[float, list[EnergySample]]) -> dict:
    """
    Approximate breakdown of ΔG by energy term using enegrp.dat fields.
    Groups pair_elec, pair_vdw, and correction terms across all windows.
    NOT true alchemical sub-legs — the complex leg is a single transformation.
    """
```

---

### Phase 8 — CLI

#### Task 17 — CLI entry points

User-facing command-line interface.

**File:** `openfep/cli.py`  
**Entry points in `pyproject.toml`:** `openfep = "openfep.cli:main"`

**Subcommands:**

```
openfep rbfe run   --input complex.mae --output ./rbfe_out --ligands lig_a,lig_b
openfep rbfe analyze --output ./rbfe_out
openfep abfe run   --input complex.mae --output ./abfe_out --ligand lig_name
openfep abfe analyze --output ./abfe_out
openfep network    --input mols.sdf --output network.png
```

**`run` flags shared:**
- `--ff OPLS4` (default; passed to MSJ writer)
- `--n-lambda 12` (RBFE) / `--n-lambda auto` (triggers a3fe adaptive)
- `--hpc slurm` / `--hpc pbs` / omit for local
- `--queue gpu --gpu 1`

**`analyze` output:**
- `results.csv` with ΔG / ΔΔG / uncertainty per pair
- `overlap.png` per lambda window
- `convergence.png` (forward vs backward)
- Prints summary table to stdout

---

## Key Correctness Constraints

These must be verified in code review before merging any MSJ generation code:

| Constraint | Source |
|---|---|
| Production stage named `lambda_hopping` | `subjob_msj.py:_set_lambda_hopping()` |
| `polarization_restraints` only when OPLS4 | `subjob_msj.py:1232-1235` |
| equil stages: `full` except last = `decay` | `subjob_msj.py:1232-1235` |
| ABFE complex leg `assign_custom_charge.mode = keep` | `absolute_binding.py:328` |
| `FepAbsoluteBindingFepPrimer` in md.msj (FINAL stage) | `absolute_binding.py:231` |
| ABFE complex.msj starts with `load_restraints_from_structure` | `absolute_binding.py:364` |
| HMR timesteps `[0.004, 0.004, 0.008]` + migration=0.024 | `msj_constants.py` |
| `energy_group` always-on (openfep deviation) | deviation from Schrodinger default |
| GCMC `scale_solvent_vdw=0.75` | `msj_constants.py` |
| Charged ligand RBFE: `make_alchemical_water=True` | `msj_constants.py` |
| Charged ligand ABFE: `add_alchemical_ions=True` | `absolute_binding.py` |
| Buffer: complex=5.0, charged_complex=8.0, solvent=10.0, vacuum=100.0 | `msj_constants.py` |
| Backbone ASL includes nucleic_acids | `absolute_binding.py:backbone_asl` |
| ABFE MD leg: RECEPTOR+SOLVENT+MEMBRANE+LIGAND | `keep_struc_tags.py:68-72` |
| ABFE complex leg: RECEPTOR+MEMBRANE+SOLVENT+COMPLEX | `keep_struc_tags.py:73-76` |

---

## Build Order

Tasks must be completed in this order (each row depends on the row above):

```
[1] constants     →  [2] mae_parser
[2] mae_parser    →  [5] rbfe_prep, [6] abfe_prep
[3] mappers       →  [4] network
[4] network       →  [5] rbfe_prep (pair selection)
[7] presets       →  [10] rbfe_writer, [11] abfe_writer
[9] stages        →  [10] rbfe_writer, [11] abfe_writer
[10,11] writers   →  [12] local_runner, [13] hpc_runner
[14] enegrp_parser→  [15] estimator
[15] estimator    →  [16] rbfe_result, abfe_result
[16] results      →  [17] cli
```

Tasks 3, 4, 8 (mappers, network, a3fe) are parallel to the prep/MSJ track and can be developed independently.

---

## Testing Strategy

- **Unit tests** for Tasks 1, 2, 9, 14 — pure parsing/data-structure work, no external deps
- **Integration tests** for Tasks 10, 11 — compare generated MSJ text against known-good fixtures derived from Schrodinger reference output (scrubbed of proprietary content)
- **Mock runner tests** for Tasks 12, 13 — mock `subprocess.Popen` to verify CLI arg construction
- **Regression tests** for Tasks 15, 16 — use synthetic `_enegrp.dat` with known ΔG values

**Fixture directory:** `tests/fixtures/` — store synthetic `.mae`, `.msj`, `_enegrp.dat` samples here.

---

## Dependencies

```toml
[project]
dependencies = [
    "rdkit",           # mol reading, Lomap
    "kartograf",       # 3D atom mapping
    "lomap2",          # perturbation network scoring
    "networkx",        # FEP graph
    "alchemlyb",       # MBAR via pymbar
    "numpy",
    "matplotlib",      # overlap/convergence plots
    "click",           # CLI
]
```

Schrodinger utilities (`multisim`, `fep_mapper`) must be on `$PATH` at runtime but are not Python dependencies.
