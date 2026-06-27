# OpenFEP-Desmond Design Spec
**Date:** 2026-06-25
**Status:** Approved for implementation

---

## 1. Overview

OpenFEP-Desmond is an open-source Python package that drives Schrodinger's Desmond MD engine for Free Energy Perturbation (FEP) calculations, replacing the Schrodinger Python stack (scisol, desmond wrappers, msj_generator) with license-free code. The Desmond engine binary is accessed via the user's existing Schrodinger or D.E. Shaw Research academic license.

**Scope:** RBFE (Relative Binding Free Energy) and ABFE (Absolute Binding Free Energy).

**What this is not:** A reimplementation of Schrodinger's proprietary code. All algorithms are reimplemented from published literature and open-source references (a3fe, OpenFE, Lomap, Kartograf).

---

## 2. Architecture

```
Single complex .mae (all CTs tagged with s_fep_struc_tag)
          │
          ▼
┌─────────────────────────┐
│  1. FEP Graph Builder   │  Kartograf (atom mapping) + Lomap (network)
│     [RBFE only]         │  → NetworkX perturbation graph + atom maps
└──────────┬──────────────┘
           │ network.json (edge list + atom mappings + LOMAP scores)
           ▼
┌─────────────────────────┐
│  2. System Preparer     │  RBFE: dual topology → annotated complex .mae
│                         │  ABFE: Boresch setup → complex .mae (MD + FEP)
└──────────┬──────────────┘
           │ {jobname}_complex.mae with all CTs tagged (RECEPTOR/LIGAND/SOLVENT/...)
           ▼
┌─────────────────────────┐
│  3. Lambda Scheduler    │  Named preset | user-defined list | adaptive (a3fe)
└──────────┬──────────────┘
           │ lambda list per leg
           ▼
┌─────────────────────────┐
│  4. MSJ Generator       │  SEA-format .msj files (no schrodinger.utils.sea dep)
└──────────┬──────────────┘
           │ {jobname}_{leg}.msj
           ▼
┌─────────────────────────┐
│  5. Job Runner          │  Local (subprocess multisim) or HPC (SLURM/PBS)
└──────────┬──────────────┘
           │ {jobname}_enegrp.dat, {jobname}.ene, {jobname}-out.mae
           ▼
┌─────────────────────────┐
│  6. Analysis Pipeline   │  Parse _enegrp.dat → u_nk → pymbar MBAR → ΔG
│                         │  RBFE: cycle closure | ABFE: forward/backward + legs
└─────────────────────────┘
```

**Package layout:**
```
openfep/
├── graph/          # Atom mapping + perturbation network
├── topology/       # Dual topology (RBFE) + ABFE system setup
├── lambda_opt/     # Lambda scheduling: presets, user-defined, adaptive
├── msj/            # MSJ text file generation
├── runner/         # Local + HPC job execution
├── analysis/       # enegrp parsing + MBAR + convergence
├── cli/            # openfep rbfe/abfe subcommands
├── constants.py    # Mirrors desmond/constants.py values (no Schrodinger import)
└── data/           # (no bundled ABFE schedules — read from $MMSHARE/data/desmond/abfep/ at runtime)
```

---

## 3. Component Specifications

### 3.1 FEP Graph Builder (RBFE)

**Inputs:** List of ligand `.mae` files (one ligand per file, before complex preparation)
**Output:** `network.json` — edge list with atom mappings and LOMAP scores

**Primary mappers (in order of preference):**

1. **Kartograf** — 3D geometry-aware per-pair atom mapping
   - Converts `.mae` → RDKit mol (via Schrodinger rdkit bridge or SDF export)
   - `KartografAtomMapper.suggest_mappings(molA, molB)` with `MappingRMSDScorer`
   - Uses shape-based alignment + rule-based filters; superior to pure MCS for ring system changes
   - Produces the atom index map stored in `network.json`

2. **Lomap** — perturbation network topology and edge scoring (runs alongside Kartograf)
   - `lomap.generate_lomap_network(ligands, mappers=[LomapAtomMapper()], scorer=default_lomap_score)`
   - LOMAP score (0–1) drives MST edge weights (higher score = preferred edge = more similar pair)
   - Determines which edges are included in the network; Kartograf provides the atom mapping for each edge

3. **RDKit MCS** — fallback when Kartograf/Lomap unavailable or no 3D coords

**Network construction:**
- Build complete weighted graph from pairwise LOMAP scores
- Extract Minimum Spanning Tree (MST) as base network
- Add `n_redundant_edges` extra edges (ranked by score) for RBFE cycle closure

**Charge detection (auto-selects lambda schedule):**
- Net charge per ligand computed from `.mae` atom properties
- If `abs(net_charge) > 0.5`: flag edge as `charged=True` → downstream uses `charge:24` schedule

**Python API:**
```python
from openfep.graph import FEPGraph

graph = FEPGraph.from_mae_files(
    ligand_files=["lig1.mae", "lig2.mae", "lig3.mae"],  # individual ligand files pre-docking
    mapper="kartograf",       # "kartograf" | "lomap" | "mcs"
    n_redundant_edges=2,
)
graph.build()
graph.to_json("network.json")
# graph.edges → list of Edge(ligA, ligB, atom_map, lomap_score, charged)
```

**Serialisation:** `network.json` stores edge list, atom index mappings, LOMAP scores, and charge flags. No Schrodinger `.fmp` format dependency.

---

### 3.2 System Preparer

**Input format — single annotated complex `.mae`:**

Schrodinger FEP takes a single `.mae` file containing multiple CTs (Coordinate Tables), each tagged with the `s_fep_struc_tag` CT-level property. The `extract_structures` stage in each subjob MSJ reads this one file and keeps only the CTs relevant to that leg.

| `s_fep_struc_tag` value | CT contains |
|-------------------------|-------------|
| `receptor` | Protein / receptor (without ligand) |
| `ligand` | Dual-topology ligand pair (RBFE) or free ligand (ABFE solvent leg) |
| `complex` | Protein-ligand bound complex (ABFE complex leg) |
| `solvent` | Pre-solvated water box CT |
| `membrane` | Membrane bilayer CT (optional) |

**Per-leg CT selection** (from `keep_struc_tags.py`):

| Leg | CTs kept |
|-----|---------|
| RBFE complex | `receptor` + `membrane` + `solvent` + `ligand` |
| RBFE solvent | `ligand` |
| RBFE vacuum | `ligand` |
| ABFE MD | `receptor` + `solvent` + `membrane` + `ligand` |
| ABFE complex | `receptor` + `membrane` + `solvent` + `complex` |
| ABFE solvent | `ligand` |

The `extract_structures` stage handles all splitting at runtime — our MSJ generator does not produce separate per-leg files.

**RBFE — Dual Topology via `fep_mapper` MSJ stage:**

`fep_mapper` is **not** a standalone CLI utility. It is a Desmond MSJ stage (`class FepMapper(cmj.StageBase)`, `desmond/stage/fep_mapper.py:490`) that runs inside `multisim`. Internally it calls `run -FROM scisol fep_mapper.py` as a subprocess. `desmond/msj_generator/common.py:406-423` configures this stage and raises `Exception("No fep_mapper stage in the msj file")` if it is absent.

**openfep approach (Option 1 — include stage in MSJ):** openfep does not produce a pre-merged dual-topology `.mae` before the run. Instead, openfep configures the `fep_mapper` stage inside the generated `main.msj`. The atom mapping from Kartograf is exported to an atom-mapping file and referenced by the stage parameter. The dual-topology merge happens inside `multisim` at runtime using Schrodinger's own `fep_mapper.py` — no reimplementation needed.

**`fep_mapper` stage parameters:**
- `graph_file` — path to the `.fmp` FEP graph file (edge list + structures)
- `atom_mapping` — path to Kartograf-generated atom-index mapping file
- `receptor` — number of environment CTs (1 + bool(solvent) + bool(membrane))
- `align_core_only` — restrict alignment to mapped core atoms (optional)

**Key atom properties set by `fep_mapper.py` on the merged ligand CT:**
- `i_fep_mapping` — atom index cross-reference A→B
- `i_des_atom_domain` — ALCHEM_A=7 (ligA atoms), ALCHEM_B=8 (ligB atoms)
- `i_fep_subst` — substitution group index
- CT-level: `s_fep_struc_tag = "ligand"`

**What openfep generates (RBFE):**
- `{jobname}.fmp` — FEP graph file (edge list, structures; required by `fep_mapper` stage)
- `{jobname}.msj` — main MSJ with `fep_mapper` stage as first stage; leg subjob MSJs follow
- Atom-mapping file per edge — Kartograf output in the format `fep_mapper.py -atom-mapping` expects

```python
from openfep.topology import RBFESystemBuilder

builder = RBFESystemBuilder()
fmp_file, msj_file = builder.prepare(
    complex_mae="receptor_lig1_complex.mae",  # pre-docked complex (receptor + ligA)
    ligB="lig2.mae",                          # second ligand for dual topology
    atom_mapping=edge.atom_map,               # from Kartograf
    output_dir=Path("./systems/edge_lig1_lig2/"),
)
# output: {jobname}.fmp (graph), {jobname}.msj (main), {edge}_atom_mapping.txt
# fep_mapper stage merges structures inside multisim at runtime
```

**ABFE — Boresch Restraints:**

Two legs: bound (protein-ligand complex) and solvent (ligand in water).

**Schrodinger ABFE is a three-phase workflow — all three phases must be reproduced:**

```
Phase 1: MD pre-equilibration (md.msj)
  → equilibrates the complex, generates trajectory for restraint geometry determination
  → backbone restraints: 50 kcal/mol/Å² on "((protein and backbone) or (nucleic_acids and nucleic_backbone)) and not a.ele H"
  → trajectory saved every 3.6 ps
  → duration: md_sim_time (default 1000 ps)

Phase 2: FEP primer (FepAbsoluteBindingFepPrimer — FINAL stage of md.msj, NOT complex.msj)
  → placed at end of md.msj by _set_fep_absolute_binding_fep_primer() (absolute_binding.py:231)
  → reads MD trajectory output to auto-select Boresch anchor atoms
    criteria: geometric stability, min_angle=45°, low RMSD variation in MD trajectory
  → sets pose_conf_restraint: {enable: True, name: "soft", sigma: 0.0, alpha: 1.0, fc: 0.0}
  → writes restraint geometry into the md.msj output structure consumed by complex.msj

Phase 3: FEP alchemical simulation (complex.msj + solvent.msj)
  → load_restraints_from_structure stage inserted after assign_forcefield
  → alchemical lambda windows (complex=68, solvent=60 by default)
```

Restraint type: Boresch (1 bond r, 2 angles θ_A/θ_B, 3 dihedrals φ_A/φ_B/φ_C).

Analytical restraint free energy correction applied at analysis stage:

```
ΔG_restraint = -RT × ln(
    8π²V₀ / (r²_aA0 × sin(θ_A0) × sin(θ_B0))
    × (σ_r × σ_θA × σ_θB × σ_φA × σ_φB × σ_φC)
)
```

where V₀ = 1661 Å³ (standard state), equilibrium values (r_aA0, θ_A0, θ_B0) are extracted from the prepared structure, and σ values are the restraint force constant widths.

```python
from openfep.topology import ABFESystemBuilder

builder = ABFESystemBuilder()
system = builder.prepare(
    complex_mae="receptor_lig1_complex.mae",  # pre-docked complex (receptor + ligand)
    restraint_type="boresch",                 # or "flat_bottom"
    box_padding=10.0,                         # Å
    md_sim_time=1000,                         # ps — Phase 1 MD duration
    min_angle=45,                             # degrees — Boresch anchor geometry filter
    ligand_asl="a.i_fep_absolute_binding_ligand 1",
    receptor_asl="protein or nucleic_acids",
    # restraint_atoms: None = auto-selected from MD trajectory
    restraint_atoms=None,
)
# output: {jobname}_abfe.mae — single file containing:
#   CT tagged s_fep_struc_tag="receptor"
#   CT tagged s_fep_struc_tag="complex"  (bound pose, for complex leg)
#   CT tagged s_fep_struc_tag="ligand"   (free ligand, for solvent leg)
#   CT tagged s_fep_struc_tag="solvent"
# + {jobname}_restraints.json: atom indices, equilibrium r/θ/φ, force constants
```

---

### 3.3 Lambda Scheduler

Three modes sharing a common interface:

**Mode 1 — Named preset** (mirrors `desmond/constants.py` exactly):

| Preset | Windows | When to use |
|--------|---------|-------------|
| `default:12` | 12 | Standard perturbations, neutral ligands |
| `flexible:16` | 16 | Scaffold hops / core-hopping |
| `charge:24` | 24 | Charged ligands (auto-selected when `charged=True`) |

**ABFE lambda window counts** (from `fep_schedule.py:282-301`):

| Protocol | Restrained? | Complex windows | Solvent windows |
|----------|-------------|-----------------|-----------------|
| Default (neutral) | No  | 68  | 60 |
| Default (neutral) | Yes | 80  | 68 |
| Charged           | No  | 108 | 60 |
| Charged           | Yes | 128 | 68 |

**ABFE schedule source (architectural decision):** Schrodinger reads full Gibbs weight vectors from MSJ files at `$MMSHARE/data/desmond/abfep/{leg}_schedule{_restrained}{_chg}.msj` (SEA format: `.schedule.weights.N.val` + `.flexible[0].A.val` fields), via `_get_abfep_schedule()` in `subjob_msj.py:176-191`. openfep must read these files from the user's Schrodinger installation at runtime — users need Desmond installed to run FEP regardless, so `$MMSHARE` is always available. These are proprietary files and must not be redistributed. The `openfep/data/abfe/` directory is therefore not needed.

**Mode 2 — User-defined list:**
```python
lambda_schedule = [0.0, 0.05, 0.10, 0.20, 0.35, 0.50, 0.65, 0.80, 0.90, 0.95, 1.0]
```

**Mode 3 — Adaptive (a3fe algorithm):**

```python
from openfep.lambda_opt import AdaptiveLambdaConfig

adaptive_config = AdaptiveLambdaConfig(
    pilot_time=500.0,           # ps per window in pilot run
    n_pilot_windows=6,          # coarse initial λ spacing
    gradient_threshold=0.5,     # kcal/mol/ns — equilibration criterion
    runtime_constant=0.001,     # C in time allocation formula
    max_windows=24,             # cap on total windows
    max_runtime=20000.0,        # ps — per-window cap
)
```

**Adaptive algorithm (three phases):**

*Phase 1 — Pilot run:*
- Run short simulations (`pilot_time` ps) at `n_pilot_windows` evenly-spaced λ values
- Collect `dG/dλ` gradient timeseries per window

*Phase 2 — Optimal window placement:*
- Compute variance of `dG/dλ` at each pilot λ value
- Integrate `√Var[dG/dλ]` across λ-space using 1D quadrature
- Place `max_windows` new λ values such that the cumulative integral is equidistributed between adjacent pairs (equal uncertainty contribution per window)
- Implemented in `AdaptiveLambdaOptimizer.get_optimal_lam_vals()`

*Phase 3 — Adaptive time allocation:*
- Per-window production time: `t_optimal_k = √(t_current_k / C) × σ(ΔF_k)`
- Windows with higher gradient variance receive longer runs
- Equilibration detection via block gradient method: window is equilibrated when rolling `|dG/dt| < gradient_threshold` (kcal/mol/ns)
- Data before equilibration is discarded before MBAR

---

### 3.4 MSJ Generator

Generates Desmond MSJ files as plain SEA-format text. No dependency on `schrodinger.utils.sea` — files are written as formatted strings and read by `multisim` at runtime.

**Stage sequence per leg (RBFE):**
```
task → extract_structures → assign_forcefield → build_geometry
     → simulate [NVT equil, T=10K,   polarization_restraints=full]
     → simulate [NPT equil, T=300K,  polarization_restraints=decay]   ← last equil
     → lambda_hopping [production, ensemble per leg]                   ← Drude fully active
     → trim
     → fep_analysis
```

Note: Schrodinger renames the production `simulate` stage to `lambda_hopping` via `_set_lambda_hopping()`. This name change activates Desmond's lambda-hopping code path. All prior stages remain as `simulate`.

**ABFE MSJ files (4 files per ligand, not 2):**
```
{jobname}.msj            — launcher that calls the other three
{jobname}_md.msj         — Phase 1+2: MD pre-equilibration + FepAbsoluteBindingFepPrimer (final stage)
{jobname}_complex.msj    — Phase 3: load_restraints_from_structure → alchemical complex leg
{jobname}_solvent.msj    — Phase 3: alchemical solvent leg
```

**Stage sequence for ABFE md.msj (phases 1+2):**
```
task → extract_structures → assign_forcefield → build_geometry
     → simulate [backbone-restrained MD, 1000 ps, NPT, traj every 3.6 ps]
     → FepAbsoluteBindingFepPrimer    ← LAST stage of md.msj; reads trajectory;
                                         auto-selects Boresch anchors; writes
                                         restraint geometry into output structure
```

**Stage sequence for ABFE complex.msj (phase 3):**
```
task → extract_structures → assign_forcefield
     → load_restraints_from_structure          ← reads geometry from md.msj output
     → build_geometry
     → simulate [NVT equil,  polarization_restraints=full]
     → simulate [NPT equil,  polarization_restraints=decay]   ← last equil
     → lambda_hopping [production, muVT, 68 windows]          ← Drude fully active
     → trim
     → fep_analysis
```

**Stage sequence for ABFE solvent MSJ:**
```
task → extract_structures → assign_forcefield → build_geometry
     → simulate [NVT equil,  polarization_restraints=full]
     → simulate [NPT equil,  polarization_restraints=decay]   ← last equil
     → lambda_hopping [production, NPT, 60 windows]
     → trim
     → fep_analysis
```

> **REST scaling subtlety (solvent leg, restrained protocol):** When Boresch restraints are enabled, `_get_abfep_solvent_restrained_lambda_windows()` counts the number of non-zero `flexible[0].A.val` entries in the schedule file and returns that count as `end_win`. `_set_solute_tempering_temperature_ladder()` is then called with `start_win=0, end_win=end_win` — this omits REST scaling on the restrained λ-windows of the solvent leg (where the ligand is still restrained to its bound pose). When `restraint_enabled=False`, `end_win=0` so no REST scaling is applied. Source: `absolute_binding.py:298-302` and `subjob_msj.py:194-200`.

**Ensemble per leg** (from `desmond/msj_generator/common.py` logic):
- Complex leg: `muVT` (default) — GCMC grand-canonical ensemble
- Solvent leg: `NPT`
- Vacuum leg: `NVT`

**muVT / GCMC configuration** (complex leg only) — set via `_set_gcmc()`:
```
gcmc {
  scale_solvent_vdw = 0.75    # softens water vdW for insertion/deletion
}
```
`muVT` allows water molecules to enter/exit the binding site during alchemical annihilation, which is important for correct thermodynamics when buried waters are displaced by the ligand.

**Required MSJ parameters** — all must be set or Desmond will error or produce wrong results:

| Parameter | Where set | Value |
|-----------|-----------|-------|
| `fep.type` | `task.set_family.simulate.fep.type` | `"small_molecule"` (RBFE) or `"absolute_binding"` (ABFE) |
| `fep_analysis.fep_type` | `fep_analysis` stage | same as above |
| `backend.is_for_fep` | every `lambda_hopping` stage | `True` |
| `forcefield` | `assign_forcefield` | `"OPLS4"` |
| `hydrogen_mass_repartition` | `assign_forcefield` | `True` (HMR on by default) |
| `polarization_restraints` | all equil stages except last | `"full"` — **OPLS4 (F17) only**; gated on `opls_version == OPLSVersion.F17` (`subjob_msj.py:1232`); omit entirely for OPLS3e (F16) or OPLS_2005 (F14) |
| `polarization_restraints` | last equil stage only | `"decay"` — **OPLS4 (F17) only**; same gate as above |
| `timestep` | `lambda_hopping` stage | `[0.004, 0.004, 0.008]` ps (4 fs outer step with HMR) |
| `backend.migration.interval` | `lambda_hopping` stage | `0.024` (with HMR) |
| `randomize_velocity.seed` | `lambda_hopping` stage | `2014` |
| `energy_group` | `lambda_hopping` stage | `{name="$JOBNAME$[_replica$REPLICA$]_enegrp.dat" first=0.0 interval=1.2}` — **deviation from Schrodinger**: Schrodinger only emits energy_group on complex legs by default; openfep always enables it on all legs (required for MBAR on all legs) |
| `assign_custom_charge.mode` | subjob MSJ | **RBFE:** `"assign"` for complex leg; `"keep"` for solvent leg. **ABFE:** `"keep"` for ALL legs (md, complex, solvent) — source: `absolute_binding.py:328` uses `CUSTOM_CHARGE_MODE.KEEP` unconditionally for ABFE |

**Buffer widths per leg** (from `msj_constants.py` — these are minimums; user can increase):

| Leg | Default buffer | Charged override |
|-----|---------------|-----------------|
| Complex | 5.0 Å | 8.0 Å |
| Solvent | **10.0 Å** | 10.0 Å |
| Vacuum | **100.0 Å** | 100.0 Å |

**Charged ligand handling** — triggers when `abs(net_charge) > 0.5`:
- `build_geometry.neutralize_system = True`
- `build_geometry.salt = {negative_ion: Cl, positive_ion: Na}`
- `build_geometry.salt.concentration = max(0.15, user_salt_conc)` M
- `build_geometry.box_shape = "cubic"` (charged systems must use cubic box)
- RBFE: `assign_forcefield.make_alchemical_water = True`; lambda schedule → `charge:24`
- ABFE: `assign_forcefield.add_alchemical_ions = True`; the `charge:24` RBFE preset does **not** apply — instead the ABFE charged Gibbs schedule is loaded from `{leg}_schedule_chg.msj` (or `{leg}_schedule_restrained_chg.msj`) at `$MMSHARE/data/desmond/abfep/`; window counts become 108/60 (unrestrained) or 128/68 (restrained) — see Lambda Scheduler table above

**Lambda syntax** (from `subjob_msj.py` `_set_task_num_windows_and_schedule`):
```
task {
  set_family = {
    simulate { fep.lambda = "default:12" }
  }
}
```
For ABFE, the full Gibbs schedule is read at runtime from `$MMSHARE/data/desmond/abfep/{leg}_schedule{_restrained}{_chg}.msj` and injected into `task.set_family.desmond.backend.force.term.gibbs` via `_set_custom_lambda_schedule()`. openfep must parse this SEA-format file without Schrodinger's `sea` library — use the plain-text SEA parser (same approach used in the MSJ generator).

For user-defined or adaptive schedules, the explicit λ list is written directly into the MSJ.

**File naming** (mirrors `get_msj_filename()`):
- Main: `{jobname}.msj`
- RBFE per-leg: `{jobname}_{leg_type}.msj`  (leg_type: `complex`, `solvent`, `vacuum`)
- ABFE: `{jobname}.msj`, `{jobname}_md.msj`, `{jobname}_complex.msj`, `{jobname}_solvent.msj`
- Output structure: `{jobname}_{leg_type}-out.mae`

**Constants mirrored from `desmond/constants.py` (no Schrodinger import):**
- Forcefield: `OPLS4` (`OPLSVersion.F17`; Schrodinger ships F16/OPLS3e as their default, but user chose OPLS4)
- Water models: `SPC` (default), `TIP3P`, `TIP4P`, `TIP4PEW`, etc.
- Simulation times: production=5000 ps, equilibration=20 ps
- Random seed: 2014
- Buffer widths: complex=5.0 Å, solvent=10.0 Å, vacuum=100.0 Å
- HMR timestep: [0.004, 0.004, 0.008] ps, migration interval=0.024
- ENERGY_GROUP: `name="$JOBNAME$[_replica$REPLICA$]_enegrp.dat"`, first=0.0, interval=1.2

**Python API:**
```python
from openfep.msj import RBFEMSJWriter, ABFEMSJWriter

writer = RBFEMSJWriter(
    jobname="lig1_lig2",
    lambda_schedule="default:12",   # str preset | list | "adaptive"
    sim_time=5000.0,                # ps
    water_model="SPC",
    forcefield="OPLS4",
    hmr=True,                       # hydrogen mass repartitioning
    rand_seed=2014,
    salt_concentration=0.0,         # M; raised to 0.15 for charged ligands
    charged=False,                  # auto-set from net_charge detection
)
writer.write_complex_msj()   # → lig1_lig2_complex.msj
writer.write_solvent_msj()   # → lig1_lig2_solvent.msj
writer.preview()             # prints MSJ without writing

# ABFE — generates all 4 files
abfe_writer = ABFEMSJWriter(
    jobname="lig1_abfe",
    md_sim_time=1000,               # ps — Phase 1 MD pre-equilibration
    sim_time=5000.0,                # ps — Phase 3 FEP production
    forcefield="OPLS4",
    hmr=True,
    rand_seed=2014,
)
abfe_writer.write_all()  # → lig1_abfe.msj + lig1_abfe_md.msj + lig1_abfe_complex.msj + lig1_abfe_solvent.msj
```

---

### 3.5 Job Runner

Invokes `$SCHRODINGER/utilities/multisim`. Two backends share a common `BaseRunner` interface.

**RBFE submission command (per-leg):**
```
$SCHRODINGER/utilities/multisim \
  -m {jobname}_{leg}.msj \
  -i {input}.mae \
  -JOBNAME {jobname}_{leg} \
  -HOST {host}:{gpus} \
  -maxjob 1 -cpu {cpus}
```

**ABFE submission — single command from main.msj (orchestrates all 3 subjobs):**
```
$SCHRODINGER/utilities/multisim \
  -m {jobname}.msj \
  -i {jobname}_abfe.mae \
  -JOBNAME {jobname} \
  -HOST {host}:{gpus} \
  -maxjob 1 -cpu {cpus}
```
`main.msj` orchestrates the sequence: MD → complex FEP → solvent FEP. `multisim` handles dependencies between the subjobs.

**LocalRunner:**
- Invokes `multisim` directly via `subprocess.run`
- Monitors job by polling log file for completion markers
- Raises `JobFailedError` on non-zero exit or missing output files

**SlurmRunner / PBSRunner:**
- Writes a `.sh` batch script with appropriate `#SBATCH` / `#PBS` headers wrapping the same `multisim` command
- Submits via `sbatch` / `qsub`; polls via `squeue -j {job_id}` / `qstat {job_id}`

**Output contract** — runner verifies these exist on completion:
- `{jobname}_{leg}_enegrp.dat` (per lambda replica) — FEP pair energy groups (MBAR input)
  - Pattern from `msj_constants.py`: `$JOBNAME$[_replica$REPLICA$]_enegrp.dat`
  - Only produced if `energy_group` block is present in `lambda_hopping` stage
- `{jobname}_{leg}.ene` — energy sequence file (QC/convergence)
- `{jobname}_{leg}-out.mae` — output structure
- Intermediate `.cms.gz` and `.tgz` files are deleted by the `trim` stage to save disk

**Checkpoint/restart:**
- Before submitting, runner checks for existing `_enegrp.dat` partial output
- If found and `allow_restart=True`, extends the run rather than restarting from scratch
- Job state serialised to `{jobname}_runner_state.json`

**Python API:**
```python
from openfep.runner import LocalRunner, SlurmRunner

runner = LocalRunner(schrodinger_path="/opt/schrodinger2024")

# RBFE: submit complex and solvent legs in parallel — both read the same tagged .mae
# multisim passes the .mae to Desmond; extract_structures splits the CTs per leg
job_complex = runner.submit(msj="lig1_lig2_complex.msj",
                            input_mae="lig1_lig2_fep.mae",   # single file, all CTs tagged
                            jobname="lig1_lig2_complex", gpus=1)
job_solvent = runner.submit(msj="lig1_lig2_solvent.msj",
                            input_mae="lig1_lig2_fep.mae",   # same input file
                            jobname="lig1_lig2_solvent", gpus=1)

# ABFE: submit via main.msj (orchestrates MD → FEP sequence)
job_abfe = runner.submit(msj="lig1_abfe.msj",
                         input_mae="lig1_abfe.mae",          # single file, all CTs tagged
                         jobname="lig1_abfe", gpus=1)

job_complex.wait()
outputs = job_complex.collect_outputs()
# outputs.enegrp (list of per-replica files), outputs.ene, outputs.mae
```

---

### 3.6 Analysis Pipeline

**Stage A — enegrp Parser:**

Reads `_enegrp.dat` from each lambda window. Reimplements `ene_utils.parse_enegrp_file()` with no Schrodinger dependency.

Energy terms included in reduced potential (from `EneGrpProp` in `ene_utils.py`):
- `pair_elec` — electrostatic alchemical cross-interaction
- `pair_vdw` — vdW alchemical cross-interaction
- `Dispersion_Correction` — long-range vdW correction
- `Self_Energy_Correction` — PME self-energy correction
- `Net_Charge_Correction` — charge correction (non-zero for charged ligands)

Reduced potential: `u(x, λ_j) = β × (pair_elec_j + pair_vdw_j + DISP_CORR_j + ENERGY_CORR_j + CHARGE_CORR_j)`

Output: `u_nk` DataFrame — index `(time, lambda)`, columns = all lambda states. Compatible with `pymbar.MBAR`.

**Stage B — Free Energy Estimator:**

```python
from openfep.analysis.estimator import FEPEstimator

est = FEPEstimator(method="MBAR")   # or "BAR"
# Decorrelation first
from alchemlyb.preprocessing import statistical_inefficiency
u_nk_sub = statistical_inefficiency(u_nk_production)
result = est.compute(u_nk_sub)
# result.dG, result.uncertainty, result.overlap (lambda overlap matrix)
```

**RBFE ΔΔG:**
```
ΔΔG = ΔG_complex - ΔG_solvent
uncertainty = √(σ²_complex + σ²_solvent)
```

**ABFE ΔG:**

Schrodinger's ABFE complex leg annihilates the ligand (elec + vdW) in the **restrained** bound state — restraints stay on for the entire complex leg. The Boresch analytical correction accounts for the free energy cost of imposing those restraints on the free ligand in solution.

```
ΔG_binding = ΔG_complex      (MBAR ΔG from complex leg:
                               ligand annihilated in bound state, restraints on)
           - ΔG_solvent       (MBAR ΔG from solvent leg:
                               ligand annihilated in water, unrestrained)
           + ΔG_restraint     (Boresch analytical correction, Boresch et al. 2003:
                               free energy of imposing restraints on free ligand)
```

> Note: the Schrodinger ABFE complex leg does NOT decompose into separate restraint-on/elec/vdW sub-legs. It is a single alchemical transformation from restrained bound ligand → nothing. The decomposition shown in some FEP+ documentation refers to older multi-stage protocols; the current single-lambda implementation uses MBAR across all windows simultaneously. Cross-check the actual FEP+ ABFE protocol paper before implementing the correction term sign convention.

**Convergence checks:**

| Check | RBFE | ABFE |
|-------|------|------|
| Lambda overlap matrix | Yes | Yes |
| Cycle closure (`\|ΔG_AB + ΔG_BC + ΔG_CA\|`) | Yes | No |
| Forward/backward time-series convergence | Optional | Yes — critical |
| Energy-term decomposition (approx.) | No | Yes — via `enegrp.dat` `pair_elec`/`pair_vdw`/correction terms, not true alchemical sub-legs |
| Block gradient equilibration detection | Yes (adaptive) | Yes |

**ABFE convergence API:**
```python
from openfep.analysis import ABFEConvergence

conv = ABFEConvergence(u_nk_bound, u_nk_solvent)
conv.overlap_matrix()       # per-leg lambda overlap heatmap
conv.forward_backward()     # ΔG vs simulation time, forward & reverse
conv.energy_decomposition()  # approximate elec/vdW/correction breakdown from enegrp.dat
                             # (NOT true alchemical sub-legs — the complex leg is a single
                             #  transformation; decomposition groups windows by dominant term)
```

---

### 3.7 CLI

Entry point: `openfep` (defined in `pyproject.toml` `[project.scripts]`).

**RBFE:**
```bash
# Step 1: build perturbation graph from individual ligand mae files
openfep rbfe graph    --ligands lig*.mae --output network.json

# Step 2: prepare dual-topology complex .mae (one per edge in the graph)
#   Input: pre-docked complex .mae (receptor + ligA) + ligB .mae per edge
openfep rbfe prepare  --graph network.json \
                      --complexes "lig1_complex.mae,lig2_complex.mae,..." \
                      --output-dir systems/
#   Output: systems/{edgename}_fep.mae — single file, all CTs tagged by s_fep_struc_tag

# Step 3: generate MSJ files (reads the tagged .mae, splits at runtime via extract_structures)
openfep rbfe msj      --graph network.json --systems-dir systems/ \
                      --lambda-schedule default:12 --sim-time 5000 --output-dir msj/

openfep rbfe run      --msj-dir msj/ --runner slurm --partition gpu --walltime 24:00:00
openfep rbfe analyze  --graph network.json --results-dir results/ --output rbfe_results.csv
```

**ABFE:**
```bash
# Step 1: prepare tagged complex .mae (receptor + complex CT + ligand CT)
#   Input: pre-docked complex .mae (receptor + ligand already docked)
openfep abfe prepare  --complex receptor_lig_complex.mae \
                      --restraint boresch --output-dir systems/
#   Output: systems/{ligname}_abfe.mae — single file with all CTs tagged

# Step 2: generate 4 MSJ files (main + md + complex + solvent)
openfep abfe msj      --systems-dir systems/ --lambda-schedule default \
                      --sim-time 5000 --output-dir msj/

openfep abfe run      --msj-dir msj/ --runner slurm --partition gpu
openfep abfe analyze  --systems-dir systems/ --results-dir results/ \
                      --output abfe_results.csv
```

**Lambda schedule options (all commands accepting `--lambda-schedule`):**
```bash
--lambda-schedule default:12
--lambda-schedule flexible:16
--lambda-schedule charge:24
--lambda-schedule "0.0,0.1,0.3,0.5,0.7,0.9,1.0"   # explicit CSV list
--lambda-schedule adaptive --pilot-time 500 --max-windows 24  # a3fe method
```

**Output CSV columns:**
- RBFE: `ligand_A, ligand_B, ddG_kcal_mol, uncertainty, cycle_closure_error`
- ABFE: `ligand, dG_kcal_mol, uncertainty, overlap_quality, forward_backward_hysteresis`

---

## 4. Key Design Decisions

### Grounded in Schrodinger source (no proprietary code copied)

| Decision | Source file | Rationale |
|----------|-------------|-----------|
| MSJ as SEA plain text | `desmond/multisim/parser.py` | Generate as strings; `multisim` parses at runtime |
| Lambda syntax `"default:12"` | `desmond/msj_generator/subjob_msj.py` | Exact syntax Desmond expects |
| Default λ counts 12/16/24 | `desmond/constants.py DefaultFEPSimulationTimes` | Matches Schrodinger defaults |
| ABFE lambda schedules (68/60, 80/68, 108/60, 128/68) | `desmond/fep_schedule.py:282-301`, `desmond/msj_generator/subjob_msj.py:176-191` | Read from `$MMSHARE/data/desmond/abfep/{leg}_schedule{_restrained}{_chg}.msj` at runtime; window counts derived from schedule file, not hardcoded |
| Energy files: `_enegrp.dat` + `.ene` | `desmond/config_utils.py`, `ene_utils.py` | Correct Desmond output format (not `.edr`) |
| File naming `{jobname}_{leg}.msj` | `desmond/util.py::get_msj_filename()` | Matches what `multisim` expects |
| ABFE 4 MSJ files (main + md + complex + solvent) | `desmond/msj_generator/fep/absolute_binding.py` | Schrodinger ABFE is a 3-phase workflow |
| Ensemble: muVT/NPT/NVT per leg | `desmond/msj_generator/common.py` L152-186 | Matches Schrodinger runtime behaviour |
| MBAR terms: pair_elec + pair_vdw + corrections | `desmond/ene_utils.py EneGrpProp` | Includes all correction terms |
| Charged protocol 24 windows | `desmond/constants.py SIMULATION_PROTOCOL.CHARGED` | Required for net-charged ligands |
| Charged: cubic box + salt 0.15M min + alchemical ions | `desmond/msj_generator/subjob_msj.py _set_net_charge_box()` | Charged systems need cubic box and neutralization |
| Forcefield: `OPLS4` (`OPLSVersion.F17` in MSJ = `forcefield = OPLS4`) | `scisol-v6.1/.../prepare_reinit.py` `forcefield = OPLS4`, `desmond/msj_generator/subjob_msj.py` F17 branch | F17=OPLS4, F16=OPLS3e (Schrodinger default), F14=OPLS_2005 |
| OPLS4 polarization restraints on equil stages | `desmond/msj_generator/subjob_msj.py` L1232-1235 | Drude FF susceptible to polarization catastrophe; all-but-last equil=`full`, last equil=`decay` |
| Input: single complex `.mae` with `s_fep_struc_tag` CT property | `desmond/constants.py FEP_STRUC_TAG`, `keep_struc_tags.py` | `extract_structures` stage splits legs at runtime; no separate per-leg input files |
| Production stage named `lambda_hopping` not `simulate` | `desmond/msj_generator/subjob_msj.py _set_lambda_hopping()` | Different Desmond code path; required for FEP |
| HMR: `hydrogen_mass_repartition=True`, timestep=[0.004, 0.004, 0.008] | `desmond/msj_constants.py _HMR_TIMESTEPS`, `_set_hmr()` | Default on; enables 4 fs outer timestep |
| HMR migration interval: 0.024 | `desmond/msj_constants.py _HMR_MIGRATION_INTERVAL` | Required when HMR enabled |
| Random seed: 2014 | `desmond/msj_generator/workflow_params.py rand_seed` | Reproducibility default |
| `energy_group` block required in `lambda_hopping` | `desmond/msj_constants.py ENERGY_GROUP` | Must be explicit; OFF by default in `compute_energy_groups` |
| `fep.type` in task + `fep_analysis.fep_type` | `desmond/msj_generator/subjob_msj.py _set_fep_type()` | Tells Desmond which alchemical protocol to run |
| `backend.is_for_fep = True` on every production stage | `desmond/msj_generator/subjob_msj.py _set_is_for_fep()` | Required for Desmond FEP licence check |
| Buffer: complex=5.0Å, solvent=10.0Å, vacuum=100.0Å | `desmond/msj_constants.py` | Leg-specific minimum widths, not a single value |
| ABFE MD pre-run: 1000 ps backbone-restrained | `desmond/msj_generator/fep/absolute_binding.py` `_DEFAULT_MD_RESTRAINT` | Phase 1: equilibration + restraint geometry |
| `load_restraints_from_structure` stage (ABFE) | `desmond/msj_generator/fep/absolute_binding.py _insert_stage_load_restraints_from_structure()` | Loads Boresch geometry from Phase 1 output |
| `assign_custom_charge.mode`: RBFE complex=`assign`, all others=`keep` (RBFE solvent, ABFE all legs) | `desmond/msj_generator/subjob_msj.py _set_assign_custom_charge()`, `absolute_binding.py:328 CUSTOM_CHARGE_MODE.KEEP` | RBFE: charges assigned on complex leg, inherited by solvent. ABFE: charges always kept (all 3 legs) |
| `trim` stage after production | `desmond/msj_constants.py _TRIM_TEMPLATE` | Removes intermediate cms.gz + tgz files |

### From open-source review

| Decision | Source | Rationale |
|----------|--------|-----------|
| Kartograf as primary mapper | OpenFE/kartograf | 3D geometry-aware, avoids bad MCS mappings |
| Lomap for network scoring | OpenFE/Lomap | Purpose-built LOMAP score for FEP perturbation networks |
| Adaptive λ algorithm | a3fe `stage.py` | Equidistributed root-variance placement + proportional time allocation |
| Block gradient equilibration | a3fe `lambda_window.py` | 0.5 kcal/mol/ns threshold, rolling average |
| Boresch restraints for ABFE | Literature (Boresch 2003) + a3fe | Standard, analytically correctable |
| PB charge correction flag | FEP-SPell-ABFE | Charged ligands need explicit handling |

---

## 5. Dependencies

```toml
[project]
name = "openfep"
requires-python = ">=3.10"

[project.dependencies]
rdkit = ">=2023.9"
networkx = ">=3.0"
numpy = ">=1.24"
pandas = ">=2.0"
pymbar = ">=4.0"
alchemlyb = ">=2.3"
kartograf = ">=0.4"
lomap3 = ">=3.0"
click = ">=8.0"

[project.optional-dependencies]
# No extra optional deps currently; adaptive λ algorithm is reimplemented internally

[project.scripts]
openfep = "openfep.cli:main"
```

**Runtime requirement (not a Python dependency):**
- `$SCHRODINGER` environment variable pointing to a Schrodinger 2023+ installation with Desmond and `multisim`

---

## 6. Known Limitations & Future Work

| Limitation | Impact | Future fix |
|-----------|--------|-----------|
| No ensemble/replica runs | Single replica per window; less robust uncertainty | v2: `n_replicas` param, Gelman-Rubin convergence |
| No core-hopping protocol | Large scaffold changes may not converge | v2: detect ring-size change, auto-select `flexible:16` + soft-bond alpha |
| No APBS PB charge correction | Charged ligands rely on counter-ion approach only | v2: optional APBS integration |
| No provenance tracking | Results not reproducible from JSON alone | v2: gufe-compatible data model |
| Adaptive λ requires pilot run GPU time | Additional cost before production | Accept as design trade-off |
| ABFE restraint atom selection defaults to auto | Auto-selection may pick suboptimal anchors for flexible loops | v2: validate anchor stability during pilot run |
