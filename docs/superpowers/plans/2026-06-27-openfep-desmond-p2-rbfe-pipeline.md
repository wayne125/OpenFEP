# openfep-Desmond — Sub-plan 2: RBFE Pipeline

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Deliver a complete RBFE calculation pipeline: lambda presets, atom mapping, perturbation network, RBFE system preparation, and MSJ file generation. After this sub-plan, `openfep` can write all MSJ files needed to run an RBFE edge on Desmond.

**Architecture:** Five tasks in dependency order. Tasks 1-3 provide scheduling and mapping primitives; Task 4 uses them to prepare the system; Task 5 assembles final MSJ files using stage builders from Sub-plan 1.

**Tech Stack:** RDKit (mol manipulation, MAE→Mol conversion), Kartograf (3D atom mapping), Lomap2 (MCS scoring), NetworkX (FEP graph), existing SEABlock stage builders from `openfep/msj/sea.py`.

**Prerequisite:** Sub-plan 1 complete — `openfep/msj/sea.py` contains `SEABlock`, `render`, `render_msj`, `task_stage`, `assign_forcefield_stage`, `build_geometry_stage`, `simulate_stage`, `lambda_hopping_stage`, `trim_stage`, `gcmc_stage`. `openfep/constants.py` and `openfep/mae_parser.py` exist and pass all tests.

## Global Constraints

- No `from schrodinger` anywhere in `openfep/` — no Schrodinger Python stack.
- Production stage named `lambda_hopping`, never `simulate`.
- `polarization_restraints` only written when `polarization_restraint` param is not None (caller gates on OPLS4).
- All equil stages except last: `polarization_restraint="full"`. Last equil stage: `polarization_restraint="decay"`.
- RBFE complex leg: `assign_custom_charge.mode = "assign"`. Solvent leg: `mode = "keep"`.
- `energy_group` always written in `lambda_hopping_stage` (openfep deviation — no condition).
- HMR always on: `timestep = [0.004 0.004 0.008]`, `migration_interval = 0.024`.
- Working directory: `C:\Users\jp18b\Desktop\SCH_FEP`
- Python ≥ 3.11. Test runner: `pytest tests/ -v`.
- Import all numeric constants from `openfep.constants` — never hardcode `0.15`, `5.0`, `12`, etc.

---

## File Map

| File | Action | Responsibility |
|------|--------|----------------|
| `openfep/lambda_scheduler/__init__.py` | Create | Empty |
| `openfep/lambda_scheduler/presets.py` | Create | `rbfe_lambdas()` — linear RBFE lambda list |
| `openfep/msj/sea.py` | Modify | Append `minimize_stage()` |
| `openfep/mapping/__init__.py` | Create | Empty |
| `openfep/mapping/base.py` | Create | `AtomMapping` dataclass, `AtomMapper` ABC |
| `openfep/mapping/kartograf_mapper.py` | Create | `KartografMapper` wrapping `kartograf` |
| `openfep/mapping/rdkit_mcs_mapper.py` | Create | `RDKitMCSMapper` using `rdkit.Chem.rdFMCS` |
| `openfep/mapping/network.py` | Create | `FEPNetwork` — build, score, prune, list pairs |
| `openfep/prep/__init__.py` | Create | Empty |
| `openfep/prep/rbfe_prep.py` | Create | `RBFESystem` dataclass, `prepare_rbfe()` |
| `openfep/msj/rbfe_writer.py` | Create | `write_rbfe_msj()` — complex + solvent + main MSJ |
| `tests/test_presets.py` | Create | Lambda preset tests |
| `tests/test_mapping.py` | Create | Mapper + network tests |
| `tests/test_rbfe_prep.py` | Create | RBFE prep tests |
| `tests/test_rbfe_writer.py` | Create | MSJ writer output tests |

---

### Task 1: Lambda Presets + minimize_stage

**Files:**
- Create: `openfep/lambda_scheduler/__init__.py`
- Create: `openfep/lambda_scheduler/presets.py`
- Modify: `openfep/msj/sea.py` — append `minimize_stage()`
- Create: `tests/test_presets.py`

**Interfaces:**
- Consumes: `openfep.constants.RBFE_DEFAULT_LAMBDAS`, `openfep.constants.BACKBONE_ASL`, `SEABlock` from `openfep.msj.sea`
- Produces:
  - `rbfe_lambdas(n: int = RBFE_DEFAULT_LAMBDAS) -> list[float]`
  - `minimize_stage(restraint_asl: str | None = None, restraint_fc: float = 50.0) -> SEABlock`

- [ ] **Step 1: Create `openfep/lambda_scheduler/__init__.py`**

Empty file.

- [ ] **Step 2: Write `openfep/lambda_scheduler/presets.py`**

```python
from __future__ import annotations

from openfep.constants import RBFE_DEFAULT_LAMBDAS


def rbfe_lambdas(n: int = RBFE_DEFAULT_LAMBDAS) -> list[float]:
    """Return n evenly-spaced lambda values from 0.0 to 1.0 inclusive."""
    if n < 2:
        raise ValueError(f"n must be >= 2, got {n}")
    return [round(i / (n - 1), 6) for i in range(n)]
```

- [ ] **Step 3: Write `tests/test_presets.py`**

```python
from openfep.lambda_scheduler.presets import rbfe_lambdas
from openfep.constants import RBFE_DEFAULT_LAMBDAS
import pytest


def test_rbfe_lambdas_default_count():
    lams = rbfe_lambdas()
    assert len(lams) == RBFE_DEFAULT_LAMBDAS


def test_rbfe_lambdas_endpoints():
    lams = rbfe_lambdas()
    assert lams[0] == 0.0
    assert lams[-1] == 1.0


def test_rbfe_lambdas_custom_n():
    lams = rbfe_lambdas(n=5)
    assert len(lams) == 5
    assert lams[0] == 0.0
    assert lams[-1] == 1.0


def test_rbfe_lambdas_spacing():
    lams = rbfe_lambdas(n=3)
    assert abs(lams[1] - 0.5) < 1e-9


def test_rbfe_lambdas_too_small():
    with pytest.raises(ValueError):
        rbfe_lambdas(n=1)
```

- [ ] **Step 4: Run test to verify it passes**

```
pytest tests/test_presets.py -v
```
Expected: 5 passed.

- [ ] **Step 5: Append `minimize_stage` to `openfep/msj/sea.py`**

Append after `gcmc_stage`:

```python
def minimize_stage(
    restraint_asl: str | None = None,
    restraint_fc: float = 50.0,
) -> SEABlock:
    """
    Desmond minimize stage.
    restraint_asl: ASL expression (without "asl:" prefix) — written only when not None.
    restraint_fc: force constant in kcal/mol/Å².
    """
    block = SEABlock("minimize")
    if restraint_asl is not None:
        r = SEABlock(
            "restraints.new",
            ("atom", f"asl:{restraint_asl}"),
            ("force_constant", restraint_fc),
        )
        block.add("restraints.new", r)
    return block
```

- [ ] **Step 6: Add tests for `minimize_stage` to `tests/test_stages.py`**

Append to `tests/test_stages.py`:

```python
from openfep.msj.sea import minimize_stage


def test_minimize_stage_name():
    assert minimize_stage().name == "minimize"


def test_minimize_stage_no_restraints():
    text = render(minimize_stage())
    assert "restraints" not in text


def test_minimize_stage_with_restraints():
    text = render(minimize_stage(restraint_asl="backbone"))
    assert "restraints.new" in text
    assert "asl:backbone" in text


def test_minimize_stage_force_constant():
    text = render(minimize_stage(restraint_asl="backbone", restraint_fc=25.0))
    assert "25.0" in text
```

- [ ] **Step 7: Run full test suite**

```
pytest tests/ -v
```
Expected: all tests pass (previous 56 + 5 preset + 4 minimize = 65 total).

- [ ] **Step 8: Commit**

```
git add openfep/lambda_scheduler/ openfep/msj/sea.py tests/test_presets.py tests/test_stages.py
git commit -m "feat: lambda presets and minimize_stage"
```

---

### Task 2: Atom Mapper ABC + Kartograf Mapper

**Files:**
- Create: `openfep/mapping/__init__.py`
- Create: `openfep/mapping/base.py`
- Create: `openfep/mapping/kartograf_mapper.py`
- Create: `tests/test_mapping.py` (partial — Kartograf tests only)

**Interfaces:**
- Consumes: `rdkit.Chem.Mol`, `kartograf.atom_mapper.KartografAtomMapper`
- Produces:
  - `AtomMapping(mol_a, mol_b, a_to_b: dict[int,int], score: float)`
  - `AtomMapper` ABC with `def map(self, mol_a: Chem.Mol, mol_b: Chem.Mol) -> AtomMapping`
  - `KartografMapper(map_hydrogens: bool = False)`

**Note on molecule creation for tests:** Use RDKit SMILES parsing (`Chem.MolFromSmiles` + `AllChem.EmbedMolecule` + `AllChem.UFFOptimizeMolecule`) to create 3D molecules. Benzene and toluene are good simple test molecules (one heavy atom difference → clear MCS).

- [ ] **Step 1: Create `openfep/mapping/__init__.py`**

Empty file.

- [ ] **Step 2: Write `openfep/mapping/base.py`**

```python
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

from rdkit.Chem import Mol


@dataclass
class AtomMapping:
    """Atom-index correspondence between two RDKit molecules."""
    mol_a: Mol
    mol_b: Mol
    a_to_b: dict[int, int]   # atom index in mol_a → atom index in mol_b
    score: float              # 0–1, higher = better mapping quality


class AtomMapper(ABC):
    @abstractmethod
    def map(self, mol_a: Mol, mol_b: Mol) -> AtomMapping:
        """Return the best atom mapping between mol_a and mol_b."""
```

- [ ] **Step 3: Write `openfep/mapping/kartograf_mapper.py`**

```python
from __future__ import annotations

from rdkit.Chem import Mol
from kartograf.atom_mapper import KartografAtomMapper as _Kartograf

from openfep.mapping.base import AtomMapper, AtomMapping


class KartografMapper(AtomMapper):
    """3D-geometry-based atom mapper using Kartograf."""

    def __init__(self, map_hydrogens: bool = False) -> None:
        self._mapper = _Kartograf(atom_map_hydrogens=map_hydrogens)

    def map(self, mol_a: Mol, mol_b: Mol) -> AtomMapping:
        mapping_result = self._mapper.suggest_mappings(mol_a, mol_b)[0]
        a_to_b: dict[int, int] = dict(mapping_result.componentA_to_componentB)
        # Score: fraction of heavy atoms mapped (Tanimoto-like)
        n_heavy_a = sum(1 for a in mol_a.GetAtoms() if a.GetAtomicNum() != 1)
        n_heavy_b = sum(1 for a in mol_b.GetAtoms() if a.GetAtomicNum() != 1)
        mapped = len(a_to_b)
        score = mapped / max(n_heavy_a, n_heavy_b)
        return AtomMapping(mol_a=mol_a, mol_b=mol_b, a_to_b=a_to_b, score=score)
```

- [ ] **Step 4: Write failing tests in `tests/test_mapping.py`**

```python
from __future__ import annotations

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from openfep.mapping.base import AtomMapper, AtomMapping
from openfep.mapping.kartograf_mapper import KartografMapper


def _mol3d(smiles: str) -> Chem.Mol:
    """Create an embedded 3D RDKit Mol from a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.UFFOptimizeMolecule(mol)
    return Chem.RemoveHs(mol)


@pytest.fixture(scope="module")
def benzene():
    return _mol3d("c1ccccc1")


@pytest.fixture(scope="module")
def toluene():
    return _mol3d("Cc1ccccc1")


# ── AtomMapper ABC ────────────────────────────────────────────────────────────

def test_atom_mapper_is_abstract():
    with pytest.raises(TypeError):
        AtomMapper()


# ── AtomMapping dataclass ─────────────────────────────────────────────────────

def test_atom_mapping_fields(benzene, toluene):
    m = AtomMapping(mol_a=benzene, mol_b=toluene, a_to_b={0: 1}, score=0.5)
    assert m.score == 0.5
    assert m.a_to_b == {0: 1}


# ── KartografMapper ───────────────────────────────────────────────────────────

def test_kartograf_returns_atom_mapping(benzene, toluene):
    mapper = KartografMapper()
    result = mapper.map(benzene, toluene)
    assert isinstance(result, AtomMapping)


def test_kartograf_mapping_score_positive(benzene, toluene):
    result = KartografMapper().map(benzene, toluene)
    assert 0.0 < result.score <= 1.0


def test_kartograf_maps_benzene_ring_atoms(benzene, toluene):
    result = KartografMapper().map(benzene, toluene)
    # Benzene has 6 atoms, toluene has 7 — all 6 benzene atoms should map
    assert len(result.a_to_b) >= 6


def test_kartograf_mapping_indices_valid(benzene, toluene):
    result = KartografMapper().map(benzene, toluene)
    n_a = benzene.GetNumAtoms()
    n_b = toluene.GetNumAtoms()
    for ia, ib in result.a_to_b.items():
        assert 0 <= ia < n_a
        assert 0 <= ib < n_b


def test_kartograf_is_atom_mapper_subclass():
    assert issubclass(KartografMapper, AtomMapper)
```

- [ ] **Step 5: Run test to verify**

```
pytest tests/test_mapping.py -v -k "not lomap and not network"
```
Expected: all Kartograf tests pass (skip lomap/network not yet implemented).

- [ ] **Step 6: Commit**

```
git add openfep/mapping/ tests/test_mapping.py
git commit -m "feat: AtomMapper ABC and KartografMapper"
```

---

### Task 3: Lomap Mapper + FEP Network

**Files:**
- Create: `openfep/mapping/lomap_mapper.py`
- Create: `openfep/mapping/network.py`
- Extend: `tests/test_mapping.py` — add Lomap + network tests

**Interfaces:**
- Consumes: `lomap.MCS`, `networkx`, `AtomMapper`, `AtomMapping` from Task 2
- Produces:
  - `LomapMapper(max3d: float = 1.0)`
  - `FEPNetwork(molecules: list[Mol], mapper: AtomMapper)` with `.build() -> nx.Graph`, `.optimal_network() -> nx.Graph`, `.rbfe_pairs() -> list[tuple[Mol, Mol, AtomMapping]]`

- [ ] **Step 1: Write `openfep/mapping/lomap_mapper.py`**

```python
from __future__ import annotations

from rdkit.Chem import Mol
import lomap

from openfep.mapping.base import AtomMapper, AtomMapping


class LomapMapper(AtomMapper):
    """MCS-based atom mapper using Lomap2."""

    def __init__(self, max3d: float = 1.0) -> None:
        self._max3d = max3d

    def map(self, mol_a: Mol, mol_b: Mol) -> AtomMapping:
        mcs = lomap.MCS.getMapping(mol_a, mol_b, hydrogens=False, max3d=self._max3d)
        a_to_b: dict[int, int] = dict(mcs.heavy_atom_mcs_map())
        score: float = mcs.score()
        return AtomMapping(mol_a=mol_a, mol_b=mol_b, a_to_b=a_to_b, score=score)
```

**Note:** Verify the exact `lomap.MCS` API against the installed lomap2 version. The method names `heavy_atom_mcs_map()` and `score()` are from lomap2's public API — confirm with `import lomap; help(lomap.MCS)` if needed.

- [ ] **Step 2: Write `openfep/mapping/network.py`**

```python
from __future__ import annotations

import itertools

import networkx as nx
from rdkit.Chem import Mol

from openfep.mapping.base import AtomMapper, AtomMapping

_MIN_SCORE = 0.2   # Lomap default pruning threshold


class FEPNetwork:
    """
    Build and score a perturbation network from a list of molecules.
    Edges are weighted by mapping score (higher = better).
    """

    def __init__(self, molecules: list[Mol], mapper: AtomMapper) -> None:
        if len(molecules) < 2:
            raise ValueError("FEPNetwork requires at least 2 molecules")
        self._molecules = molecules
        self._mapper = mapper
        self._mappings: dict[tuple[int, int], AtomMapping] = {}

    def build(self) -> nx.Graph:
        """All-pairs mapping. Returns graph with nodes=mol-indices, edges=score."""
        g = nx.Graph()
        g.add_nodes_from(range(len(self._molecules)))
        for i, j in itertools.combinations(range(len(self._molecules)), 2):
            mapping = self._mapper.map(self._molecules[i], self._molecules[j])
            self._mappings[(i, j)] = mapping
            if mapping.score >= _MIN_SCORE:
                g.add_edge(i, j, weight=mapping.score, mapping=mapping)
        return g

    def optimal_network(self) -> nx.Graph:
        """MST on inverted weights (maximize mapping quality)."""
        g = self.build()
        # nx.maximum_spanning_tree maximises weight directly
        return nx.maximum_spanning_tree(g, weight="weight")

    def rbfe_pairs(self) -> list[tuple[Mol, Mol, AtomMapping]]:
        """Return (mol_a, mol_b, mapping) for each edge in the optimal network."""
        tree = self.optimal_network()
        pairs = []
        for i, j in tree.edges():
            key = (min(i, j), max(i, j))
            mapping = self._mappings[key]
            pairs.append((self._molecules[i], self._molecules[j], mapping))
        return pairs
```

- [ ] **Step 3: Add Lomap + network tests to `tests/test_mapping.py`**

Append to `tests/test_mapping.py`:

```python
from openfep.mapping.lomap_mapper import LomapMapper
from openfep.mapping.network import FEPNetwork
import networkx as nx


# ── LomapMapper ───────────────────────────────────────────────────────────────

def test_lomap_returns_atom_mapping(benzene, toluene):
    result = LomapMapper().map(benzene, toluene)
    assert isinstance(result, AtomMapping)


def test_lomap_score_positive(benzene, toluene):
    result = LomapMapper().map(benzene, toluene)
    assert 0.0 < result.score <= 1.0


def test_lomap_is_atom_mapper_subclass():
    assert issubclass(LomapMapper, AtomMapper)


# ── FEPNetwork ────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def three_mols(benzene, toluene):
    ethylbenzene = _mol3d("CCc1ccccc1")
    return [benzene, toluene, ethylbenzene]


def test_fepnetwork_build_returns_graph(three_mols):
    net = FEPNetwork(three_mols, KartografMapper())
    g = net.build()
    assert isinstance(g, nx.Graph)


def test_fepnetwork_build_has_all_nodes(three_mols):
    g = FEPNetwork(three_mols, KartografMapper()).build()
    assert g.number_of_nodes() == 3


def test_fepnetwork_optimal_network_is_spanning_tree(three_mols):
    net = FEPNetwork(three_mols, KartografMapper())
    tree = net.optimal_network()
    assert nx.is_tree(tree)


def test_fepnetwork_rbfe_pairs_length(three_mols):
    pairs = FEPNetwork(three_mols, KartografMapper()).rbfe_pairs()
    # MST of 3 nodes has exactly 2 edges
    assert len(pairs) == 2


def test_fepnetwork_rbfe_pairs_contain_mappings(three_mols):
    pairs = FEPNetwork(three_mols, KartografMapper()).rbfe_pairs()
    for mol_a, mol_b, mapping in pairs:
        assert isinstance(mapping, AtomMapping)
        assert mapping.score > 0.0


def test_fepnetwork_requires_two_mols(benzene):
    with pytest.raises(ValueError):
        FEPNetwork([benzene], KartografMapper())
```

- [ ] **Step 4: Run full mapping tests**

```
pytest tests/test_mapping.py -v
```
Expected: all tests pass.

- [ ] **Step 5: Run full suite for regressions**

```
pytest tests/ -v
```
Expected: all tests pass.

- [ ] **Step 6: Commit**

```
git add openfep/mapping/lomap_mapper.py openfep/mapping/network.py tests/test_mapping.py
git commit -m "feat: LomapMapper and FEPNetwork"
```

---

### Task 4: RBFE System Preparer

**Files:**
- Create: `openfep/prep/__init__.py`
- Create: `openfep/prep/rbfe_prep.py`
- Create: `tests/fixtures/rbfe_complex.mae` — synthetic 2-CT fixture (receptor + ligand)
- Create: `tests/test_rbfe_prep.py`

**Interfaces:**
- Consumes: `openfep.mae_parser.parse_mae`, `openfep.constants` (buffers, StrucTag), `AtomMapping` from Task 2
- Produces:
  ```python
  @dataclass
  class RBFESystem:
      fmp_path: Path                        # written .fmp file
      atom_mapping_paths: dict[str, Path]   # edge_key -> atom_mapping.txt
      receptor_count: int                   # 1 + bool(solvent_ct) + bool(membrane_ct)
      is_charged: bool
      buffer_width: float                   # COMPLEX_BUFFER_WIDTH or NET_CHARGE_COMPLEX_BUFFER_WIDTH
      has_membrane: bool

  def prepare_rbfe(
      mae_path: Path,
      output_dir: Path,
      mappings: list[tuple[str, str, AtomMapping]],  # [(lig_a_name, lig_b_name, mapping)]
      jobname: str = "rbfe",
  ) -> RBFESystem
  ```

**FMP file format — CRITICAL:** Before implementing `write_fmp_file()`, read `C:\Users\jp18b\Desktop\SCH_FEP\desmond\stage\fep_mapper.py` to find the exact `.fmp` format the `fep_mapper` MSJ stage expects. Look for `open(…".fmp"…)` or `fmp_file` references. The `.fmp` file likely encodes the FEP perturbation graph as an edge list. Implement exactly what the reference expects. If the format is a SEA block, use `render()` from `openfep.msj.sea`. If it is plain text, write it directly.

**Net charge detection:** The ligand CT's `s_fep_net_charge` property (integer CT property) tells you if the ligand is charged. If the property is absent, assume neutral (not charged). Extend `CTBlock` in `mae_parser.py` if needed to carry this property, OR parse it directly in `prepare_rbfe` from `ct_block.raw_block` using a simple regex.

**Receptor count logic:**
```python
receptor_count = 1  # always: the protein
if any(ct.struc_tag == StrucTag.SOLVENT for ct in ct_blocks):
    receptor_count += 1
if any(ct.struc_tag == StrucTag.MEMBRANE for ct in ct_blocks):
    receptor_count += 1
```

**Atom mapping file format** (one pair per file, tab-separated):
```
# atom_a_idx  atom_b_idx
0    2
1    3
2    4
```

- [ ] **Step 1: Create `openfep/prep/__init__.py`**

Empty file.

- [ ] **Step 2: Create fixture `tests/fixtures/rbfe_complex.mae`**

Minimal synthetic MAE with one receptor CT and one ligand CT. Copy the structure of `tests/fixtures/two_ct.mae` but ensure the ligand has `s_fep_net_charge` set to 0:

```
 f_m_ct {
  s_m_title
  s_fep_struc_tag
  i_fep_net_charge
  r_m_mass
 :::
  "Receptor"
  receptor
  0
  12300.0
 :::
 }
 f_m_ct {
  s_m_title
  s_fep_struc_tag
  i_fep_net_charge
  r_m_mass
 :::
  "LigandA"
  ligand
  0
  285.4
 :::
 }
```

- [ ] **Step 3: Read the .fmp format from the reference**

Read `C:\Users\jp18b\Desktop\SCH_FEP\desmond\stage\fep_mapper.py` (read-only reference). Find the section that reads or writes `.fmp` files. Note the exact format.

- [ ] **Step 4: Write `openfep/prep/rbfe_prep.py`**

```python
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from openfep.constants import (
    COMPLEX_BUFFER_WIDTH,
    NET_CHARGE_COMPLEX_BUFFER_WIDTH,
    StrucTag,
)
from openfep.mae_parser import parse_mae
from openfep.mapping.base import AtomMapping


@dataclass
class RBFESystem:
    fmp_path: Path
    atom_mapping_paths: dict[str, Path] = field(default_factory=dict)
    receptor_count: int = 1
    is_charged: bool = False
    buffer_width: float = COMPLEX_BUFFER_WIDTH
    has_membrane: bool = False


def prepare_rbfe(
    mae_path: Path,
    output_dir: Path,
    mappings: list[tuple[str, str, AtomMapping]],
    jobname: str = "rbfe",
) -> RBFESystem:
    """
    Prepare RBFE input files.
    mappings: list of (lig_a_name, lig_b_name, AtomMapping)
    Writes: {jobname}.fmp and {edge_key}_atom_mapping.txt files in output_dir.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    ct_blocks = parse_mae(mae_path)

    # Detect net charge from any ligand CT
    is_charged = _detect_net_charge(ct_blocks)
    has_membrane = any(ct.struc_tag == StrucTag.MEMBRANE for ct in ct_blocks)
    receptor_count = _compute_receptor_count(ct_blocks)
    buffer_width = (
        NET_CHARGE_COMPLEX_BUFFER_WIDTH if is_charged else COMPLEX_BUFFER_WIDTH
    )

    # Write atom mapping files
    mapping_paths: dict[str, Path] = {}
    for lig_a, lig_b, mapping in mappings:
        edge_key = f"{lig_a}_to_{lig_b}"
        p = output_dir / f"{edge_key}_atom_mapping.txt"
        _write_atom_mapping(p, mapping)
        mapping_paths[edge_key] = p

    # Write .fmp file
    fmp_path = output_dir / f"{jobname}.fmp"
    _write_fmp(fmp_path, mappings, mae_path)

    return RBFESystem(
        fmp_path=fmp_path,
        atom_mapping_paths=mapping_paths,
        receptor_count=receptor_count,
        is_charged=is_charged,
        buffer_width=buffer_width,
        has_membrane=has_membrane,
    )


def _detect_net_charge(ct_blocks) -> bool:
    """Return True if any ligand CT has non-zero i_fep_net_charge."""
    import re
    for ct in ct_blocks:
        if ct.struc_tag != StrucTag.LIGAND:
            continue
        # Find i_fep_net_charge in header and read its value
        lines = ct.raw_block.splitlines()
        prop_names: list[str] = []
        values: list[str] = []
        sep_count = 0
        state = "header"
        for line in lines[1:]:
            stripped = line.strip()
            if stripped == ":::":
                sep_count += 1
                state = "values" if sep_count == 1 else "done"
                if state == "done":
                    break
                continue
            if state == "header":
                if stripped and not stripped.startswith("#"):
                    prop_names.append(stripped)
            elif state == "values":
                if stripped:
                    values.append(stripped)
        if "i_fep_net_charge" in prop_names:
            idx = prop_names.index("i_fep_net_charge")
            if idx < len(values):
                try:
                    if int(values[idx]) != 0:
                        return True
                except ValueError:
                    pass
    return False


def _compute_receptor_count(ct_blocks) -> int:
    count = 1
    if any(ct.struc_tag == StrucTag.SOLVENT for ct in ct_blocks):
        count += 1
    if any(ct.struc_tag == StrucTag.MEMBRANE for ct in ct_blocks):
        count += 1
    return count


def _write_atom_mapping(path: Path, mapping: AtomMapping) -> None:
    lines = ["# atom_a_idx\tatom_b_idx\n"]
    for ia, ib in sorted(mapping.a_to_b.items()):
        lines.append(f"{ia}\t{ib}\n")
    path.write_text("".join(lines), encoding="utf-8")


def _write_fmp(
    path: Path,
    mappings: list[tuple[str, str, AtomMapping]],
    mae_path: Path,
) -> None:
    """
    Write the FEP Map (.fmp) file.
    FORMAT: determined by reading desmond/stage/fep_mapper.py.
    Implement the exact format expected by the fep_mapper MSJ stage.
    The placeholder below is a reasonable guess — verify and replace.
    """
    # TODO: replace with exact format from desmond/stage/fep_mapper.py
    lines = [f"# FEP Map — generated by openfep\n"]
    lines.append(f"structure_file = {mae_path}\n")
    for lig_a, lig_b, _ in mappings:
        lines.append(f"edge {lig_a} {lig_b}\n")
    path.write_text("".join(lines), encoding="utf-8")
```

**IMPORTANT:** The `_write_fmp` function above is a placeholder. After reading `desmond/stage/fep_mapper.py`, replace the implementation with the exact format the `fep_mapper` stage expects.

- [ ] **Step 5: Write `tests/test_rbfe_prep.py`**

```python
from pathlib import Path
import pytest
from openfep.prep.rbfe_prep import prepare_rbfe, RBFESystem
from openfep.mapping.base import AtomMapping
from rdkit import Chem
from rdkit.Chem import AllChem

FIXTURE_MAE = Path(__file__).parent / "fixtures" / "rbfe_complex.mae"


def _dummy_mapping() -> AtomMapping:
    mol_a = Chem.MolFromSmiles("c1ccccc1")
    mol_b = Chem.MolFromSmiles("Cc1ccccc1")
    return AtomMapping(mol_a=mol_a, mol_b=mol_b, a_to_b={0: 1, 1: 2}, score=0.9)


def test_prepare_rbfe_returns_system(tmp_path):
    system = prepare_rbfe(
        mae_path=FIXTURE_MAE,
        output_dir=tmp_path,
        mappings=[("LigandA", "LigandB", _dummy_mapping())],
    )
    assert isinstance(system, RBFESystem)


def test_prepare_rbfe_fmp_created(tmp_path):
    system = prepare_rbfe(FIXTURE_MAE, tmp_path, [("A", "B", _dummy_mapping())])
    assert system.fmp_path.exists()


def test_prepare_rbfe_atom_mapping_created(tmp_path):
    system = prepare_rbfe(FIXTURE_MAE, tmp_path, [("A", "B", _dummy_mapping())])
    assert "A_to_B" in system.atom_mapping_paths
    assert system.atom_mapping_paths["A_to_B"].exists()


def test_prepare_rbfe_atom_mapping_content(tmp_path):
    system = prepare_rbfe(FIXTURE_MAE, tmp_path, [("A", "B", _dummy_mapping())])
    content = system.atom_mapping_paths["A_to_B"].read_text()
    assert "0\t1" in content


def test_prepare_rbfe_neutral_buffer(tmp_path):
    from openfep.constants import COMPLEX_BUFFER_WIDTH
    system = prepare_rbfe(FIXTURE_MAE, tmp_path, [("A", "B", _dummy_mapping())])
    assert system.buffer_width == COMPLEX_BUFFER_WIDTH


def test_prepare_rbfe_receptor_count(tmp_path):
    system = prepare_rbfe(FIXTURE_MAE, tmp_path, [("A", "B", _dummy_mapping())])
    assert system.receptor_count == 1  # fixture has only receptor + ligand, no solvent/membrane


def test_prepare_rbfe_not_charged(tmp_path):
    system = prepare_rbfe(FIXTURE_MAE, tmp_path, [("A", "B", _dummy_mapping())])
    assert system.is_charged is False
```

- [ ] **Step 6: Run tests**

```
pytest tests/test_rbfe_prep.py -v
```
Expected: all 7 tests pass.

- [ ] **Step 7: Run full suite**

```
pytest tests/ -v
```
Expected: all tests pass.

- [ ] **Step 8: Commit**

```
git add openfep/prep/ tests/test_rbfe_prep.py tests/fixtures/rbfe_complex.mae
git commit -m "feat: RBFE system preparer"
```

---

### Task 5: RBFE MSJ Writer

**Files:**
- Create: `openfep/msj/rbfe_writer.py`
- Create: `tests/test_rbfe_writer.py`

**Interfaces:**
- Consumes: all stage builders from `openfep.msj.sea`, `RBFESystem` from Task 4, `rbfe_lambdas` from Task 1, `openfep.constants`
- Produces:
  ```python
  def write_rbfe_msj(
      system: RBFESystem,
      output_dir: Path,
      lambdas: list[float] | None = None,    # None → use rbfe_lambdas()
      temperature: float = 300.0,
      n_equil_stages: int = 4,               # total equilibration simulate stages
      time_equil_ns: list[float] | None = None,  # per-stage times; None → defaults
  ) -> dict[str, Path]:
      """
      Write main.msj + complex.msj + solvent.msj.
      Returns {"main": Path, "complex": Path, "solvent": Path}.
      """
  ```

**Equilibration stage sequence:**  
Before implementing, read `C:\Users\jp18b\Desktop\SCH_FEP\desmond\msj_generator\fep\small_molecule.py` to find the exact number of equilibration stages and their times. The common pattern is 4 equilibration stages. Default times (verify against reference):
- Stage 1: 0.12 ns (NPT, `polarization_restraints="full"`)
- Stage 2: 0.12 ns (NPT, `polarization_restraints="full"`)
- Stage 3: 0.24 ns (NPT, `polarization_restraints="full"`)
- Stage 4: 0.24 ns (NPT, `polarization_restraints="decay"`) ← LAST equil stage

Default production time: 5.0 ns per leg.

**Stage sequence for `complex.msj`:**
```
task_stage(fep_type="small_molecule", schedule="default", n_lambda=len(lambdas))
assign_forcefield_stage(custom_charge_mode="assign")
build_geometry_stage(buffer=system.buffer_width, neutralize=system.is_charged,
                     make_alchemical_water=system.is_charged)
[gcmc_stage()]          # only if system.has_membrane
minimize_stage(restraint_asl=BACKBONE_ASL)
minimize_stage()        # no restraints
simulate_stage(0.12, temperature, "NPT", polarization_restraint="full")
simulate_stage(0.12, temperature, "NPT", polarization_restraint="full")
simulate_stage(0.24, temperature, "NPT", polarization_restraint="full")
simulate_stage(0.24, temperature, "NPT", polarization_restraint="decay")
lambda_hopping_stage(time_ns=5.0, ensemble="NPT")
trim_stage()
```

**Stage sequence for `solvent.msj`:**
```
task_stage(fep_type="small_molecule", schedule="default", n_lambda=len(lambdas))
assign_forcefield_stage(custom_charge_mode="keep")
build_geometry_stage(buffer=SOLVENT_BUFFER_WIDTH, neutralize=system.is_charged,
                     make_alchemical_water=system.is_charged)
minimize_stage()
simulate_stage(0.12, temperature, "NPT", polarization_restraint="full")
simulate_stage(0.24, temperature, "NPT", polarization_restraint="decay")
lambda_hopping_stage(time_ns=5.0, ensemble="NPT")
trim_stage()
```

**`main.msj` content:**
```
task_stage(fep_type="small_molecule", schedule="default", n_lambda=len(lambdas))
fep_mapper {
  graph_file = "{jobname}.fmp"
  atom_mapping = "{edge_key}_atom_mapping.txt"
  receptor = {system.receptor_count}
}
```

The `fep_mapper` block is a plain SEABlock (name="fep_mapper") with keys for `graph_file`, `atom_mapping`, and `receptor`. Write it using `SEABlock("fep_mapper", ...)`.

- [ ] **Step 1: Write `openfep/msj/rbfe_writer.py`**

```python
from __future__ import annotations

from pathlib import Path

from openfep.constants import BACKBONE_ASL, SOLVENT_BUFFER_WIDTH
from openfep.lambda_scheduler.presets import rbfe_lambdas
from openfep.msj.sea import (
    SEABlock,
    render_msj,
    task_stage,
    assign_forcefield_stage,
    build_geometry_stage,
    simulate_stage,
    lambda_hopping_stage,
    minimize_stage,
    trim_stage,
    gcmc_stage,
)
from openfep.prep.rbfe_prep import RBFESystem

_DEFAULT_EQUIL_TIMES_NS = [0.12, 0.12, 0.24, 0.24]
_DEFAULT_PRODUCTION_NS = 5.0


def write_rbfe_msj(
    system: RBFESystem,
    output_dir: Path,
    lambdas: list[float] | None = None,
    temperature: float = 300.0,
    equil_times_ns: list[float] | None = None,
    production_ns: float = _DEFAULT_PRODUCTION_NS,
    jobname: str = "rbfe",
    edge_key: str = "",
) -> dict[str, Path]:
    """
    Write main.msj, complex.msj, solvent.msj to output_dir.
    Returns {"main": Path, "complex": Path, "solvent": Path}.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if lambdas is None:
        lambdas = rbfe_lambdas()
    if equil_times_ns is None:
        equil_times_ns = _DEFAULT_EQUIL_TIMES_NS

    n_lambda = len(lambdas)

    complex_msj = _build_complex_msj(system, n_lambda, temperature, equil_times_ns, production_ns)
    solvent_msj = _build_solvent_msj(system, n_lambda, temperature, equil_times_ns, production_ns)
    main_msj = _build_main_msj(system, n_lambda, jobname, edge_key)

    paths = {}
    for name, stages in [("complex", complex_msj), ("solvent", solvent_msj), ("main", main_msj)]:
        p = output_dir / f"{name}.msj"
        p.write_text(render_msj(stages), encoding="utf-8")
        paths[name] = p
    return paths


def _build_complex_msj(
    system: RBFESystem,
    n_lambda: int,
    temperature: float,
    equil_times_ns: list[float],
    production_ns: float,
) -> list[SEABlock]:
    stages: list[SEABlock] = [
        task_stage("small_molecule", "default", n_lambda),
        assign_forcefield_stage(custom_charge_mode="assign"),
        build_geometry_stage(
            buffer=system.buffer_width,
            neutralize=system.is_charged,
            make_alchemical_water=system.is_charged,
        ),
    ]
    if system.has_membrane:
        stages.append(gcmc_stage())
    stages += [
        minimize_stage(restraint_asl=BACKBONE_ASL),
        minimize_stage(),
    ]
    for i, t in enumerate(equil_times_ns):
        is_last = i == len(equil_times_ns) - 1
        pol = "decay" if is_last else "full"
        stages.append(simulate_stage(t, temperature, "NPT", polarization_restraint=pol))
    stages += [
        lambda_hopping_stage(production_ns, "NPT"),
        trim_stage(),
    ]
    return stages


def _build_solvent_msj(
    system: RBFESystem,
    n_lambda: int,
    temperature: float,
    equil_times_ns: list[float],
    production_ns: float,
) -> list[SEABlock]:
    # Solvent leg uses only the last two equil stages (lighter protocol)
    solvent_equil = equil_times_ns[-2:]
    stages: list[SEABlock] = [
        task_stage("small_molecule", "default", n_lambda),
        assign_forcefield_stage(custom_charge_mode="keep"),
        build_geometry_stage(
            buffer=SOLVENT_BUFFER_WIDTH,
            neutralize=system.is_charged,
            make_alchemical_water=system.is_charged,
        ),
        minimize_stage(),
    ]
    for i, t in enumerate(solvent_equil):
        is_last = i == len(solvent_equil) - 1
        pol = "decay" if is_last else "full"
        stages.append(simulate_stage(t, temperature, "NPT", polarization_restraint=pol))
    stages += [
        lambda_hopping_stage(production_ns, "NPT"),
        trim_stage(),
    ]
    return stages


def _build_main_msj(
    system: RBFESystem,
    n_lambda: int,
    jobname: str,
    edge_key: str,
) -> list[SEABlock]:
    fep_mapper = SEABlock(
        "fep_mapper",
        ("graph_file", f"{jobname}.fmp"),
        ("atom_mapping", f"{edge_key}_atom_mapping.txt"),
        ("receptor", system.receptor_count),
    )
    return [
        task_stage("small_molecule", "default", n_lambda),
        fep_mapper,
    ]
```

**Note on solvent equil stages:** The choice of "last 2 equil stages" for the solvent leg is an approximation. After reading `desmond/msj_generator/fep/small_molecule.py`, adjust the solvent equil protocol to match the reference exactly.

- [ ] **Step 2: Write `tests/test_rbfe_writer.py`**

```python
from pathlib import Path
import pytest
from openfep.msj.rbfe_writer import write_rbfe_msj
from openfep.prep.rbfe_prep import RBFESystem
from openfep.constants import (
    COMPLEX_BUFFER_WIDTH,
    SOLVENT_BUFFER_WIDTH,
    BACKBONE_ASL,
)

FIXTURE_MAE = Path(__file__).parent / "fixtures" / "rbfe_complex.mae"


@pytest.fixture
def neutral_system(tmp_path):
    return RBFESystem(
        fmp_path=tmp_path / "rbfe.fmp",
        atom_mapping_paths={"A_to_B": tmp_path / "A_to_B_atom_mapping.txt"},
        receptor_count=1,
        is_charged=False,
        buffer_width=COMPLEX_BUFFER_WIDTH,
        has_membrane=False,
    )


@pytest.fixture
def charged_system(tmp_path):
    from openfep.constants import NET_CHARGE_COMPLEX_BUFFER_WIDTH
    return RBFESystem(
        fmp_path=tmp_path / "rbfe.fmp",
        atom_mapping_paths={},
        receptor_count=1,
        is_charged=True,
        buffer_width=NET_CHARGE_COMPLEX_BUFFER_WIDTH,
        has_membrane=False,
    )


def test_write_rbfe_msj_returns_three_files(neutral_system, tmp_path):
    paths = write_rbfe_msj(neutral_system, tmp_path)
    assert "main" in paths
    assert "complex" in paths
    assert "solvent" in paths


def test_complex_msj_file_exists(neutral_system, tmp_path):
    paths = write_rbfe_msj(neutral_system, tmp_path)
    assert paths["complex"].exists()


def test_solvent_msj_file_exists(neutral_system, tmp_path):
    paths = write_rbfe_msj(neutral_system, tmp_path)
    assert paths["solvent"].exists()


def test_main_msj_file_exists(neutral_system, tmp_path):
    paths = write_rbfe_msj(neutral_system, tmp_path)
    assert paths["main"].exists()


# ── complex.msj content ───────────────────────────────────────────────────────

def test_complex_msj_starts_with_task(neutral_system, tmp_path):
    paths = write_rbfe_msj(neutral_system, tmp_path)
    content = paths["complex"].read_text()
    assert content.strip().startswith("task {")


def test_complex_msj_assign_forcefield_assign(neutral_system, tmp_path):
    content = write_rbfe_msj(neutral_system, tmp_path)["complex"].read_text()
    assert '"assign"' in content


def test_complex_msj_lambda_hopping_present(neutral_system, tmp_path):
    content = write_rbfe_msj(neutral_system, tmp_path)["complex"].read_text()
    assert "lambda_hopping {" in content


def test_complex_msj_trim_present(neutral_system, tmp_path):
    content = write_rbfe_msj(neutral_system, tmp_path)["complex"].read_text()
    assert "trim {" in content


def test_complex_msj_polarization_full_and_decay(neutral_system, tmp_path):
    content = write_rbfe_msj(neutral_system, tmp_path)["complex"].read_text()
    assert '"full"' in content
    assert '"decay"' in content


def test_complex_msj_no_polarization_in_lambda_hopping(neutral_system, tmp_path):
    content = write_rbfe_msj(neutral_system, tmp_path)["complex"].read_text()
    # polarization_restraints must NOT appear inside lambda_hopping block
    lh_start = content.index("lambda_hopping {")
    lh_end = content.index("\n\n", lh_start)
    lh_block = content[lh_start:lh_end]
    assert "polarization_restraints" not in lh_block


def test_complex_msj_backbone_asl_in_minimize(neutral_system, tmp_path):
    content = write_rbfe_msj(neutral_system, tmp_path)["complex"].read_text()
    assert BACKBONE_ASL in content


def test_complex_msj_buffer_width(neutral_system, tmp_path):
    content = write_rbfe_msj(neutral_system, tmp_path)["complex"].read_text()
    assert str(COMPLEX_BUFFER_WIDTH) in content


# ── solvent.msj content ───────────────────────────────────────────────────────

def test_solvent_msj_assign_keep(neutral_system, tmp_path):
    content = write_rbfe_msj(neutral_system, tmp_path)["solvent"].read_text()
    assert '"keep"' in content


def test_solvent_msj_solvent_buffer(neutral_system, tmp_path):
    content = write_rbfe_msj(neutral_system, tmp_path)["solvent"].read_text()
    assert str(SOLVENT_BUFFER_WIDTH) in content


def test_solvent_msj_lambda_hopping_present(neutral_system, tmp_path):
    content = write_rbfe_msj(neutral_system, tmp_path)["solvent"].read_text()
    assert "lambda_hopping {" in content


# ── main.msj content ──────────────────────────────────────────────────────────

def test_main_msj_has_fep_mapper(neutral_system, tmp_path):
    content = write_rbfe_msj(neutral_system, tmp_path)["main"].read_text()
    assert "fep_mapper {" in content


def test_main_msj_graph_file(neutral_system, tmp_path):
    content = write_rbfe_msj(neutral_system, tmp_path, jobname="myjob")["main"].read_text()
    assert "myjob.fmp" in content


def test_main_msj_receptor_count(neutral_system, tmp_path):
    content = write_rbfe_msj(neutral_system, tmp_path)["main"].read_text()
    assert "receptor = 1" in content


# ── charged system ────────────────────────────────────────────────────────────

def test_charged_complex_msj_neutralize(charged_system, tmp_path):
    content = write_rbfe_msj(charged_system, tmp_path)["complex"].read_text()
    assert "neutralize_system = true" in content


def test_charged_complex_msj_alchemical_water(charged_system, tmp_path):
    content = write_rbfe_msj(charged_system, tmp_path)["complex"].read_text()
    assert "make_alchemical_water = true" in content
```

- [ ] **Step 3: Run tests**

```
pytest tests/test_rbfe_writer.py -v
```
Expected: all tests pass.

- [ ] **Step 4: Run full suite for regressions**

```
pytest tests/ -v
```
Expected: all tests pass (56 base + 5 presets + 4 minimize + mapping + prep + writer).

- [ ] **Step 5: Commit**

```
git add openfep/msj/rbfe_writer.py tests/test_rbfe_writer.py
git commit -m "feat: RBFE MSJ writer"
```

---

## Key Correctness Constraints (sub-plan 2 scope)

| Constraint | Where enforced |
|---|---|
| RBFE complex: `assign_custom_charge.mode = "assign"` | `rbfe_writer.py:_build_complex_msj` |
| RBFE solvent: `assign_custom_charge.mode = "keep"` | `rbfe_writer.py:_build_solvent_msj` |
| Equil stages except last: `polarization_restraints = "full"` | `rbfe_writer.py` loop |
| Last equil stage: `polarization_restraints = "decay"` | `rbfe_writer.py` loop |
| `lambda_hopping` NOT `simulate` | `lambda_hopping_stage()` from sub-plan 1 |
| `polarization_restraints` absent from `lambda_hopping` | `lambda_hopping_stage()` from sub-plan 1 |
| `energy_group` always in `lambda_hopping` | `lambda_hopping_stage()` from sub-plan 1 |
| Charged RBFE: `neutralize_system=True`, `make_alchemical_water=True` | `build_geometry_stage()` |
| Buffer: complex=5.0, charged_complex=8.0, solvent=10.0 | constants imported in writer |
| FMP format matches `fep_mapper` expectation | implementer verifies from desmond/ reference |

## Build Order

```
[Task 1] presets + minimize_stage  →  [Task 5] rbfe_writer (uses both)
[Task 2] base + KartografMapper    →  [Task 3] FEPNetwork, [Task 4] rbfe_prep
[Task 3] FEPNetwork                →  [Task 4] rbfe_prep (mappings input)
[Task 4] RBFESystem                →  [Task 5] rbfe_writer (system input)
```
