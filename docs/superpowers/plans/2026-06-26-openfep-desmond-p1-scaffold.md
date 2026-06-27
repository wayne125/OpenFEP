# openfep-Desmond — Sub-plan 1: Scaffold + Constants + MAE Parser + SEA Stage Builder

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stand up the `openfep` package with all shared constants, a MAE CT parser, and the SEA-format stage-builder primitives that every subsequent sub-plan depends on.

**Architecture:** Pure Python library with no Schrodinger imports. MSJ stages are rendered as SEA-format text strings via `SEABlock`. The MAE parser is a plain-text state machine. Constants mirror `desmond/constants.py` values but carry no Schrodinger dependency.

**Tech Stack:** Python 3.11+, pytest, no Schrodinger imports in `openfep/` itself.

## Global Constraints

- No `from schrodinger` anywhere in `openfep/` — all Schrodinger values are reimplemented as plain Python constants.
- Forcefield string: `"OPLS4"` (OPLSVersion.F17). Not `"opls4"`, not `"OPLS3e"`.
- `polarization_restraints` is only written to MSJ when `forcefield == "OPLS4"`.
- Production stage is named `lambda_hopping`, never `simulate`.
- All SEA bool values: `true` / `false` (lowercase), not Python `True`/`False`.
- Python minimum: 3.11 (uses `list[str]` etc. without `from __future__ import annotations`).
- Test runner: `pytest`. All tests live under `tests/`. Run with `pytest tests/ -v`.

---

## File Map

| File | Action | Responsibility |
|------|--------|---------------|
| `pyproject.toml` | Create | Package metadata + dependencies |
| `openfep/__init__.py` | Create | Version export |
| `openfep/constants.py` | Create | All numeric/string constants, StrucTag, KEEP_STRUC_TAGS |
| `openfep/mae_parser.py` | Create | Parse `.mae` CT blocks, extract `s_fep_struc_tag` |
| `openfep/msj/__init__.py` | Create | Empty |
| `openfep/msj/sea.py` | Create | SEABlock renderer + all MSJ stage constructors |
| `tests/__init__.py` | Create | Empty |
| `tests/fixtures/two_ct.mae` | Create | Synthetic MAE fixture (receptor CT + ligand CT) |
| `tests/test_constants.py` | Create | Smoke-test that all expected symbols exist and have correct values |
| `tests/test_mae_parser.py` | Create | Parse fixture; verify struc_tag extraction; missing-tag skip |
| `tests/test_sea.py` | Create | SEABlock rendering + every stage constructor output |

---

### Task 1: Project Scaffold

**Files:**
- Create: `pyproject.toml`
- Create: `openfep/__init__.py`
- Create: `openfep/msj/__init__.py`
- Create: `tests/__init__.py`

**Interfaces:**
- Produces: `openfep` importable as a package; `pytest tests/ -v` runs and finds no tests yet (0 collected = pass)

- [ ] **Step 1: Write `pyproject.toml`**

```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "openfep"
version = "0.1.0"
description = "Open-source FEP toolkit driving Schrodinger Desmond without proprietary Python stack"
requires-python = ">=3.11"
dependencies = [
    "rdkit>=2023.9",
    "kartograf>=1.1",
    "lomap2>=2.2",
    "networkx>=3.2",
    "alchemlyb>=2.3",
    "numpy>=1.26",
    "matplotlib>=3.8",
    "click>=8.1",
]

[project.optional-dependencies]
dev = ["pytest>=8.0", "pytest-cov"]

[project.scripts]
openfep = "openfep.cli:main"
```

- [ ] **Step 2: Create package files**

`openfep/__init__.py`:
```python
__version__ = "0.1.0"
```

`openfep/msj/__init__.py`:
```python
```

`tests/__init__.py`:
```python
```

- [ ] **Step 3: Install in editable mode**

```bash
pip install -e ".[dev]"
```

Expected: no errors; `python -c "import openfep; print(openfep.__version__)"` prints `0.1.0`.

- [ ] **Step 4: Verify test runner finds nothing (clean baseline)**

```bash
pytest tests/ -v
```

Expected: `no tests ran` or `0 passed`. No import errors.

- [ ] **Step 5: Commit**

```bash
git add pyproject.toml openfep/__init__.py openfep/msj/__init__.py tests/__init__.py
git commit -m "feat: scaffold openfep package with pyproject.toml"
```

---

### Task 2: Constants

**Files:**
- Create: `openfep/constants.py`
- Create: `tests/test_constants.py`

**Interfaces:**
- Produces:
  - `openfep.constants.FORCEFIELD: str` → `"OPLS4"`
  - `openfep.constants.StrucTag` — class with `RECEPTOR`, `LIGAND`, `COMPLEX`, `SOLVENT`, `MEMBRANE: str`
  - `openfep.constants.KEEP_STRUC_TAGS: dict[str, dict[str, list[str]]]`
  - `openfep.constants.BACKBONE_ASL: str`
  - All buffer widths, HMR values, ABFE lambda counts (see implementation step)

- [ ] **Step 1: Write the failing test**

`tests/test_constants.py`:
```python
from openfep.constants import (
    FORCEFIELD,
    StrucTag,
    KEEP_STRUC_TAGS,
    COMPLEX_BUFFER_WIDTH,
    NET_CHARGE_COMPLEX_BUFFER_WIDTH,
    SOLVENT_BUFFER_WIDTH,
    VACUUM_BUFFER_WIDTH,
    HMR_TIMESTEPS,
    HMR_MIGRATION_INTERVAL,
    MIN_CHARGED_SALT_CONC,
    GCMC_SOLVENT_VDW_SCALE_FACTOR,
    ENEGRP_NAME_TEMPLATE,
    ENEGRP_FIRST,
    ENEGRP_INTERVAL,
    ABFE_COMPLEX_LAMBDAS,
    ABFE_SOLVENT_LAMBDAS,
    ABFE_RESTRAINED_COMPLEX_LAMBDAS,
    ABFE_RESTRAINED_SOLVENT_LAMBDAS,
    ABFE_CHARGED_COMPLEX_LAMBDAS,
    ABFE_CHARGED_SOLVENT_LAMBDAS,
    ABFE_CHARGED_RESTRAINED_COMPLEX_LAMBDAS,
    ABFE_CHARGED_RESTRAINED_SOLVENT_LAMBDAS,
    RBFE_DEFAULT_LAMBDAS,
    BACKBONE_ASL,
    RECEPTOR_ASL,
)


def test_forcefield_is_opls4():
    assert FORCEFIELD == "OPLS4"


def test_struc_tag_values():
    assert StrucTag.RECEPTOR == "receptor"
    assert StrucTag.LIGAND == "ligand"
    assert StrucTag.COMPLEX == "complex"
    assert StrucTag.SOLVENT == "solvent"
    assert StrucTag.MEMBRANE == "membrane"


def test_keep_struc_tags_rbfe_complex_includes_membrane():
    tags = KEEP_STRUC_TAGS["rbfe"]["complex"]
    assert StrucTag.MEMBRANE in tags
    assert StrucTag.RECEPTOR in tags
    assert StrucTag.LIGAND in tags


def test_keep_struc_tags_abfe_md_includes_membrane():
    tags = KEEP_STRUC_TAGS["abfe"]["md"]
    assert StrucTag.MEMBRANE in tags
    assert StrucTag.RECEPTOR in tags
    assert StrucTag.LIGAND in tags


def test_abfe_complex_leg_keeps_complex_not_ligand():
    tags = KEEP_STRUC_TAGS["abfe"]["complex"]
    assert StrucTag.COMPLEX in tags
    assert StrucTag.LIGAND not in tags


def test_buffer_widths():
    assert COMPLEX_BUFFER_WIDTH == 5.0
    assert NET_CHARGE_COMPLEX_BUFFER_WIDTH == 8.0
    assert SOLVENT_BUFFER_WIDTH == 10.0
    assert VACUUM_BUFFER_WIDTH == 100.0


def test_hmr():
    assert HMR_TIMESTEPS == [0.004, 0.004, 0.008]
    assert HMR_MIGRATION_INTERVAL == 0.024


def test_abfe_lambda_counts_default():
    assert ABFE_COMPLEX_LAMBDAS == 68
    assert ABFE_SOLVENT_LAMBDAS == 60
    assert ABFE_RESTRAINED_COMPLEX_LAMBDAS == 80
    assert ABFE_RESTRAINED_SOLVENT_LAMBDAS == 68


def test_abfe_lambda_counts_charged():
    assert ABFE_CHARGED_COMPLEX_LAMBDAS == 108
    assert ABFE_CHARGED_SOLVENT_LAMBDAS == 60
    assert ABFE_CHARGED_RESTRAINED_COMPLEX_LAMBDAS == 128
    assert ABFE_CHARGED_RESTRAINED_SOLVENT_LAMBDAS == 68


def test_backbone_asl_includes_nucleic():
    assert "nucleic_acids" in BACKBONE_ASL
    assert "protein" in BACKBONE_ASL
    assert "not a.ele H" in BACKBONE_ASL


def test_enegrp_template():
    assert "$JOBNAME$" in ENEGRP_NAME_TEMPLATE
    assert "_enegrp.dat" in ENEGRP_NAME_TEMPLATE
```

- [ ] **Step 2: Run test to confirm it fails**

```bash
pytest tests/test_constants.py -v
```

Expected: `ImportError: cannot import name 'FORCEFIELD' from 'openfep.constants'` (module doesn't exist yet).

- [ ] **Step 3: Write `openfep/constants.py`**

```python
# Forcefield
FORCEFIELD = "OPLS4"  # OPLSVersion.F17 string in MSJ. Not OPLS3e (F16).

# Buffer widths (Å) — from desmond/msj_constants.py
COMPLEX_BUFFER_WIDTH = 5.0
NET_CHARGE_COMPLEX_BUFFER_WIDTH = 8.0
SOLVENT_BUFFER_WIDTH = 10.0
VACUUM_BUFFER_WIDTH = 100.0

# HMR — from desmond/msj_constants.py _HMR_TIMESTEPS / _HMR_MIGRATION_INTERVAL
HMR_TIMESTEPS = [0.004, 0.004, 0.008]   # ps
HMR_MIGRATION_INTERVAL = 0.024           # ps

# Salt
MIN_CHARGED_SALT_CONC = 0.15            # M

# GCMC — from desmond/msj_constants.py _GCMC_SOLVENT_VDW_SCALE_FACTOR
GCMC_SOLVENT_VDW_SCALE_FACTOR = 0.75

# ENERGY_GROUP output — from desmond/msj_constants.py ENERGY_GROUP
ENEGRP_NAME_TEMPLATE = "$JOBNAME$[_replica$REPLICA$]_enegrp.dat"
ENEGRP_FIRST = 0.0    # ps
ENEGRP_INTERVAL = 1.2  # ps

# RBFE
RBFE_DEFAULT_LAMBDAS = 12

# ABFE lambda window counts — from desmond/fep_schedule.py:282-301
# Default (neutral) protocol
ABFE_COMPLEX_LAMBDAS = 68
ABFE_SOLVENT_LAMBDAS = 60
ABFE_RESTRAINED_COMPLEX_LAMBDAS = 80
ABFE_RESTRAINED_SOLVENT_LAMBDAS = 68
# Charged protocol
ABFE_CHARGED_COMPLEX_LAMBDAS = 108
ABFE_CHARGED_SOLVENT_LAMBDAS = 60
ABFE_CHARGED_RESTRAINED_COMPLEX_LAMBDAS = 128
ABFE_CHARGED_RESTRAINED_SOLVENT_LAMBDAS = 68

# Backbone ASL for receptor restraints — from absolute_binding.py
BACKBONE_ASL = (
    "((protein and backbone) or "
    "(nucleic_acids and nucleic_backbone)) and not a.ele H"
)
RECEPTOR_ASL = "protein or nucleic_acids"


class StrucTag:
    """Values of the s_fep_struc_tag CT property — from desmond/constants.py FEP_STRUC_TAG."""
    RECEPTOR = "receptor"
    LIGAND   = "ligand"
    COMPLEX  = "complex"
    SOLVENT  = "solvent"
    MEMBRANE = "membrane"


# Per-leg CT selection — mirrors desmond/msj_generator/keep_struc_tags.py
KEEP_STRUC_TAGS: dict[str, dict[str, list[str]]] = {
    "rbfe": {
        "complex": [
            StrucTag.RECEPTOR, StrucTag.MEMBRANE,
            StrucTag.SOLVENT,  StrucTag.LIGAND,
        ],
        "solvent": [StrucTag.LIGAND],
        "vacuum":  [StrucTag.LIGAND],
    },
    "abfe": {
        "md": [
            StrucTag.RECEPTOR, StrucTag.SOLVENT,
            StrucTag.MEMBRANE, StrucTag.LIGAND,
        ],
        "complex": [
            StrucTag.RECEPTOR, StrucTag.MEMBRANE,
            StrucTag.SOLVENT,  StrucTag.COMPLEX,
        ],
        "solvent": [StrucTag.LIGAND],
    },
}
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/test_constants.py -v
```

Expected: all tests pass. `14 passed`.

- [ ] **Step 5: Commit**

```bash
git add openfep/constants.py tests/test_constants.py
git commit -m "feat: add openfep constants (StrucTag, buffer widths, ABFE lambda counts, HMR)"
```

---

### Task 3: MAE CT Parser

**Files:**
- Create: `openfep/mae_parser.py`
- Create: `tests/fixtures/two_ct.mae`
- Create: `tests/test_mae_parser.py`

**Interfaces:**
- Consumes: nothing
- Produces:
  - `openfep.mae_parser.CTBlock` — dataclass with `struc_tag: str`, `raw_block: str`
  - `openfep.mae_parser.parse_mae(path: str | Path) -> list[CTBlock]`
  - CT blocks without `s_fep_struc_tag` are silently skipped (not returned)

**MAE format primer:** A `.mae` file is plain text. CT blocks start with ` f_m_ct {` (may have leading space). Each CT has a header section (property names, one per line) then `:::` then a values section (one value per line, in the same order as names) then `:::`. We only need CT-level properties — they appear before `m_atom[...]` or `m_bond[...]` sub-blocks.

- [ ] **Step 1: Create the synthetic fixture**

`tests/fixtures/two_ct.mae`:
```
 f_m_ct {
  s_m_title
  s_fep_struc_tag
  r_m_mass
 :::
  "Receptor protein"
  receptor
  12300.0
 :::
  m_atom[2] {
   # First column is atom index #
   i_m_mmod_type
   :::
   1 1
   2 1
   :::
  }
 }
 f_m_ct {
  s_m_title
  s_fep_struc_tag
  r_m_mass
 :::
  "Ligand A"
  ligand
  285.4
 :::
  m_atom[1] {
   # First column is atom index #
   i_m_mmod_type
   :::
   1 6
   :::
  }
 }
 f_m_ct {
  s_m_title
  r_m_mass
 :::
  "No tag CT"
  500.0
 :::
 }
```

- [ ] **Step 2: Write the failing tests**

`tests/test_mae_parser.py`:
```python
from pathlib import Path
import pytest
from openfep.mae_parser import parse_mae, CTBlock

FIXTURE = Path(__file__).parent / "fixtures" / "two_ct.mae"


def test_parse_returns_only_tagged_cts():
    blocks = parse_mae(FIXTURE)
    # Third CT has no s_fep_struc_tag, must be skipped
    assert len(blocks) == 2


def test_first_ct_is_receptor():
    blocks = parse_mae(FIXTURE)
    assert blocks[0].struc_tag == "receptor"


def test_second_ct_is_ligand():
    blocks = parse_mae(FIXTURE)
    assert blocks[1].struc_tag == "ligand"


def test_raw_block_contains_ct_header():
    blocks = parse_mae(FIXTURE)
    assert "f_m_ct" in blocks[0].raw_block


def test_ct_block_dataclass():
    b = CTBlock(struc_tag="receptor", raw_block="f_m_ct {}")
    assert b.struc_tag == "receptor"
    assert b.raw_block == "f_m_ct {}"


def test_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        parse_mae("/nonexistent/path.mae")
```

- [ ] **Step 3: Run tests to confirm failure**

```bash
pytest tests/test_mae_parser.py -v
```

Expected: `ImportError: cannot import name 'parse_mae' from 'openfep.mae_parser'`.

- [ ] **Step 4: Write `openfep/mae_parser.py`**

```python
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass
class CTBlock:
    struc_tag: str    # value of s_fep_struc_tag CT property
    raw_block: str    # full text of the f_m_ct { ... } block


def parse_mae(path: str | Path) -> list[CTBlock]:
    """
    Parse a .mae file and return CT blocks that have an s_fep_struc_tag property.
    CTs without this property are silently skipped.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"MAE file not found: {path}")
    text = path.read_text(encoding="latin-1")
    return _extract_tagged_cts(text)


def _extract_tagged_cts(text: str) -> list[CTBlock]:
    blocks = []
    # Split on CT start markers; first chunk is pre-header, skip it
    raw_blocks = _split_ct_blocks(text)
    for raw in raw_blocks:
        tag = _extract_struc_tag(raw)
        if tag is not None:
            blocks.append(CTBlock(struc_tag=tag, raw_block=raw))
    return blocks


def _split_ct_blocks(text: str) -> list[str]:
    """Return list of raw f_m_ct block texts (including the opening line)."""
    blocks = []
    lines = text.splitlines(keepends=True)
    current: list[str] = []
    in_block = False
    depth = 0

    for line in lines:
        stripped = line.strip()
        if stripped == "f_m_ct {":
            if in_block and current:
                blocks.append("".join(current))
            current = [line]
            in_block = True
            depth = 1
        elif in_block:
            current.append(line)
            depth += stripped.count("{") - stripped.count("}")
            if depth <= 0:
                blocks.append("".join(current))
                current = []
                in_block = False
                depth = 0

    if current:
        blocks.append("".join(current))
    return blocks


def _extract_struc_tag(raw_block: str) -> str | None:
    """
    Extract s_fep_struc_tag value from a CT block.

    MAE CT-level properties follow this layout:
        f_m_ct {
          prop_name_1         ← header: one property name per line
          prop_name_2
         :::
          value_1            ← values: same order as names
          value_2
         :::
    We find s_fep_struc_tag in the header, then read the value at the same index.
    """
    lines = raw_block.splitlines()

    # Find the header section: between the opening { and the first :::
    prop_names: list[str] = []
    sep_count = 0
    values: list[str] = []

    state = "header"
    for line in lines[1:]:   # skip "f_m_ct {"
        stripped = line.strip()
        if stripped == ":::":
            sep_count += 1
            if sep_count == 1:
                state = "values"
            elif sep_count == 2:
                break
            continue

        if state == "header":
            if stripped and not stripped.startswith("#") and not stripped.startswith("m_"):
                prop_names.append(stripped)
        elif state == "values":
            if stripped:
                values.append(stripped.strip('"'))

    if "s_fep_struc_tag" not in prop_names:
        return None
    idx = prop_names.index("s_fep_struc_tag")
    if idx >= len(values):
        return None
    return values[idx]
```

- [ ] **Step 5: Run tests**

```bash
pytest tests/test_mae_parser.py -v
```

Expected: `6 passed`.

- [ ] **Step 6: Commit**

```bash
git add openfep/mae_parser.py tests/mae_parser_test.py tests/fixtures/two_ct.mae
git commit -m "feat: add MAE CT parser (s_fep_struc_tag extraction)"
```

---

### Task 4: SEA Block Renderer

**Files:**
- Create: `openfep/msj/sea.py`
- Create: `tests/test_sea.py`

**Interfaces:**
- Produces:
  - `openfep.msj.sea.SEABlock(name: str, *children: tuple[str, Any])` — a named SEA block
  - `openfep.msj.sea.render(block: SEABlock, indent: int = 0) -> str` — SEA text output
  - `openfep.msj.sea.sea_val(v: Any) -> str` — serialize a scalar/list value
  - `openfep.msj.sea.render_msj(stages: list[SEABlock]) -> str` — join stages with blank lines

**SEA format rules** (from `desmond/multisim/parser.py` behavior):
- Block: `name {\n  key = value\n}` — blocks open with `name {`, close with `}`
- String values: `"quoted"` with double quotes
- Bool: `true` / `false` (lowercase)
- Number: bare (no quotes)
- Nested block as value: `key = name { ... }` OR bare `name { ... }` (unnamed nested map)
- List: `[ item1 item2 ... ]` space-separated

- [ ] **Step 1: Write the failing tests**

`tests/test_sea.py`:
```python
from openfep.msj.sea import SEABlock, render, sea_val, render_msj


def test_sea_val_string():
    assert sea_val("OPLS4") == '"OPLS4"'


def test_sea_val_bool_true():
    assert sea_val(True) == "true"


def test_sea_val_bool_false():
    assert sea_val(False) == "false"


def test_sea_val_int():
    assert sea_val(42) == "42"


def test_sea_val_float():
    assert sea_val(3.14) == "3.14"


def test_sea_val_list_of_numbers():
    assert sea_val([0.004, 0.004, 0.008]) == "[0.004 0.004 0.008]"


def test_simple_block_renders():
    block = SEABlock("assign_forcefield", ("forcefield", "OPLS4"))
    text = render(block)
    assert text.startswith("assign_forcefield {")
    assert 'forcefield = "OPLS4"' in text
    assert text.endswith("}")


def test_block_bool_field_lowercase():
    block = SEABlock("assign_forcefield", ("hydrogen_mass_repartition", True))
    text = render(block)
    assert "hydrogen_mass_repartition = true" in text


def test_nested_block_renders():
    inner = SEABlock("simulate", ("fep.lambda", "default:12"))
    outer = SEABlock("set_family", ("simulate", inner))
    text = render(outer)
    assert "set_family {" in text
    assert "simulate {" in text
    assert 'fep.lambda = "default:12"' in text


def test_render_msj_joins_with_blank_line():
    a = SEABlock("task", ("task", "generic"))
    b = SEABlock("simulate", ("time", 120.0))
    text = render_msj([a, b])
    assert "task {" in text
    assert "simulate {" in text
    assert "\n\n" in text   # blank line between stages


def test_indentation():
    block = SEABlock("simulate", ("time", 120.0))
    text = render(block, indent=1)
    assert text.startswith("  simulate {")   # 2 spaces per indent level
    assert "    time = 120.0" in text         # child indented by 4
```

- [ ] **Step 2: Run tests to confirm failure**

```bash
pytest tests/test_sea.py -v
```

Expected: `ImportError`.

- [ ] **Step 3: Write `openfep/msj/sea.py`**

```python
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

_INDENT = "  "  # 2 spaces


@dataclass
class SEABlock:
    """
    A named SEA block:  name { key = val ... }
    children is a list of (key, value) pairs in insertion order.
    value can be: str, int, float, bool, SEABlock, or list of the above.
    """
    name: str
    children: list[tuple[str, Any]] = field(default_factory=list)

    def __init__(self, name: str, *children: tuple[str, Any]):
        self.name = name
        self.children = list(children)

    def add(self, key: str, value: Any) -> "SEABlock":
        """Append a (key, value) pair and return self for chaining."""
        self.children.append((key, value))
        return self


def sea_val(v: Any) -> str:
    """Serialize a scalar Python value to SEA text."""
    if isinstance(v, bool):           # must test bool before int (bool subclasses int)
        return "true" if v else "false"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        # Avoid unnecessary trailing zeros but keep precision
        return str(v)
    if isinstance(v, str):
        return f'"{v}"'
    if isinstance(v, list):
        return "[" + " ".join(sea_val(item) for item in v) + "]"
    if isinstance(v, SEABlock):
        raise TypeError(
            "Use render() for SEABlock values, not sea_val(). "
            "Embed via render(block, indent=parent_indent+1)."
        )
    raise TypeError(f"Unsupported SEA value type: {type(v).__name__}")


def render(block: SEABlock, indent: int = 0) -> str:
    """Render a SEABlock to SEA/MSJ text."""
    pad = _INDENT * indent
    child_pad = _INDENT * (indent + 1)
    lines = [f"{pad}{block.name} {{"]
    for key, val in block.children:
        if isinstance(val, SEABlock):
            nested = render(val, indent + 1)
            lines.append(f"{child_pad}{key} = {nested.lstrip()}")
        else:
            lines.append(f"{child_pad}{key} = {sea_val(val)}")
    lines.append(f"{pad}}}")
    return "\n".join(lines)


def render_msj(stages: list[SEABlock]) -> str:
    """Render a sequence of stages as a complete MSJ file."""
    return "\n\n".join(render(s) for s in stages) + "\n"
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/test_sea.py -v
```

Expected: `10 passed`.

- [ ] **Step 5: Commit**

```bash
git add openfep/msj/sea.py tests/test_sea.py
git commit -m "feat: add SEA block renderer for MSJ generation"
```

---

### Task 5: MSJ Stage Builders

**Files:**
- Modify: `openfep/msj/sea.py` — append stage builder functions (no new file; stage builders use `SEABlock` directly)
- Create: `tests/test_stages.py`

**Interfaces:**
- Consumes:
  - `openfep.msj.sea.SEABlock`, `render`
  - `openfep.constants.FORCEFIELD`, `HMR_TIMESTEPS`, `HMR_MIGRATION_INTERVAL`, `GCMC_SOLVENT_VDW_SCALE_FACTOR`, `ENEGRP_NAME_TEMPLATE`, `ENEGRP_FIRST`, `ENEGRP_INTERVAL`
- Produces (all in `openfep.msj.sea`):
  - `task_stage(fep_type, schedule, n_lambda) -> SEABlock`
  - `assign_forcefield_stage(forcefield, hmr, custom_charge_mode) -> SEABlock`
  - `build_geometry_stage(buffer, neutralize, salt_conc, box_shape, make_alchemical_water, add_alchemical_ions) -> SEABlock`
  - `simulate_stage(time_ns, temperature, ensemble, polarization_restraint) -> SEABlock`
  - `lambda_hopping_stage(time_ns, ensemble, energy_groups) -> SEABlock`
  - `trim_stage() -> SEABlock`
  - `load_restraints_from_structure_stage() -> SEABlock`
  - `gcmc_stage() -> SEABlock`

**Critical constraints:**
- `polarization_restraint` parameter on `simulate_stage`: only written to the block when the forcefield is OPLS4 and the caller provides a non-None value. The caller (MSJ writer) is responsible for passing the right value based on forcefield; this function writes it when given.
- `lambda_hopping_stage` always writes `energy_group` block (deviation from Schrodinger default).
- Stage name must be `lambda_hopping`, not `simulate`.

- [ ] **Step 1: Write the failing tests**

`tests/test_stages.py`:
```python
from openfep.msj.sea import (
    render,
    task_stage,
    assign_forcefield_stage,
    build_geometry_stage,
    simulate_stage,
    lambda_hopping_stage,
    trim_stage,
    load_restraints_from_structure_stage,
    gcmc_stage,
)
from openfep.constants import (
    FORCEFIELD,
    ENEGRP_NAME_TEMPLATE,
    HMR_TIMESTEPS,
    HMR_MIGRATION_INTERVAL,
    GCMC_SOLVENT_VDW_SCALE_FACTOR,
)


# ── task_stage ─────────────────────────────────────────────────────────────────

def test_task_stage_name():
    block = task_stage(fep_type="small_molecule", schedule="default", n_lambda=12)
    assert block.name == "task"


def test_task_stage_contains_lambda_string():
    block = task_stage(fep_type="small_molecule", schedule="default", n_lambda=12)
    text = render(block)
    assert '"default:12"' in text


def test_task_stage_fep_type():
    text = render(task_stage(fep_type="absolute_binding", schedule="default", n_lambda=68))
    assert '"absolute_binding"' in text


# ── assign_forcefield_stage ────────────────────────────────────────────────────

def test_assign_forcefield_name():
    block = assign_forcefield_stage()
    assert block.name == "assign_forcefield"


def test_assign_forcefield_opls4():
    text = render(assign_forcefield_stage())
    assert f'"OPLS4"' in text


def test_assign_forcefield_hmr_on():
    text = render(assign_forcefield_stage(hmr=True))
    assert "hydrogen_mass_repartition = true" in text


def test_assign_forcefield_hmr_off():
    text = render(assign_forcefield_stage(hmr=False))
    assert "hydrogen_mass_repartition = false" in text


def test_assign_custom_charge_mode_assign():
    text = render(assign_forcefield_stage(custom_charge_mode="assign"))
    assert '"assign"' in text


def test_assign_custom_charge_mode_keep():
    text = render(assign_forcefield_stage(custom_charge_mode="keep"))
    assert '"keep"' in text


# ── build_geometry_stage ───────────────────────────────────────────────────────

def test_build_geometry_buffer():
    text = render(build_geometry_stage(buffer=5.0))
    assert "5.0" in text


def test_build_geometry_neutralize():
    text = render(build_geometry_stage(buffer=8.0, neutralize=True))
    assert "neutralize_system = true" in text


def test_build_geometry_cubic_when_charged():
    text = render(build_geometry_stage(buffer=8.0, box_shape="cubic"))
    assert '"cubic"' in text


def test_build_geometry_alchemical_water():
    text = render(build_geometry_stage(buffer=5.0, make_alchemical_water=True))
    assert "make_alchemical_water = true" in text


# ── simulate_stage ─────────────────────────────────────────────────────────────

def test_simulate_stage_name():
    block = simulate_stage(time_ns=0.12, temperature=300.0, ensemble="NPT")
    assert block.name == "simulate"


def test_simulate_stage_ensemble():
    text = render(simulate_stage(time_ns=0.12, temperature=300.0, ensemble="NPT"))
    assert '"NPT"' in text


def test_simulate_stage_polarization_full():
    text = render(simulate_stage(
        time_ns=0.12, temperature=300.0, ensemble="NPT",
        polarization_restraint="full"
    ))
    assert '"full"' in text
    assert "polarization_restraints" in text


def test_simulate_stage_polarization_none_not_written():
    text = render(simulate_stage(
        time_ns=0.12, temperature=300.0, ensemble="NPT",
        polarization_restraint=None
    ))
    assert "polarization_restraints" not in text


def test_simulate_stage_polarization_decay():
    text = render(simulate_stage(
        time_ns=0.12, temperature=300.0, ensemble="NPT",
        polarization_restraint="decay"
    ))
    assert '"decay"' in text


# ── lambda_hopping_stage ───────────────────────────────────────────────────────

def test_lambda_hopping_name():
    block = lambda_hopping_stage(time_ns=5.0, ensemble="NPT")
    assert block.name == "lambda_hopping"


def test_lambda_hopping_energy_group_always_present():
    text = render(lambda_hopping_stage(time_ns=5.0, ensemble="NPT"))
    assert "energy_group" in text
    assert ENEGRP_NAME_TEMPLATE in text


def test_lambda_hopping_ensemble():
    text = render(lambda_hopping_stage(time_ns=5.0, ensemble="muVT"))
    assert '"muVT"' in text


def test_lambda_hopping_hmr_timesteps():
    text = render(lambda_hopping_stage(time_ns=5.0, ensemble="NPT"))
    # HMR_TIMESTEPS = [0.004, 0.004, 0.008]
    assert "0.004" in text
    assert "0.008" in text


def test_lambda_hopping_migration_interval():
    text = render(lambda_hopping_stage(time_ns=5.0, ensemble="NPT"))
    assert str(HMR_MIGRATION_INTERVAL) in text


# ── utility stages ─────────────────────────────────────────────────────────────

def test_trim_stage_name():
    assert trim_stage().name == "trim"


def test_load_restraints_stage_name():
    assert load_restraints_from_structure_stage().name == "load_restraints_from_structure"


def test_gcmc_stage_name():
    assert gcmc_stage().name == "gcmc"


def test_gcmc_scale_solvent_vdw():
    text = render(gcmc_stage())
    assert str(GCMC_SOLVENT_VDW_SCALE_FACTOR) in text
```

- [ ] **Step 2: Run tests to confirm failure**

```bash
pytest tests/test_stages.py -v
```

Expected: `ImportError: cannot import name 'task_stage'`.

- [ ] **Step 3: Append stage builders to `openfep/msj/sea.py`**

Add after `render_msj()`:

```python
# ── Stage builder functions ───────────────────────────────────────────────────
# Each returns a SEABlock. Callers compose them into a stage list and call render_msj().

from openfep.constants import (
    FORCEFIELD,
    HMR_TIMESTEPS,
    HMR_MIGRATION_INTERVAL,
    GCMC_SOLVENT_VDW_SCALE_FACTOR,
    ENEGRP_NAME_TEMPLATE,
    ENEGRP_FIRST,
    ENEGRP_INTERVAL,
)


def task_stage(fep_type: str, schedule: str, n_lambda: int) -> SEABlock:
    """
    The first stage of every subjob MSJ.
    schedule: "default" | "flexible" | "charge"
    n_lambda: number of lambda windows
    fep_type: "small_molecule" (RBFE) | "absolute_binding" (ABFE)
    """
    simulate_inner = SEABlock(
        "simulate",
        ("fep.type", fep_type),
        ("fep.lambda", f"{schedule}:{n_lambda}"),
    )
    set_family = SEABlock("set_family", ("simulate", simulate_inner))
    return SEABlock("task", ("task", "generic"), ("set_family", set_family))


def assign_forcefield_stage(
    forcefield: str = FORCEFIELD,
    hmr: bool = True,
    custom_charge_mode: str = "assign",
) -> SEABlock:
    """
    Assigns forcefield and partial charges. Always first substantive stage.
    custom_charge_mode: "assign" (RBFE complex) or "keep" (all other legs).
    """
    block = SEABlock("assign_forcefield")
    block.add("forcefield", forcefield)
    block.add("hydrogen_mass_repartition", hmr)
    charge_block = SEABlock("assign_custom_charge", ("mode", custom_charge_mode))
    block.add("assign_custom_charge", charge_block)
    return block


def build_geometry_stage(
    buffer: float,
    neutralize: bool = False,
    salt_conc: float = 0.15,
    box_shape: str = "orthorhombic",
    make_alchemical_water: bool = False,
    add_alchemical_ions: bool = False,
) -> SEABlock:
    """Solvate + box setup. buffer in Å."""
    block = SEABlock("build_geometry")
    block.add("buffer_width", buffer)
    if neutralize:
        block.add("neutralize_system", True)
        salt_block = SEABlock(
            "salt",
            ("negative_ion", "Cl"),
            ("positive_ion", "Na"),
            ("concentration", salt_conc),
        )
        block.add("salt", salt_block)
        block.add("box_shape", box_shape)
    if make_alchemical_water:
        block.add("make_alchemical_water", True)
    if add_alchemical_ions:
        block.add("add_alchemical_ions", True)
    return block


def simulate_stage(
    time_ns: float,
    temperature: float,
    ensemble: str,
    polarization_restraint: str | None = None,
) -> SEABlock:
    """
    A single equilibration simulate stage.
    ensemble: "NPT" | "NVT" | "muVT"
    polarization_restraint: "full" | "decay" | None
      - Only written when not None (caller gates on forcefield == OPLS4)
    """
    block = SEABlock("simulate")
    block.add("time", time_ns * 1000)   # convert ns → ps for Desmond
    block.add("temperature", temperature)
    block.add("ensemble", ensemble)
    if polarization_restraint is not None:
        block.add("polarization_restraints", polarization_restraint)
    return block


def lambda_hopping_stage(
    time_ns: float,
    ensemble: str,
) -> SEABlock:
    """
    Production stage. Must be named 'lambda_hopping' (not 'simulate').
    energy_group block is always written (openfep deviation from Schrodinger default).
    HMR timesteps always applied.
    ensemble: "NPT" (solvent/vacuum) | "muVT" (complex with membrane/GCMC)
    """
    block = SEABlock("lambda_hopping")
    block.add("time", time_ns * 1000)
    block.add("ensemble", ensemble)
    # HMR
    block.add("timestep", HMR_TIMESTEPS)
    block.add("migration_interval", HMR_MIGRATION_INTERVAL)
    # energy_group — always on (openfep deviation: Schrodinger only enables on complex leg)
    eg = SEABlock(
        "energy_group",
        ("name", ENEGRP_NAME_TEMPLATE),
        ("first", ENEGRP_FIRST),
        ("interval", ENEGRP_INTERVAL),
    )
    block.add("energy_group", eg)
    return block


def trim_stage() -> SEABlock:
    """Erases intermediate *-in.cms.gz, *.cms.gz, *-out.tgz after production."""
    return SEABlock("trim")


def load_restraints_from_structure_stage() -> SEABlock:
    """ABFE complex.msj: reads Boresch geometry from md.msj output structure."""
    return SEABlock("load_restraints_from_structure")


def gcmc_stage() -> SEABlock:
    """Complex leg with membrane: GCMC grand-canonical water insertion/deletion."""
    return SEABlock("gcmc", ("scale_solvent_vdw", GCMC_SOLVENT_VDW_SCALE_FACTOR))
```

- [ ] **Step 4: Run all tests**

```bash
pytest tests/ -v
```

Expected: all tests pass (`30+ passed`). Check that earlier tasks' tests still pass.

- [ ] **Step 5: Commit**

```bash
git add openfep/msj/sea.py tests/test_stages.py
git commit -m "feat: add MSJ stage builder functions (assign_forcefield, simulate, lambda_hopping, etc.)"
```

---

## Self-Review Checklist

- [x] **Spec coverage:** Constants (§2 buffer widths, §3.3 lambda counts, §3.4 HMR/GCMC/enegrp), MAE parser (§3.2 s_fep_struc_tag), SEA renderer (§3.4 MSJ generator), all stage constructors (§3.4 stage sequences). No spec sections are skipped.
- [x] **Placeholder scan:** No TBD/TODO. Every step has complete code.
- [x] **Type consistency:** `CTBlock.struc_tag: str` used in Task 3; `KEEP_STRUC_TAGS` maps `str → list[str]` — consistent with `struc_tag`. `SEABlock` defined in Task 4, consumed in Task 5 — same class, same module.
- [x] **polarization_restraint gating:** `simulate_stage(polarization_restraint=None)` skips the field; caller (sub-plan 2 MSJ writer) checks `forcefield == OPLS4` before passing a value. Tests cover both `None` and `"full"` / `"decay"`.
- [x] **lambda_hopping naming:** Stage name is `"lambda_hopping"`, not `"simulate"`. Test asserts `block.name == "lambda_hopping"`.
- [x] **ABFE membrane in KEEP_STRUC_TAGS:** Both `abfe.md` and `rbfe.complex` include `StrucTag.MEMBRANE`. Tests assert this explicitly.

---

## Execution Options

Plan complete and saved to `docs/superpowers/plans/2026-06-26-openfep-desmond-p1-scaffold.md`.

**Two execution options:**

**1. Subagent-Driven (recommended)** — dispatch a fresh subagent per task, review between tasks.

**2. Inline Execution** — execute tasks in this session using executing-plans.

Which approach?
