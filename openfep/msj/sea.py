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
            lines.append(render(val, indent + 1))
        else:
            lines.append(f"{child_pad}{key} = {sea_val(val)}")
    lines.append(f"{pad}}}")
    return "\n".join(lines)


def render_msj(stages: list[SEABlock]) -> str:
    """Render a sequence of stages as a complete MSJ file."""
    return "\n\n".join(render(s) for s in stages) + "\n"


# ── Stage builder functions ────────────────────────────────────────────────────
# Each returns a SEABlock. Callers compose them into a stage list and call render_msj().

from openfep.constants import (
    FORCEFIELD,
    HMR_TIMESTEPS,
    HMR_MIGRATION_INTERVAL,
    GCMC_SOLVENT_VDW_SCALE_FACTOR,
    ENEGRP_NAME_TEMPLATE,
    ENEGRP_FIRST,
    ENEGRP_INTERVAL,
    MIN_CHARGED_SALT_CONC,
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
    salt_conc: float = MIN_CHARGED_SALT_CONC,
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
      - Only written to block when not None (caller gates on forcefield == OPLS4)
    time_ns is converted to ps internally (multiply by 1000).
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
    Production stage. MUST be named 'lambda_hopping' (not 'simulate').
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
