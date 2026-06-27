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
    Write main.msj, complex.msj, and solvent.msj to output_dir.

    Returns {"main": Path, "complex": Path, "solvent": Path}.

    Parameters
    ----------
    system:
        RBFESystem produced by prepare_rbfe() or constructed directly.
    output_dir:
        Directory where MSJ files are written (created if absent).
    lambdas:
        Lambda window values. None → use rbfe_lambdas() defaults.
    temperature:
        Simulation temperature in Kelvin.
    equil_times_ns:
        Per-stage equilibration times (ns). None → [0.12, 0.12, 0.24, 0.24].
    production_ns:
        Lambda-hopping production time (ns) for both legs.
    jobname:
        Base name used in the main.msj fep_mapper graph_file entry.
    edge_key:
        Stem used in the atom_mapping filename in main.msj.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if lambdas is None:
        lambdas = rbfe_lambdas()
    if equil_times_ns is None:
        equil_times_ns = _DEFAULT_EQUIL_TIMES_NS

    n_lambda = len(lambdas)

    complex_stages = _build_complex_msj(
        system, n_lambda, temperature, equil_times_ns, production_ns
    )
    solvent_stages = _build_solvent_msj(
        system, n_lambda, temperature, equil_times_ns, production_ns
    )
    main_stages = _build_main_msj(system, n_lambda, jobname, edge_key)

    paths: dict[str, Path] = {}
    for name, stages in [
        ("complex", complex_stages),
        ("solvent", solvent_stages),
        ("main", main_stages),
    ]:
        p = output_dir / f"{name}.msj"
        p.write_text(render_msj(stages), encoding="utf-8")
        paths[name] = p

    return paths


# ── Private stage-list builders ───────────────────────────────────────────────


def _build_complex_msj(
    system: RBFESystem,
    n_lambda: int,
    temperature: float,
    equil_times_ns: list[float],
    production_ns: float,
) -> list[SEABlock]:
    """Build the complex-leg stage list."""
    stages: list[SEABlock] = [
        task_stage("small_molecule", "default", n_lambda),
        # Complex leg MUST use custom_charge_mode="assign"
        assign_forcefield_stage(custom_charge_mode="assign"),
        build_geometry_stage(
            buffer=system.buffer_width,
            neutralize=system.is_charged,
            make_alchemical_water=system.is_charged,
        ),
    ]

    if system.has_membrane:
        stages.append(gcmc_stage())

    # Two minimization stages: first with backbone restraints, second free
    stages.append(minimize_stage(restraint_asl=BACKBONE_ASL))
    stages.append(minimize_stage())

    # Equilibration: all except the last use polarization_restraint="full";
    # the last uses "decay"
    for i, t in enumerate(equil_times_ns):
        is_last = (i == len(equil_times_ns) - 1)
        pol = "decay" if is_last else "full"
        stages.append(simulate_stage(t, temperature, "NPT", polarization_restraint=pol))

    # Production + cleanup
    stages.append(lambda_hopping_stage(production_ns, "NPT"))
    stages.append(trim_stage())

    return stages


def _build_solvent_msj(
    system: RBFESystem,
    n_lambda: int,
    temperature: float,
    equil_times_ns: list[float],
    production_ns: float,
) -> list[SEABlock]:
    """Build the solvent-leg stage list (lighter equilibration protocol)."""
    # Solvent leg uses only the last two equil stages
    solvent_equil = equil_times_ns[-2:]

    stages: list[SEABlock] = [
        task_stage("small_molecule", "default", n_lambda),
        # Solvent leg MUST use custom_charge_mode="keep"
        assign_forcefield_stage(custom_charge_mode="keep"),
        build_geometry_stage(
            buffer=SOLVENT_BUFFER_WIDTH,
            neutralize=system.is_charged,
            make_alchemical_water=system.is_charged,
        ),
        minimize_stage(),
    ]

    for i, t in enumerate(solvent_equil):
        is_last = (i == len(solvent_equil) - 1)
        pol = "decay" if is_last else "full"
        stages.append(simulate_stage(t, temperature, "NPT", polarization_restraint=pol))

    stages.append(lambda_hopping_stage(production_ns, "NPT"))
    stages.append(trim_stage())

    return stages


def _build_main_msj(
    system: RBFESystem,
    n_lambda: int,
    jobname: str,
    edge_key: str,
) -> list[SEABlock]:
    """Build the main-leg stage list (task + fep_mapper orchestrator)."""
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
