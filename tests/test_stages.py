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
    minimize_stage,
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
    assert '"OPLS4"' in text


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
    assert '"cubic"' not in text   # box_shape only written when neutralize=True


def test_build_geometry_cubic_with_neutralize():
    text = render(build_geometry_stage(buffer=8.0, neutralize=True, box_shape="cubic"))
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


# ── minimize_stage ────────────────────────────────────────────────────────────────

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


# ── ns→ps time conversion ──────────────────────────────────────────────────────

def test_simulate_stage_time_converted_to_ps():
    from openfep.msj.sea import simulate_stage, render
    block = simulate_stage(time_ns=0.12, temperature=300.0, ensemble="NPT")
    rendered = render(block)
    assert "time = 120.0" in rendered  # 0.12 ns → 120.0 ps


def test_lambda_hopping_time_converted_to_ps():
    from openfep.msj.sea import lambda_hopping_stage, render
    block = lambda_hopping_stage(time_ns=5.0, ensemble="NPT")
    rendered = render(block)
    assert "time = 5000.0" in rendered  # 5.0 ns → 5000.0 ps
