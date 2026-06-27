"""
Tests for openfep.msj.rbfe_writer.write_rbfe_msj.

Correctness constraints verified here:
- complex.msj uses custom_charge_mode="assign"
- solvent.msj uses custom_charge_mode="keep"
- Equil stages: all-but-last → polarization_restraint="full", last → "decay"
- lambda_hopping block does NOT contain polarization_restraints
- backbone_asl appears in the first minimize_stage call (complex leg)
- Charged system: neutralize_system=true and make_alchemical_water=true
- Buffer widths come from constants (not hardcoded)
- main.msj has fep_mapper block with graph_file, atom_mapping, receptor keys
"""
from pathlib import Path
import pytest

from openfep.msj.rbfe_writer import write_rbfe_msj
from openfep.prep.rbfe_prep import RBFESystem
from openfep.constants import (
    COMPLEX_BUFFER_WIDTH,
    NET_CHARGE_COMPLEX_BUFFER_WIDTH,
    SOLVENT_BUFFER_WIDTH,
    BACKBONE_ASL,
)


# ── Fixtures ──────────────────────────────────────────────────────────────────


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
    return RBFESystem(
        fmp_path=tmp_path / "rbfe.fmp",
        atom_mapping_paths={},
        receptor_count=1,
        is_charged=True,
        buffer_width=NET_CHARGE_COMPLEX_BUFFER_WIDTH,
        has_membrane=False,
    )


@pytest.fixture
def membrane_system(tmp_path):
    return RBFESystem(
        fmp_path=tmp_path / "rbfe.fmp",
        atom_mapping_paths={},
        receptor_count=2,
        is_charged=False,
        buffer_width=COMPLEX_BUFFER_WIDTH,
        has_membrane=True,
    )


# ── File-level tests ──────────────────────────────────────────────────────────


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


def test_output_dir_created_if_absent(neutral_system, tmp_path):
    out = tmp_path / "new_subdir"
    assert not out.exists()
    write_rbfe_msj(neutral_system, out)
    assert out.exists()


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


def test_complex_msj_has_two_minimize_stages(neutral_system, tmp_path):
    content = write_rbfe_msj(neutral_system, tmp_path)["complex"].read_text()
    assert content.count("minimize {") == 2


def test_complex_msj_has_build_geometry(neutral_system, tmp_path):
    content = write_rbfe_msj(neutral_system, tmp_path)["complex"].read_text()
    assert "build_geometry {" in content


def test_complex_msj_has_assign_forcefield(neutral_system, tmp_path):
    content = write_rbfe_msj(neutral_system, tmp_path)["complex"].read_text()
    assert "assign_forcefield {" in content


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


def test_solvent_msj_trim_present(neutral_system, tmp_path):
    content = write_rbfe_msj(neutral_system, tmp_path)["solvent"].read_text()
    assert "trim {" in content


def test_solvent_msj_starts_with_task(neutral_system, tmp_path):
    content = write_rbfe_msj(neutral_system, tmp_path)["solvent"].read_text()
    assert content.strip().startswith("task {")


def test_solvent_msj_no_backbone_asl(neutral_system, tmp_path):
    """Solvent leg has no protein; backbone ASL must not appear."""
    content = write_rbfe_msj(neutral_system, tmp_path)["solvent"].read_text()
    assert BACKBONE_ASL not in content


def test_solvent_msj_no_gcmc(neutral_system, tmp_path):
    content = write_rbfe_msj(neutral_system, tmp_path)["solvent"].read_text()
    assert "gcmc {" not in content


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


def test_main_msj_has_atom_mapping_key(neutral_system, tmp_path):
    content = write_rbfe_msj(neutral_system, tmp_path)["main"].read_text()
    assert "atom_mapping" in content


def test_main_msj_starts_with_task(neutral_system, tmp_path):
    content = write_rbfe_msj(neutral_system, tmp_path)["main"].read_text()
    assert content.strip().startswith("task {")


# ── charged system ────────────────────────────────────────────────────────────


def test_charged_complex_msj_neutralize(charged_system, tmp_path):
    content = write_rbfe_msj(charged_system, tmp_path)["complex"].read_text()
    assert "neutralize_system = true" in content


def test_charged_complex_msj_alchemical_water(charged_system, tmp_path):
    content = write_rbfe_msj(charged_system, tmp_path)["complex"].read_text()
    assert "make_alchemical_water = true" in content


def test_charged_solvent_msj_neutralize(charged_system, tmp_path):
    content = write_rbfe_msj(charged_system, tmp_path)["solvent"].read_text()
    assert "neutralize_system = true" in content


def test_charged_solvent_msj_alchemical_water(charged_system, tmp_path):
    content = write_rbfe_msj(charged_system, tmp_path)["solvent"].read_text()
    assert "make_alchemical_water = true" in content


def test_neutral_complex_msj_no_neutralize(neutral_system, tmp_path):
    content = write_rbfe_msj(neutral_system, tmp_path)["complex"].read_text()
    assert "neutralize_system" not in content


# ── membrane system ───────────────────────────────────────────────────────────


def test_membrane_complex_msj_has_gcmc(membrane_system, tmp_path):
    content = write_rbfe_msj(membrane_system, tmp_path)["complex"].read_text()
    assert "gcmc {" in content


def test_no_membrane_complex_msj_no_gcmc(neutral_system, tmp_path):
    content = write_rbfe_msj(neutral_system, tmp_path)["complex"].read_text()
    assert "gcmc {" not in content


def test_membrane_main_msj_receptor_count(membrane_system, tmp_path):
    content = write_rbfe_msj(membrane_system, tmp_path)["main"].read_text()
    assert "receptor = 2" in content


# ── lambda customization ──────────────────────────────────────────────────────


def test_custom_lambda_count_reflected(neutral_system, tmp_path):
    """Custom lambda list length appears in the task stage."""
    custom_lambdas = [0.0, 0.25, 0.5, 0.75, 1.0]
    content = write_rbfe_msj(neutral_system, tmp_path, lambdas=custom_lambdas)["complex"].read_text()
    # n_lambda = 5 → "default:5" in fep.lambda value
    assert "default:5" in content


def test_default_lambda_count(neutral_system, tmp_path):
    """Default rbfe_lambdas() uses RBFE_DEFAULT_LAMBDAS = 12 windows."""
    from openfep.constants import RBFE_DEFAULT_LAMBDAS
    content = write_rbfe_msj(neutral_system, tmp_path)["complex"].read_text()
    assert f"default:{RBFE_DEFAULT_LAMBDAS}" in content


# ── polarization ordering (I4) ─────────────────────────────────────────────────


def test_complex_msj_polarization_ordering(neutral_system, tmp_path):
    """Last equil stage must be decay; all prior must be full."""
    content = write_rbfe_msj(neutral_system, tmp_path)["complex"].read_text()
    # Split into stage blocks by double-newline
    blocks = [b.strip() for b in content.split("\n\n") if b.strip()]
    simulate_blocks = [b for b in blocks if b.startswith("simulate {")]
    assert len(simulate_blocks) >= 2  # need at least 2 equil stages
    for block in simulate_blocks[:-1]:
        assert '"full"' in block, f"Non-last simulate should have full: {block}"
    assert '"decay"' in simulate_blocks[-1], "Last simulate should have decay"
    assert '"full"' not in simulate_blocks[-1], "Last simulate should NOT have full"


# ── backbone ASL ordering (I5) ─────────────────────────────────────────────────


def test_complex_msj_backbone_asl_only_in_first_minimize(neutral_system, tmp_path):
    """First minimize has backbone restraint; second is free."""
    content = write_rbfe_msj(neutral_system, tmp_path)["complex"].read_text()
    blocks = [b.strip() for b in content.split("\n\n") if b.strip()]
    minimize_blocks = [b for b in blocks if b.startswith("minimize {")]
    assert len(minimize_blocks) == 2
    assert BACKBONE_ASL in minimize_blocks[0], "First minimize should have backbone ASL"
    assert BACKBONE_ASL not in minimize_blocks[1], "Second minimize should be free"
