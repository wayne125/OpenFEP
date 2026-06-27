from pathlib import Path
import pytest
from openfep.prep.rbfe_prep import prepare_rbfe, RBFESystem
from openfep.mapping.base import AtomMapping
from rdkit import Chem

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
