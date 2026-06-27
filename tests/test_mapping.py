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


# ── LomapMapper ───────────────────────────────────────────────────────────────

from openfep.mapping.lomap_mapper import LomapMapper
from openfep.mapping.network import FEPNetwork
import networkx as nx


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
    net = FEPNetwork(three_mols, LomapMapper())
    g = net.build()
    assert isinstance(g, nx.Graph)


def test_fepnetwork_build_has_all_nodes(three_mols):
    g = FEPNetwork(three_mols, LomapMapper()).build()
    assert g.number_of_nodes() == 3


def test_fepnetwork_optimal_network_is_spanning_tree(three_mols):
    net = FEPNetwork(three_mols, LomapMapper())
    tree = net.optimal_network()
    assert nx.is_tree(tree)


def test_fepnetwork_rbfe_pairs_length(three_mols):
    pairs = FEPNetwork(three_mols, LomapMapper()).rbfe_pairs()
    # MST of 3 nodes has exactly 2 edges
    assert len(pairs) == 2


def test_fepnetwork_rbfe_pairs_contain_mappings(three_mols):
    pairs = FEPNetwork(three_mols, LomapMapper()).rbfe_pairs()
    for mol_a, mol_b, mapping in pairs:
        assert isinstance(mapping, AtomMapping)
        assert mapping.score > 0.0


def test_fepnetwork_requires_two_mols(benzene):
    with pytest.raises(ValueError):
        FEPNetwork([benzene], LomapMapper())


def test_fepnetwork_no_edges_raises_value_error():
    """Dissimilar molecules with no common scaffold should raise ValueError."""
    benzene = _mol3d("c1ccccc1")
    hexane = _mol3d("CCCCCC")
    with pytest.raises(ValueError):
        FEPNetwork([benzene, hexane], LomapMapper()).optimal_network()
