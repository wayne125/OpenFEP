from __future__ import annotations

from rdkit.Chem import Mol
from kartograf.atom_mapper import KartografAtomMapper as _Kartograf

from openfep.mapping.base import AtomMapper, AtomMapping


class KartografMapper(AtomMapper):
    """3D-geometry-based atom mapper using Kartograf."""

    def __init__(self, map_hydrogens: bool = False) -> None:
        self._mapper = _Kartograf(atom_map_hydrogens=map_hydrogens)

    def map(self, mol_a: Mol, mol_b: Mol) -> AtomMapping:
        a_to_b: dict[int, int] = self._mapper.suggest_mapping_from_rdmols(mol_a, mol_b)
        # Score: fraction of heavy atoms mapped (Tanimoto-like)
        n_heavy_a = sum(1 for a in mol_a.GetAtoms() if a.GetAtomicNum() != 1)
        n_heavy_b = sum(1 for a in mol_b.GetAtoms() if a.GetAtomicNum() != 1)
        mapped = len(a_to_b)
        denom = max(n_heavy_a, n_heavy_b)
        score = mapped / denom if denom > 0 else 0.0
        return AtomMapping(mol_a=mol_a, mol_b=mol_b, a_to_b=a_to_b, score=score)
