from __future__ import annotations

import lomap
from rdkit.Chem import Mol

from openfep.mapping.base import AtomMapper, AtomMapping


class LomapMapper(AtomMapper):
    """MCS-based atom mapper using Lomap2."""

    def __init__(self, max3d: float = 1.0) -> None:
        self._max3d = max3d

    def map(self, mol_a: Mol, mol_b: Mol) -> AtomMapping:
        pairs: list[tuple[int, int]] = list(
            lomap.MCS.getMapping(mol_a, mol_b, hydrogens=False)
        )
        a_to_b: dict[int, int] = dict(pairs)
        n_heavy_a = sum(1 for a in mol_a.GetAtoms() if a.GetAtomicNum() != 1)
        n_heavy_b = sum(1 for a in mol_b.GetAtoms() if a.GetAtomicNum() != 1)
        denom = max(n_heavy_a, n_heavy_b)
        score = len(a_to_b) / denom if denom > 0 else 0.0
        return AtomMapping(mol_a=mol_a, mol_b=mol_b, a_to_b=a_to_b, score=score)
