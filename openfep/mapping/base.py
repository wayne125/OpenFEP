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
