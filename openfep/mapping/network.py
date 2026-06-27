from __future__ import annotations

import itertools

import networkx as nx
from rdkit.Chem import Mol

from openfep.mapping.base import AtomMapper, AtomMapping

_MIN_SCORE = 0.2  # Lomap default pruning threshold


class FEPNetwork:
    """
    Build and score a perturbation network from a list of molecules.
    Edges are weighted by mapping score (higher = better).
    """

    def __init__(self, molecules: list[Mol], mapper: AtomMapper) -> None:
        if len(molecules) < 2:
            raise ValueError("FEPNetwork requires at least 2 molecules")
        self._molecules = molecules
        self._mapper = mapper
        self._mappings: dict[tuple[int, int], AtomMapping] = {}

    def build(self) -> nx.Graph:
        """All-pairs mapping. Returns graph with nodes=mol-indices, edges=score."""
        g = nx.Graph()
        g.add_nodes_from(range(len(self._molecules)))
        for i, j in itertools.combinations(range(len(self._molecules)), 2):
            mapping = self._mapper.map(self._molecules[i], self._molecules[j])
            self._mappings[(i, j)] = mapping
            if mapping.score >= _MIN_SCORE:
                g.add_edge(i, j, weight=mapping.score, mapping=mapping)
        return g

    def optimal_network(self) -> nx.Graph:
        """MST on weights (maximise mapping quality)."""
        g = self.build()
        if g.number_of_edges() == 0:
            raise ValueError(
                "FEP network has no edges: all pairwise mapping scores are below "
                f"the minimum threshold ({_MIN_SCORE}). Try a less strict mapper "
                "or add more structurally similar molecules."
            )
        # nx.maximum_spanning_tree maximises weight directly
        return nx.maximum_spanning_tree(g, weight="weight")

    def rbfe_pairs(self) -> list[tuple[Mol, Mol, AtomMapping]]:
        """Return (mol_a, mol_b, mapping) for each edge in the optimal network."""
        tree = self.optimal_network()
        pairs = []
        for i, j in tree.edges():
            key = (min(i, j), max(i, j))
            mapping = self._mappings[key]
            pairs.append((self._molecules[i], self._molecules[j], mapping))
        return pairs
