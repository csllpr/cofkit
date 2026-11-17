"""Ring discovery must distinguish a chemical closure from a lattice winding."""
import pytest
from rdkit import Chem

from cofkit import periodic_aromaticity as pa


def _phenylene():
    molecule = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    Chem.Kekulize(molecule, clearAromaticFlags=True)
    # Replace the para H atoms with the inter-cell phenylene link.
    editable = Chem.RWMol(molecule)
    hydrogens = [next(a.GetIdx() for a in editable.GetAtomWithIdx(i).GetNeighbors()
                      if a.GetSymbol() == "H") for i in [0, 3]]
    for i in sorted(hydrogens, reverse=True):
        editable.RemoveAtom(i)
    molecule = editable.GetMol()
    edges = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), (0, 0, 0), b.GetBondTypeAsDouble())
             for b in molecule.GetBonds()]
    return [a.GetSymbol() for a in molecule.GetAtoms()], edges + [(0, 3, (1, 0, 0), 1.)]


def test_true_ring_is_recovered_when_quotient_sssr_only_contains_winding_shortcuts():
    symbols, edges = _phenylene()
    result = pa.normalize_periodic_bond_orders(symbols, edges)
    assert result[:6] == (1.5,) * 6
    assert result[-1] == 1.
    # Reverse both atom and edge listings, including gain orientation.
    n = len(symbols)
    reordered = [(n - 1 - v, n - 1 - u, tuple(-x for x in g), order)
                 for u, v, g, order in reversed(edges)]
    assert pa.normalize_periodic_bond_orders(list(reversed(symbols)), reordered) == tuple(reversed(result))


def test_resource_exhaustion_does_not_return_partial_aromatic_assignments(monkeypatch):
    symbols, edges = _phenylene()
    monkeypatch.setattr(pa, "PERIODIC_RING_MAX_SEARCH_STATES", 1)
    with pytest.raises(pa.PeriodicAromaticityError, match="state budget"):
        pa.normalize_periodic_bond_orders(symbols, edges)


def test_exocyclic_multiple_bond_is_not_aromatized_by_winding_cycles():
    symbols, edges = _phenylene()
    ring_orders = (1., 2., 1., 1., 2., 1.)
    quinoid = [(u, v, gain, ring_orders[i] if i < 6 else (2. if i == len(edges) - 1 else order))
               for i, (u, v, gain, order) in enumerate(edges)]
    aromatic = pa.normalize_periodic_bond_orders(symbols, edges)
    normalized = pa.normalize_periodic_bond_orders(symbols, quinoid)
    assert aromatic[:6] == (1.5,) * 6
    # RDKit can aromatize the internal carbocycle in both representations,
    # but the exocyclic multiple bond must remain explicit, not a false ring.
    assert aromatic[-1] == 1.
    assert normalized[-1] == 2.
    assert normalized != aromatic


@pytest.mark.parametrize("edges", [
    [(0, 0, (1, 0, 0), 1.)],
    [(0, 1, (0, 0, 0), 1.), (0, 1, (1, 0, 0), 1.)],
])
def test_unrepresentable_self_or_parallel_edges_are_explicit(edges):
    with pytest.raises(pa.PeriodicAromaticityError, match="larger supercell"):
        pa.normalize_periodic_bond_orders(["C", "C"], edges)


def _small_ctf_cell():
    # CTF-0 reference case 0183: small enough that a genuine pore cycle visits
    # distinct images of the same atom. Retain the original explicit H/order
    # evidence, with zero-based indices in CIF row order.
    symbols = ["C", "C", "C", "H", "N", "C", "C", "C", "H", "N", "C", "C", "C", "H", "N"]
    rows = [(0, 6, 2), (0, 1, 1), (0, 12, 1), (1, 3, 1), (1, 10, 2),
            (2, 5, 1), (2, 14, 2), (2, 4, 1), (4, 7, 2), (5, 6, 1),
            (5, 11, 2), (6, 8, 1), (7, 9, 1), (7, 10, 1), (9, 12, 2),
            (10, 11, 1), (11, 13, 1), (12, 14, 1)]
    gains = {(0, 12): (-1, 0, 0), (7, 10): (0, -1, 0)}
    return symbols, [(u, v, gains.get((u, v), (0, 0, 0)), float(order)) for u, v, order in rows]


def test_repeated_periodic_ring_images_are_lifted_and_folded_consistently():
    symbols, edges = _small_ctf_cell()
    orders = pa.normalize_periodic_bond_orders(symbols, edges)
    assert all(order == (1. if symbols[u] == "H" or symbols[v] == "H" else 1.5)
               for (u, v, _, _), order in zip(edges, orders))
    aromatic = [(u, v, gain, order) for (u, v, gain, _), order in zip(edges, orders)]
    assert pa.normalize_periodic_bond_orders(symbols, aromatic) == orders
    # Reverse input ordering and make a nontrivial unimodular basis change.
    changed = [(v, u, (-g[0] - g[1], -g[1], -g[2]), order) for u, v, g, order in reversed(edges)]
    assert pa.normalize_periodic_bond_orders(symbols, changed) == tuple(reversed(orders))


def test_adaptive_lift_budget_is_explicit(monkeypatch):
    symbols, edges = _small_ctf_cell()
    monkeypatch.setattr(pa, "PERIODIC_AROMATICITY_MAX_LIFT_ATOMS", len(symbols))
    with pytest.raises(pa.PeriodicAromaticityError, match="lifted-atom budget"):
        pa.normalize_periodic_bond_orders(symbols, edges)
