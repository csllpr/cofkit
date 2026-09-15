"""Conservative valence-preserving repair of quinoid imine bond assignments.

A finite quotient graph admits double-bond matchings that run through several
periodic precursor copies. An unconstrained bond-order converter can therefore
make short imine links single and the longer aryl-N bonds double. Local linkage
detection must see a consistent assignment before cutting the graph.
"""

from __future__ import annotations

import math
from collections import Counter
from typing import Mapping

import networkx as nx
from rdkit import Chem

# Contact lengths are ranked for the tie-break in integer micro-angstroms
# inside the [1.0, 1.7] angstrom window, matching the CIF's own precision;
# longer or unknown contacts score zero.
_LENGTH_REFERENCE_MICRO = 1_700_000
_MAX_LENGTH_SCORE = 700_000


def _length_score(distance: float | None) -> int:
    """Rank contact lengths so that shorter bonds attract the double bonds."""

    if distance is None or not math.isfinite(distance):
        return 0
    score = _LENGTH_REFERENCE_MICRO - round(distance * 1_000_000)
    return max(0, min(_MAX_LENGTH_SCORE, score))


def _match_residual(residual: nx.Graph) -> set[frozenset[int]]:
    """Keep incoming doubles first, then put doubles on the shortest contacts."""

    # One preserved double bond must outweigh every length score in the same
    # residual, so the primary objective stays lexicographically first.
    preserve_step = (_MAX_LENGTH_SCORE + 1) * (residual.number_of_nodes() // 2)
    weighted = nx.Graph()
    weighted.add_nodes_from(sorted(residual))
    for first, second, data in residual.edges(data=True):
        weighted.add_edge(
            first,
            second,
            weight=preserve_step * data["preserve"] + data["length"],
        )
    return {
        frozenset(edge) for edge in nx.max_weight_matching(weighted, maxcardinality=True)
    }


def normalize_imine_bond_orders(
    mol,
    distances: Mapping[frozenset[int], float],
    *,
    edge_multiplicities: Mapping[frozenset[int], int] | None = None,
) -> dict[str, object]:
    """Repair only neutral, H-explicit imines with decisive length evidence.

    Every participating C/N atom already has exactly one formal double bond.
    A perfect matching, constrained by short CH=N / longer N-C contacts, moves
    those bonds without changing valences, H counts, charges, or connectivity.
    The matching re-solves the whole conjugated component, so other C=C and
    C-N double bonds inside that component can move as well, and
    ``changed_bonds`` is typically larger than ``restored_imine_bonds``. Among
    equally conservative assignments the one that places the double bonds on
    the shorter contacts wins, so the outcome follows the contact geometry
    instead of the atom row order; only assignments whose contact-length sums
    are exactly equal stay interchangeable. Bonds outside the matching graph
    (aromatic, triple and carbonyl bonds) stay fixed. Unsatisfiable components
    are left intact and reported. Periodic edge multiplicities are included in
    the valence check. If a proposal changes periodic valences, retry with all
    parallel-image bond orders fixed before rejecting that component.
    No filenames, precursor identities, expected
    species counts, or topology hints are used.
    """
    multiplicities = edge_multiplicities or {}
    eligible = set()
    for atom in mol.GetAtoms():
        orders = [bond.GetBondTypeAsDouble() for bond in atom.GetBonds()]
        if (
            atom.GetAtomicNum() in {6, 7}
            and atom.GetFormalCharge() == 0
            and not atom.GetIsAromatic()
            and orders.count(2.0) == 1
            and all(order in {1.0, 2.0} for order in orders)
            and all(bond.GetStereo() == Chem.BondStereo.STEREONONE for bond in atom.GetBonds())
        ):
            eligible.add(atom.GetIdx())

    graph = nx.Graph()
    graph.add_nodes_from(sorted(eligible))
    original = set()
    for bond in mol.GetBonds():
        first, second = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if first in eligible and second in eligible:
            double = bond.GetBondTypeAsDouble() == 2.0
            graph.add_edge(
                first,
                second,
                preserve=int(double),
                length=_length_score(distances.get(frozenset((first, second)))),
            )
            if double:
                original.add(frozenset((first, second)))

    # An atom whose double-bond partner is outside this graph cannot be
    # rematched. This includes C=O and charged/unsupported partner environments.
    covered = {node for edge in original for node in edge}
    graph.remove_nodes_from(eligible - covered)
    Chem.GetSymmSSSR(mol)
    targets = set()
    for first, second in graph.edges:
        carbon, nitrogen = mol.GetAtomWithIdx(first), mol.GetAtomWithIdx(second)
        if carbon.GetAtomicNum() == 7:
            carbon, nitrogen = nitrogen, carbon
        if (carbon.GetAtomicNum(), nitrogen.GetAtomicNum()) != (6, 7):
            continue
        carbon_idx, nitrogen_idx = carbon.GetIdx(), nitrogen.GetIdx()
        bond = mol.GetBondBetweenAtoms(carbon_idx, nitrogen_idx)
        if any(bond.IsInRingSize(size) for size in range(3, 7)):
            continue
        # Explicit H is required: missing-H CIFs cannot justify this repair.
        if Counter(atom.GetAtomicNum() for atom in carbon.GetNeighbors()) != Counter((1, 6, 7)):
            continue
        if Counter(atom.GetAtomicNum() for atom in nitrogen.GetNeighbors()) != Counter((6, 6)):
            continue
        first_instance = carbon.GetProp("instance_id") if carbon.HasProp("instance_id") else ""
        second_instance = nitrogen.GetProp("instance_id") if nitrogen.HasProp("instance_id") else ""
        if first_instance and first_instance == second_instance:
            continue
        anchor = next(atom for atom in nitrogen.GetNeighbors() if atom.GetIdx() != carbon_idx)
        pair = frozenset((carbon_idx, nitrogen_idx))
        short = distances.get(pair, float("nan"))
        long = distances.get(frozenset((nitrogen_idx, anchor.GetIdx())), float("nan"))
        if 1.20 <= short <= 1.34 and long >= short + 0.06:
            targets.add(pair)

    report: dict[str, object] = {
        "strategy": "geometry_constrained_valence_matching",
        "candidate_imine_bonds": len(targets),
        "restored_imine_bonds": 0,
        "changed_bonds": [],
        "unresolved_components": 0,
    }
    if not targets - original:
        return report

    selected = set(original)
    for component in nx.connected_components(graph):
        forced = {edge for edge in targets if edge <= component}
        if not forced - original:
            continue
        incoming = {edge for edge in original if edge <= component}
        component_graph = graph.subgraph(component)
        proposal = _match_component(component_graph, forced)
        if proposal is not None and not _preserves_periodic_valences(incoming, proposal, multiplicities):
            parallel = {
                frozenset((first, second))
                for first, second in component_graph.edges
                if multiplicities.get(frozenset((first, second)), 1) > 1
            }
            retry_graph = component_graph.copy()
            retry_graph.remove_edges_from(tuple(edge) for edge in parallel - incoming)
            proposal = _match_component(retry_graph, forced | (parallel & incoming))
        if proposal is None or not _preserves_periodic_valences(incoming, proposal, multiplicities):
            report["unresolved_components"] += 1
            continue
        selected.difference_update(incoming)
        selected.update(proposal)

    changes = []
    for edge in sorted(original ^ selected, key=lambda edge: tuple(sorted(edge))):
        first, second = sorted(edge)
        bond = mol.GetBondBetweenAtoms(first, second)
        new_order = 2 if edge in selected else 1
        changes.append({"atoms": [first, second], "before": int(bond.GetBondTypeAsDouble()), "after": new_order})
        bond.SetBondType(Chem.BondType.DOUBLE if new_order == 2 else Chem.BondType.SINGLE)
    mol.UpdatePropertyCache(strict=False)
    report["changed_bonds"] = changes
    report["restored_imine_bonds"] = len((targets - original) & selected)
    return report


def _match_component(graph: nx.Graph, forced: set[frozenset[int]]) -> set[frozenset[int]] | None:
    endpoints = {node for edge in forced for node in edge}
    if len(endpoints) != 2 * len(forced) or any(not graph.has_edge(*edge) for edge in forced):
        return None
    residual = graph.subgraph(set(graph) - endpoints)
    matching = set()
    # Pinning linkage bonds usually separates the small precursor cores.
    for nodes in nx.connected_components(residual):
        matching.update(_match_residual(residual.subgraph(nodes)))
    if len(matching) * 2 != len(residual):
        return None
    return matching | forced


def _preserves_periodic_valences(original, selected, multiplicities) -> bool:
    # Unchanged orders cancel. Each reassigned quotient edge changes every
    # distinct periodic-image bond represented by that edge.
    delta = Counter()
    for edge in original ^ selected:
        change = multiplicities.get(edge, 1) * (1 if edge in selected else -1)
        for node in edge:
            delta[node] += change
    return all(change == 0 for change in delta.values())
