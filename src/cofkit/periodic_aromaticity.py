"""Aromatic bond normalization with lattice-aware molecular ring discovery."""
from __future__ import annotations

from collections import deque
from collections.abc import Iterable
from numbers import Integral

from rdkit import Chem

from .bond_types import normalize_bond_order

# Cited: RDKit's default fused-ring aromaticity candidate bound,
# Code/GraphMol/Aromaticity.cpp, maxFusedAromaticRingSize. We apply this bound
# to all candidate rings in the lift; larger-ring aromaticity is outside this
# normalization policy (it is not inferred from a winding quotient cycle).
# https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/Aromaticity.cpp
AROMATIC_RING_MAX_ATOMS = 24
# Heuristic — pending calibration: resource ceilings for pathological periodic
# graphs. Exhaustion is explicit; partial ring searches never produce types.
PERIODIC_RING_MAX_SEARCH_STATES = 100_000
PERIODIC_RING_MAX_SHORTEST_PATHS = 4096
PERIODIC_AROMATICITY_MAX_LIFT_ATOMS = 100_000

Gain = tuple[int, int, int]
PeriodicBond = tuple[int, int, Gain, float | None]


class PeriodicAromaticityError(ValueError):
    """The supplied periodic model cannot be normalized without guessing."""


class _RingNeedsLift(Exception):
    def __init__(self, axis, separation):
        self.axis = axis
        self.separation = separation


def _add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def _rings(symbols, edges):
    """Find all shortest zero-gain cycles through each edge in the lift.

    States include the atom AND its lattice image. All tied shortest paths are
    retained, avoiding a ring choice based on atom order or integer basis.
    """
    adjacency = [[] for _ in symbols]
    for index, (u, v, gain, _) in enumerate(edges):
        adjacency[u].append((v, gain, index))
        adjacency[v].append((u, tuple(-x for x in gain), index))
    rings = {}
    zero = (0, 0, 0)
    for edge_index, (first, second, gain, _) in enumerate(edges):
        # Terminal edges (including explicit H) cannot belong to a ring.
        if len(adjacency[first]) < 2 or len(adjacency[second]) < 2:
            continue
        root, target = (second, gain), (first, zero)
        distances = {root: 0}
        parents = {}
        queue = deque([root])
        found = None
        while queue:
            state = queue.popleft()
            distance = distances[state]
            if (found is not None and distance >= found) or distance >= AROMATIC_RING_MAX_ATOMS - 1:
                continue
            atom, image = state
            for neighbor, shift, index in adjacency[atom]:
                next_state = (neighbor, _add(image, shift))
                if index == edge_index and {state, next_state} == {root, target}:
                    continue
                next_distance = distance + 1
                if next_state not in distances:
                    distances[next_state] = next_distance
                    parents[next_state] = []
                    queue.append(next_state)
                    if len(distances) > PERIODIC_RING_MAX_SEARCH_STATES:
                        raise PeriodicAromaticityError("Periodic ring search exceeded its state budget.")
                if distances[next_state] == next_distance:
                    parents[next_state].append((state, index))
                if next_state == target:
                    found = next_distance
        if found is None:
            continue
        stack = [(target, [target], [])]
        paths = 0
        while stack:
            state, states, bonds = stack.pop()
            if state == root:
                paths += 1
                if paths > PERIODIC_RING_MAX_SHORTEST_PATHS:
                    raise PeriodicAromaticityError("Periodic ring search exceeded its shortest-path budget.")
                atoms = [item[0] for item in states]
                images = {}
                for atom, image in states:
                    if atom in images:
                        difference = tuple(x - y for x, y in zip(image, images[atom]))
                        axis = next(i for i, value in enumerate(difference) if value)
                        raise _RingNeedsLift(axis, difference[axis])
                    images[atom] = image
                ring_bonds = bonds + [edge_index]
                rings.setdefault(tuple(sorted(ring_bonds)), (atoms, ring_bonds))
                continue
            for parent, index in parents[state]:
                stack.append((parent, states + [parent], bonds + [index]))
    return tuple(rings.values())


def normalize_periodic_bond_orders(
    symbols: Iterable[str], edges: Iterable[PeriodicBond], *, formal_charges: Iterable[int] | None = None,
) -> tuple[float, ...]:
    """Return aromatic-normalized orders in input edge order.

    ``edges`` contain zero-based endpoints, integer relative lattice gain and
    explicit bond order. Coordinates, H placement and the input records are
    never changed. Self/parallel-image bonds require a larger explicit cell,
    as on the LAMMPS input path. Candidate rings have at most
    ``AROMATIC_RING_MAX_ATOMS`` atoms. No missing bond orders are inferred.
    All supplied H must be explicit; charges default to neutral.
    """
    symbols = tuple(symbols)
    edges = tuple((u, v, tuple(gain), order) for u, v, gain, order in edges)
    charges = tuple(formal_charges) if formal_charges is not None else (0,) * len(symbols)
    if len(charges) != len(symbols):
        raise PeriodicAromaticityError("Formal charges must be aligned with atom symbols.")
    editable = Chem.RWMol()
    for symbol, charge in zip(symbols, charges):
        atom = Chem.Atom(symbol)
        atom.SetNoImplicit(True)
        atom.SetFormalCharge(charge)
        editable.AddAtom(atom)
    seen = set()
    for u, v, gain, order in edges:
        if not (0 <= u < len(symbols) and 0 <= v < len(symbols)):
            raise PeriodicAromaticityError("Periodic bond endpoints are outside the atom list.")
        if len(gain) != 3 or not all(isinstance(x, Integral) for x in gain):
            raise PeriodicAromaticityError("Periodic bond translations must contain three integers.")
        if u == v or frozenset((u, v)) in seen:
            raise PeriodicAromaticityError("Self/parallel-image bonds require an explicit larger supercell.")
        seen.add(frozenset((u, v)))
        if order is None:
            raise PeriodicAromaticityError("Force-field typing requires explicit bond orders for every bond.")
        kind = {1.: Chem.BondType.SINGLE, 1.5: Chem.BondType.AROMATIC,
                2.: Chem.BondType.DOUBLE, 3.: Chem.BondType.TRIPLE}.get(normalize_bond_order(order))
        if kind is None:
            raise PeriodicAromaticityError(f"Unsupported explicit bond order for force-field typing: {order}.")
        editable.AddBond(u, v, kind)
    mol = editable.GetMol()
    try:
        rings = _rings(symbols, edges)
    except _RingNeedsLift as need:
        # Choose a cover that separates the actual repeated ring images. This
        # is topology-driven context growth, not a fixed supercell assumption.
        factor = 2
        while need.separation % factor == 0:
            factor += 1
        if len(symbols) * factor > PERIODIC_AROMATICITY_MAX_LIFT_ATOMS:
            raise PeriodicAromaticityError("Periodic aromaticity exceeded its lifted-atom budget.")
        expanded = []
        n = len(symbols)
        for cell in range(factor):
            for u, v, gain, order in edges:
                period, remainder = divmod(cell + gain[need.axis], factor)
                shift = list(gain)
                shift[need.axis] = period
                expanded.append((cell * n + u, remainder * n + v, tuple(shift), order))
        normalized = normalize_periodic_bond_orders(symbols * factor, expanded, formal_charges=charges * factor)
        result = normalized[:len(edges)]
        if any(normalized[i] != result[i % len(edges)] for i in range(len(normalized))):
            raise PeriodicAromaticityError("Aromatic assignments disagree between equivalent periodic ring images.")
        return result
    try:
        mol.UpdatePropertyCache(strict=False)
        # Populate only rings closed in the infinite lift, before either
        # kekulization or aromaticity. Never run quotient SSSR afterwards.
        info = mol.GetRingInfo()
        for atoms, bonds in rings:
            info.AddRing(atoms, bonds)
        if not rings:
            # A no-ring graph with explicit aromatic notation is unresolved.
            if any(b.GetIsAromatic() for b in mol.GetBonds()):
                raise PeriodicAromaticityError("Explicit aromatic bonds have no closed molecular ring.")
            return tuple(float(normalize_bond_order(edge[3])) for edge in edges)
        Chem.Kekulize(mol, clearAromaticFlags=True)
        mol.UpdatePropertyCache(strict=False)
        Chem.SetAromaticity(mol)
    except (ValueError, RuntimeError) as exc:
        raise PeriodicAromaticityError(f"Cannot normalize aromatic bond notation: {exc}") from exc
    return tuple(b.GetBondTypeAsDouble() for b in mol.GetBonds())
