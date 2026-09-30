"""Spatial node-shape classification for topology compatibility.

A 4-connecting monomer's node shape is classified from the geometry of its
embedded conformer: the four connector (reactive-motif) positions, measured
as vectors from their centroid. The sorted pairwise-angle fingerprint of
those vectors separates the supported families:

- planar centered quadrilaterals (square / rectangular): both opposite
  connector pairs are antiparallel (~180 degrees apart) and the remaining
  adjacent angles come in a beta / 180 - beta pairing; beta near 90 degrees
  is square, otherwise rectangular.
- tetrahedral: all six pairwise angles sit near the ideal tetrahedral
  angle.

Any insufficiency — a motif count other than four, missing or non-finite
coordinates, coincident or irregularly spaced connectors, or a fingerprint
matching no family — yields ``unknown``, and ``unknown`` never excludes a
topology (``shapes_compatible`` returns ``None``).

``classify_monomer_graph_shape`` is retained as a graph-automorphism family
label for diagnostics only: graph distances cannot distinguish a
square-planar node from a tetrahedral one, so it is not used for topology
filtering.
"""
from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache

from .geometry import dot, norm
from .model import MonomerSpec

SHAPE_SQUARE = "square"
SHAPE_RECTANGULAR = "rectangular"
SHAPE_TETRAHEDRAL = "tetrahedral"
SHAPE_UNKNOWN = "unknown"
SHAPE_C4_TD = "c4_td"
SHAPE_D2H = "d2h"
KNOWN_SHAPE_LABELS = frozenset({SHAPE_SQUARE, SHAPE_RECTANGULAR, SHAPE_TETRAHEDRAL})

# Derived: ideal tetrahedral angle, acos(-1/3) in degrees.
TETRAHEDRAL_NODE_ANGLE_DEGREES = math.degrees(math.acos(-1.0 / 3.0))

# Heuristic — pending calibration: relative (max - min) / mean spread of the
# centroid-to-connector radii tolerated for a symmetric node. Minimized
# tetraaminobenzene / tetraphenylmethane nodes measure <= 0.01; distorted
# conformers of nominally square nodes (e.g. a phenyl arm flipped out of the
# porphyrin plane) exceed this and honestly classify as unknown.
NODE_SHAPE_RADIAL_SPREAD_TOLERANCE = 0.15
# Heuristic — pending calibration: deviation from 180 degrees tolerated for
# an opposite connector pair of a planar centered node (embedded
# 1,2,4,5-tetraaminobenzene measures ~177 degrees).
NODE_SHAPE_OPPOSITE_ANGLE_TOLERANCE_DEGREES = 15.0
# Heuristic — pending calibration: consistency window for the beta /
# 180 - beta adjacent-angle pairing of a planar quadrilateral node.
NODE_SHAPE_PLANAR_PAIRING_TOLERANCE_DEGREES = 12.0
# Heuristic — pending calibration: beta within this of 90 degrees classifies
# square rather than rectangular; the D2h reference node (1,2,4,5-substituted
# benzene) measures beta ~ 60 degrees.
NODE_SHAPE_SQUARE_ANGLE_TOLERANCE_DEGREES = 15.0
# Heuristic — pending calibration: window around the ideal tetrahedral
# angle; minimized tetraphenylmethane nodes measure 97-116 degrees.
NODE_SHAPE_TETRAHEDRAL_ANGLE_TOLERANCE_DEGREES = 15.0
# Heuristic — pending calibration: planar adjacent angles below this mean
# nearly collinear or coincident connectors, not a meaningful quadrilateral.
NODE_SHAPE_MIN_PLANAR_BETA_DEGREES = 30.0
# Derived numerical guard: a mean connector radius at or below this carries
# no directional evidence.
_NODE_SHAPE_MIN_MEAN_RADIUS_ANGSTROM = 1e-6


@dataclass(frozen=True)
class NodeShapeSignature:
    label: str


def unknown_signature() -> NodeShapeSignature:
    return NodeShapeSignature(SHAPE_UNKNOWN)


def classify_monomer_graph_shape(spec: MonomerSpec) -> str:
    """Graph-automorphism family label (``c4_td`` / ``d2h`` / ``unknown``).

    Diagnostic only: graph distances and color refinement cannot distinguish
    a square-planar node from a tetrahedral one, so topology filtering uses
    the geometric ``classify_monomer_node_shape`` instead.
    """
    if len(spec.motifs) != 4 or not spec.bonds:
        return SHAPE_UNKNOWN
    adj = {}
    for a, b, _ in spec.bonds:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    sites = []
    for motif in spec.motifs:
        members = set(motif.atom_ids)
        candidates = [a for a in members if a < len(spec.atom_symbols) and spec.atom_symbols[a] != "H"]
        if not candidates:
            return SHAPE_UNKNOWN
        sites.append(max(candidates, key=lambda a: sum(n not in members for n in adj.get(a, ()))))
    if len(set(sites)) != 4:
        return SHAPE_UNKNOWN
    distances = []
    for i, source in enumerate(sites):
        seen = {source: 0}
        queue = [source]
        for node in queue:
            for n in adj.get(node, ()):
                if n not in seen:
                    seen[n] = seen[node] + 1
                    queue.append(n)
        for target in sites[i + 1:]:
            if target not in seen:
                return SHAPE_UNKNOWN
            distances.append(seen[target])
    colors = {n: (spec.atom_symbols[n], len(adj.get(n, ()))) for n in adj}
    for _ in range(max(1, len(adj))):
        keys = {n: (colors[n], tuple(sorted(colors.get(x, ()) for x in adj.get(n, ())))) for n in adj}
        palette = {k: i for i, k in enumerate(sorted(set(keys.values()), key=repr))}
        colors = {n: palette[k] for n, k in keys.items()}
    orbit_counts = sorted(Counter(colors.get(s, -1) for s in sites).values())
    counts = sorted(Counter(distances).values())
    if counts == [6] or (orbit_counts == [4] and len(set(distances)) <= 2):
        return SHAPE_C4_TD
    if counts in ([2, 4], [2, 2, 2]):
        return SHAPE_D2H
    return SHAPE_UNKNOWN


def _connector_evidence(spec: MonomerSpec) -> tuple[float, ...] | None:
    """Sorted pairwise connector angles in degrees, or None when insufficient.

    Evidence is the four motif-frame origins (the connector/reactive-atom
    positions of the embedded conformer) measured from their centroid; the
    pairwise-angle fingerprint is invariant under rigid motion of the
    monomer. Connectors with wildly unequal radii are not a centered node
    and yield no evidence.
    """
    if len(spec.motifs) != 4:
        return None
    origins = []
    for motif in spec.motifs:
        origin = motif.frame.origin
        if len(origin) != 3 or not all(math.isfinite(component) for component in origin):
            return None
        origins.append(origin)
    center = tuple(sum(origin[axis] for origin in origins) / 4.0 for axis in range(3))
    vectors = [tuple(origin[axis] - center[axis] for axis in range(3)) for origin in origins]
    radii = [norm(vector) for vector in vectors]
    mean_radius = sum(radii) / 4.0
    if mean_radius <= _NODE_SHAPE_MIN_MEAN_RADIUS_ANGSTROM:
        return None
    if (max(radii) - min(radii)) / mean_radius > NODE_SHAPE_RADIAL_SPREAD_TOLERANCE:
        return None
    angles = []
    for i in range(4):
        for j in range(i + 1, 4):
            cosine = max(-1.0, min(1.0, dot(vectors[i], vectors[j]) / (radii[i] * radii[j])))
            angles.append(math.degrees(math.acos(cosine)))
    return tuple(sorted(angles))


def _classify_planar_angles(angles: tuple[float, ...]) -> str | None:
    """Square / rectangular verdict from a centered-quadrilateral fingerprint."""
    opposite_a, opposite_b = angles[4], angles[5]
    if min(opposite_a, opposite_b) < 180.0 - NODE_SHAPE_OPPOSITE_ANGLE_TOLERANCE_DEGREES:
        return None
    if abs(opposite_a - opposite_b) > NODE_SHAPE_OPPOSITE_ANGLE_TOLERANCE_DEGREES:
        return None
    beta_pair, complement_pair = angles[0:2], angles[2:4]
    if abs(beta_pair[0] - beta_pair[1]) > NODE_SHAPE_PLANAR_PAIRING_TOLERANCE_DEGREES:
        return None
    if abs(complement_pair[0] - complement_pair[1]) > NODE_SHAPE_PLANAR_PAIRING_TOLERANCE_DEGREES:
        return None
    beta = (beta_pair[0] + beta_pair[1]) / 2.0
    complement = (complement_pair[0] + complement_pair[1]) / 2.0
    if abs(beta + complement - 180.0) > NODE_SHAPE_PLANAR_PAIRING_TOLERANCE_DEGREES:
        return None
    if beta < NODE_SHAPE_MIN_PLANAR_BETA_DEGREES:
        return None
    if abs(beta - 90.0) <= NODE_SHAPE_SQUARE_ANGLE_TOLERANCE_DEGREES:
        return SHAPE_SQUARE
    return SHAPE_RECTANGULAR


def classify_monomer_node_shape(spec) -> NodeShapeSignature:
    """Spatial connector-shape signature from conformer geometry.

    Returns ``unknown`` whenever the geometric evidence is insufficient;
    unknown excludes no topology.
    """
    angles = _connector_evidence(spec)
    if angles is None:
        return unknown_signature()
    planar_label = _classify_planar_angles(angles)
    if planar_label is not None:
        return NodeShapeSignature(planar_label)
    if all(abs(angle - TETRAHEDRAL_NODE_ANGLE_DEGREES) <= NODE_SHAPE_TETRAHEDRAL_ANGLE_TOLERANCE_DEGREES for angle in angles):
        return NodeShapeSignature(SHAPE_TETRAHEDRAL)
    return unknown_signature()


@lru_cache(maxsize=None)
def classify_topology_node_shape(topology_id) -> NodeShapeSignature:
    return NodeShapeSignature(
        {"sql": SHAPE_SQUARE, "kgm": SHAPE_RECTANGULAR, "dia": SHAPE_TETRAHEDRAL}.get(topology_id, SHAPE_UNKNOWN)
    )


def shapes_compatible(a, b):
    """True/False for confident labels, None (no exclusion) when either side is unknown."""
    if a.label not in KNOWN_SHAPE_LABELS or b.label not in KNOWN_SHAPE_LABELS:
        return None
    return a.label == b.label
