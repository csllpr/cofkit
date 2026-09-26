from __future__ import annotations

import sys
from collections import deque
from dataclasses import dataclass
from functools import lru_cache
from math import atan2, cos, pi, sin

import gemmi

from .geometry import Vec3, add, classify_2d_cell_parameters, norm, normalize, scale, sub
from .topology_index import TopologyDefinition
from .topologies import load_topology


# Topology ids whose unsupported-expansion warning has already been printed.
# list_supported_single_node_topology_ids runs inside per-pair batch loops, so
# the omission warning must fire once per topology id, not once per call.
_OMITTED_TOPOLOGY_WARNING_SEEN: set[str] = set()

_ANGLE_TOLERANCE = 1e-4


@dataclass(frozen=True)
class ExpandedSingleNodeNodeSite:
    id: str
    fractional_position: tuple[float, float, float]
    cartesian_position: Vec3
    directions: tuple[Vec3, ...]
    edge_ids: tuple[str, ...]
    sublattice: int | None = None


@dataclass(frozen=True)
class ExpandedSingleNodeEdge:
    id: str
    start_node_id: str
    end_node_id: str
    end_image: tuple[int, int, int]
    base_vector: Vec3
    center_fractional: tuple[float, float, float]


@dataclass(frozen=True)
class ExpandedSingleNodeTopology:
    topology_id: str
    connectivity: int
    cell: tuple[Vec3, Vec3, Vec3]
    metric_family: str
    node_sites: tuple[ExpandedSingleNodeNodeSite, ...]
    edge_sites: tuple[ExpandedSingleNodeEdge, ...]
    is_bipartite: bool


@dataclass(frozen=True)
class SingleNodeTopologyLayout:
    topology_id: str
    connectivity: int
    directions: tuple[Vec3, ...]
    metric_family: str
    inference_mode: str
    symmetry_orbit_size: int
    placement_model: str
    supports_current_builder: bool
    unsupported_reason: str | None
    supports_node_node: bool
    supports_node_linker: bool


_SUPPORTED_SINGLE_NODE_TOPOLOGIES: dict[str, dict[str, object]] = {
    "hcb": {
        "inference_mode": "expanded",
    },
    "hca": {
        "inference_mode": "expanded",
    },
    "fes": {
        "inference_mode": "expanded",
    },
    "fxt": {
        "inference_mode": "explicit",
    },
    "sql": {
        "inference_mode": "expanded",
    },
    "kgm": {
        "inference_mode": "expanded",
    },
    "htb": {
        "inference_mode": "expanded",
    },
    "hxl": {
        "inference_mode": "expanded",
    },
}


def list_supported_single_node_topology_ids(
    connectivity: int,
    *,
    pair_mode: str,
) -> tuple[str, ...]:
    ids: list[str] = []
    for topology_id in _SUPPORTED_SINGLE_NODE_TOPOLOGIES:
        try:
            layout = resolve_single_node_topology_layout(topology_id)
        except (KeyError, ValueError) as exc:
            if topology_id not in _OMITTED_TOPOLOGY_WARNING_SEEN:
                _OMITTED_TOPOLOGY_WARNING_SEEN.add(topology_id)
                print(
                    f"warning: single-node topology {topology_id!r} omitted from the supported list: {exc}",
                    file=sys.stderr,
                )
            continue
        if layout.connectivity != connectivity:
            continue
        if not layout.supports_current_builder:
            continue
        if pair_mode == "node-node" and not layout.supports_node_node:
            continue
        if pair_mode == "node-linker" and not layout.supports_node_linker:
            continue
        ids.append(topology_id)
    return tuple(ids)


@lru_cache(maxsize=None)
def resolve_single_node_topology_layout(topology_id: str) -> SingleNodeTopologyLayout:
    support = _SUPPORTED_SINGLE_NODE_TOPOLOGIES.get(topology_id)
    if support is None:
        raise KeyError(f"unsupported single-node topology {topology_id!r}")

    definition = load_topology(topology_id)
    if definition.dimensionality != "2D":
        raise ValueError(f"single-node topology {topology_id!r} must be 2D")
    if len(definition.node_definitions) != 1:
        raise ValueError(f"single-node topology {topology_id!r} must have exactly one node definition")

    node = definition.node_definitions[0]
    connectivity = int(node.connectivity)
    metric_family = _metric_family(tuple(definition.metadata.get("cell_parameters", ())))
    inference_mode = str(support["inference_mode"])
    expanded = expand_single_node_topology(topology_id)
    if inference_mode == "expanded":
        directions = _layout_directions_from_expanded_topology(expanded, connectivity)
    elif inference_mode == "explicit":
        explicit_directions = _explicit_edge_directions(definition)
        if len(explicit_directions) != connectivity:
            raise ValueError("explicit single-node topology layout does not provide the full coordination shell")
        directions = _sort_directions(explicit_directions)
    else:
        raise ValueError(f"unsupported single-node inference mode {inference_mode!r}")
    placement_model = "p1-two-node" if len(expanded.node_sites) == 2 else "p1-expanded"

    return SingleNodeTopologyLayout(
        topology_id=topology_id,
        connectivity=connectivity,
        directions=directions,
        metric_family=metric_family,
        inference_mode=inference_mode,
        symmetry_orbit_size=len(expanded.node_sites),
        placement_model=placement_model,
        supports_current_builder=True,
        unsupported_reason=None,
        supports_node_node=expanded.is_bipartite,
        supports_node_linker=True,
    )


@lru_cache(maxsize=None)
def expand_single_node_topology(topology_id: str) -> ExpandedSingleNodeTopology:
    definition = load_topology(topology_id)
    return expand_single_node_topology_definition(definition)


def expand_single_node_topology_definition(definition: TopologyDefinition) -> ExpandedSingleNodeTopology:
    if definition.dimensionality != "2D":
        raise ValueError(f"single-node topology {definition.id!r} must be 2D")
    if len(definition.node_definitions) != 1:
        raise ValueError(f"single-node topology {definition.id!r} must have exactly one node definition")

    cell = _cell_vectors_2d(tuple(definition.metadata.get("cell_parameters", ())))
    metric_family = _metric_family(tuple(definition.metadata.get("cell_parameters", ())))
    connectivity = int(definition.node_definitions[0].connectivity)
    operations = _point_operations_2d(str(definition.metadata.get("space_group", "")))
    if not operations:
        raise ValueError(
            f"topology {definition.id!r} uses a space group that is not yet expanded into a P1 single-node cell"
        )

    asymmetric_node = _fractional_point(definition.node_definitions[0].position)
    node_positions: list[tuple[float, float, float]] = []
    for operation in operations:
        wrapped = _wrap_fractional_point(_apply_point_operation(operation, asymmetric_node))
        if any(_same_fractional_position(wrapped, existing) for existing in node_positions):
            continue
        node_positions.append(wrapped)
    node_positions = sorted(node_positions, key=lambda position: (position[0], position[1], position[2]))
    node_ids = tuple(f"n{index + 1}" for index in range(len(node_positions)))

    raw_edges: dict[
        tuple[int, int, tuple[int, int, int]],
        tuple[Vec3, tuple[float, float, float]],
    ] = {}
    for definition_edge in definition.edge_definitions:
        start_point = _fractional_point(definition_edge.start)
        end_point = _fractional_point(definition_edge.end)
        for operation in operations:
            start_fractional = _apply_point_operation(operation, start_point)
            end_fractional = _apply_point_operation(operation, end_point)
            start_wrapped = _wrap_fractional_point(start_fractional)
            end_wrapped = _wrap_fractional_point(end_fractional)
            start_index = _match_fractional_position(start_wrapped, tuple(node_positions))
            end_index = _match_fractional_position(end_wrapped, tuple(node_positions))
            if start_index is None or end_index is None:
                continue

            start_image = _fractional_image(start_fractional, start_wrapped)
            end_image = _fractional_image(end_fractional, end_wrapped)
            relative_image = (
                end_image[0] - start_image[0],
                end_image[1] - start_image[1],
                0,
            )
            base_vector = _edge_vector_from_indices(cell, tuple(node_positions), start_index, end_index, relative_image)
            canonical = _canonical_periodic_edge(start_index, end_index, relative_image, base_vector)
            key, canonical_vector = canonical
            center_fractional = _edge_center_fractional(
                cell,
                tuple(node_positions),
                key[0],
                key[1],
                key[2],
                canonical_vector,
            )
            raw_edges[key] = (canonical_vector, center_fractional)

    edge_keys = sorted(raw_edges)
    edge_ids = tuple(f"e{index + 1}" for index in range(len(edge_keys)))
    edge_sites = tuple(
        ExpandedSingleNodeEdge(
            id=edge_id,
            start_node_id=node_ids[key[0]],
            end_node_id=node_ids[key[1]],
            end_image=key[2],
            base_vector=raw_edges[key][0],
            center_fractional=raw_edges[key][1],
        )
        for edge_id, key in zip(edge_ids, edge_keys)
    )
    edge_by_id = {edge.id: edge for edge in edge_sites}

    is_bipartite, sublattice_by_id = _bipartite_coloring(tuple(edge_sites), node_ids)
    node_sites = []
    for node_id, fractional_position in zip(node_ids, node_positions):
        incident: list[tuple[float, Vec3, str]] = []
        for edge in edge_sites:
            if edge.start_node_id == node_id:
                vector = edge.base_vector
                incident.append((atan2(vector[1], vector[0]), normalize(vector), _edge_endpoint_id(edge.id, "start")))
            if edge.end_node_id == node_id:
                vector = scale(edge.base_vector, -1.0)
                incident.append((atan2(vector[1], vector[0]), normalize(vector), _edge_endpoint_id(edge.id, "end")))
        incident.sort(key=lambda item: item[0])
        node_sites.append(
            ExpandedSingleNodeNodeSite(
                id=node_id,
                fractional_position=fractional_position,
                cartesian_position=_fractional_to_cartesian(fractional_position, cell),
                directions=tuple(direction for _, direction, _ in incident),
                edge_ids=tuple(edge_id for _, _, edge_id in incident),
                sublattice=sublattice_by_id.get(node_id),
            )
        )

    return ExpandedSingleNodeTopology(
        topology_id=definition.id,
        connectivity=connectivity,
        cell=cell,
        metric_family=metric_family,
        node_sites=tuple(node_sites),
        edge_sites=tuple(edge_sites),
        is_bipartite=is_bipartite,
    )


def _explicit_edge_directions(definition) -> tuple[Vec3, ...]:
    cell = _cell_vectors_2d(tuple(definition.metadata.get("cell_parameters", ())))
    node_position = _fractional_point(definition.node_definitions[0].position)
    directions: list[Vec3] = []
    for edge in definition.edge_definitions:
        start = _fractional_point(edge.start)
        end = _fractional_point(edge.end)
        if _same_fractional_position(start, node_position):
            other = end
        elif _same_fractional_position(end, node_position):
            other = start
        else:
            continue
        delta = _fractional_to_cartesian(
            (
                other[0] - node_position[0],
                other[1] - node_position[1],
                other[2] - node_position[2],
            ),
            cell,
        )
        if norm(delta) < 1e-8:
            continue
        directions.append(normalize(delta))
    return _deduplicate_directions(tuple(directions))


def _metric_family(cell_parameters: tuple[float, ...]) -> str:
    # Delegates to the shared geometry classifier via the parameter adapter
    # (audit A8).  Definition cells use the γ = 120° convention; the setting
    # string is dropped here because callers only consume the family.
    kind, _setting = classify_2d_cell_parameters(cell_parameters)
    return kind


def _sort_directions(directions: tuple[Vec3, ...]) -> tuple[Vec3, ...]:
    deduplicated = _deduplicate_directions(directions)
    return tuple(sorted(deduplicated, key=lambda direction: atan2(direction[1], direction[0])))


def _deduplicate_directions(directions: tuple[Vec3, ...]) -> tuple[Vec3, ...]:
    unique: list[Vec3] = []
    for direction in directions:
        unit = normalize(direction)
        if any(_direction_distance(unit, existing) < _ANGLE_TOLERANCE for existing in unique):
            continue
        unique.append(unit)
    return tuple(unique)


def _direction_distance(first: Vec3, second: Vec3) -> float:
    return abs(_wrap_angle(atan2(first[1], first[0]) - atan2(second[1], second[0])))


def _wrap_angle(value: float) -> float:
    while value <= -pi:
        value += 2.0 * pi
    while value > pi:
        value -= 2.0 * pi
    return value


def _fractional_point(values: tuple[float, ...]) -> tuple[float, float, float]:
    if len(values) == 3:
        return values
    if len(values) == 2:
        return (values[0], values[1], 0.0)
    raise ValueError(f"unsupported fractional point {values!r}")


def _same_fractional_position(
    first: tuple[float, float, float],
    second: tuple[float, float, float],
) -> bool:
    return all(abs(_canonical_fractional(a) - _canonical_fractional(b)) <= 1e-4 for a, b in zip(first, second))


def _canonical_fractional(value: float) -> float:
    wrapped = value % 1.0
    if abs(wrapped - 1.0) <= 1e-6:
        return 0.0
    return wrapped


def _wrap_fractional_point(values: tuple[float, float, float]) -> tuple[float, float, float]:
    return (
        _canonical_fractional(values[0]),
        _canonical_fractional(values[1]),
        _canonical_fractional(values[2]),
    )


def _cell_vectors_2d(cell_parameters: tuple[float, ...]) -> tuple[Vec3, Vec3, Vec3]:
    if len(cell_parameters) >= 6:
        a, b, c, _alpha, _beta, gamma = cell_parameters[:6]
    elif len(cell_parameters) == 3:
        a, b, gamma = cell_parameters
        c = 1.0
    else:
        raise ValueError("single-node topology definition is missing cell parameters")
    gamma_radians = gamma * pi / 180.0
    return (
        (a, 0.0, 0.0),
        (b * cos(gamma_radians), b * sin(gamma_radians), 0.0),
        (0.0, 0.0, c),
    )


def _fractional_to_cartesian(
    fractional: tuple[float, float, float],
    cell: tuple[Vec3, Vec3, Vec3],
) -> Vec3:
    return (
        fractional[0] * cell[0][0] + fractional[1] * cell[1][0] + fractional[2] * cell[2][0],
        fractional[0] * cell[0][1] + fractional[1] * cell[1][1] + fractional[2] * cell[2][1],
        fractional[0] * cell[0][2] + fractional[1] * cell[1][2] + fractional[2] * cell[2][2],
    )


# A 2D point operation in fractional coordinates: a 2x2 in-plane matrix plus
# an in-plane fractional translation.
FractionalPointOperation2D = tuple[
    tuple[tuple[float, float], tuple[float, float]],
    tuple[float, float],
]


def _point_operations_2d(space_group: str) -> tuple[FractionalPointOperation2D, ...]:
    # General 2D plane-group support via gemmi: take the 3D space-group
    # operations and keep those that map the layer plane onto itself
    # (z -> +z or z -> -z with no fractional z translation), projected to
    # their in-plane fractional matrix and translation. This covers every
    # plane group gemmi can parse (P6/mmm, P4/mmm, Cmmm, P4/mbm, Pbam, ...)
    # instead of a hand-written per-group table.
    if not space_group:
        return ()
    group = gemmi.find_spacegroup_by_name(space_group)
    if group is None:
        return ()
    denominator = gemmi.Op.DEN
    operations: list[FractionalPointOperation2D] = []
    for operation in group.operations():
        rot = operation.rot
        if rot[0][2] != 0 or rot[1][2] != 0 or rot[2][0] != 0 or rot[2][1] != 0:
            continue
        if abs(rot[2][2]) != denominator:
            continue
        z_shift = operation.tran[2] / denominator
        if abs(z_shift - round(z_shift)) > 1e-9:
            continue
        matrix = (
            (rot[0][0] / denominator, rot[0][1] / denominator),
            (rot[1][0] / denominator, rot[1][1] / denominator),
        )
        shift = (
            _canonical_fractional(operation.tran[0] / denominator),
            _canonical_fractional(operation.tran[1] / denominator),
        )
        operations.append((matrix, shift))
    return _deduplicate_point_operations(tuple(operations))


def _deduplicate_point_operations(
    operations: tuple[FractionalPointOperation2D, ...],
) -> tuple[FractionalPointOperation2D, ...]:
    unique: list[FractionalPointOperation2D] = []
    for matrix, shift in operations:
        if any(
            all(abs(matrix[row][col] - existing_matrix[row][col]) < 1e-8 for row in range(2) for col in range(2))
            and all(abs(shift[axis] - existing_shift[axis]) < 1e-8 for axis in range(2))
            for existing_matrix, existing_shift in unique
        ):
            continue
        unique.append((matrix, shift))
    return tuple(unique)


def _apply_point_operation(
    operation: FractionalPointOperation2D,
    point: tuple[float, float, float],
) -> tuple[float, float, float]:
    matrix, shift = operation
    return (
        matrix[0][0] * point[0] + matrix[0][1] * point[1] + shift[0],
        matrix[1][0] * point[0] + matrix[1][1] * point[1] + shift[1],
        point[2],
    )


def _cartesian_to_fractional_2d(
    point: Vec3,
    cell: tuple[Vec3, Vec3, Vec3],
) -> tuple[float, float, float]:
    first, second, _ = cell
    determinant = first[0] * second[1] - first[1] * second[0]
    if abs(determinant) < 1e-8:
        raise ValueError("single-node topology cell is singular")
    return (
        (point[0] * second[1] - point[1] * second[0]) / determinant,
        (first[0] * point[1] - first[1] * point[0]) / determinant,
        0.0,
    )


def _match_fractional_position(
    position: tuple[float, float, float],
    positions: tuple[tuple[float, float, float], ...],
) -> int | None:
    for index, existing in enumerate(positions):
        if _same_fractional_position(position, existing):
            return index
    return None


def _fractional_image(
    unwrapped: tuple[float, float, float],
    wrapped: tuple[float, float, float],
) -> tuple[int, int, int]:
    return (
        int(round(unwrapped[0] - wrapped[0])),
        int(round(unwrapped[1] - wrapped[1])),
        int(round(unwrapped[2] - wrapped[2])),
    )


def _edge_vector_from_indices(
    cell: tuple[Vec3, Vec3, Vec3],
    node_positions: tuple[tuple[float, float, float], ...],
    start_index: int,
    end_index: int,
    end_image: tuple[int, int, int],
) -> Vec3:
    start = _fractional_to_cartesian(node_positions[start_index], cell)
    end = add(
        _fractional_to_cartesian(node_positions[end_index], cell),
        _periodic_offset(cell, end_image),
    )
    return sub(end, start)


def _canonical_periodic_edge(
    start_index: int,
    end_index: int,
    end_image: tuple[int, int, int],
    base_vector: Vec3,
) -> tuple[tuple[int, int, tuple[int, int, int]], Vec3]:
    forward = (start_index, end_index, end_image)
    reverse = (end_index, start_index, (-end_image[0], -end_image[1], -end_image[2]))
    if reverse < forward:
        return reverse, scale(base_vector, -1.0)
    return forward, base_vector


def _edge_center_fractional(
    cell: tuple[Vec3, Vec3, Vec3],
    node_positions: tuple[tuple[float, float, float], ...],
    start_index: int,
    end_index: int,
    end_image: tuple[int, int, int],
    base_vector: Vec3,
) -> tuple[float, float, float]:
    start = _fractional_to_cartesian(node_positions[start_index], cell)
    midpoint = add(start, scale(base_vector, 0.5))
    return _wrap_fractional_point(_cartesian_to_fractional_2d(midpoint, cell))


def _periodic_offset(
    cell: tuple[Vec3, Vec3, Vec3],
    image: tuple[int, int, int],
) -> Vec3:
    return add(add(scale(cell[0], image[0]), scale(cell[1], image[1])), scale(cell[2], image[2]))


def _bipartite_coloring(
    edge_sites: tuple[ExpandedSingleNodeEdge, ...],
    node_ids: tuple[str, ...],
) -> tuple[bool, dict[str, int]]:
    # Exact quotient-graph test: each edge carries its integer image shift,
    # and a periodic 2-coloring has the form
    #     color(node, image) = sublattice(node) XOR q . image   (mod 2)
    # for some parity wavevector q in (Z2)^3. Such a coloring exists iff the
    # parity-expanded graph — states (node, image mod 2), with an edge from
    # (start, p) to (end, p + end_image mod 2) — is bipartite, which this BFS
    # checks directly. Cycles of any length are handled exactly, including
    # self-loops with nonzero image (an odd-parity self-loop constrains two
    # distinct parity states and stays bipartite; an even-parity self-loop is
    # an immediate conflict). The test can report non-bipartite for a
    # disconnected lifted net whose components are individually bipartite but
    # admit no coloring of the realizable form above; that conservative
    # direction is safe for sublattice role assignment.
    if not edge_sites:
        return True, {node_id: 0 for node_id in node_ids}

    colors: dict[tuple[str, int, int, int], int] = {}
    for node_id in node_ids:
        seed = (node_id, 0, 0, 0)
        if seed in colors:
            continue
        colors[seed] = 0
        queue = deque((seed,))
        while queue:
            current_id, parity_x, parity_y, parity_z = queue.popleft()
            current_color = colors[(current_id, parity_x, parity_y, parity_z)]
            for edge in edge_sites:
                neighbors: list[tuple[str, int, int, int]] = []
                if edge.start_node_id == current_id:
                    neighbors.append(
                        (
                            edge.end_node_id,
                            (parity_x + edge.end_image[0]) % 2,
                            (parity_y + edge.end_image[1]) % 2,
                            (parity_z + edge.end_image[2]) % 2,
                        )
                    )
                if edge.end_node_id == current_id:
                    neighbors.append(
                        (
                            edge.start_node_id,
                            (parity_x - edge.end_image[0]) % 2,
                            (parity_y - edge.end_image[1]) % 2,
                            (parity_z - edge.end_image[2]) % 2,
                        )
                    )
                for other_key in neighbors:
                    if other_key not in colors:
                        colors[other_key] = 1 - current_color
                        queue.append(other_key)
                        continue
                    if colors[other_key] == current_color:
                        return False, {}

    return True, {node_id: colors[(node_id, 0, 0, 0)] for node_id in node_ids}


def _layout_directions_from_expanded_topology(
    expanded: ExpandedSingleNodeTopology,
    connectivity: int,
) -> tuple[Vec3, ...]:
    for node_site in expanded.node_sites:
        if len(node_site.directions) == connectivity:
            return _sort_directions(node_site.directions)
    raise ValueError(
        f"expanded single-node topology {expanded.topology_id!r} does not provide a full {connectivity}-connected star"
    )


def _edge_endpoint_id(edge_id: str, endpoint: str) -> str:
    return f"{edge_id}:{endpoint}"


__all__ = [
    "ExpandedSingleNodeEdge",
    "ExpandedSingleNodeNodeSite",
    "ExpandedSingleNodeTopology",
    "SingleNodeTopologyLayout",
    "expand_single_node_topology",
    "expand_single_node_topology_definition",
    "list_supported_single_node_topology_ids",
    "resolve_single_node_topology_layout",
]
