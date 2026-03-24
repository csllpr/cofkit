from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import atan2, cos, pi, sin, sqrt

from .geometry import Vec3, add, norm, normalize, scale, sub
from .topologies import load_topology
from .topology_index import TopologyDefinition
from .topology_symmetry import expand_topology_definition


@dataclass(frozen=True)
class ExpandedIndexedNodeSite:
    id: str
    fractional_position: tuple[float, float, float]
    cartesian_position: Vec3
    directions: tuple[Vec3, ...]
    edge_ids: tuple[str, ...]
    connectivity: int
    source_node_index: int
    source_label: str
    sublattice: int | None = None


@dataclass(frozen=True)
class ExpandedIndexedEdge:
    id: str
    start_node_id: str
    end_node_id: str
    end_image: tuple[int, int, int]
    base_vector: Vec3
    center_fractional: tuple[float, float, float]


@dataclass(frozen=True)
class ExpandedIndexedTopology:
    topology_id: str
    dimensionality: str
    cell: tuple[Vec3, Vec3, Vec3]
    metric_family: str
    node_sites: tuple[ExpandedIndexedNodeSite, ...]
    edge_sites: tuple[ExpandedIndexedEdge, ...]
    is_bipartite: bool


@lru_cache(maxsize=None)
def expand_indexed_topology(topology_id: str) -> ExpandedIndexedTopology:
    return expand_indexed_topology_definition(load_topology(topology_id))


def expand_indexed_topology_definition(definition: TopologyDefinition) -> ExpandedIndexedTopology:
    expanded = expand_topology_definition(definition)
    cell = _cell_vectors(definition)
    metric_family = _metric_family(definition)
    node_positions = {
        node_site.id: _normalize_fractional_position(node_site.fractional_position)
        for node_site in expanded.node_sites
    }
    edge_sites = tuple(
        ExpandedIndexedEdge(
            id=edge.id,
            start_node_id=edge.start_node_id,
            end_node_id=edge.end_node_id,
            end_image=_normalize_image(edge.end_image),
            base_vector=_periodic_edge_vector(
                cell,
                node_positions[edge.start_node_id],
                node_positions[edge.end_node_id],
                _normalize_image(edge.end_image),
            ),
            center_fractional=_edge_center_fractional(
                node_positions[edge.start_node_id],
                node_positions[edge.end_node_id],
                _normalize_image(edge.end_image),
            ),
        )
        for edge in expanded.edge_sites
    )

    incident_by_node_id: dict[str, list[tuple[object, Vec3, str]]] = {
        node_site.id: [] for node_site in expanded.node_sites
    }
    for edge in edge_sites:
        forward = _safe_normalize(edge.base_vector)
        reverse = _safe_normalize(scale(edge.base_vector, -1.0))
        incident_by_node_id[edge.start_node_id].append(
            (_incident_sort_key(forward, definition.dimensionality), forward, _edge_endpoint_key(edge.id, "start"))
        )
        incident_by_node_id[edge.end_node_id].append(
            (_incident_sort_key(reverse, definition.dimensionality), reverse, _edge_endpoint_key(edge.id, "end"))
        )

    sublattice_by_source = _source_orbit_sublattice_map(expanded, edge_sites)
    node_sites = tuple(
        ExpandedIndexedNodeSite(
            id=node_site.id,
            fractional_position=node_positions[node_site.id],
            cartesian_position=_fractional_to_cartesian(node_positions[node_site.id], cell),
            directions=tuple(
                direction
                for _, direction, _ in sorted(incident_by_node_id[node_site.id], key=lambda item: item[0])
            ),
            edge_ids=tuple(
                edge_id
                for _, _, edge_id in sorted(incident_by_node_id[node_site.id], key=lambda item: item[0])
            ),
            connectivity=int(node_site.connectivity),
            source_node_index=int(node_site.source_node_index),
            source_label=str(node_site.source_label),
            sublattice=sublattice_by_source.get(int(node_site.source_node_index)),
        )
        for node_site in expanded.node_sites
    )
    return ExpandedIndexedTopology(
        topology_id=definition.id,
        dimensionality=definition.dimensionality,
        cell=cell,
        metric_family=metric_family,
        node_sites=node_sites,
        edge_sites=edge_sites,
        is_bipartite=bool(definition.metadata.get("zero_linker_bipartite")),
    )


def _cell_vectors(definition: TopologyDefinition) -> tuple[Vec3, Vec3, Vec3]:
    cell_parameters = tuple(float(value) for value in definition.metadata.get("cell_parameters", ()))
    if definition.dimensionality == "2D":
        if len(cell_parameters) >= 6:
            a, b, c, _alpha, _beta, gamma = cell_parameters[:6]
        elif len(cell_parameters) == 3:
            a, b, gamma = cell_parameters
            c = 1.0
        else:
            raise ValueError(f"topology {definition.id!r} is missing 2D cell parameters")
        gamma_radians = gamma * pi / 180.0
        return (
            (a, 0.0, 0.0),
            (b * cos(gamma_radians), b * sin(gamma_radians), 0.0),
            (0.0, 0.0, c),
        )

    if len(cell_parameters) >= 6:
        a, b, c, alpha, beta, gamma = cell_parameters[:6]
    elif len(cell_parameters) == 3:
        a, b, c = cell_parameters
        alpha = beta = gamma = 90.0
    else:
        raise ValueError(f"topology {definition.id!r} is missing 3D cell parameters")

    alpha_radians = alpha * pi / 180.0
    beta_radians = beta * pi / 180.0
    gamma_radians = gamma * pi / 180.0
    cos_alpha = cos(alpha_radians)
    cos_beta = cos(beta_radians)
    cos_gamma = cos(gamma_radians)
    sin_gamma = sin(gamma_radians)
    if abs(sin_gamma) < 1e-8:
        raise ValueError(f"topology {definition.id!r} has a degenerate gamma angle")

    cell_a = (a, 0.0, 0.0)
    cell_b = (b * cos_gamma, b * sin_gamma, 0.0)
    c_x = c * cos_beta
    c_y = c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma
    c_z_squared = c * c - c_x * c_x - c_y * c_y
    c_z = sqrt(max(c_z_squared, 0.0))
    cell_c = (c_x, c_y, c_z)
    return (cell_a, cell_b, cell_c)


def _metric_family(definition: TopologyDefinition) -> str:
    cell_parameters = tuple(float(value) for value in definition.metadata.get("cell_parameters", ()))
    if definition.dimensionality == "2D":
        if len(cell_parameters) >= 6:
            a, b, _c, _alpha, _beta, gamma = cell_parameters[:6]
        elif len(cell_parameters) == 3:
            a, b, gamma = cell_parameters
        else:
            return "unknown"
        if abs(a - b) < 1e-3 and abs(gamma - 120.0) < 1e-2:
            return "hexagonal"
        if abs(a - b) < 1e-3 and abs(gamma - 90.0) < 1e-2:
            return "square"
        if abs(gamma - 90.0) < 1e-2:
            return "orthogonal"
        return "oblique"

    if len(cell_parameters) >= 6:
        a, b, c, alpha, beta, gamma = cell_parameters[:6]
    elif len(cell_parameters) == 3:
        a, b, c = cell_parameters
        alpha = beta = gamma = 90.0
    else:
        return "unknown"
    if all(abs(angle - 90.0) < 1e-2 for angle in (alpha, beta, gamma)):
        if max(a, b, c) - min(a, b, c) < 1e-3:
            return "cubic"
        return "orthorhombic"
    return "triclinic"


def _normalize_fractional_position(position: tuple[float, ...]) -> tuple[float, float, float]:
    if len(position) >= 3:
        return (float(position[0]), float(position[1]), float(position[2]))
    if len(position) == 2:
        return (float(position[0]), float(position[1]), 0.0)
    raise ValueError(f"expected a 2D or 3D fractional position, got {position!r}")


def _normalize_image(image: tuple[int, ...]) -> tuple[int, int, int]:
    if len(image) >= 3:
        return (int(image[0]), int(image[1]), int(image[2]))
    if len(image) == 2:
        return (int(image[0]), int(image[1]), 0)
    raise ValueError(f"expected a 2D or 3D periodic image, got {image!r}")


def _fractional_to_cartesian(
    fractional: tuple[float, float, float],
    cell: tuple[Vec3, Vec3, Vec3],
) -> Vec3:
    return add(
        add(scale(cell[0], fractional[0]), scale(cell[1], fractional[1])),
        scale(cell[2], fractional[2]),
    )


def _periodic_edge_vector(
    cell: tuple[Vec3, Vec3, Vec3],
    start_fractional: tuple[float, float, float],
    end_fractional: tuple[float, float, float],
    end_image: tuple[int, int, int],
) -> Vec3:
    start = _fractional_to_cartesian(start_fractional, cell)
    end = add(
        _fractional_to_cartesian(end_fractional, cell),
        add(add(scale(cell[0], end_image[0]), scale(cell[1], end_image[1])), scale(cell[2], end_image[2])),
    )
    return sub(end, start)


def _edge_center_fractional(
    start_fractional: tuple[float, float, float],
    end_fractional: tuple[float, float, float],
    end_image: tuple[int, int, int],
) -> tuple[float, float, float]:
    return (
        (start_fractional[0] + end_fractional[0] + end_image[0]) * 0.5 % 1.0,
        (start_fractional[1] + end_fractional[1] + end_image[1]) * 0.5 % 1.0,
        (start_fractional[2] + end_fractional[2] + end_image[2]) * 0.5 % 1.0,
    )


def _source_orbit_sublattice_map(
    expanded,
    edge_sites: tuple[ExpandedIndexedEdge, ...],
) -> dict[int, int]:
    node_by_id = {node_site.id: node_site for node_site in expanded.node_sites}
    source_indices = tuple(sorted({int(node_site.source_node_index) for node_site in expanded.node_sites}))
    if len(source_indices) != 2:
        return {}
    for edge in edge_sites:
        start_source = int(node_by_id[edge.start_node_id].source_node_index)
        end_source = int(node_by_id[edge.end_node_id].source_node_index)
        if start_source == end_source:
            return {}
    return {source_indices[0]: 0, source_indices[1]: 1}


def _incident_sort_key(direction: Vec3, dimensionality: str) -> object:
    if dimensionality == "2D":
        return atan2(direction[1], direction[0])
    return tuple(round(value, 6) for value in direction)


def _safe_normalize(vector: Vec3) -> Vec3:
    if norm(vector) < 1e-8:
        return (0.0, 0.0, 1.0)
    return normalize(vector)


def _edge_endpoint_key(edge_id: str, endpoint: str) -> str:
    return f"{edge_id}:{endpoint}"


__all__ = [
    "ExpandedIndexedEdge",
    "ExpandedIndexedNodeSite",
    "ExpandedIndexedTopology",
    "expand_indexed_topology",
    "expand_indexed_topology_definition",
]
