from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

from .geometry import Vec3, add, norm, normalize, scale
from .single_node_topologies import (
    ExpandedSingleNodeEdge,
    ExpandedSingleNodeNodeSite,
    ExpandedSingleNodeTopology,
    SingleNodeTopologyLayout,
)
from .topology_index import TopologyDefinition
from .topologies import load_topology


@dataclass(frozen=True)
class _SupportedThreeDSingleNodeTopology:
    connectivity: int
    metric_family: str
    placement_model: str
    supports_node_node: bool
    supports_node_linker: bool


_SUPPORTED_3D_SINGLE_NODE_TOPOLOGIES: dict[str, _SupportedThreeDSingleNodeTopology] = {
    "dia": _SupportedThreeDSingleNodeTopology(
        connectivity=4,
        metric_family="cubic",
        placement_model="p1-two-node-3d",
        supports_node_node=True,
        supports_node_linker=True,
    ),
    "pcu": _SupportedThreeDSingleNodeTopology(
        connectivity=6,
        metric_family="cubic",
        placement_model="p1-self-edge-3d",
        supports_node_node=False,
        supports_node_linker=True,
    ),
}


def list_supported_3d_single_node_topology_ids(
    connectivity: int,
    *,
    pair_mode: str,
) -> tuple[str, ...]:
    ids: list[str] = []
    for topology_id in _SUPPORTED_3D_SINGLE_NODE_TOPOLOGIES:
        try:
            layout = resolve_three_d_single_node_topology_layout(topology_id)
        except (KeyError, ValueError):
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
def resolve_three_d_single_node_topology_layout(topology_id: str) -> SingleNodeTopologyLayout:
    support = _SUPPORTED_3D_SINGLE_NODE_TOPOLOGIES.get(topology_id)
    if support is None:
        raise KeyError(f"unsupported 3D single-node topology {topology_id!r}")

    definition = load_topology(topology_id)
    if definition.dimensionality != "3D":
        raise ValueError(f"single-node 3D topology {topology_id!r} must be 3D")
    if len(definition.node_definitions) != 1:
        raise ValueError(f"single-node 3D topology {topology_id!r} must have exactly one node definition")

    expanded = expand_three_d_single_node_topology(topology_id)
    directions = _layout_directions_from_expanded_topology(expanded, support.connectivity)
    return SingleNodeTopologyLayout(
        topology_id=topology_id,
        connectivity=support.connectivity,
        directions=directions,
        metric_family=support.metric_family,
        inference_mode="explicit-3d",
        symmetry_orbit_size=len(expanded.node_sites),
        placement_model=support.placement_model,
        supports_current_builder=True,
        unsupported_reason=None,
        supports_node_node=support.supports_node_node,
        supports_node_linker=support.supports_node_linker,
    )


@lru_cache(maxsize=None)
def expand_three_d_single_node_topology(topology_id: str) -> ExpandedSingleNodeTopology:
    definition = load_topology(topology_id)
    return expand_three_d_single_node_topology_definition(definition)


def expand_three_d_single_node_topology_definition(definition: TopologyDefinition) -> ExpandedSingleNodeTopology:
    if definition.dimensionality != "3D":
        raise ValueError(f"single-node 3D topology {definition.id!r} must be 3D")
    if len(definition.node_definitions) != 1:
        raise ValueError(f"single-node 3D topology {definition.id!r} must have exactly one node definition")

    if definition.id == "dia":
        return _build_dia_topology(definition)
    if definition.id == "pcu":
        return _build_pcu_topology(definition)
    raise KeyError(f"unsupported 3D single-node topology {definition.id!r}")


def _build_dia_topology(definition) -> ExpandedSingleNodeTopology:
    conventional_a = float(definition.metadata.get("cell_parameters", (1.0,))[0])
    cell = (
        (0.0, 0.5 * conventional_a, 0.5 * conventional_a),
        (0.5 * conventional_a, 0.0, 0.5 * conventional_a),
        (0.5 * conventional_a, 0.5 * conventional_a, 0.0),
    )
    node_positions = {
        "n1": (0.0, 0.0, 0.0),
        "n2": (0.25, 0.25, 0.25),
    }
    raw_edges = (
        ("e1", "n1", "n2", (0, 0, 0)),
        ("e2", "n1", "n2", (-1, 0, 0)),
        ("e3", "n1", "n2", (0, -1, 0)),
        ("e4", "n1", "n2", (0, 0, -1)),
    )
    edge_sites = tuple(
        ExpandedSingleNodeEdge(
            id=edge_id,
            start_node_id=start_node_id,
            end_node_id=end_node_id,
            end_image=end_image,
            base_vector=_periodic_edge_vector(cell, node_positions[start_node_id], node_positions[end_node_id], end_image),
            center_fractional=_edge_center_fractional(node_positions[start_node_id], node_positions[end_node_id], end_image),
        )
        for edge_id, start_node_id, end_node_id, end_image in raw_edges
    )
    node_sites = (
        _build_node_site(
            node_id="n1",
            fractional_position=node_positions["n1"],
            cell=cell,
            edge_sites=edge_sites,
            sublattice=0,
        ),
        _build_node_site(
            node_id="n2",
            fractional_position=node_positions["n2"],
            cell=cell,
            edge_sites=edge_sites,
            sublattice=1,
        ),
    )
    return ExpandedSingleNodeTopology(
        topology_id=definition.id,
        connectivity=4,
        cell=cell,
        metric_family="cubic",
        node_sites=node_sites,
        edge_sites=edge_sites,
        is_bipartite=True,
    )


def _build_pcu_topology(definition) -> ExpandedSingleNodeTopology:
    a, b, c = tuple(float(value) for value in definition.metadata.get("cell_parameters", (1.0, 1.0, 1.0))[:3])
    cell = (
        (a, 0.0, 0.0),
        (0.0, b, 0.0),
        (0.0, 0.0, c),
    )
    node_position = (0.0, 0.0, 0.0)
    edge_sites = (
        ExpandedSingleNodeEdge(
            id="e1",
            start_node_id="n1",
            end_node_id="n1",
            end_image=(1, 0, 0),
            base_vector=(a, 0.0, 0.0),
            center_fractional=(0.5, 0.0, 0.0),
        ),
        ExpandedSingleNodeEdge(
            id="e2",
            start_node_id="n1",
            end_node_id="n1",
            end_image=(0, 1, 0),
            base_vector=(0.0, b, 0.0),
            center_fractional=(0.0, 0.5, 0.0),
        ),
        ExpandedSingleNodeEdge(
            id="e3",
            start_node_id="n1",
            end_node_id="n1",
            end_image=(0, 0, 1),
            base_vector=(0.0, 0.0, c),
            center_fractional=(0.0, 0.0, 0.5),
        ),
    )
    node_sites = (
        _build_node_site(
            node_id="n1",
            fractional_position=node_position,
            cell=cell,
            edge_sites=edge_sites,
            sublattice=0,
        ),
    )
    return ExpandedSingleNodeTopology(
        topology_id=definition.id,
        connectivity=6,
        cell=cell,
        metric_family="cubic",
        node_sites=node_sites,
        edge_sites=edge_sites,
        is_bipartite=True,
    )


def _build_node_site(
    *,
    node_id: str,
    fractional_position: tuple[float, float, float],
    cell: tuple[Vec3, Vec3, Vec3],
    edge_sites: tuple[ExpandedSingleNodeEdge, ...],
    sublattice: int,
) -> ExpandedSingleNodeNodeSite:
    incident: list[tuple[str, Vec3]] = []
    for edge in edge_sites:
        if edge.start_node_id == node_id:
            incident.append((_edge_endpoint_key(edge.id, "start"), normalize(edge.base_vector)))
        if edge.end_node_id == node_id:
            incident.append((_edge_endpoint_key(edge.id, "end"), normalize(scale(edge.base_vector, -1.0))))
    directions = tuple(direction for _, direction in incident)
    edge_ids = tuple(edge_id for edge_id, _ in incident)
    return ExpandedSingleNodeNodeSite(
        id=node_id,
        fractional_position=fractional_position,
        cartesian_position=_fractional_to_cartesian(fractional_position, cell),
        directions=directions,
        edge_ids=edge_ids,
        sublattice=sublattice,
    )


def _layout_directions_from_expanded_topology(
    expanded: ExpandedSingleNodeTopology,
    connectivity: int,
) -> tuple[Vec3, ...]:
    for node_site in expanded.node_sites:
        if len(node_site.directions) == connectivity:
            return node_site.directions
    raise ValueError(
        f"expanded 3D single-node topology {expanded.topology_id!r} does not provide a full {connectivity}-connected star"
    )


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
    return (end[0] - start[0], end[1] - start[1], end[2] - start[2])


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


def _edge_endpoint_key(edge_id: str, endpoint: str) -> str:
    return f"{edge_id}:{endpoint}"


__all__ = [
    "expand_three_d_single_node_topology",
    "expand_three_d_single_node_topology_definition",
    "list_supported_3d_single_node_topology_ids",
    "resolve_three_d_single_node_topology_layout",
]
