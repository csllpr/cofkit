"""Frozen linear symmetry-expansion oracle before topology preparation caching."""
from cofkit.topology_symmetry import (
    ExpandedPeriodicNodeSite, ExpandedPeriodicEdge, ExpandedPeriodicTopology,
    _gemmi, _normalize_position, _apply_operation, _project_position,
    _wrap_fractional_position, _same_fractional_position, _match_position,
    _fractional_image, _canonical_periodic_edge,
)
from cofkit.topology_index import TopologyDefinition


def exhaustive_topology_expansion(definition: TopologyDefinition) -> ExpandedPeriodicTopology:
    if definition.dimensionality not in {"2D", "3D"}:
        raise ValueError(f"topology {definition.id!r} must be 2D or 3D")
    if not definition.node_definitions:
        raise ValueError(f"topology {definition.id!r} does not define any nodes")
    if not definition.edge_definitions:
        raise ValueError(f"topology {definition.id!r} does not define any edges")
    dimension_count = 2 if definition.dimensionality == "2D" else 3

    gemmi = _gemmi()
    space_group_name = str(definition.metadata.get("space_group", "")).strip()
    if not space_group_name:
        raise ValueError(f"topology {definition.id!r} does not declare a space group")
    space_group = gemmi.find_spacegroup_by_name(space_group_name)
    if space_group is None:
        raise ValueError(f"space group {space_group_name!r} is not recognized by gemmi")

    operations = tuple(space_group.operations())
    if not operations:
        raise ValueError(f"space group {space_group_name!r} does not provide symmetry operations")

    raw_node_sites: list[tuple[tuple[float, ...], int, str, int]] = []
    for source_index, node in enumerate(definition.node_definitions):
        position = _normalize_position(node.position)
        for operation in operations:
            transformed = _apply_operation(operation, position)
            wrapped = _wrap_fractional_position(_project_position(transformed, dimension_count))
            if any(_same_fractional_position(wrapped, existing[0]) for existing in raw_node_sites):
                continue
            raw_node_sites.append((wrapped, source_index, node.label, int(node.connectivity)))

    raw_node_sites.sort(key=lambda item: (item[0], item[1], item[2]))
    node_sites = tuple(
        ExpandedPeriodicNodeSite(
            id=f"n{index + 1}",
            fractional_position=position,
            source_node_index=source_index,
            source_label=source_label,
            connectivity=connectivity,
        )
        for index, (position, source_index, source_label, connectivity) in enumerate(raw_node_sites)
    )
    raw_edge_keys: set[tuple[int, int, tuple[int, ...]]] = set()
    for edge in definition.edge_definitions:
        start = _normalize_position(edge.start)
        end = _normalize_position(edge.end)
        for operation in operations:
            transformed_start = _apply_operation(operation, start)
            transformed_end = _apply_operation(operation, end)
            start_projected = _project_position(transformed_start, dimension_count)
            end_projected = _project_position(transformed_end, dimension_count)
            start_wrapped = _wrap_fractional_position(start_projected)
            end_wrapped = _wrap_fractional_position(end_projected)
            start_index = _match_position(start_wrapped, node_sites)
            end_index = _match_position(end_wrapped, node_sites)
            if start_index is None or end_index is None:
                raise ValueError(
                    f"symmetry-expanded edge endpoints for topology {definition.id!r} do not resolve to node sites"
                )
            start_image = _fractional_image(start_projected, start_wrapped)
            end_image = _fractional_image(end_projected, end_wrapped)
            relative_image = tuple(
                end_component - start_component
                for start_component, end_component in zip(start_image, end_image)
            )
            raw_edge_keys.add(_canonical_periodic_edge(start_index, end_index, relative_image))

    edge_keys = sorted(raw_edge_keys)
    edge_sites = tuple(
        ExpandedPeriodicEdge(
            id=f"e{index + 1}",
            start_node_id=node_sites[start_index].id,
            end_node_id=node_sites[end_index].id,
            end_image=end_image,
        )
        for index, (start_index, end_index, end_image) in enumerate(edge_keys)
    )
    return ExpandedPeriodicTopology(
        topology_id=definition.id,
        dimensionality=definition.dimensionality,
        node_sites=node_sites,
        edge_sites=edge_sites,
    )

