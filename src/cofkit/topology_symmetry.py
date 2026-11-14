from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import product
import math

from .topology_index import TopologyDefinition, TopologyNodeDefinition, TopologyEdgeDefinition

# Fractional-coordinate tolerance for treating a value as sitting on a cell
# boundary (wrapped to 0/1) during periodic wrap. Heuristic — pending
# calibration. Single owner; topology_analysis imports this definition.
_FRACTIONAL_WRAP_TOLERANCE = 1e-4

# Heuristic — pending calibration: bounded working set for immutable expanded
# topologies and their derived decomposition graphs. Keys include every input
# used by expansion, so editing a definition cannot return an obsolete graph.
TOPOLOGY_PREPARATION_CACHE_SIZE = 256


@dataclass(frozen=True)
class ExpandedPeriodicNodeSite:
    id: str
    fractional_position: tuple[float, ...]
    source_node_index: int
    source_label: str
    connectivity: int


@dataclass(frozen=True)
class ExpandedPeriodicEdge:
    id: str
    start_node_id: str
    end_node_id: str
    end_image: tuple[int, ...]


@dataclass(frozen=True)
class ExpandedPeriodicTopology:
    topology_id: str
    dimensionality: str
    node_sites: tuple[ExpandedPeriodicNodeSite, ...]
    edge_sites: tuple[ExpandedPeriodicEdge, ...]


def expand_topology_definition(definition: TopologyDefinition) -> ExpandedPeriodicTopology:
    return _cached_topology_expansion(
        definition.id, definition.dimensionality,
        str(definition.metadata.get("space_group", "")).strip(),
        tuple((node.label, tuple(node.position), node.connectivity) for node in definition.node_definitions),
        tuple((tuple(edge.start), tuple(edge.end)) for edge in definition.edge_definitions),
    )


@lru_cache(maxsize=TOPOLOGY_PREPARATION_CACHE_SIZE)
def _cached_topology_expansion(topology_id, dimensionality, space_group, nodes, edges):
    # Only these fields enter the expansion algorithm. Other mutable metadata
    # (including annotations added during repository loading) cannot affect it.
    definition = TopologyDefinition(
        id=topology_id, name=topology_id, dimensionality=dimensionality,
        node_definitions=tuple(TopologyNodeDefinition(label, position, connectivity) for label, position, connectivity in nodes),
        edge_definitions=tuple(TopologyEdgeDefinition(start, end) for start, end in edges),
        metadata={"space_group": space_group},
    )
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
    raw_position_index = _PositionIndex()
    for source_index, node in enumerate(definition.node_definitions):
        position = _normalize_position(node.position)
        for operation in operations:
            transformed = _apply_operation(operation, position)
            wrapped = _wrap_fractional_position(_project_position(transformed, dimension_count))
            if raw_position_index.match(wrapped) is not None:
                continue
            raw_position_index.add(wrapped)
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
    node_position_index = _PositionIndex()
    for node_site in node_sites:
        node_position_index.add(node_site.fractional_position)
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
            start_index = node_position_index.match(start_wrapped)
            end_index = node_position_index.match(end_wrapped)
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


def expand_two_d_topology_definition(definition: TopologyDefinition) -> ExpandedPeriodicTopology:
    if definition.dimensionality != "2D":
        raise ValueError(f"topology {definition.id!r} must be 2D")
    return expand_topology_definition(definition)


def expand_three_d_topology_definition(definition: TopologyDefinition) -> ExpandedPeriodicTopology:
    if definition.dimensionality != "3D":
        raise ValueError(f"topology {definition.id!r} must be 3D")
    return expand_topology_definition(definition)


def has_gemmi() -> bool:
    try:
        _gemmi()
    except ImportError:
        return False
    return True


def _gemmi():
    try:
        import gemmi
    except ImportError as exc:
        raise ImportError("gemmi is required for generic 3D topology symmetry expansion") from exc
    return gemmi


def _apply_operation(operation, position: tuple[float, float, float]) -> tuple[float, float, float]:
    transformed = operation.apply_to_xyz(list(position))
    return (float(transformed[0]), float(transformed[1]), float(transformed[2]))


def _normalize_position(position: tuple[float, ...]) -> tuple[float, float, float]:
    if len(position) >= 3:
        return (float(position[0]), float(position[1]), float(position[2]))
    if len(position) == 2:
        return (float(position[0]), float(position[1]), 0.0)
    raise ValueError(f"expected 2D or 3D fractional coordinates, got {position!r}")


def _project_position(position: tuple[float, float, float], dimension_count: int) -> tuple[float, ...]:
    return tuple(float(position[index]) for index in range(dimension_count))


def _wrap_fractional_position(position: tuple[float, ...]) -> tuple[float, ...]:
    return tuple(_canonical_fractional(value) for value in position)


def _canonical_fractional(value: float) -> float:
    if abs(value - round(value)) <= _FRACTIONAL_WRAP_TOLERANCE:
        return 0.0
    wrapped = value % 1.0
    if abs(wrapped - 1.0) <= _FRACTIONAL_WRAP_TOLERANCE:
        return 0.0
    return wrapped


def _same_fractional_position(
    first: tuple[float, ...],
    second: tuple[float, ...],
    *,
    tolerance: float = _FRACTIONAL_WRAP_TOLERANCE,
) -> bool:
    if len(first) != len(second):
        return False
    return all(abs(_canonical_fractional(a) - _canonical_fractional(b)) <= tolerance for a, b in zip(first, second))


class _PositionIndex:
    """Spatial candidates for the existing coordinate-equality predicate.

    Preserve the first matching position in insertion order, including ties;
    bins only filter candidates and never replace the exact predicate.
    """

    def __init__(self):
        self.positions: list[tuple[float, ...]] = []
        self.bins: dict[tuple[int, ...], list[int]] = {}

    def _key(self, position):
        return tuple(math.floor(_canonical_fractional(value) / _FRACTIONAL_WRAP_TOLERANCE) for value in position)

    def add(self, position):
        self.bins.setdefault(self._key(position), []).append(len(self.positions))
        self.positions.append(position)

    def match(self, position):
        # Derived: |a-b| <= tolerance spans at most one tolerance-wide bin.
        # Include an extra bin per side for rounded bin-boundary division.
        keys = product(*(range(value - 2, value + 3) for value in self._key(position)))
        candidates = sorted(index for key in keys for index in self.bins.get(key, ()))
        for index in candidates:
            if _same_fractional_position(position, self.positions[index]):
                return index
        return None


def _match_position(
    position: tuple[float, ...],
    node_sites: tuple[ExpandedPeriodicNodeSite, ...],
) -> int | None:
    for index, node_site in enumerate(node_sites):
        if _same_fractional_position(position, node_site.fractional_position):
            return index
    return None


def _fractional_image(
    unwrapped: tuple[float, ...],
    wrapped: tuple[float, ...],
) -> tuple[int, ...]:
    return tuple(int(round(unwrapped[index] - wrapped[index])) for index in range(len(unwrapped)))


def _canonical_periodic_edge(
    start_index: int,
    end_index: int,
    end_image: tuple[int, ...],
) -> tuple[int, int, tuple[int, ...]]:
    forward = (start_index, end_index, end_image)
    reverse = (end_index, start_index, tuple(-component for component in end_image))
    return reverse if reverse < forward else forward


__all__ = [
    "ExpandedPeriodicEdge",
    "ExpandedPeriodicNodeSite",
    "ExpandedPeriodicTopology",
    "expand_topology_definition",
    "expand_two_d_topology_definition",
    "expand_three_d_topology_definition",
    "has_gemmi",
]
