from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Mapping

from .topology_index import TopologyDefinition, TopologyEdgeDefinition, TopologyIndexEntry


_ZERO_LINKER_METADATA_KEYS = (
    "zero_linker_compatible",
    "zero_linker_bipartite",
    "zero_linker_role_count_lower_bound",
    "zero_linker_reason",
    "zero_linker_scan_method",
    "zero_linker_scan_version",
    "zero_linker_builder_supported",
    "two_monomer_compatible",
    "two_monomer_modes",
    "two_monomer_node_node_modes",
    "two_monomer_node_linker_modes",
    "two_monomer_reason",
    "two_monomer_node_node_reason",
    "two_monomer_node_linker_reason",
)
_ZERO_LINKER_SCAN_VERSION = 3
_FRACTIONAL_WRAP_TOLERANCE = 1e-4
_SUPPORTED_2D_SINGLE_NODE_IDS = {"hcb", "hca", "fes", "fxt", "sql", "kgm", "htb", "hxl"}
_SUPPORTED_3D_SINGLE_NODE_NODE_NODE_SUPPORT = {"dia": True, "pcu": False}


@dataclass(frozen=True)
class _PeriodicNodeGraph:
    node_connectivities: tuple[int, ...]
    source_node_indices: tuple[int, ...]
    gain_edges: tuple[tuple[int, int, tuple[int, ...]], ...]
    dimension_count: int


@dataclass(frozen=True)
class ZeroLinkerScanResult:
    compatible: bool | None
    bipartite: bool | None
    role_count_lower_bound: int | None
    reason: str
    scan_method: str
    builder_supported: bool
    two_monomer_compatible: bool | None
    two_monomer_modes: tuple[str, ...]
    two_monomer_node_node_modes: tuple[str, ...]
    two_monomer_node_linker_modes: tuple[str, ...]
    two_monomer_reason: str
    two_monomer_node_node_reason: str
    two_monomer_node_linker_reason: str

    def to_metadata(self) -> dict[str, object]:
        return {
            "zero_linker_compatible": self.compatible,
            "zero_linker_bipartite": self.bipartite,
            "zero_linker_role_count_lower_bound": self.role_count_lower_bound,
            "zero_linker_reason": self.reason,
            "zero_linker_scan_method": self.scan_method,
            "zero_linker_scan_version": _ZERO_LINKER_SCAN_VERSION,
            "zero_linker_builder_supported": self.builder_supported,
            "two_monomer_compatible": self.two_monomer_compatible,
            "two_monomer_modes": list(self.two_monomer_modes),
            "two_monomer_node_node_modes": list(self.two_monomer_node_node_modes),
            "two_monomer_node_linker_modes": list(self.two_monomer_node_linker_modes),
            "two_monomer_reason": self.two_monomer_reason,
            "two_monomer_node_node_reason": self.two_monomer_node_node_reason,
            "two_monomer_node_linker_reason": self.two_monomer_node_linker_reason,
        }


def has_zero_linker_metadata(metadata: object) -> bool:
    if not isinstance(metadata, Mapping):
        return False
    return metadata.get("zero_linker_scan_version") == _ZERO_LINKER_SCAN_VERSION and all(
        key in metadata for key in _ZERO_LINKER_METADATA_KEYS
    )


def annotate_topology_definition_zero_linker_metadata(definition: TopologyDefinition) -> TopologyDefinition:
    if has_zero_linker_metadata(definition.metadata):
        return definition
    metadata = dict(definition.metadata)
    metadata.update(scan_zero_linker_compatibility(definition).to_metadata())
    return replace(definition, metadata=metadata)


def annotate_topology_index_entry_zero_linker_metadata(
    entry: TopologyIndexEntry,
    definition: TopologyDefinition | None = None,
) -> TopologyIndexEntry:
    if has_zero_linker_metadata(entry.metadata):
        return entry
    scan = scan_zero_linker_compatibility(definition) if definition is not None else _entry_only_scan(entry)
    metadata = dict(entry.metadata)
    metadata.update(scan.to_metadata())
    return replace(entry, metadata=metadata)


def scan_zero_linker_compatibility(definition: TopologyDefinition) -> ZeroLinkerScanResult:
    special_case = _scan_supported_single_node_topology(definition)
    if special_case is not None:
        return special_case

    explicit_scan = _scan_explicit_quotient_graph(definition)
    if explicit_scan is not None:
        return explicit_scan

    symmetry_expanded_scan = _scan_symmetry_expanded_topology(definition)
    if symmetry_expanded_scan is not None:
        return symmetry_expanded_scan

    return _unavailable_scan_result(
        builder_supported=_builder_support_from_definition(definition),
        reason=(
            "asymmetric-unit edge endpoints do not resolve to an explicit periodic node quotient graph without "
            "additional space-group expansion"
        ),
    )


def _unavailable_scan_result(
    *,
    builder_supported: bool,
    reason: str,
) -> ZeroLinkerScanResult:
    return ZeroLinkerScanResult(
        compatible=None,
        bipartite=None,
        role_count_lower_bound=None,
        reason=reason,
        scan_method="unavailable",
        builder_supported=builder_supported,
        two_monomer_compatible=None,
        two_monomer_modes=(),
        two_monomer_node_node_modes=(),
        two_monomer_node_linker_modes=(),
        two_monomer_reason=reason,
        two_monomer_node_node_reason=reason,
        two_monomer_node_linker_reason=reason,
    )


def _entry_only_scan(entry: TopologyIndexEntry) -> ZeroLinkerScanResult:
    if entry.n_node_definitions > 1:
        return _unavailable_scan_result(
            builder_supported=True,
            reason="bulk index entry has not been hydrated with a full topology definition yet",
        )
    if entry.id in _SUPPORTED_3D_SINGLE_NODE_NODE_NODE_SUPPORT:
        return _unavailable_scan_result(
            builder_supported=_SUPPORTED_3D_SINGLE_NODE_NODE_NODE_SUPPORT[entry.id],
            reason="bulk index entry has not been hydrated with a full topology definition yet",
        )
    return _unavailable_scan_result(
        builder_supported=False,
        reason="bulk index entry has not been hydrated with a full topology definition yet",
    )


def _scan_result_from_periodic_graph(
    *,
    graph: _PeriodicNodeGraph,
    scan_method: str,
    builder_supported: bool,
) -> ZeroLinkerScanResult:
    bipartite = _is_periodic_gain_graph_bipartite(
        node_count=len(graph.node_connectivities),
        gain_edges=graph.gain_edges,
        dimension_count=graph.dimension_count,
    )
    (
        two_monomer_compatible,
        two_monomer_modes,
        node_node_modes,
        node_linker_modes,
        two_monomer_reason,
        node_node_reason,
        node_linker_reason,
    ) = _two_monomer_modes_from_periodic_graph(graph, bipartite)
    return ZeroLinkerScanResult(
        compatible=bipartite,
        bipartite=bipartite,
        role_count_lower_bound=2 if bipartite else 3,
        reason=(
            "periodic node graph is bipartite"
            if bipartite
            else "periodic node graph is not bipartite"
        ),
        scan_method=scan_method,
        builder_supported=builder_supported,
        two_monomer_compatible=two_monomer_compatible,
        two_monomer_modes=two_monomer_modes,
        two_monomer_node_node_modes=node_node_modes,
        two_monomer_node_linker_modes=node_linker_modes,
        two_monomer_reason=two_monomer_reason,
        two_monomer_node_node_reason=node_node_reason,
        two_monomer_node_linker_reason=node_linker_reason,
    )


def _scan_supported_single_node_topology(definition: TopologyDefinition) -> ZeroLinkerScanResult | None:
    if len(definition.node_definitions) != 1:
        return None

    topology_id = definition.id
    if definition.dimensionality == "2D" and topology_id in _SUPPORTED_2D_SINGLE_NODE_IDS:
        from .single_node_topologies import expand_single_node_topology_definition

        expanded = expand_single_node_topology_definition(definition)
        return _scan_result_from_periodic_graph(
            graph=_periodic_graph_from_single_node_expanded(
                node_count=len(expanded.node_sites),
                gain_edges=tuple(
                    (
                        _node_index_by_id(expanded.node_sites)[edge.start_node_id],
                        _node_index_by_id(expanded.node_sites)[edge.end_node_id],
                        edge.end_image[:2],
                    )
                    for edge in expanded.edge_sites
                ),
                connectivity=int(definition.node_definitions[0].connectivity),
                dimension_count=2,
            ),
            scan_method="expanded-p1-periodic-bipartite-scan",
            builder_supported=expanded.is_bipartite,
        )

    if definition.dimensionality == "3D" and topology_id in _SUPPORTED_3D_SINGLE_NODE_NODE_NODE_SUPPORT:
        from .single_node_topologies_3d import expand_three_d_single_node_topology_definition

        expanded = expand_three_d_single_node_topology_definition(definition)
        return _scan_result_from_periodic_graph(
            graph=_periodic_graph_from_single_node_expanded(
                node_count=len(expanded.node_sites),
                gain_edges=tuple(
                    (
                        _node_index_by_id(expanded.node_sites)[edge.start_node_id],
                        _node_index_by_id(expanded.node_sites)[edge.end_node_id],
                        edge.end_image,
                    )
                    for edge in expanded.edge_sites
                ),
                connectivity=int(definition.node_definitions[0].connectivity),
                dimension_count=3,
            ),
            scan_method="expanded-p1-periodic-bipartite-scan",
            builder_supported=_SUPPORTED_3D_SINGLE_NODE_NODE_NODE_SUPPORT[topology_id],
        )
    return None


def _periodic_graph_from_single_node_expanded(
    *,
    node_count: int,
    gain_edges: tuple[tuple[int, int, tuple[int, ...]], ...],
    connectivity: int,
    dimension_count: int,
) -> _PeriodicNodeGraph:
    return _PeriodicNodeGraph(
        node_connectivities=(connectivity,) * node_count,
        source_node_indices=(0,) * node_count,
        gain_edges=gain_edges,
        dimension_count=dimension_count,
    )


def _scan_explicit_quotient_graph(definition: TopologyDefinition) -> ZeroLinkerScanResult | None:
    if not definition.node_definitions:
        return _unavailable_scan_result(
            builder_supported=False,
            reason="topology definition does not include explicit node definitions",
        )
    if not definition.edge_definitions:
        return _scan_result_from_periodic_graph(
            graph=_PeriodicNodeGraph(
                node_connectivities=tuple(int(node.connectivity) for node in definition.node_definitions),
                source_node_indices=tuple(range(len(definition.node_definitions))),
                gain_edges=(),
                dimension_count=2 if definition.dimensionality == "2D" else 3,
            ),
            scan_method="quotient-graph-parity-scan",
            builder_supported=False,
        )

    dimension_count = 2 if definition.dimensionality == "2D" else 3
    node_positions = tuple(_normalize_position(node.position, dimension_count) for node in definition.node_definitions)
    label_to_indices: dict[str, list[int]] = {}
    for index, node in enumerate(definition.node_definitions):
        label_to_indices.setdefault(node.label, []).append(index)

    gain_edges: list[tuple[int, int, tuple[int, ...]]] = []
    for edge in definition.edge_definitions:
        if not edge.start or not edge.end:
            return None
        start_index = _match_explicit_endpoint(edge.start_label, edge.start, label_to_indices, node_positions, dimension_count)
        end_index = _match_explicit_endpoint(edge.end_label, edge.end, label_to_indices, node_positions, dimension_count)
        if start_index is None or end_index is None:
            return None
        gain_edges.append(
            (
                start_index,
                end_index,
                _relative_image(edge, node_positions[start_index], node_positions[end_index], dimension_count),
            )
        )
    return _scan_result_from_periodic_graph(
        graph=_PeriodicNodeGraph(
            node_connectivities=tuple(int(node.connectivity) for node in definition.node_definitions),
            source_node_indices=tuple(range(len(definition.node_definitions))),
            gain_edges=tuple(gain_edges),
            dimension_count=dimension_count,
        ),
        scan_method="quotient-graph-parity-scan",
        builder_supported=_builder_support_from_definition(definition),
    )


def _scan_symmetry_expanded_topology(definition: TopologyDefinition) -> ZeroLinkerScanResult | None:
    try:
        from .topology_symmetry import expand_topology_definition
    except ImportError:
        return None

    try:
        expanded = expand_topology_definition(definition)
    except (ImportError, KeyError, ValueError):
        return None

    node_index_by_id = {
        node_site.id: index
        for index, node_site in enumerate(expanded.node_sites)
    }
    return _scan_result_from_periodic_graph(
        graph=_PeriodicNodeGraph(
            node_connectivities=tuple(int(node_site.connectivity) for node_site in expanded.node_sites),
            source_node_indices=tuple(int(node_site.source_node_index) for node_site in expanded.node_sites),
            gain_edges=tuple(
                (
                    node_index_by_id[edge.start_node_id],
                    node_index_by_id[edge.end_node_id],
                    edge.end_image,
                )
                for edge in expanded.edge_sites
            ),
            dimension_count=2 if definition.dimensionality == "2D" else 3,
        ),
        scan_method="gemmi-space-group-expanded-p1-scan",
        builder_supported=_builder_support_from_definition(definition),
    )


def _builder_support_from_definition(definition: TopologyDefinition) -> bool:
    if len(definition.node_definitions) > 1:
        return True
    return _SUPPORTED_3D_SINGLE_NODE_NODE_NODE_SUPPORT.get(definition.id, False)


def _node_index_by_id(node_sites) -> dict[str, int]:
    return {
        node_site.id: index
        for index, node_site in enumerate(node_sites)
    }


def _two_monomer_modes_from_periodic_graph(
    graph: _PeriodicNodeGraph,
    bipartite: bool,
) -> tuple[bool, tuple[str, ...], tuple[str, ...], tuple[str, ...], str, str, str]:
    if not graph.gain_edges:
        node_node_reason = "topology does not include any edges, so no node-node pairing can form a network"
        node_linker_reason = "topology does not include any edges, so no node-linker pairing can form a network"
        summary = f"node-node: {node_node_reason}; node-linker: {node_linker_reason}"
        return False, (), (), (), summary, node_node_reason, node_linker_reason

    distinct_connectivities = tuple(sorted(set(int(value) for value in graph.node_connectivities)))

    if len(distinct_connectivities) == 1:
        connectivity = distinct_connectivities[0]
        node_linker_modes = (f"{connectivity}+2",)
        node_linker_reason = (
            f"all node sites have connectivity {connectivity}, so a {connectivity}+2 node-linker pairing is "
            "chemically consistent"
        )
    else:
        node_linker_modes = ()
        node_linker_reason = (
            "node-linker pairing requires one uniform node connectivity, but the topology contains "
            f"{_format_connectivity_set(distinct_connectivities)}"
        )

    if len(distinct_connectivities) == 1:
        connectivity = distinct_connectivities[0]
        if bipartite:
            node_node_modes = (f"{connectivity}+{connectivity}",)
            node_node_reason = (
                f"the periodic graph supports an alternating {connectivity}+{connectivity} node-node pairing"
            )
        else:
            node_node_modes = ()
            node_node_reason = (
                f"the periodic graph is not bipartite, so an alternating {connectivity}+{connectivity} node-node "
                "pairing is not possible"
            )
    elif len(distinct_connectivities) == 2:
        pair = (distinct_connectivities[0], distinct_connectivities[1])
        if _supports_fixed_connectivity_node_node_pair(graph, pair):
            node_node_modes = (f"{pair[0]}+{pair[1]}",)
            node_node_reason = (
                f"the periodic graph admits a fixed-connectivity {pair[0]}+{pair[1]} node-node alternation"
            )
        elif bipartite:
            node_node_modes = ()
            node_node_reason = (
                "the periodic graph is bipartite only through a coloring that mixes node connectivities across the "
                "two monomer roles"
            )
        else:
            node_node_modes = ()
            node_node_reason = (
                f"the periodic graph is not bipartite, so a fixed-connectivity {pair[0]}+{pair[1]} node-node "
                "pairing is not possible"
            )
    else:
        node_node_modes = ()
        node_node_reason = (
            "node-node pairing requires at most two monomer connectivities, but the topology contains "
            f"{_format_connectivity_set(distinct_connectivities)}"
        )

    two_monomer_modes = tuple(
        [f"node-node:{mode}" for mode in node_node_modes]
        + [f"node-linker:{mode}" for mode in node_linker_modes]
    )
    two_monomer_compatible = bool(two_monomer_modes)
    if node_node_modes and node_linker_modes:
        summary = (
            f"supports node-node {', '.join(node_node_modes)} and node-linker {', '.join(node_linker_modes)} "
            "two-monomer pairings"
        )
    elif node_node_modes:
        summary = f"supports node-node {', '.join(node_node_modes)} two-monomer pairings"
    elif node_linker_modes:
        summary = f"supports node-linker {', '.join(node_linker_modes)} two-monomer pairings"
    else:
        summary = f"node-node: {node_node_reason}; node-linker: {node_linker_reason}"
    return (
        two_monomer_compatible,
        two_monomer_modes,
        node_node_modes,
        node_linker_modes,
        summary,
        node_node_reason,
        node_linker_reason,
    )


def _supports_fixed_connectivity_node_node_pair(
    graph: _PeriodicNodeGraph,
    connectivity_pair: tuple[int, int],
) -> bool:
    first, second = connectivity_pair
    if first == second:
        return _is_periodic_gain_graph_bipartite(
            node_count=len(graph.node_connectivities),
            gain_edges=graph.gain_edges,
            dimension_count=graph.dimension_count,
        )

    for assignment in ({first: 0, second: 1}, {first: 1, second: 0}):
        rows: list[int] = []
        rhs: list[int] = []
        for start_index, end_index, shift in graph.gain_edges:
            start_connectivity = graph.node_connectivities[start_index]
            end_connectivity = graph.node_connectivities[end_index]
            start_color = assignment.get(start_connectivity)
            end_color = assignment.get(end_connectivity)
            if start_color is None or end_color is None:
                return False
            row = 0
            for axis, component in enumerate(shift):
                if component % 2:
                    row ^= 1 << axis
            rows.append(row)
            rhs.append(1 ^ start_color ^ end_color)
        if _solve_binary_linear_system(rows, rhs, graph.dimension_count):
            return True
    return False


def _format_connectivity_set(connectivities: tuple[int, ...]) -> str:
    return ", ".join(str(value) for value in connectivities)


def _normalize_position(position: tuple[float, ...], dimension_count: int) -> tuple[float, ...]:
    if len(position) >= dimension_count:
        return tuple(float(position[index]) for index in range(dimension_count))
    return tuple(float(value) for value in position) + (0.0,) * (dimension_count - len(position))


def _match_explicit_endpoint(
    label: str | None,
    position: tuple[float, ...],
    label_to_indices: dict[str, list[int]],
    node_positions: tuple[tuple[float, ...], ...],
    dimension_count: int,
) -> int | None:
    normalized_position = _normalize_position(position, dimension_count)
    wrapped = _wrap_fractional_position(normalized_position)

    if label:
        label_matches = label_to_indices.get(label, [])
        if len(label_matches) == 1:
            return label_matches[0]
        if len(label_matches) > 1:
            matched = [index for index in label_matches if _same_fractional_position(node_positions[index], wrapped)]
            if len(matched) == 1:
                return matched[0]
            return None

    matches = [index for index, node_position in enumerate(node_positions) if _same_fractional_position(node_position, wrapped)]
    if len(matches) == 1:
        return matches[0]
    return None


def _relative_image(
    edge: TopologyEdgeDefinition,
    start_position: tuple[float, ...],
    end_position: tuple[float, ...],
    dimension_count: int,
) -> tuple[int, ...]:
    start = _normalize_position(edge.start, dimension_count)
    end = _normalize_position(edge.end, dimension_count)
    start_image = _fractional_image(start, start_position)
    end_image = _fractional_image(end, end_position)
    return tuple(end_component - start_component for start_component, end_component in zip(start_image, end_image))


def _fractional_image(
    unwrapped: tuple[float, ...],
    node_position: tuple[float, ...],
) -> tuple[int, ...]:
    return tuple(int(round(component - anchor)) for component, anchor in zip(unwrapped, node_position))


def _wrap_fractional_position(position: tuple[float, ...]) -> tuple[float, ...]:
    return tuple(_canonical_fractional(value) for value in position)


def _same_fractional_position(
    first: tuple[float, ...],
    second: tuple[float, ...],
    *,
    tolerance: float = 1e-4,
) -> bool:
    if len(first) != len(second):
        return False
    return all(abs(_canonical_fractional(a) - _canonical_fractional(b)) <= tolerance for a, b in zip(first, second))


def _canonical_fractional(value: float) -> float:
    if abs(value - round(value)) <= _FRACTIONAL_WRAP_TOLERANCE:
        return 0.0
    wrapped = value % 1.0
    if abs(wrapped - 1.0) <= _FRACTIONAL_WRAP_TOLERANCE:
        return 0.0
    return wrapped


def _is_periodic_gain_graph_bipartite(
    *,
    node_count: int,
    gain_edges: tuple[tuple[int, int, tuple[int, ...]], ...],
    dimension_count: int,
) -> bool:
    if not gain_edges:
        return True

    rows: list[int] = []
    rhs: list[int] = []
    variable_count = node_count + dimension_count
    for start_index, end_index, shift in gain_edges:
        row = 0
        row ^= 1 << start_index
        row ^= 1 << end_index
        for axis, component in enumerate(shift):
            if component % 2:
                row ^= 1 << (node_count + axis)
        rows.append(row)
        rhs.append(1)
    return _solve_binary_linear_system(rows, rhs, variable_count)


def _solve_binary_linear_system(
    rows: list[int],
    rhs: list[int],
    variable_count: int,
) -> bool:
    pivot_row = 0
    equation_count = len(rows)

    for column in range(variable_count):
        pivot_index = None
        for row_index in range(pivot_row, equation_count):
            if (rows[row_index] >> column) & 1:
                pivot_index = row_index
                break
        if pivot_index is None:
            continue
        rows[pivot_row], rows[pivot_index] = rows[pivot_index], rows[pivot_row]
        rhs[pivot_row], rhs[pivot_index] = rhs[pivot_index], rhs[pivot_row]

        for row_index in range(equation_count):
            if row_index == pivot_row:
                continue
            if (rows[row_index] >> column) & 1:
                rows[row_index] ^= rows[pivot_row]
                rhs[row_index] ^= rhs[pivot_row]
        pivot_row += 1

    for row, value in zip(rows, rhs):
        if row == 0 and value:
            return False
    return True


__all__ = [
    "ZeroLinkerScanResult",
    "annotate_topology_definition_zero_linker_metadata",
    "annotate_topology_index_entry_zero_linker_metadata",
    "has_zero_linker_metadata",
    "scan_zero_linker_compatibility",
]
