from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping

from .planner import TopologyHint


@dataclass(frozen=True)
class TopologyNodeDefinition:
    label: str
    position: tuple[float, ...]
    connectivity: int
    metadata: Mapping[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class TopologyEdgeDefinition:
    start: tuple[float, ...]
    end: tuple[float, ...]
    start_label: str | None = None
    end_label: str | None = None
    metadata: Mapping[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class TopologyDefinition:
    id: str
    name: str
    dimensionality: str
    node_definitions: tuple[TopologyNodeDefinition, ...]
    edge_definitions: tuple[TopologyEdgeDefinition, ...]
    metadata: Mapping[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class TopologyIndexEntry:
    id: str
    name: str
    dimensionality: str
    n_node_definitions: int
    n_edge_definitions: int
    node_connectivities: tuple[int, ...]
    space_group: str | None = None
    cell_parameters: tuple[float, ...] = ()
    source_archive: str | None = None
    source_member: str | None = None
    metadata: Mapping[str, object] = field(default_factory=dict)

    def to_hint(self) -> TopologyHint:
        hint_metadata = {
            "n_node_definitions": self.n_node_definitions,
            "n_edge_definitions": self.n_edge_definitions,
            "node_connectivities": self.node_connectivities,
            "space_group": self.space_group,
            "cell_parameters": self.cell_parameters,
            "source_archive": self.source_archive,
            "source_member": self.source_member,
        }
        hint_metadata.update(self.metadata)
        return TopologyHint(
            id=self.id,
            dimensionality=self.dimensionality,
            node_coordination=self.node_connectivities,
            metadata=hint_metadata,
        )
