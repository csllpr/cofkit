from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Mapping

from .planner import TopologyHint
from .topology_data import load_workspace_topology_definition, load_workspace_topology_index
from .topology_analysis import (
    annotate_topology_definition_zero_linker_metadata,
    annotate_topology_index_entry_zero_linker_metadata,
    has_zero_linker_metadata,
)
from .topology_importers import RCSRArchiveImporter, discover_rcsr_archives
from .topology_index import TopologyDefinition, TopologyIndexEntry


@dataclass(frozen=True)
class BuiltinTopologyFallback:
    hint: TopologyHint
    description: str = ""


_BUILTIN_FALLBACKS: dict[str, BuiltinTopologyFallback] = {
    "hcb": BuiltinTopologyFallback(
        hint=TopologyHint(id="hcb", dimensionality="2D", node_coordination=(3, 3)),
        description="2D honeycomb net for paired 3-connected monomers",
    ),
    "sql": BuiltinTopologyFallback(
        hint=TopologyHint(id="sql", dimensionality="2D", node_coordination=(4, 4)),
        description="2D square lattice for paired 4-connected monomers",
    ),
    "kgm": BuiltinTopologyFallback(
        hint=TopologyHint(id="kgm", dimensionality="2D", node_coordination=(3, 6)),
        description="2D kagome-like mixed 3/6-connected net",
    ),
}


class TopologyRepository:
    def __init__(
        self,
        entries: tuple[TopologyIndexEntry, ...] = (),
        definitions: Mapping[str, TopologyDefinition] | None = None,
        builtin_fallbacks: Mapping[str, BuiltinTopologyFallback] | None = None,
        definition_loader: Callable[[str], TopologyDefinition] | None = None,
    ):
        self._entries = {entry.id: entry for entry in entries}
        self._definitions = dict(definitions or {})
        self._builtin_fallbacks = dict(builtin_fallbacks or _BUILTIN_FALLBACKS)
        self._definition_loader = definition_loader

    @classmethod
    def from_rcsr_archives(
        cls,
        archive_paths: Mapping[str, str | Path] | None = None,
        importer: RCSRArchiveImporter | None = None,
        include_builtin_fallbacks: bool = True,
    ) -> "TopologyRepository":
        loader = importer or RCSRArchiveImporter()
        discovered = discover_rcsr_archives() if archive_paths is None else archive_paths
        entries: list[TopologyIndexEntry] = []
        definitions: dict[str, TopologyDefinition] = {}

        for dimensionality, archive_path in discovered.items():
            try:
                for definition in loader.iter_definitions(archive_path=archive_path, dimensionality=dimensionality):
                    definitions[definition.id] = definition
                    entries.append(loader.definition_to_entry(definition))
            except OSError:
                continue
            except ValueError:
                continue
        fallback_map = _BUILTIN_FALLBACKS if include_builtin_fallbacks else {}
        return cls(entries=tuple(entries), definitions=definitions, builtin_fallbacks=fallback_map)

    @classmethod
    def from_workspace_data(
        cls,
        data_dir: str | Path | None = None,
        importer: RCSRArchiveImporter | None = None,
        include_builtin_fallbacks: bool = True,
    ) -> "TopologyRepository":
        entries = load_workspace_topology_index(data_dir=data_dir)
        parser = importer or RCSRArchiveImporter()
        fallback_map = _BUILTIN_FALLBACKS if include_builtin_fallbacks else {}
        return cls(
            entries=entries,
            definitions={},
            builtin_fallbacks=fallback_map,
            definition_loader=lambda topology_id: load_workspace_topology_definition(
                topology_id=topology_id,
                data_dir=data_dir,
                importer=parser,
            ),
        )

    def list_index(
        self,
        dimensionality: str | None = None,
        node_connectivities: tuple[int, ...] = (),
    ) -> tuple[TopologyIndexEntry, ...]:
        entries = tuple(self._entries.values())
        if dimensionality is not None:
            entries = tuple(entry for entry in entries if entry.dimensionality == dimensionality)
        if node_connectivities:
            wanted = tuple(sorted(node_connectivities))
            entries = tuple(entry for entry in entries if tuple(sorted(entry.node_connectivities)) == wanted)
        if entries:
            return tuple(sorted(entries, key=lambda entry: entry.id))

        if dimensionality is not None:
            fallback_entries = tuple(
                self._fallback_entry(topology_id)
                for topology_id, fallback in sorted(self._builtin_fallbacks.items())
                if fallback.hint.dimensionality == dimensionality
                and (not node_connectivities or tuple(sorted(fallback.hint.node_coordination)) == tuple(sorted(node_connectivities)))
            )
        else:
            fallback_entries = tuple(
                self._fallback_entry(topology_id) for topology_id in sorted(self._builtin_fallbacks)
            )
        return fallback_entries

    def list_hints(
        self,
        dimensionality: str | None = None,
        node_connectivities: tuple[int, ...] = (),
    ) -> tuple[TopologyHint, ...]:
        return tuple(entry.to_hint() for entry in self.list_index(dimensionality, node_connectivities))

    def get_index_entry(self, topology_id: str) -> TopologyIndexEntry:
        if topology_id in self._entries:
            entry = self._entries[topology_id]
            if not has_zero_linker_metadata(entry.metadata):
                if topology_id in self._definitions:
                    entry = annotate_topology_index_entry_zero_linker_metadata(entry, self._definitions[topology_id])
                elif self._definition_loader is not None:
                    definition = self.load(topology_id)
                    entry = annotate_topology_index_entry_zero_linker_metadata(entry, definition)
                self._entries[topology_id] = entry
            return entry
        if topology_id in self._builtin_fallbacks:
            return self._fallback_entry(topology_id)
        raise KeyError(f"unknown topology {topology_id!r}")

    def load(self, topology_id: str) -> TopologyDefinition:
        if topology_id in self._definitions:
            return self._definitions[topology_id]
        if topology_id in self._entries and self._definition_loader is not None:
            definition = annotate_topology_definition_zero_linker_metadata(self._definition_loader(topology_id))
            self._definitions[topology_id] = definition
            self._entries[topology_id] = annotate_topology_index_entry_zero_linker_metadata(
                self._entries[topology_id],
                definition,
            )
            return definition
        if topology_id in self._builtin_fallbacks:
            fallback = self._builtin_fallbacks[topology_id]
            return annotate_topology_definition_zero_linker_metadata(
                TopologyDefinition(
                    id=fallback.hint.id,
                    name=fallback.hint.id,
                    dimensionality=fallback.hint.dimensionality,
                    node_definitions=tuple(),
                    edge_definitions=tuple(),
                    metadata={"description": fallback.description, "source": "builtin-fallback"},
                )
            )
        raise KeyError(f"unknown topology {topology_id!r}")

    def get_hint(self, topology_id: str) -> TopologyHint:
        return self.get_index_entry(topology_id).to_hint()

    def _fallback_entry(self, topology_id: str) -> TopologyIndexEntry:
        fallback = self._builtin_fallbacks[topology_id]
        return TopologyIndexEntry(
            id=fallback.hint.id,
            name=fallback.hint.id,
            dimensionality=fallback.hint.dimensionality,
            n_node_definitions=len(fallback.hint.node_coordination),
            n_edge_definitions=0,
            node_connectivities=fallback.hint.node_coordination,
            metadata={"description": fallback.description, "source": "builtin-fallback"},
        )
