from __future__ import annotations

from .planner import TopologyHint
from .topology_data import (
    load_workspace_topology_import_metadata,
    write_workspace_topology_index,
    workspace_topology_data_exists,
    workspace_topology_data_dir,
)
from .topology_importers import RCSRArchiveImporter, configured_rcsr_archives, discover_rcsr_archives
from .topology_index import TopologyDefinition, TopologyEdgeDefinition, TopologyIndexEntry, TopologyNodeDefinition
from .topology_repository import BuiltinTopologyFallback, TopologyRepository

_DEFAULT_REPOSITORY: TopologyRepository | None = None


def default_topology_repository() -> TopologyRepository:
    global _DEFAULT_REPOSITORY
    if _DEFAULT_REPOSITORY is None:
        if workspace_topology_data_exists():
            _DEFAULT_REPOSITORY = TopologyRepository.from_workspace_data()
        else:
            _DEFAULT_REPOSITORY = TopologyRepository.from_rcsr_archives(
                archive_paths=configured_rcsr_archives()
            )
    return _DEFAULT_REPOSITORY


def imported_topology_data_dir() -> str:
    return str(workspace_topology_data_dir())


def imported_topology_data_metadata() -> dict[str, object]:
    return dict(load_workspace_topology_import_metadata())


def rebuild_imported_topology_index() -> str:
    return str(write_workspace_topology_index())


def get_topology_hint(topology_id: str) -> TopologyHint:
    return default_topology_repository().get_hint(topology_id)


def all_topology_hints() -> tuple[TopologyHint, ...]:
    return default_topology_repository().list_hints()


def list_topology_index(
    dimensionality: str | None = None,
    node_connectivities: tuple[int, ...] = (),
) -> tuple[TopologyIndexEntry, ...]:
    return default_topology_repository().list_index(
        dimensionality=dimensionality,
        node_connectivities=node_connectivities,
    )


def load_topology(topology_id: str) -> TopologyDefinition:
    return default_topology_repository().load(topology_id)
