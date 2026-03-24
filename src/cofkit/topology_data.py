from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
from typing import Mapping

from .topology_importers import RCSRArchiveImporter
from .topology_index import TopologyDefinition, TopologyIndexEntry


_DEFAULT_TOPOLOGY_DATA_DIR = Path(__file__).resolve().parent / "data" / "topologies"
_INDEX_FILENAME = "index.json"
_IMPORT_METADATA_FILENAME = "import-metadata.json"
_BUNDLED_SOURCE_SCHEME = "bundled-topology://"


def workspace_topology_data_dir(data_dir: str | Path | None = None) -> Path:
    if data_dir is None:
        return _DEFAULT_TOPOLOGY_DATA_DIR
    return Path(data_dir).expanduser()


def workspace_topology_index_path(data_dir: str | Path | None = None) -> Path:
    return workspace_topology_data_dir(data_dir) / _INDEX_FILENAME


def workspace_topology_import_metadata_path(data_dir: str | Path | None = None) -> Path:
    return workspace_topology_data_dir(data_dir) / _IMPORT_METADATA_FILENAME


def workspace_topology_data_exists(data_dir: str | Path | None = None) -> bool:
    return workspace_topology_index_path(data_dir).is_file()


def load_workspace_topology_index(data_dir: str | Path | None = None) -> tuple[TopologyIndexEntry, ...]:
    payload = _read_json(workspace_topology_index_path(data_dir))
    data_root = workspace_topology_data_dir(data_dir)
    entries = []
    for entry in payload.get("entries", ()):
        item = dict(entry)
        item["node_connectivities"] = tuple(item.get("node_connectivities", ()))
        item["cell_parameters"] = tuple(item.get("cell_parameters", ()))
        item["source_archive"] = _normalize_workspace_source_archive(
            source_archive=item.get("source_archive"),
            workspace_file=item.get("metadata", {}).get("workspace_file"),
            data_root=data_root,
        )
        entries.append(TopologyIndexEntry(**item))
    return tuple(sorted(entries, key=lambda entry: entry.id))


def load_workspace_topology_import_metadata(data_dir: str | Path | None = None) -> Mapping[str, object]:
    path = workspace_topology_import_metadata_path(data_dir)
    if not path.is_file():
        return {}
    return _sanitize_workspace_import_metadata(_read_json(path))


def load_workspace_topology_definition(
    topology_id: str,
    data_dir: str | Path | None = None,
    importer: RCSRArchiveImporter | None = None,
) -> TopologyDefinition:
    index_payload = _read_json(workspace_topology_index_path(data_dir))
    entries = {
        entry["id"]: entry
        for entry in index_payload.get("entries", ())
        if isinstance(entry, dict) and "id" in entry
    }
    if topology_id not in entries:
        raise KeyError(f"unknown topology {topology_id!r}")

    entry = entries[topology_id]
    relative_path = entry.get("metadata", {}).get("workspace_file")
    if not relative_path:
        raise KeyError(f"workspace file metadata missing for topology {topology_id!r}")

    data_root = workspace_topology_data_dir(data_dir)
    topology_path = data_root / str(relative_path)
    parser = importer or RCSRArchiveImporter()
    for definition in parser.iter_path_definitions(
        topology_path=topology_path,
        dimensionality=entry["dimensionality"],
        source_archive=_bundled_source_reference(str(relative_path)),
        source_member=topology_path.name,
    ):
        if definition.id == topology_id:
            return definition
    raise KeyError(f"workspace topology file did not contain topology {topology_id!r}")


def build_workspace_topology_index_entries(
    data_dir: str | Path | None = None,
    importer: RCSRArchiveImporter | None = None,
) -> tuple[TopologyIndexEntry, ...]:
    data_root = workspace_topology_data_dir(data_dir)
    parser = importer or RCSRArchiveImporter()
    entries: list[TopologyIndexEntry] = []

    for dimensionality, relative_dir in (("2D", "2d"), ("3D", "3d")):
        topology_dir = data_root / relative_dir
        if not topology_dir.is_dir():
            continue
        for topology_path in sorted(topology_dir.glob("*.cgd")):
            relative_path = topology_path.relative_to(data_root).as_posix()
            try:
                definitions = tuple(
                    parser.iter_path_definitions(
                        topology_path=topology_path,
                        dimensionality=dimensionality,
                        source_archive=_bundled_source_reference(relative_path),
                        source_member=topology_path.name,
                    )
                )
            except (OSError, ValueError):
                continue
            for definition in definitions:
                entry = parser.definition_to_entry(definition)
                metadata = dict(entry.metadata)
                metadata["workspace_file"] = relative_path
                entries.append(
                    replace(
                        entry,
                        metadata=metadata,
                        source_archive=_bundled_source_reference(relative_path),
                        source_member=topology_path.name,
                    )
                )
    return tuple(sorted(entries, key=lambda entry: entry.id))


def build_workspace_topology_index_payload(
    data_dir: str | Path | None = None,
    importer: RCSRArchiveImporter | None = None,
) -> dict[str, object]:
    entries = build_workspace_topology_index_entries(data_dir=data_dir, importer=importer)
    return {
        "data_root": ".",
        "entries": [_topology_index_entry_to_payload(entry) for entry in entries],
        "format_version": 1,
    }


def write_workspace_topology_index(
    data_dir: str | Path | None = None,
    importer: RCSRArchiveImporter | None = None,
) -> Path:
    payload = build_workspace_topology_index_payload(data_dir=data_dir, importer=importer)
    path = workspace_topology_index_path(data_dir)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def _read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _bundled_source_reference(relative_path: str) -> str:
    return f"{_BUNDLED_SOURCE_SCHEME}{Path(relative_path).as_posix()}"


def _normalize_workspace_source_archive(
    source_archive: object,
    workspace_file: object,
    data_root: Path,
) -> str | None:
    if isinstance(workspace_file, str) and workspace_file:
        return _bundled_source_reference(workspace_file)
    if not isinstance(source_archive, str) or not source_archive:
        return None
    archive_path = Path(source_archive)
    try:
        relative = archive_path.resolve().relative_to(data_root.resolve())
    except Exception:
        return source_archive
    return _bundled_source_reference(relative.as_posix())


def _sanitize_workspace_import_metadata(payload: Mapping[str, object]) -> dict[str, object]:
    sanitized = dict(payload)

    archives_found = payload.get("archives_found")
    if isinstance(archives_found, Mapping):
        sanitized["archives_found"] = {
            str(key): _sanitize_archive_reference(value)
            for key, value in archives_found.items()
        }

    import_summary = payload.get("import_summary")
    if isinstance(import_summary, Mapping):
        sanitized_summary: dict[str, object] = {}
        for key, value in import_summary.items():
            if not isinstance(value, Mapping):
                sanitized_summary[str(key)] = value
                continue
            item = dict(value)
            if "archive" in item:
                item["archive"] = _sanitize_archive_reference(item.get("archive"))
            sanitized_summary[str(key)] = item
        sanitized["import_summary"] = sanitized_summary

    return sanitized


def _sanitize_archive_reference(value: object) -> object:
    if not isinstance(value, str) or not value:
        return value
    return Path(value).name


def _topology_index_entry_to_payload(entry: TopologyIndexEntry) -> dict[str, object]:
    return {
        "cell_parameters": list(entry.cell_parameters),
        "dimensionality": entry.dimensionality,
        "id": entry.id,
        "metadata": dict(entry.metadata),
        "n_edge_definitions": entry.n_edge_definitions,
        "n_node_definitions": entry.n_node_definitions,
        "name": entry.name,
        "node_connectivities": list(entry.node_connectivities),
        "source_archive": entry.source_archive,
        "source_member": entry.source_member,
        "space_group": entry.space_group,
    }
