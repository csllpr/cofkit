from __future__ import annotations

import math
import os
import shlex
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Mapping, Sequence

from .topology_index import (
    TopologyDefinition,
    TopologyEdgeDefinition,
    TopologyIndexEntry,
    TopologyNodeDefinition,
)
from .topology_analysis import annotate_topology_definition_zero_linker_metadata


_SUPPORTED_SUFFIXES = (".cgd", ".txt", ".net", ".arc", ".data")


@dataclass(frozen=True)
class _RawTopology:
    metadata: Mapping[str, object] = field(default_factory=dict)
    nodes: tuple[tuple[str, tuple[str, ...]], ...] = ()
    edges: tuple[tuple[str, tuple[str, ...]], ...] = ()


class RCSRArchiveImporter:
    """Parses RCSR topology bundles that contain Systre/CGD-style net definitions."""

    def build_index(
        self,
        archive_path: str | Path,
        dimensionality: str,
    ) -> tuple[TopologyIndexEntry, ...]:
        entries = [
            self._definition_to_entry(definition)
            for definition in self.iter_definitions(archive_path=archive_path, dimensionality=dimensionality)
        ]
        return tuple(sorted(entries, key=lambda entry: entry.id))

    def iter_definitions(
        self,
        archive_path: str | Path,
        dimensionality: str,
    ) -> Iterator[TopologyDefinition]:
        archive = Path(archive_path).expanduser()
        with zipfile.ZipFile(archive) as bundle:
            for member_name in bundle.namelist():
                if member_name.endswith("/"):
                    continue
                if not self._looks_like_supported_member(member_name):
                    continue
                payload = bundle.read(member_name)
                text = self._decode_payload(payload)
                try:
                    yield from self.iter_text_definitions(
                        text=text,
                        dimensionality=dimensionality,
                        source_archive=str(archive),
                        source_member=member_name,
                    )
                except ValueError:
                    continue

    def iter_text_definitions(
        self,
        text: str,
        dimensionality: str,
        source_archive: str,
        source_member: str,
    ) -> Iterator[TopologyDefinition]:
        for raw in self._parse_cgd_text(text):
            yield self._raw_to_definition(
                raw=raw,
                dimensionality=dimensionality,
                source_archive=source_archive,
                source_member=source_member,
            )

    def iter_path_definitions(
        self,
        topology_path: str | Path,
        dimensionality: str,
        source_archive: str | None = None,
        source_member: str | None = None,
    ) -> Iterator[TopologyDefinition]:
        path = Path(topology_path).expanduser()
        text = path.read_text(encoding="utf-8", errors="replace")
        yield from self.iter_text_definitions(
            text=text,
            dimensionality=dimensionality,
            source_archive=source_archive or str(path),
            source_member=source_member or path.name,
        )

    def definition_to_entry(self, definition: TopologyDefinition) -> TopologyIndexEntry:
        return self._definition_to_entry(definition)

    def _definition_to_entry(self, definition: TopologyDefinition) -> TopologyIndexEntry:
        return TopologyIndexEntry(
            id=definition.id,
            name=definition.name,
            dimensionality=definition.dimensionality,
            n_node_definitions=len(definition.node_definitions),
            n_edge_definitions=len(definition.edge_definitions),
            node_connectivities=tuple(node.connectivity for node in definition.node_definitions),
            space_group=self._as_optional_str(definition.metadata.get("space_group")),
            cell_parameters=tuple(definition.metadata.get("cell_parameters", ())),
            source_archive=self._as_optional_str(definition.metadata.get("source_archive")),
            source_member=self._as_optional_str(definition.metadata.get("source_member")),
            metadata={
                key: value
                for key, value in definition.metadata.items()
                if key not in {"space_group", "cell_parameters", "source_archive", "source_member"}
            },
        )

    def _raw_to_definition(
        self,
        raw: _RawTopology,
        dimensionality: str,
        source_archive: str,
        source_member: str,
    ) -> TopologyDefinition:
        dim = 2 if dimensionality.upper() == "2D" else 3
        topology_id = str(raw.metadata.get("name", raw.metadata.get("id", source_member))).strip().lower()
        name = str(raw.metadata.get("name", topology_id))

        node_records = [self._parse_node_record(record, dim) for record in raw.nodes]
        edge_records = [self._parse_edge_record(record, dim) for record in raw.edges]
        connectivities = self._infer_node_connectivities(node_records, edge_records)
        node_definitions = tuple(
            TopologyNodeDefinition(
                label=node["label"],
                position=node["position"],
                connectivity=connectivities[index],
                metadata=node["metadata"],
            )
            for index, node in enumerate(node_records)
        )
        edge_definitions = tuple(
            TopologyEdgeDefinition(
                start=edge["start"],
                end=edge["end"],
                start_label=edge.get("start_label"),
                end_label=edge.get("end_label"),
                metadata=edge["metadata"],
            )
            for edge in edge_records
        )

        metadata = dict(raw.metadata)
        metadata["source_archive"] = source_archive
        metadata["source_member"] = source_member
        definition = TopologyDefinition(
            id=topology_id,
            name=name,
            dimensionality=dimensionality,
            node_definitions=node_definitions,
            edge_definitions=edge_definitions,
            metadata=metadata,
        )
        return annotate_topology_definition_zero_linker_metadata(definition)

    def _parse_cgd_text(self, text: str) -> tuple[_RawTopology, ...]:
        topologies: list[_RawTopology] = []
        current_metadata: dict[str, object] = {}
        current_nodes: list[tuple[str, tuple[str, ...]]] = []
        current_edges: list[tuple[str, tuple[str, ...]]] = []
        inside_block = False

        for raw_line in text.splitlines():
            line = self._strip_comments(raw_line).strip()
            if not line:
                continue
            tokens = tuple(shlex.split(line, comments=False, posix=True))
            if not tokens:
                continue

            keyword = tokens[0].upper()
            values = tokens[1:]
            if keyword == "CRYSTAL":
                if current_metadata or current_nodes or current_edges:
                    topologies.append(
                        _RawTopology(
                            metadata=dict(current_metadata),
                            nodes=tuple(current_nodes),
                            edges=tuple(current_edges),
                        )
                    )
                    current_metadata = {}
                    current_nodes = []
                    current_edges = []
                inside_block = True
                current_metadata.setdefault("record_type", "CRYSTAL")
                continue
            if keyword in {"END", "ENDCRYSTAL"} and inside_block:
                topologies.append(
                    _RawTopology(
                        metadata=dict(current_metadata),
                        nodes=tuple(current_nodes),
                        edges=tuple(current_edges),
                    )
                )
                current_metadata = {}
                current_nodes = []
                current_edges = []
                inside_block = False
                continue

            if not inside_block:
                continue

            if keyword == "NAME" and values:
                current_metadata["name"] = values[0]
                continue
            if keyword == "ID" and values:
                current_metadata["id"] = values[0]
                continue
            if keyword == "GROUP" and values:
                current_metadata["space_group"] = " ".join(values)
                continue
            if keyword == "CELL":
                current_metadata["cell_parameters"] = tuple(self._parse_float(value) for value in values)
                continue
            if keyword in {"NODE", "ATOM"}:
                current_nodes.append((keyword, values))
                continue
            if keyword.startswith("EDGE"):
                current_edges.append((keyword, values))
                continue
            current_metadata.setdefault("extra_records", []).append((keyword, values))

        if current_metadata or current_nodes or current_edges:
            topologies.append(
                _RawTopology(
                    metadata=dict(current_metadata),
                    nodes=tuple(current_nodes),
                    edges=tuple(current_edges),
                )
            )
        return tuple(topologies)

    def _parse_node_record(
        self,
        record: tuple[str, tuple[str, ...]],
        dimension_count: int,
    ) -> dict[str, object]:
        _, values = record
        floats = self._float_indices(values)
        coords, coordinate_slice = self._extract_node_coordinates(values, floats, dimension_count)
        label_tokens = [value for index, value in enumerate(values) if index not in coordinate_slice]
        label = label_tokens[0] if label_tokens else f"node-{abs(hash(values))}"
        explicit_connectivity = None
        for token in label_tokens[1:]:
            if token.isdigit():
                explicit_connectivity = int(token)
                break
        return {
            "label": label,
            "position": coords,
            "metadata": {
                "explicit_connectivity": explicit_connectivity,
                "raw_record": values,
            },
        }

    def _parse_edge_record(
        self,
        record: tuple[str, tuple[str, ...]],
        dimension_count: int,
    ) -> dict[str, object]:
        keyword, values = record
        floats = self._float_indices(values)
        if len(floats) >= dimension_count * 2:
            start, end, coordinate_slice = self._extract_edge_coordinates(values, floats, dimension_count)
            labels = [value for index, value in enumerate(values) if index not in coordinate_slice]
            start_label = labels[0] if len(labels) >= 1 else None
            end_label = labels[1] if len(labels) >= 2 else None
        elif len(values) >= 2:
            start = ()
            end = ()
            start_label = values[0]
            end_label = values[1]
        else:
            raise ValueError(f"edge record is missing endpoints: {values!r}")
        return {
            "start": start,
            "end": end,
            "start_label": start_label,
            "end_label": end_label,
            "metadata": {"record_type": keyword, "raw_record": values},
        }

    def _infer_node_connectivities(
        self,
        node_records: Sequence[Mapping[str, object]],
        edge_records: Sequence[Mapping[str, object]],
    ) -> tuple[int, ...]:
        counts = [0 for _ in node_records]
        label_to_indices: dict[str, list[int]] = {}
        position_lookup = {
            index: tuple(node["position"]) for index, node in enumerate(node_records)
        }
        explicit_nodes: set[int] = set()
        for index, node in enumerate(node_records):
            label = str(node["label"])
            label_to_indices.setdefault(label, []).append(index)
            explicit = node["metadata"].get("explicit_connectivity")
            if isinstance(explicit, int) and explicit > 0:
                counts[index] = explicit
                explicit_nodes.add(index)

        for edge in edge_records:
            start_index = self._match_edge_endpoint(edge.get("start_label"), tuple(edge["start"]), label_to_indices, position_lookup)
            end_index = self._match_edge_endpoint(edge.get("end_label"), tuple(edge["end"]), label_to_indices, position_lookup)
            if start_index is not None and start_index not in explicit_nodes:
                counts[start_index] += 1
            if end_index is not None and end_index not in explicit_nodes:
                counts[end_index] += 1
        return tuple(counts)

    def _match_edge_endpoint(
        self,
        label: str | None,
        position: tuple[float, ...],
        label_to_indices: Mapping[str, list[int]],
        position_lookup: Mapping[int, tuple[float, ...]],
    ) -> int | None:
        if label:
            indices = label_to_indices.get(label)
            if indices:
                return indices[0]
        if not position:
            return None
        canonical_position = tuple(self._canonical_fractional(value) for value in position)
        for index, node_position in position_lookup.items():
            if self._same_position(canonical_position, node_position):
                return index
        return None

    def _same_position(self, first: tuple[float, ...], second: tuple[float, ...], tolerance: float = 1e-4) -> bool:
        if len(first) != len(second):
            return False
        return all(abs(self._canonical_fractional(a) - self._canonical_fractional(b)) <= tolerance for a, b in zip(first, second))

    def _canonical_fractional(self, value: float) -> float:
        if math.isclose(value, round(value), abs_tol=1e-6):
            value = 0.0
        wrapped = value % 1.0
        if math.isclose(wrapped, 1.0, abs_tol=1e-6):
            return 0.0
        return wrapped

    def _looks_like_supported_member(self, member_name: str) -> bool:
        path = Path(member_name)
        return path.suffix.lower() in _SUPPORTED_SUFFIXES or path.suffix == ""

    def _decode_payload(self, payload: bytes) -> str:
        for encoding in ("utf-8", "latin-1"):
            try:
                return payload.decode(encoding)
            except UnicodeDecodeError:
                continue
        return payload.decode("utf-8", errors="replace")

    def _float_indices(self, values: Sequence[str]) -> tuple[tuple[int, float], ...]:
        parsed: list[tuple[int, float]] = []
        for index, value in enumerate(values):
            try:
                parsed.append((index, self._parse_float(value)))
            except ValueError:
                continue
        return tuple(parsed)

    def _parse_float(self, value: str) -> float:
        return float(value.replace("D", "E").replace("d", "e"))

    def _extract_node_coordinates(
        self,
        values: Sequence[str],
        floats: Sequence[tuple[int, float]],
        dimension_count: int,
    ) -> tuple[tuple[float, ...], tuple[int, ...]]:
        if len(floats) < dimension_count:
            raise ValueError(f"node record does not contain {dimension_count} coordinates: {values!r}")
        if dimension_count == 2 and len(values) >= 5 and len(floats) >= 3:
            coordinate_indices = tuple(index for index, _ in floats[-3:])
            x_index, y_index = coordinate_indices[:2]
            return (float(values[x_index]), float(values[y_index])), coordinate_indices
        coordinate_indices = tuple(index for index, _ in floats[-dimension_count:])
        return tuple(float(values[index]) for index in coordinate_indices), coordinate_indices

    def _extract_edge_coordinates(
        self,
        values: Sequence[str],
        floats: Sequence[tuple[int, float]],
        dimension_count: int,
    ) -> tuple[tuple[float, ...], tuple[float, ...], tuple[int, ...]]:
        if dimension_count == 2 and len(floats) >= 6:
            coordinate_indices = tuple(index for index, _ in floats[-6:])
            numeric = [float(values[index]) for index in coordinate_indices]
            return (numeric[0], numeric[1]), (numeric[3], numeric[4]), coordinate_indices
        coordinate_indices = tuple(index for index, _ in floats[-(dimension_count * 2) :])
        numeric = [float(values[index]) for index in coordinate_indices]
        return tuple(numeric[:dimension_count]), tuple(numeric[dimension_count:]), coordinate_indices

    def _strip_comments(self, line: str) -> str:
        for marker in ("#", "!"):
            if marker in line:
                line = line.split(marker, 1)[0]
        return line

    def _as_optional_str(self, value: object) -> str | None:
        if value is None:
            return None
        return str(value)


def discover_rcsr_archives(downloads_dir: str | Path | None = None) -> dict[str, str]:
    root = Path(downloads_dir).expanduser() if downloads_dir is not None else Path.home() / "Downloads"
    discovered: dict[str, str] = {}

    topo_2d = root / "topo_2d.zip"
    if topo_2d.is_file():
        discovered["2D"] = str(topo_2d)

    preferred_3d = root / "topo_3d.zip"
    if preferred_3d.is_file():
        discovered["3D"] = str(preferred_3d)
    else:
        candidates = [path for path in root.rglob("*.zip") if _looks_like_3d_topology_archive(path)]
        if candidates:
            discovered["3D"] = str(sorted(candidates, key=_topology_archive_sort_key)[0])
    return discovered


def configured_rcsr_archives() -> dict[str, str]:
    archives: dict[str, str] = {}
    archive_2d = os.environ.get("COFKIT_TOPOLOGY_2D_ARCHIVE")
    archive_3d = os.environ.get("COFKIT_TOPOLOGY_3D_ARCHIVE")
    downloads_dir = os.environ.get("COFKIT_TOPOLOGY_DOWNLOADS_DIR")

    if archive_2d:
        archives["2D"] = str(Path(archive_2d).expanduser())
    if archive_3d:
        archives["3D"] = str(Path(archive_3d).expanduser())
    if not archives and downloads_dir:
        archives.update(discover_rcsr_archives(downloads_dir))
    return archives


def _looks_like_3d_topology_archive(path: Path) -> bool:
    name = path.name.lower()
    return "3d" in name and any(token in name for token in ("topo", "topology", "rcsr", "net"))


def _topology_archive_sort_key(path: Path) -> tuple[int, int, str]:
    name = path.name.lower()
    priority = 0 if name == "topo_3d.zip" else 1
    penalty = 0
    for token in ("topo", "topology", "rcsr", "bundle"):
        if token not in name:
            penalty += 1
    return (priority, penalty, name)
