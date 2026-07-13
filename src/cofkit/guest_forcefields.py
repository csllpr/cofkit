from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path


SUPPORTED_GUEST_FORCEFIELDS = ("dreiding", "uff")


class GuestForceFieldMetadataError(ValueError):
    """Raised when packaged guest force-field metadata is invalid."""


@dataclass(frozen=True)
class GuestForceFieldMetadata:
    name: str
    parameter_family: str
    parameter_source: str
    compatible_framework_forcefields: tuple[str, ...]
    compatibility_note: str | None = None
    provenance_note: str | None = None

    def supports(self, forcefield: str) -> bool:
        return forcefield.strip().lower().replace("-", "_") in self.compatible_framework_forcefields

    def to_dict(self) -> dict[str, object]:
        result: dict[str, object] = {
            "name": self.name,
            "parameter_family": self.parameter_family,
            "parameter_source": self.parameter_source,
            "compatible_framework_forcefields": list(self.compatible_framework_forcefields),
        }
        if self.compatibility_note is not None:
            result["compatibility_note"] = self.compatibility_note
        if self.provenance_note is not None:
            result["provenance_note"] = self.provenance_note
        return result


def packaged_guest_forcefield_metadata_path() -> Path:
    return Path(__file__).resolve().parent / "data" / "graspa" / "guest_forcefields.json"


@lru_cache(maxsize=1)
def load_packaged_guest_forcefield_metadata() -> tuple[GuestForceFieldMetadata, ...]:
    path = packaged_guest_forcefield_metadata_path()
    if not path.is_file():
        raise GuestForceFieldMetadataError(f"Packaged guest force-field metadata is missing: {path}")
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise GuestForceFieldMetadataError(f"Invalid packaged guest force-field metadata: {path}") from exc
    if not isinstance(raw, dict) or raw.get("schema_version") != 1 or not isinstance(raw.get("guests"), list):
        raise GuestForceFieldMetadataError(
            f"Packaged guest force-field metadata must use schema_version 1 and contain a guests list: {path}"
        )

    entries: list[GuestForceFieldMetadata] = []
    seen_names: set[str] = set()
    for index, item in enumerate(raw["guests"]):
        context = f"{path}:guests[{index}]"
        if not isinstance(item, dict):
            raise GuestForceFieldMetadataError(f"Guest metadata entry must be an object: {context}")
        name = _required_nonblank_string(item, "name", context=context)
        name_key = name.casefold()
        if name_key in seen_names:
            raise GuestForceFieldMetadataError(f"Duplicate packaged guest metadata name {name!r}: {path}")
        seen_names.add(name_key)

        parameter_family = _required_forcefield(item, "parameter_family", context=context)
        parameter_source = _required_nonblank_string(item, "parameter_source", context=context)
        compatible_raw = item.get("compatible_framework_forcefields")
        if not isinstance(compatible_raw, list) or not compatible_raw:
            raise GuestForceFieldMetadataError(
                f"compatible_framework_forcefields must be a non-empty list: {context}"
            )
        compatible_framework_forcefields = tuple(
            _normalize_forcefield(value, context=f"{context}:compatible_framework_forcefields")
            for value in compatible_raw
        )
        if len(set(compatible_framework_forcefields)) != len(compatible_framework_forcefields):
            raise GuestForceFieldMetadataError(f"compatible_framework_forcefields contains duplicates: {context}")
        if parameter_family not in compatible_framework_forcefields:
            raise GuestForceFieldMetadataError(
                f"Parameter family {parameter_family!r} must be listed as compatible: {context}"
            )
        note_raw = item.get("compatibility_note")
        if note_raw is not None and (not isinstance(note_raw, str) or not note_raw.strip()):
            raise GuestForceFieldMetadataError(f"compatibility_note must be a non-blank string: {context}")
        provenance_note_raw = item.get("provenance_note")
        if provenance_note_raw is not None and (
            not isinstance(provenance_note_raw, str) or not provenance_note_raw.strip()
        ):
            raise GuestForceFieldMetadataError(f"provenance_note must be a non-blank string: {context}")
        entries.append(
            GuestForceFieldMetadata(
                name=name,
                parameter_family=parameter_family,
                parameter_source=parameter_source,
                compatible_framework_forcefields=compatible_framework_forcefields,
                compatibility_note=None if note_raw is None else note_raw.strip(),
                provenance_note=None if provenance_note_raw is None else provenance_note_raw.strip(),
            )
        )
    definition_dir = path.parent / "widom_template"
    definition_names = {
        definition_path.stem
        for definition_path in definition_dir.glob("*.def")
        if definition_path.name not in {"force_field.def", "force_field_mixing_rules.def", "pseudo_atoms.def"}
    }
    metadata_names = {entry.name for entry in entries}
    if metadata_names != definition_names:
        missing_metadata = sorted(definition_names - metadata_names)
        missing_definitions = sorted(metadata_names - definition_names)
        details: list[str] = []
        if missing_metadata:
            details.append(f"definitions without metadata: {', '.join(missing_metadata)}")
        if missing_definitions:
            details.append(f"metadata without definitions: {', '.join(missing_definitions)}")
        raise GuestForceFieldMetadataError(
            f"Packaged guest metadata and molecule definitions do not match in {definition_dir}: "
            + "; ".join(details)
        )
    return tuple(entries)


def packaged_guest_forcefield_catalog() -> dict[str, GuestForceFieldMetadata]:
    return {entry.name: entry for entry in load_packaged_guest_forcefield_metadata()}


def _required_nonblank_string(raw: dict[str, object], key: str, *, context: str) -> str:
    value = raw.get(key)
    if not isinstance(value, str) or not value.strip():
        raise GuestForceFieldMetadataError(f"{key} must be a non-blank string: {context}")
    return value.strip()


def _required_forcefield(raw: dict[str, object], key: str, *, context: str) -> str:
    return _normalize_forcefield(raw.get(key), context=f"{context}:{key}")


def _normalize_forcefield(value: object, *, context: str) -> str:
    if not isinstance(value, str):
        raise GuestForceFieldMetadataError(f"Force field must be a string: {context}")
    normalized = value.strip().lower().replace("-", "_")
    if normalized not in SUPPORTED_GUEST_FORCEFIELDS:
        raise GuestForceFieldMetadataError(
            f"Unsupported packaged guest force field {value!r} at {context}; "
            f"expected one of: {', '.join(SUPPORTED_GUEST_FORCEFIELDS)}"
        )
    return normalized
