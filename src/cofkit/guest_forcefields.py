from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from .forcefields import supported_forcefield_families


SUPPORTED_GUEST_FORCEFIELDS = supported_forcefield_families()
SUPPORTED_GUEST_PARAMETER_FAMILIES = (
    *SUPPORTED_GUEST_FORCEFIELDS,
    "genericmofs",
    "raspa_example_molecule",
)
SUPPORTED_GUEST_MODEL_TAGS = ("DREIDING", "GENERICMOFS", "RASPA")
SUPPORTED_GUEST_VDW_TREATMENTS = ("truncated", "shifted")
SUPPORTED_GUEST_MIXING_RULES = ("lorentz_berthelot",)
PACKAGED_GUEST_FORCEFIELD_SCHEMA_VERSION = 3


class GuestForceFieldMetadataError(ValueError):
    """Raised when packaged guest force-field metadata is invalid."""


@dataclass(frozen=True)
class GuestForceFieldMetadata:
    name: str
    species: str
    model_tag: str
    parameter_family: str
    parameter_source: str
    rotatable: bool
    default_widom: bool
    vdw_treatment: str
    tail_corrections: bool
    mixing_rule: str
    supports_hybrid_guest_restart: bool
    pseudo_atom_rows: tuple[str, ...] = ()
    mixing_rule_rows: tuple[str, ...] = ()
    compatibility_note: str | None = None
    provenance_note: str | None = None
    source_url: str | None = None
    source_commit: str | None = None
    license: str | None = None

    def to_dict(self) -> dict[str, object]:
        result: dict[str, object] = {
            "name": self.name,
            "species": self.species,
            "model_tag": self.model_tag,
            "parameter_family": self.parameter_family,
            "parameter_source": self.parameter_source,
            "rotatable": self.rotatable,
            "default_widom": self.default_widom,
            "vdw_treatment": self.vdw_treatment,
            "tail_corrections": self.tail_corrections,
            "mixing_rule": self.mixing_rule,
            "supports_hybrid_guest_restart": self.supports_hybrid_guest_restart,
        }
        if self.compatibility_note is not None:
            result["compatibility_note"] = self.compatibility_note
        if self.provenance_note is not None:
            result["provenance_note"] = self.provenance_note
        if self.source_url is not None:
            result["source_url"] = self.source_url
        if self.source_commit is not None:
            result["source_commit"] = self.source_commit
        if self.license is not None:
            result["license"] = self.license
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
    if (
        not isinstance(raw, dict)
        or raw.get("schema_version") != PACKAGED_GUEST_FORCEFIELD_SCHEMA_VERSION
        or not isinstance(raw.get("guests"), list)
    ):
        raise GuestForceFieldMetadataError(
            "Packaged guest force-field metadata must use schema_version "
            f"{PACKAGED_GUEST_FORCEFIELD_SCHEMA_VERSION} and contain a guests list: {path}"
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

        species = _required_nonblank_string(item, "species", context=context)
        model_tag = _required_model_tag(item, "model_tag", context=context)
        if not name.casefold().endswith(f"_{model_tag}".casefold()):
            raise GuestForceFieldMetadataError(
                f"Packaged guest name {name!r} must end with its model_tag _{model_tag}: {context}"
            )
        parameter_family = _required_parameter_family(item, "parameter_family", context=context)
        parameter_source = _required_nonblank_string(item, "parameter_source", context=context)
        rotatable = _required_bool(item, "rotatable", context=context)
        default_widom = _required_bool(item, "default_widom", context=context)
        vdw_treatment = _required_choice(
            item,
            "vdw_treatment",
            choices=SUPPORTED_GUEST_VDW_TREATMENTS,
            context=context,
        )
        tail_corrections = _required_bool(item, "tail_corrections", context=context)
        mixing_rule = _required_choice(
            item,
            "mixing_rule",
            choices=SUPPORTED_GUEST_MIXING_RULES,
            context=context,
        )
        supports_hybrid_guest_restart = _required_bool(
            item,
            "supports_hybrid_guest_restart",
            context=context,
        )
        pseudo_atom_rows = _optional_rows(item, "pseudo_atom_rows", context=context)
        mixing_rule_rows = _optional_rows(item, "mixing_rule_rows", context=context)
        if bool(pseudo_atom_rows) != bool(mixing_rule_rows):
            raise GuestForceFieldMetadataError(
                f"pseudo_atom_rows and mixing_rule_rows must either both be provided or both be empty: {context}"
            )
        note_raw = item.get("compatibility_note")
        if note_raw is not None and (not isinstance(note_raw, str) or not note_raw.strip()):
            raise GuestForceFieldMetadataError(f"compatibility_note must be a non-blank string: {context}")
        provenance_note_raw = item.get("provenance_note")
        if provenance_note_raw is not None and (
            not isinstance(provenance_note_raw, str) or not provenance_note_raw.strip()
        ):
            raise GuestForceFieldMetadataError(f"provenance_note must be a non-blank string: {context}")
        source_url_raw = _optional_nonblank_string(item, "source_url", context=context)
        source_commit_raw = _optional_nonblank_string(item, "source_commit", context=context)
        license_raw = _optional_nonblank_string(item, "license", context=context)
        entries.append(
            GuestForceFieldMetadata(
                name=name,
                species=species,
                model_tag=model_tag,
                parameter_family=parameter_family,
                parameter_source=parameter_source,
                rotatable=rotatable,
                default_widom=default_widom,
                vdw_treatment=vdw_treatment,
                tail_corrections=tail_corrections,
                mixing_rule=mixing_rule,
                supports_hybrid_guest_restart=supports_hybrid_guest_restart,
                pseudo_atom_rows=pseudo_atom_rows,
                mixing_rule_rows=mixing_rule_rows,
                compatibility_note=None if note_raw is None else note_raw.strip(),
                provenance_note=None if provenance_note_raw is None else provenance_note_raw.strip(),
                source_url=source_url_raw,
                source_commit=source_commit_raw,
                license=license_raw,
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


def _required_parameter_family(raw: dict[str, object], key: str, *, context: str) -> str:
    value = raw.get(key)
    if not isinstance(value, str):
        raise GuestForceFieldMetadataError(f"Parameter family must be a string: {context}:{key}")
    normalized = value.strip().lower().replace("-", "_")
    if normalized not in SUPPORTED_GUEST_PARAMETER_FAMILIES:
        raise GuestForceFieldMetadataError(
            f"Unsupported packaged guest parameter family {value!r} at {context}:{key}; "
            f"expected one of: {', '.join(SUPPORTED_GUEST_PARAMETER_FAMILIES)}"
        )
    return normalized


def _required_model_tag(raw: dict[str, object], key: str, *, context: str) -> str:
    value = raw.get(key)
    if not isinstance(value, str):
        raise GuestForceFieldMetadataError(f"model_tag must be a string: {context}:{key}")
    normalized = value.strip().upper().replace("-", "_")
    if normalized not in SUPPORTED_GUEST_MODEL_TAGS:
        raise GuestForceFieldMetadataError(
            f"Unsupported packaged guest model_tag {value!r} at {context}:{key}; "
            f"expected one of: {', '.join(SUPPORTED_GUEST_MODEL_TAGS)}"
        )
    return normalized


def _required_bool(raw: dict[str, object], key: str, *, context: str) -> bool:
    value = raw.get(key)
    if not isinstance(value, bool):
        raise GuestForceFieldMetadataError(f"{key} must be a boolean: {context}")
    return value


def _required_choice(
    raw: dict[str, object],
    key: str,
    *,
    choices: tuple[str, ...],
    context: str,
) -> str:
    value = raw.get(key)
    if not isinstance(value, str):
        raise GuestForceFieldMetadataError(f"{key} must be a string: {context}")
    normalized = value.strip().lower().replace("-", "_")
    if normalized not in choices:
        raise GuestForceFieldMetadataError(
            f"Unsupported {key} {value!r} at {context}; expected one of: {', '.join(choices)}"
        )
    return normalized


def _optional_nonblank_string(raw: dict[str, object], key: str, *, context: str) -> str | None:
    value = raw.get(key)
    if value is None:
        return None
    if not isinstance(value, str) or not value.strip():
        raise GuestForceFieldMetadataError(f"{key} must be a non-blank string: {context}")
    return value.strip()


def _optional_rows(raw: dict[str, object], key: str, *, context: str) -> tuple[str, ...]:
    if key not in raw:
        return ()
    value = raw[key]
    if not isinstance(value, list):
        raise GuestForceFieldMetadataError(f"{key} must be a list of row strings: {context}")
    rows: list[str] = []
    for index, row in enumerate(value):
        if not isinstance(row, str) or not row.strip():
            raise GuestForceFieldMetadataError(f"{key}[{index}] must be a non-blank string: {context}")
        rows.append(row.strip())
    return tuple(rows)
