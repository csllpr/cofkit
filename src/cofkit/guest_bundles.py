from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .guest_forcefields import SUPPORTED_GUEST_FORCEFIELDS


class GuestBundleError(ValueError):
    """Raised when a parameterized guest bundle cannot be loaded."""


@dataclass(frozen=True)
class RaspaGuestBundle:
    molecule_definition_text: str
    pseudo_atom_rows: tuple[str, ...]
    mixing_rule_rows: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "pseudo_atom_rows": list(self.pseudo_atom_rows),
            "mixing_rule_rows": list(self.mixing_rule_rows),
            "molecule_definition_lines": len(self.molecule_definition_text.splitlines()),
        }


@dataclass(frozen=True)
class GuestBundle:
    name: str
    source_path: str
    aliases: tuple[str, ...]
    rotatable: bool
    parameter_family: str
    parameter_source: str
    compatible_framework_forcefields: tuple[str, ...]
    raspa: RaspaGuestBundle
    lammps: Mapping[str, Any]
    notes: tuple[str, ...] = ()

    def names(self) -> tuple[str, ...]:
        return (self.name, *self.aliases)

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "source_path": self.source_path,
            "aliases": list(self.aliases),
            "rotatable": self.rotatable,
            "parameter_family": self.parameter_family,
            "parameter_source": self.parameter_source,
            "compatible_framework_forcefields": list(self.compatible_framework_forcefields),
            "raspa": self.raspa.to_dict(),
            "lammps_keys": sorted(str(key) for key in self.lammps),
            "notes": list(self.notes),
        }


def load_guest_bundle(path: str | Path) -> GuestBundle:
    bundle_path = Path(path).expanduser().resolve()
    if not bundle_path.is_file():
        raise GuestBundleError(f"Guest bundle file does not exist: {bundle_path}")
    try:
        raw = json.loads(bundle_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise GuestBundleError(f"Guest bundle is not valid JSON: {bundle_path}") from exc
    if not isinstance(raw, Mapping):
        raise GuestBundleError(f"Guest bundle root must be a JSON object: {bundle_path}")

    version = raw.get("version", 1)
    if version != 2:
        raise GuestBundleError(
            f"Unsupported guest bundle version {version!r} in {bundle_path}; expected 2. "
            "Version 2 requires force-field provenance and compatibility metadata."
        )

    name = _required_token(raw, "name", context=str(bundle_path))
    aliases = tuple(_normalize_token(value, context=f"{bundle_path}:aliases") for value in raw.get("aliases", ()))
    if len(set(alias.casefold() for alias in aliases)) != len(aliases):
        raise GuestBundleError(f"Guest bundle aliases must be unique case-insensitively: {bundle_path}")
    if name.casefold() in {alias.casefold() for alias in aliases}:
        raise GuestBundleError(f"Guest bundle alias list repeats the canonical name {name!r}: {bundle_path}")

    rotatable = raw.get("rotatable", True)
    if not isinstance(rotatable, bool):
        raise GuestBundleError(f"Guest bundle rotatable must be a boolean: {bundle_path}")

    parameter_family = _required_forcefield(raw, "parameter_family", context=str(bundle_path))
    parameter_source = _required_nonblank_string(raw, "parameter_source", context=str(bundle_path))
    compatible_raw = raw.get("compatible_framework_forcefields")
    if not isinstance(compatible_raw, Sequence) or isinstance(compatible_raw, (str, bytes)) or not compatible_raw:
        raise GuestBundleError(
            f"Guest bundle compatible_framework_forcefields must be a non-empty list: {bundle_path}"
        )
    compatible_framework_forcefields = tuple(
        _normalize_forcefield(value, context=f"{bundle_path}:compatible_framework_forcefields")
        for value in compatible_raw
    )
    if len(set(compatible_framework_forcefields)) != len(compatible_framework_forcefields):
        raise GuestBundleError(
            f"Guest bundle compatible_framework_forcefields must not contain duplicates: {bundle_path}"
        )
    if parameter_family not in compatible_framework_forcefields:
        raise GuestBundleError(
            f"Guest bundle parameter family {parameter_family!r} must be listed in "
            f"compatible_framework_forcefields: {bundle_path}"
        )

    raspa_raw = raw.get("raspa")
    if not isinstance(raspa_raw, Mapping):
        raise GuestBundleError(f"Guest bundle {name!r} must include a raspa object: {bundle_path}")
    raspa = RaspaGuestBundle(
        molecule_definition_text=_load_inline_or_path_text(
            raspa_raw,
            inline_key="molecule_definition",
            path_key="molecule_definition_path",
            base_dir=bundle_path.parent,
            context=f"{bundle_path}:raspa",
        ),
        pseudo_atom_rows=_rows_from_mapping(
            raspa_raw,
            rows_key="pseudo_atom_rows",
            path_key="pseudo_atoms_path",
            base_dir=bundle_path.parent,
            context=f"{bundle_path}:raspa",
        ),
        mixing_rule_rows=_rows_from_mapping(
            raspa_raw,
            rows_key="mixing_rule_rows",
            path_key="mixing_rules_path",
            base_dir=bundle_path.parent,
            context=f"{bundle_path}:raspa",
        ),
    )
    if not raspa.molecule_definition_text.strip():
        raise GuestBundleError(f"Guest bundle {name!r} has an empty RASPA molecule definition: {bundle_path}")
    if not raspa.pseudo_atom_rows:
        raise GuestBundleError(f"Guest bundle {name!r} must provide at least one RASPA pseudo atom row: {bundle_path}")
    if not raspa.mixing_rule_rows:
        raise GuestBundleError(f"Guest bundle {name!r} must provide at least one RASPA mixing-rule row: {bundle_path}")

    lammps_raw = raw.get("lammps")
    if not isinstance(lammps_raw, Mapping) or not lammps_raw:
        raise GuestBundleError(
            f"Guest bundle {name!r} must include a non-empty lammps object so hybrid runs have one "
            "synchronized force-field source."
        )

    notes_raw = raw.get("notes", ())
    if isinstance(notes_raw, str):
        notes = (notes_raw,)
    elif isinstance(notes_raw, Sequence):
        notes = tuple(str(value) for value in notes_raw)
    else:
        raise GuestBundleError(f"Guest bundle notes must be a string or list of strings: {bundle_path}")

    return GuestBundle(
        name=name,
        source_path=str(bundle_path),
        aliases=aliases,
        rotatable=rotatable,
        parameter_family=parameter_family,
        parameter_source=parameter_source,
        compatible_framework_forcefields=compatible_framework_forcefields,
        raspa=raspa,
        lammps=dict(lammps_raw),
        notes=notes,
    )


def load_guest_bundles(paths: Sequence[str | Path]) -> tuple[GuestBundle, ...]:
    bundles = tuple(load_guest_bundle(path) for path in paths)
    seen: dict[str, str] = {}
    for bundle in bundles:
        for name in bundle.names():
            key = name.casefold()
            previous = seen.get(key)
            if previous is not None:
                raise GuestBundleError(
                    f"Guest bundle component name or alias {name!r} is defined by both {previous} and {bundle.source_path}."
                )
            seen[key] = bundle.source_path
    return bundles


def _required_token(raw: Mapping[str, Any], key: str, *, context: str) -> str:
    if key not in raw:
        raise GuestBundleError(f"Missing required guest bundle field {key!r}: {context}")
    return _normalize_token(raw[key], context=f"{context}:{key}")


def _normalize_token(value: Any, *, context: str) -> str:
    if not isinstance(value, str):
        raise GuestBundleError(f"Expected a string token at {context}.")
    token = value.strip()
    if not token:
        raise GuestBundleError(f"String token must not be blank at {context}.")
    if any(character.isspace() for character in token):
        raise GuestBundleError(f"String token must not contain whitespace at {context}: {value!r}")
    return token


def _required_forcefield(raw: Mapping[str, Any], key: str, *, context: str) -> str:
    if key not in raw:
        raise GuestBundleError(f"Missing required guest bundle field {key!r}: {context}")
    return _normalize_forcefield(raw[key], context=f"{context}:{key}")


def _required_nonblank_string(raw: Mapping[str, Any], key: str, *, context: str) -> str:
    value = raw.get(key)
    if not isinstance(value, str) or not value.strip():
        raise GuestBundleError(f"Guest bundle {key} must be a non-blank string: {context}")
    return value.strip()


def _normalize_forcefield(value: Any, *, context: str) -> str:
    token = _normalize_token(value, context=context).lower().replace("-", "_")
    if token not in SUPPORTED_GUEST_FORCEFIELDS:
        raise GuestBundleError(
            f"Unsupported guest bundle force field {value!r} at {context}; "
            f"expected one of: {', '.join(SUPPORTED_GUEST_FORCEFIELDS)}"
        )
    return token


def _load_inline_or_path_text(
    raw: Mapping[str, Any],
    *,
    inline_key: str,
    path_key: str,
    base_dir: Path,
    context: str,
) -> str:
    inline = raw.get(inline_key)
    path_value = raw.get(path_key)
    if inline is not None and path_value is not None:
        raise GuestBundleError(f"Provide only one of {inline_key!r} or {path_key!r} at {context}.")
    if inline is not None:
        if not isinstance(inline, str):
            raise GuestBundleError(f"{inline_key!r} must be a string at {context}.")
        return inline
    if path_value is not None:
        relative_path = _normalize_path_token(path_value, context=f"{context}:{path_key}")
        path = (base_dir / relative_path).resolve()
        if not path.is_file():
            raise GuestBundleError(f"Referenced guest bundle file does not exist: {path}")
        return path.read_text(encoding="utf-8")
    raise GuestBundleError(f"Missing required {inline_key!r} or {path_key!r} at {context}.")


def _rows_from_mapping(
    raw: Mapping[str, Any],
    *,
    rows_key: str,
    path_key: str,
    base_dir: Path,
    context: str,
) -> tuple[str, ...]:
    rows_value = raw.get(rows_key)
    path_value = raw.get(path_key)
    if rows_value is not None and path_value is not None:
        raise GuestBundleError(f"Provide only one of {rows_key!r} or {path_key!r} at {context}.")
    if rows_value is not None:
        if not isinstance(rows_value, Sequence) or isinstance(rows_value, (str, bytes)):
            raise GuestBundleError(f"{rows_key!r} must be a list of row strings at {context}.")
        return _normalize_rows(rows_value, context=f"{context}:{rows_key}")
    if path_value is not None:
        relative_path = _normalize_path_token(path_value, context=f"{context}:{path_key}")
        path = (base_dir / relative_path).resolve()
        if not path.is_file():
            raise GuestBundleError(f"Referenced guest bundle file does not exist: {path}")
        return _extract_data_rows(path.read_text(encoding="utf-8"), context=str(path))
    return ()


def _normalize_path_token(value: Any, *, context: str) -> Path:
    if not isinstance(value, str):
        raise GuestBundleError(f"Path value must be a string at {context}.")
    token = value.strip()
    if not token:
        raise GuestBundleError(f"Path value must not be blank at {context}.")
    return Path(token)


def _normalize_rows(rows: Sequence[Any], *, context: str) -> tuple[str, ...]:
    normalized: list[str] = []
    for index, row in enumerate(rows):
        if not isinstance(row, str):
            raise GuestBundleError(f"Row {index} must be a string at {context}.")
        stripped = row.strip()
        if not stripped or stripped.startswith("#"):
            continue
        normalized.append(stripped)
    return tuple(normalized)


def _extract_data_rows(text: str, *, context: str) -> tuple[str, ...]:
    rows: list[str] = []
    for raw_line in text.splitlines():
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if _is_integer_token(stripped):
            continue
        rows.append(stripped)
    if not rows:
        raise GuestBundleError(f"No data rows found in guest bundle row file: {context}")
    return tuple(rows)


def _is_integer_token(value: str) -> bool:
    try:
        int(value)
    except ValueError:
        return False
    return True
