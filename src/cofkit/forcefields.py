from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from types import MappingProxyType
from typing import Any, ClassVar


FORCEFIELD_METADATA_SCHEMA_VERSION = 1


class ForceFieldMetadataError(ValueError):
    """Raised when force-field metadata is missing, inconsistent, or stale."""


@dataclass(frozen=True)
class ForceFieldCitation:
    role: str
    title: str
    authors: str
    year: int
    doi: str | None = None
    url: str | None = None

    def to_dict(self) -> dict[str, object]:
        result: dict[str, object] = {
            "role": self.role,
            "title": self.title,
            "authors": self.authors,
            "year": self.year,
        }
        if self.doi is not None:
            result["doi"] = self.doi
        if self.url is not None:
            result["url"] = self.url
        return result


@dataclass(frozen=True)
class ForceFieldArtifact:
    role: str
    path: str
    sha256: str
    provenance: str

    def to_dict(self) -> dict[str, object]:
        return {
            "role": self.role,
            "path": self.path,
            "sha256": self.sha256,
            "provenance": self.provenance,
        }


@dataclass(frozen=True)
class ForceFieldParameterSource:
    description: str
    citations: tuple[ForceFieldCitation, ...]
    artifacts: tuple[ForceFieldArtifact, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "description": self.description,
            "citations": [citation.to_dict() for citation in self.citations],
            "artifacts": [artifact.to_dict() for artifact in self.artifacts],
        }


@dataclass(frozen=True)
class ForceFieldAtomTyping:
    backends: tuple[str, ...]
    method: str
    implementation: str
    implementation_version: str
    requires_explicit_bond_orders: bool
    notes: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "backends": list(self.backends),
            "method": self.method,
            "implementation": self.implementation,
            "implementation_version": self.implementation_version,
            "requires_explicit_bond_orders": self.requires_explicit_bond_orders,
            "notes": list(self.notes),
        }


@dataclass(frozen=True)
class ForceFieldElectrostatics:
    native_charge_model: str
    supported_charge_sources: tuple[str, ...]
    notes: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "native_charge_model": self.native_charge_model,
            "supported_charge_sources": list(self.supported_charge_sources),
            "notes": list(self.notes),
        }


@dataclass(frozen=True)
class ForceFieldNonbonded:
    potential: str
    epsilon_mixing_rule: str
    sigma_mixing_rule: str
    pair_coefficient_mode_by_backend: Mapping[str, str]
    notes: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "potential": self.potential,
            "epsilon_mixing_rule": self.epsilon_mixing_rule,
            "sigma_mixing_rule": self.sigma_mixing_rule,
            "pair_coefficient_mode_by_backend": dict(self.pair_coefficient_mode_by_backend),
            "notes": list(self.notes),
        }


@dataclass(frozen=True)
class ForceFieldIntramolecularScaling:
    backends: tuple[str, ...]
    lj_12: float
    lj_13: float
    lj_14: float
    coulomb_12: float
    coulomb_13: float
    coulomb_14: float

    def to_dict(self) -> dict[str, object]:
        return {
            "backends": list(self.backends),
            "lj_12": self.lj_12,
            "lj_13": self.lj_13,
            "lj_14": self.lj_14,
            "coulomb_12": self.coulomb_12,
            "coulomb_13": self.coulomb_13,
            "coulomb_14": self.coulomb_14,
        }


@dataclass(frozen=True)
class ForceFieldCoverage:
    parameterized_elements: tuple[str, ...]
    validated_linkages: tuple[str, ...]
    linkage_validation_status: str
    backend_elements: Mapping[str, tuple[str, ...]]
    notes: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "parameterized_elements": list(self.parameterized_elements),
            "validated_linkages": list(self.validated_linkages),
            "linkage_validation_status": self.linkage_validation_status,
            "backend_elements": {
                backend: list(elements) for backend, elements in self.backend_elements.items()
            },
            "notes": list(self.notes),
        }


@dataclass(frozen=True)
class ForceFieldLicense:
    status: str
    identifier: str | None
    source_url: str | None
    notes: str

    def to_dict(self) -> dict[str, object]:
        return {
            "status": self.status,
            "identifier": self.identifier,
            "source_url": self.source_url,
            "notes": self.notes,
        }


@dataclass(frozen=True)
class ForceFieldMetadata:
    schema_version: ClassVar[int] = FORCEFIELD_METADATA_SCHEMA_VERSION

    id: str
    family: str
    variant: str
    version: str
    display_name: str
    default_for_family: bool
    aliases: tuple[str, ...]
    scope: str
    parameter_source: ForceFieldParameterSource
    atom_typing: tuple[ForceFieldAtomTyping, ...]
    electrostatics: ForceFieldElectrostatics
    nonbonded: ForceFieldNonbonded
    intramolecular_scaling: ForceFieldIntramolecularScaling
    coverage: ForceFieldCoverage
    license: ForceFieldLicense
    notes: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": self.schema_version,
            "id": self.id,
            "family": self.family,
            "variant": self.variant,
            "version": self.version,
            "display_name": self.display_name,
            "default_for_family": self.default_for_family,
            "aliases": list(self.aliases),
            "scope": self.scope,
            "parameter_source": self.parameter_source.to_dict(),
            "atom_typing": [typing.to_dict() for typing in self.atom_typing],
            "electrostatics": self.electrostatics.to_dict(),
            "nonbonded": self.nonbonded.to_dict(),
            "intramolecular_scaling": self.intramolecular_scaling.to_dict(),
            "coverage": self.coverage.to_dict(),
            "license": self.license.to_dict(),
            "notes": list(self.notes),
        }


def packaged_forcefield_metadata_path() -> Path:
    return Path(__file__).resolve().parent / "data" / "forcefields" / "registry.json"


def load_forcefield_metadata(
    path: str | Path,
    *,
    artifact_root: str | Path | None = None,
    verify_artifacts: bool = True,
) -> tuple[ForceFieldMetadata, ...]:
    metadata_path = Path(path).expanduser().resolve()
    if not metadata_path.is_file():
        raise ForceFieldMetadataError(f"Force-field metadata file does not exist: {metadata_path}")
    try:
        raw = json.loads(metadata_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ForceFieldMetadataError(f"Invalid force-field metadata JSON: {metadata_path}") from exc
    if not isinstance(raw, Mapping) or raw.get("schema_version") != FORCEFIELD_METADATA_SCHEMA_VERSION:
        raise ForceFieldMetadataError(
            f"Force-field metadata must use schema_version {FORCEFIELD_METADATA_SCHEMA_VERSION}: "
            f"{metadata_path}"
        )
    entries_raw = raw.get("forcefields")
    if not isinstance(entries_raw, list) or not entries_raw:
        raise ForceFieldMetadataError(f"Force-field metadata must contain a non-empty forcefields list: {metadata_path}")

    entries = tuple(
        _parse_forcefield_entry(item, context=f"{metadata_path}:forcefields[{index}]")
        for index, item in enumerate(entries_raw)
    )
    _validate_registry_uniqueness(entries, metadata_path)
    if verify_artifacts:
        root = Path(artifact_root).expanduser().resolve() if artifact_root is not None else metadata_path.parent
        for entry in entries:
            verify_forcefield_artifacts(entry, artifact_root=root)
    return entries


@lru_cache(maxsize=1)
def load_packaged_forcefield_metadata() -> tuple[ForceFieldMetadata, ...]:
    return load_forcefield_metadata(
        packaged_forcefield_metadata_path(),
        artifact_root=Path(__file__).resolve().parent,
    )


def packaged_forcefield_catalog() -> dict[str, ForceFieldMetadata]:
    return {entry.id: entry for entry in load_packaged_forcefield_metadata()}


def supported_forcefield_families() -> tuple[str, ...]:
    return tuple(entry.family for entry in load_packaged_forcefield_metadata() if entry.default_for_family)


def supported_forcefield_selectors() -> tuple[str, ...]:
    selectors: list[str] = []
    for entry in load_packaged_forcefield_metadata():
        if entry.default_for_family:
            selectors.append(entry.family)
        selectors.extend((entry.id, *entry.aliases))
    return tuple(dict.fromkeys(selectors))


def resolve_forcefield_metadata(selector: str) -> ForceFieldMetadata:
    normalized = _normalize_selector(selector)
    matches: list[ForceFieldMetadata] = []
    for entry in load_packaged_forcefield_metadata():
        selectors = {_normalize_selector(entry.id), *(_normalize_selector(alias) for alias in entry.aliases)}
        if entry.default_for_family:
            selectors.add(_normalize_selector(entry.family))
        if normalized in selectors:
            matches.append(entry)
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise ForceFieldMetadataError(f"Ambiguous force-field selector {selector!r}.")
    available = ", ".join(entry.id for entry in load_packaged_forcefield_metadata())
    raise ForceFieldMetadataError(f"Unknown force-field selector {selector!r}; available parameter sets: {available}.")


def packaged_forcefield_artifact_path(selector: str, role: str) -> Path:
    metadata = resolve_forcefield_metadata(selector)
    matches = [artifact for artifact in metadata.parameter_source.artifacts if artifact.role == role]
    if len(matches) != 1:
        raise ForceFieldMetadataError(
            f"Force-field {metadata.id!r} must define exactly one artifact with role {role!r}."
        )
    path = _resolve_artifact_path(matches[0], Path(__file__).resolve().parent)
    _verify_artifact(matches[0], path, metadata.id)
    return path


def verify_forcefield_artifacts(
    metadata: ForceFieldMetadata,
    *,
    artifact_root: str | Path,
) -> dict[str, str]:
    root = Path(artifact_root).expanduser().resolve()
    verified: dict[str, str] = {}
    for artifact in metadata.parameter_source.artifacts:
        artifact_path = _resolve_artifact_path(artifact, root)
        _verify_artifact(artifact, artifact_path, metadata.id)
        verified[artifact.role] = artifact.sha256
    return verified


def _parse_forcefield_entry(raw: object, *, context: str) -> ForceFieldMetadata:
    entry = _mapping(raw, context=context)
    parameter_source_raw = _required_mapping(entry, "parameter_source", context=context)
    citations_raw = _required_sequence(parameter_source_raw, "citations", context=context, allow_empty=False)
    artifacts_raw = _required_sequence(parameter_source_raw, "artifacts", context=context, allow_empty=False)
    artifacts = tuple(_parse_artifact(item, context=f"{context}:artifact") for item in artifacts_raw)
    artifact_roles = [artifact.role for artifact in artifacts]
    if len(artifact_roles) != len(set(artifact_roles)):
        raise ForceFieldMetadataError(f"Artifact roles must be unique: {context}")
    parameter_source = ForceFieldParameterSource(
        description=_required_string(parameter_source_raw, "description", context=context),
        citations=tuple(_parse_citation(item, context=f"{context}:citation") for item in citations_raw),
        artifacts=artifacts,
    )

    typing_entries_raw = _required_sequence(entry, "atom_typing", context=context, allow_empty=False)
    atom_typing = tuple(
        _parse_atom_typing(item, context=f"{context}:atom_typing[{index}]")
        for index, item in enumerate(typing_entries_raw)
    )
    typing_backends = [backend for typing in atom_typing for backend in typing.backends]
    if len(typing_backends) != len(set(typing_backends)):
        raise ForceFieldMetadataError(f"Atom-typing backends must not overlap: {context}")
    electrostatics_raw = _required_mapping(entry, "electrostatics", context=context)
    nonbonded_raw = _required_mapping(entry, "nonbonded", context=context)
    pair_modes_raw = _required_mapping(
        nonbonded_raw, "pair_coefficient_mode_by_backend", context=context
    )
    pair_modes = MappingProxyType(
        {
            _normalize_token(
                _mapping_key(backend, context=f"{context}:pair_coefficient_mode_by_backend"),
                context=f"{context}:pair_coefficient_mode_by_backend",
            ): _normalize_token(
                _mapping_value(mode, context=f"{context}:pair_coefficient_mode_by_backend:{backend}"),
                context=f"{context}:pair_coefficient_mode_by_backend:{backend}",
            )
            for backend, mode in pair_modes_raw.items()
        }
    )
    if not pair_modes:
        raise ForceFieldMetadataError(f"pair_coefficient_mode_by_backend must not be empty: {context}")
    scaling_raw = _required_mapping(entry, "intramolecular_scaling", context=context)
    scaling_backends = _required_token_tuple(
        scaling_raw, "backends", context=context, allow_empty=False
    )
    coverage_raw = _required_mapping(entry, "coverage", context=context)
    backend_elements_raw = _required_mapping(coverage_raw, "backend_elements", context=context)
    backend_elements = MappingProxyType(
        {
            _normalize_token(
                _mapping_key(backend, context=f"{context}:backend_elements"),
                context=f"{context}:backend_elements",
            ): _string_tuple(
                elements,
                context=f"{context}:backend_elements:{backend}",
                allow_empty=False,
            )
            for backend, elements in backend_elements_raw.items()
        }
    )
    parameterized_elements = _required_string_tuple(
        coverage_raw, "parameterized_elements", context=context, allow_empty=False
    )
    parameterized_set = set(parameterized_elements)
    for backend, elements in backend_elements.items():
        missing = set(elements) - parameterized_set
        if missing:
            raise ForceFieldMetadataError(
                f"Backend {backend!r} lists elements outside parameterized_elements at {context}: "
                + ", ".join(sorted(missing))
            )
    registered_backends = set(backend_elements)
    if set(typing_backends) != registered_backends:
        raise ForceFieldMetadataError(
            f"Atom-typing backends must match coverage.backend_elements at {context}"
        )
    if set(pair_modes) != registered_backends:
        raise ForceFieldMetadataError(
            f"Nonbonded backends must match coverage.backend_elements at {context}"
        )
    if not set(scaling_backends).issubset(registered_backends):
        raise ForceFieldMetadataError(
            f"Intramolecular-scaling backends must be registered in coverage.backend_elements at {context}"
        )

    license_raw = _required_mapping(entry, "license", context=context)
    license_status = _required_string(license_raw, "status", context=context)
    if license_status not in {"verified", "unreviewed", "restricted", "unknown"}:
        raise ForceFieldMetadataError(f"Unsupported license status {license_status!r}: {context}")

    return ForceFieldMetadata(
        id=_metadata_identifier(_required_string(entry, "id", context=context), context=f"{context}:id"),
        family=_normalize_family(_required_string(entry, "family", context=context), context=f"{context}:family"),
        variant=_metadata_identifier(
            _required_string(entry, "variant", context=context), context=f"{context}:variant"
        ),
        version=_required_string(entry, "version", context=context),
        display_name=_required_string(entry, "display_name", context=context),
        default_for_family=_required_bool(entry, "default_for_family", context=context),
        aliases=_required_string_tuple(entry, "aliases", context=context, allow_empty=True),
        scope=_normalize_token(_required_string(entry, "scope", context=context), context=f"{context}:scope"),
        parameter_source=parameter_source,
        atom_typing=atom_typing,
        electrostatics=ForceFieldElectrostatics(
            native_charge_model=_normalize_token(
                _required_string(electrostatics_raw, "native_charge_model", context=context), context=context
            ),
            supported_charge_sources=_required_string_tuple(
                electrostatics_raw, "supported_charge_sources", context=context, allow_empty=False
            ),
            notes=_required_string_tuple(electrostatics_raw, "notes", context=context, allow_empty=True),
        ),
        nonbonded=ForceFieldNonbonded(
            potential=_normalize_token(_required_string(nonbonded_raw, "potential", context=context), context=context),
            epsilon_mixing_rule=_normalize_token(
                _required_string(nonbonded_raw, "epsilon_mixing_rule", context=context), context=context
            ),
            sigma_mixing_rule=_normalize_token(
                _required_string(nonbonded_raw, "sigma_mixing_rule", context=context), context=context
            ),
            pair_coefficient_mode_by_backend=pair_modes,
            notes=_required_string_tuple(nonbonded_raw, "notes", context=context, allow_empty=True),
        ),
        intramolecular_scaling=ForceFieldIntramolecularScaling(
            backends=scaling_backends,
            **{
                key: _required_scale(scaling_raw, key, context=context)
                for key in ("lj_12", "lj_13", "lj_14", "coulomb_12", "coulomb_13", "coulomb_14")
            }
        ),
        coverage=ForceFieldCoverage(
            parameterized_elements=parameterized_elements,
            validated_linkages=_required_string_tuple(
                coverage_raw, "validated_linkages", context=context, allow_empty=True
            ),
            linkage_validation_status=_normalize_token(
                _required_string(coverage_raw, "linkage_validation_status", context=context), context=context
            ),
            backend_elements=backend_elements,
            notes=_required_string_tuple(coverage_raw, "notes", context=context, allow_empty=True),
        ),
        license=ForceFieldLicense(
            status=license_status,
            identifier=_optional_string(license_raw, "identifier", context=context),
            source_url=_optional_string(license_raw, "source_url", context=context),
            notes=_required_string(license_raw, "notes", context=context),
        ),
        notes=_required_string_tuple(entry, "notes", context=context, allow_empty=True),
    )


def _parse_atom_typing(raw: object, *, context: str) -> ForceFieldAtomTyping:
    typing = _mapping(raw, context=context)
    return ForceFieldAtomTyping(
        backends=_required_token_tuple(typing, "backends", context=context, allow_empty=False),
        method=_normalize_token(_required_string(typing, "method", context=context), context=context),
        implementation=_required_string(typing, "implementation", context=context),
        implementation_version=_required_string(typing, "implementation_version", context=context),
        requires_explicit_bond_orders=_required_bool(
            typing, "requires_explicit_bond_orders", context=context
        ),
        notes=_required_string_tuple(typing, "notes", context=context, allow_empty=True),
    )


def _parse_citation(raw: object, *, context: str) -> ForceFieldCitation:
    citation = _mapping(raw, context=context)
    year = citation.get("year")
    if not isinstance(year, int) or isinstance(year, bool) or year <= 0:
        raise ForceFieldMetadataError(f"Citation year must be a positive integer: {context}")
    return ForceFieldCitation(
        role=_normalize_token(_required_string(citation, "role", context=context), context=context),
        title=_required_string(citation, "title", context=context),
        authors=_required_string(citation, "authors", context=context),
        year=year,
        doi=_optional_string(citation, "doi", context=context),
        url=_optional_string(citation, "url", context=context),
    )


def _parse_artifact(raw: object, *, context: str) -> ForceFieldArtifact:
    artifact = _mapping(raw, context=context)
    digest = _required_string(artifact, "sha256", context=context).lower()
    if len(digest) != 64 or any(character not in "0123456789abcdef" for character in digest):
        raise ForceFieldMetadataError(f"Artifact sha256 must contain 64 hexadecimal characters: {context}")
    return ForceFieldArtifact(
        role=_normalize_token(_required_string(artifact, "role", context=context), context=context),
        path=_required_string(artifact, "path", context=context),
        sha256=digest,
        provenance=_required_string(artifact, "provenance", context=context),
    )


def _validate_registry_uniqueness(entries: tuple[ForceFieldMetadata, ...], path: Path) -> None:
    ids = [_normalize_selector(entry.id) for entry in entries]
    if len(ids) != len(set(ids)):
        raise ForceFieldMetadataError(f"Force-field IDs must be unique: {path}")
    defaults: dict[str, int] = {}
    for entry in entries:
        if entry.default_for_family:
            defaults[entry.family] = defaults.get(entry.family, 0) + 1
    invalid_defaults = {family: count for family, count in defaults.items() if count != 1}
    families = {entry.family for entry in entries}
    missing_defaults = families - set(defaults)
    if invalid_defaults or missing_defaults:
        raise ForceFieldMetadataError(f"Each force-field family must define exactly one default variant: {path}")

    selectors: dict[str, str] = {}
    for entry in entries:
        entry_selectors = [entry.id, *entry.aliases]
        if entry.default_for_family:
            entry_selectors.append(entry.family)
        for selector in entry_selectors:
            key = _normalize_selector(selector)
            previous = selectors.get(key)
            if previous is not None:
                raise ForceFieldMetadataError(
                    f"Force-field selector {selector!r} is shared by {previous!r} and {entry.id!r}: {path}"
                )
            selectors[key] = entry.id


def _resolve_artifact_path(artifact: ForceFieldArtifact, root: Path) -> Path:
    root = root.resolve()
    candidate = (root / artifact.path).resolve()
    try:
        candidate.relative_to(root)
    except ValueError as exc:
        raise ForceFieldMetadataError(
            f"Force-field artifact escapes its configured root: {artifact.path!r}"
        ) from exc
    return candidate


def _verify_artifact(artifact: ForceFieldArtifact, path: Path, forcefield_id: str) -> None:
    if not path.is_file():
        raise ForceFieldMetadataError(
            f"Force-field artifact {artifact.role!r} for {forcefield_id!r} is missing: {path}"
        )
    observed = hashlib.sha256(path.read_bytes()).hexdigest()
    if observed != artifact.sha256:
        raise ForceFieldMetadataError(
            f"Force-field artifact {artifact.role!r} for {forcefield_id!r} failed SHA-256 verification: "
            f"expected {artifact.sha256}, observed {observed}."
        )


def _mapping(raw: object, *, context: str) -> Mapping[str, Any]:
    if not isinstance(raw, Mapping):
        raise ForceFieldMetadataError(f"Expected an object: {context}")
    return raw


def _mapping_key(raw: object, *, context: str) -> str:
    if not isinstance(raw, str) or not raw.strip():
        raise ForceFieldMetadataError(f"Object keys must be non-blank strings: {context}")
    return raw.strip()


def _mapping_value(raw: object, *, context: str) -> str:
    if not isinstance(raw, str) or not raw.strip():
        raise ForceFieldMetadataError(f"Object values must be non-blank strings: {context}")
    return raw.strip()


def _required_mapping(raw: Mapping[str, Any], key: str, *, context: str) -> Mapping[str, Any]:
    return _mapping(raw.get(key), context=f"{context}:{key}")


def _required_sequence(
    raw: Mapping[str, Any], key: str, *, context: str, allow_empty: bool
) -> Sequence[object]:
    value = raw.get(key)
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or (not allow_empty and not value):
        qualifier = "" if allow_empty else " non-empty"
        raise ForceFieldMetadataError(f"{key} must be a{qualifier} list: {context}")
    return value


def _required_string(raw: Mapping[str, Any], key: str, *, context: str) -> str:
    value = raw.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ForceFieldMetadataError(f"{key} must be a non-blank string: {context}")
    return value.strip()


def _optional_string(raw: Mapping[str, Any], key: str, *, context: str) -> str | None:
    value = raw.get(key)
    if value is None:
        return None
    if not isinstance(value, str) or not value.strip():
        raise ForceFieldMetadataError(f"{key} must be null or a non-blank string: {context}")
    return value.strip()


def _required_bool(raw: Mapping[str, Any], key: str, *, context: str) -> bool:
    value = raw.get(key)
    if not isinstance(value, bool):
        raise ForceFieldMetadataError(f"{key} must be a boolean: {context}")
    return value


def _required_scale(raw: Mapping[str, Any], key: str, *, context: str) -> float:
    value = raw.get(key)
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise ForceFieldMetadataError(f"{key} must be numeric: {context}")
    scale = float(value)
    if not 0.0 <= scale <= 1.0:
        raise ForceFieldMetadataError(f"{key} must be between 0 and 1: {context}")
    return scale


def _required_string_tuple(
    raw: Mapping[str, Any], key: str, *, context: str, allow_empty: bool
) -> tuple[str, ...]:
    return _string_tuple(
        _required_sequence(raw, key, context=context, allow_empty=allow_empty),
        context=f"{context}:{key}",
        allow_empty=allow_empty,
    )


def _required_token_tuple(
    raw: Mapping[str, Any], key: str, *, context: str, allow_empty: bool
) -> tuple[str, ...]:
    values = _required_string_tuple(raw, key, context=context, allow_empty=allow_empty)
    normalized = tuple(_normalize_token(value, context=f"{context}:{key}") for value in values)
    if len(normalized) != len(set(normalized)):
        raise ForceFieldMetadataError(f"Token list contains normalized duplicates: {context}:{key}")
    return normalized


def _string_tuple(raw: object, *, context: str, allow_empty: bool) -> tuple[str, ...]:
    if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes)) or (not allow_empty and not raw):
        raise ForceFieldMetadataError(f"Expected a list of strings: {context}")
    values = tuple(_required_sequence_item(value, context=context) for value in raw)
    if len(values) != len(set(values)):
        raise ForceFieldMetadataError(f"String list contains duplicates: {context}")
    return values


def _required_sequence_item(value: object, *, context: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ForceFieldMetadataError(f"List entries must be non-blank strings: {context}")
    return value.strip()


def _normalize_token(value: str, *, context: str) -> str:
    normalized = value.strip().lower().replace("-", "_")
    if not normalized or any(character.isspace() for character in normalized):
        raise ForceFieldMetadataError(f"Invalid metadata token {value!r}: {context}")
    return normalized


def _metadata_identifier(value: str, *, context: str) -> str:
    identifier = value.strip()
    allowed = frozenset("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_.")
    if not identifier or any(character not in allowed for character in identifier):
        raise ForceFieldMetadataError(f"Invalid metadata identifier {value!r}: {context}")
    return identifier


def _normalize_family(value: str, *, context: str) -> str:
    return _normalize_token(value, context=context)


def _normalize_selector(value: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ForceFieldMetadataError("Force-field selector must be a non-blank string.")
    return value.strip().lower().replace("-", "_")
