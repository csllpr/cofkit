from __future__ import annotations

from collections import Counter
from collections.abc import Mapping as ABCMapping
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping

from .model import Candidate, MonomerSpec
from .reactions import ReactionLibrary
from .topologies import get_topology_hint

try:
    from rdkit import Chem
except ImportError:  # pragma: no cover - environment-dependent
    Chem = None


@dataclass(frozen=True)
class COFidMonomer:
    connectivity: int
    reactive_group: str
    canonical_smiles: str


@dataclass(frozen=True)
class ParsedCOFid:
    monomers: tuple[COFidMonomer, ...]
    topology: str
    linkage: str


@dataclass(frozen=True)
class COFidBuildMonomer:
    connectivity: int
    reactive_group: str
    motif_kind: str
    canonical_smiles: str


@dataclass(frozen=True)
class COFidBuildRequest:
    cofid: str
    monomers: tuple[COFidBuildMonomer, ...]
    topology: str
    target_dimensionality: str
    linkage_code: str
    template_id: str


@dataclass(frozen=True)
class COFidComment:
    cofid: str
    suffix: str | None = None


_LINKAGE_CODE_TO_TEMPLATE_ID: Mapping[str, str] = {
    "azine": "azine_bridge",
    "azine_bridge": "azine_bridge",
    "beta-ketoenamine": "keto_enamine_bridge",
    "beta_ketoenamine": "keto_enamine_bridge",
    "bken": "keto_enamine_bridge",
    "boest": "boronate_ester_bridge",
    "boronate_ester": "boronate_ester_bridge",
    "boronate_ester_bridge": "boronate_ester_bridge",
    "hydrazone": "hydrazone_bridge",
    "hydrazone_bridge": "hydrazone_bridge",
    "imine": "imine_bridge",
    "imine_bridge": "imine_bridge",
    "keto_enamine_bridge": "keto_enamine_bridge",
    "vinylene": "vinylene_bridge",
    "vinylene_bridge": "vinylene_bridge",
}

_TEMPLATE_ID_TO_LINKAGE_CODE: Mapping[str, str] = {
    "azine_bridge": "azine",
    "boronate_ester_bridge": "boest",
    "hydrazone_bridge": "hydrazone",
    "imine_bridge": "imine",
    "keto_enamine_bridge": "bken",
    "vinylene_bridge": "vinylene",
}

COFID_COMMENT_PREFIX = "# COFid: "


def canonicalize_smiles(smiles: str) -> str:
    if Chem is None:
        raise RuntimeError("RDKit is required for COFid parsing and generation.")
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"invalid SMILES in COFid: {smiles!r}")
    return str(Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False))


def parse_cofid(cofid: str) -> ParsedCOFid:
    if not cofid:
        raise ValueError("COFid must not be blank.")
    if any(character.isspace() for character in cofid):
        raise ValueError("COFid must not contain whitespace.")

    blocks = cofid.split("&&")
    if len(blocks) != 3:
        raise ValueError("COFid must have exactly 3 blocks separated by &&.")

    monomer_block, topology, linkage = blocks
    if not topology or topology != topology.lower():
        raise ValueError("COFid topology must be a non-empty lowercase token.")
    if not linkage:
        raise ValueError("COFid linkage must be a non-empty token.")
    if any(separator in topology for separator in (":", ".", "&&")):
        raise ValueError("COFid topology contains a reserved separator.")
    if any(separator in linkage for separator in (":", ".", "&&")):
        raise ValueError("COFid linkage contains a reserved separator.")

    monomers: list[COFidMonomer] = []
    for monomer_text in monomer_block.split("."):
        parts = monomer_text.split(":", 2)
        if len(parts) != 3:
            raise ValueError(f"invalid COFid monomer block: {monomer_text!r}")
        connectivity_text, reactive_group, smiles = parts
        try:
            connectivity = int(connectivity_text)
        except ValueError as exc:
            raise ValueError(f"invalid COFid connectivity: {connectivity_text!r}") from exc
        if connectivity <= 0:
            raise ValueError(f"COFid connectivity must be positive: {connectivity!r}")
        if not reactive_group:
            raise ValueError("COFid reactive group must not be blank.")
        if any(separator in reactive_group for separator in (":", ".", "&&")):
            raise ValueError(f"COFid reactive group contains a reserved separator: {reactive_group!r}")
        canonical_smiles = canonicalize_smiles(smiles)
        if canonical_smiles != smiles:
            raise ValueError(
                "COFid SMILES must already be RDKit-canonical. "
                f"Expected {canonical_smiles!r}, received {smiles!r}."
            )
        monomers.append(
            COFidMonomer(
                connectivity=connectivity,
                reactive_group=reactive_group,
                canonical_smiles=canonical_smiles,
            )
        )

    expected_order = tuple(
        sorted(
            monomers,
            key=lambda monomer: (-monomer.connectivity, monomer.canonical_smiles),
        )
    )
    if tuple(monomers) != expected_order:
        raise ValueError("COFid monomers must be sorted by connectivity descending, then canonical SMILES.")

    return ParsedCOFid(
        monomers=tuple(monomers),
        topology=topology,
        linkage=linkage,
    )


def cofid_to_build_request(
    cofid: str,
    *,
    reaction_library: ReactionLibrary | None = None,
) -> COFidBuildRequest:
    parsed = parse_cofid(cofid)
    if len(parsed.monomers) != 2:
        raise ValueError("Current cofkit COFid input supports exactly 2 monomers.")

    template_id = _template_id_for_linkage(parsed.linkage)
    reaction_library = reaction_library or ReactionLibrary.builtin()
    template = reaction_library.get(template_id)
    profile = reaction_library.linkage_profile(template_id)
    if profile is None or len(profile.binary_bridge_roles) != 2:
        raise ValueError(f"COFid linkage {parsed.linkage!r} does not map to a supported binary-bridge build workflow.")

    build_monomers = tuple(
        COFidBuildMonomer(
            connectivity=monomer.connectivity,
            reactive_group=monomer.reactive_group,
            motif_kind=_input_motif_kind(monomer.reactive_group, template_id),
            canonical_smiles=monomer.canonical_smiles,
        )
        for monomer in parsed.monomers
    )

    expected_motif_kinds = Counter(role.motif_kind for role in profile.binary_bridge_roles)
    actual_motif_kinds = Counter(monomer.motif_kind for monomer in build_monomers)
    if actual_motif_kinds != expected_motif_kinds:
        raise ValueError(
            f"COFid monomer groups {tuple(monomer.reactive_group for monomer in parsed.monomers)!r} "
            f"do not match linkage {parsed.linkage!r} / template {template_id!r}."
        )

    topology_hint = get_topology_hint(parsed.topology)
    if not template.allows_dimensionality(topology_hint.dimensionality):
        raise ValueError(
            f"COFid linkage {parsed.linkage!r} / template {template_id!r} does not support "
            f"topology dimensionality {topology_hint.dimensionality!r}."
        )

    return COFidBuildRequest(
        cofid=cofid,
        monomers=build_monomers,
        topology=parsed.topology,
        target_dimensionality=topology_hint.dimensionality,
        linkage_code=parsed.linkage,
        template_id=template_id,
    )


def generate_cofid(
    candidate: Candidate,
    monomer_specs: Mapping[str, MonomerSpec] | Iterable[MonomerSpec],
) -> str:
    topology = _candidate_topology(candidate)
    template_id = _candidate_template_id(candidate)
    spec_map = _coerce_spec_map(monomer_specs)
    participating_monomer_ids = _participating_monomer_ids(candidate, spec_map)

    monomers: list[COFidMonomer] = []
    for monomer_id in participating_monomer_ids:
        spec = spec_map[monomer_id]
        motif_kind = _single_motif_kind(spec)
        source_smiles = _source_smiles(spec)
        monomers.append(
            COFidMonomer(
                connectivity=len(spec.motifs),
                reactive_group=_output_reactive_group(motif_kind),
                canonical_smiles=canonicalize_smiles(source_smiles),
            )
        )

    sorted_monomers = tuple(
        sorted(
            monomers,
            key=lambda monomer: (-monomer.connectivity, monomer.canonical_smiles),
        )
    )
    return serialize_cofid(
        monomers=sorted_monomers,
        topology=topology,
        linkage=_preferred_linkage_code(template_id),
    )


def try_generate_cofid(
    candidate: Candidate,
    monomer_specs: Mapping[str, MonomerSpec] | Iterable[MonomerSpec],
) -> str | None:
    try:
        return generate_cofid(candidate, monomer_specs)
    except (KeyError, RuntimeError, ValueError):
        return None


def serialize_cofid(
    *,
    monomers: Iterable[COFidMonomer],
    topology: str,
    linkage: str,
) -> str:
    monomer_block = ".".join(
        f"{monomer.connectivity}:{monomer.reactive_group}:{monomer.canonical_smiles}"
        for monomer in monomers
    )
    return f"{monomer_block}&&{topology}&&{linkage}"


def cofid_comment_line(cofid: str, *, suffix: str | None = None) -> str:
    if not cofid:
        raise ValueError("COFid comment value must not be blank.")
    if any(character in cofid for character in ("\r", "\n")):
        raise ValueError("COFid comment value must be a single line.")
    if suffix is None:
        return f"{COFID_COMMENT_PREFIX}{cofid}"
    normalized_suffix = str(suffix).strip()
    if not normalized_suffix:
        return f"{COFID_COMMENT_PREFIX}{cofid}"
    if any(character in normalized_suffix for character in ("\r", "\n")):
        raise ValueError("COFid comment suffix must be a single line.")
    return f"{COFID_COMMENT_PREFIX}{cofid} {normalized_suffix}"


def read_cofid_comment_from_cif(cif_path: str | Path) -> COFidComment | None:
    input_path = Path(cif_path)
    try:
        with input_path.open("r", encoding="utf-8") as handle:
            first_line = handle.readline()
    except FileNotFoundError:
        raise
    return _cofid_comment_from_line(first_line)


def read_cofid_from_cif(cif_path: str | Path) -> str | None:
    comment = read_cofid_comment_from_cif(cif_path)
    return None if comment is None else comment.cofid


def ensure_cif_has_cofid_comment(
    cif_path: str | Path,
    cofid: str | None,
    *,
    suffix: str | None = None,
) -> None:
    if not cofid:
        return
    path = Path(cif_path)
    expected_line = cofid_comment_line(cofid, suffix=suffix)
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()
    if lines and lines[0] == expected_line:
        return
    if lines and lines[0].startswith(COFID_COMMENT_PREFIX):
        lines[0] = expected_line
    else:
        lines.insert(0, expected_line)
    updated = "\n".join(lines)
    if text.endswith("\n"):
        updated += "\n"
    path.write_text(updated, encoding="utf-8")


def _cofid_comment_from_line(line: str) -> COFidComment | None:
    if not line.startswith(COFID_COMMENT_PREFIX):
        return None
    payload = line[len(COFID_COMMENT_PREFIX) :].rstrip("\r\n").strip()
    if not payload:
        return None
    tokens = payload.split(maxsplit=1)
    cofid = tokens[0]
    suffix = tokens[1].strip() if len(tokens) > 1 else None
    return COFidComment(cofid=cofid, suffix=suffix or None)


def _template_id_for_linkage(linkage: str) -> str:
    try:
        return _LINKAGE_CODE_TO_TEMPLATE_ID[linkage]
    except KeyError as exc:
        raise ValueError(
            f"unsupported COFid linkage code {linkage!r}; current cofkit COFid input only supports "
            f"{tuple(sorted(_LINKAGE_CODE_TO_TEMPLATE_ID))!r}"
        ) from exc


def _preferred_linkage_code(template_id: str) -> str:
    return _TEMPLATE_ID_TO_LINKAGE_CODE.get(template_id, template_id)


def _input_motif_kind(reactive_group: str, template_id: str) -> str:
    del template_id
    return reactive_group


def _output_reactive_group(motif_kind: str) -> str:
    return motif_kind


def _candidate_topology(candidate: Candidate) -> str:
    net_plan = candidate.metadata.get("net_plan")
    if isinstance(net_plan, ABCMapping):
        topology = net_plan.get("topology")
        if isinstance(topology, str) and topology:
            return topology
    raise ValueError("candidate metadata does not expose a topology id for COFid generation.")


def _candidate_template_id(candidate: Candidate) -> str:
    template_ids = {str(event.template_id) for event in candidate.events}
    if len(template_ids) != 1:
        raise ValueError(f"COFid generation requires exactly one reaction template, got {tuple(sorted(template_ids))!r}")
    return next(iter(template_ids))


def _coerce_spec_map(monomer_specs: Mapping[str, MonomerSpec] | Iterable[MonomerSpec]) -> dict[str, MonomerSpec]:
    if isinstance(monomer_specs, ABCMapping):
        return {str(key): value for key, value in monomer_specs.items()}
    return {spec.id: spec for spec in monomer_specs}


def _participating_monomer_ids(candidate: Candidate, spec_map: Mapping[str, MonomerSpec]) -> tuple[str, ...]:
    monomer_ids = {
        str(participant.monomer_id)
        for event in candidate.events
        for participant in event.participants
    }
    if monomer_ids:
        return tuple(sorted(monomer_ids))
    instance_to_monomer = candidate.metadata.get("instance_to_monomer")
    if isinstance(instance_to_monomer, ABCMapping):
        return tuple(sorted({str(monomer_id) for monomer_id in instance_to_monomer.values()}))
    return tuple(sorted(spec_map))


def _single_motif_kind(spec: MonomerSpec) -> str:
    motif_kinds = {motif.kind for motif in spec.motifs}
    if len(motif_kinds) != 1:
        raise ValueError(
            f"COFid generation requires one motif kind per monomer, got {tuple(sorted(motif_kinds))!r} for {spec.id!r}."
        )
    return next(iter(motif_kinds))


def _source_smiles(spec: MonomerSpec) -> str:
    source_smiles = spec.metadata.get("source_smiles")
    if isinstance(source_smiles, str) and source_smiles:
        return source_smiles
    canonical_smiles = spec.metadata.get("canonical_smiles")
    if isinstance(canonical_smiles, str) and canonical_smiles:
        return canonical_smiles
    raise ValueError(f"monomer {spec.id!r} does not expose source_smiles metadata for COFid generation.")
