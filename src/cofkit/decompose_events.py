"""Experimental event-based CIF decomposition.

This module deliberately lives beside :mod:`cofkit.decompose` instead of
replacing it.  It shares the established CIF graph construction, precursor
validation, topology matching, and COFid serialization primitives, but changes
the decision model:

1. detect immutable local linkage events;
2. resolve only chemically local overlaps;
3. enumerate bounded reconstruction hypotheses;
4. validate each complete hypothesis globally;
5. select a result only after chemical and topological validation.

The public entry point is reached with
``decompose_cif_to_cofid(..., decomposition_mode="event")``.  Legacy mode
remains the default until the two implementations have been benchmarked on an
independent labelled corpus.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field, replace
from itertools import islice, product
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from . import decompose as legacy
from .cofid import cofid_to_build_request, parse_cofid, serialize_cofid
from .decompose_cif import PeriodicCifAtoms, read_periodic_cif_atoms

try:
    from rdkit import Chem
except ImportError:  # pragma: no cover - project dependency guard
    Chem = None


EVENT_STATUS_COMPLETE = "SUCCESS_COMPLETE"
EVENT_STATUS_AMBIGUOUS = "AMBIGUOUS_MULTIPLE_DECOMPOSITIONS"
EVENT_STATUS_CHEMICAL = "FAILED_CHEMICAL_VALIDATION"
EVENT_STATUS_ENDPOINT = "FAILED_ENDPOINT_ACCOUNTING"
EVENT_STATUS_TOPOLOGY = "FAILED_TOPOLOGY_VALIDATION"
EVENT_STATUS_UNEXPLAINED = "FAILED_UNEXPLAINED_FRAMEWORK"
EVENT_STATUS_UNSUPPORTED = "UNSUPPORTED_LINKAGE"
EVENT_STATUS_SUPPRESSED = "SUPPRESSED_LOCAL_OVERLAP"
EVENT_STATUS_TRIAZINE_MOTIF = "SUPPRESSED_TRIAZINE_MOTIF"

_CANONICAL_FAMILIES = (
    "azine",
    "hydrazone",
    "bken",
    "imine",
    "boest",
    "vinylene",
    "boroxine",
    "triazine",
)
_CONFIDENCE_SCORE = {"low": 1.0, "medium": 2.0, "high": 3.0}
_MAX_FAMILY_HYPOTHESES = 256
_ORIGINAL_ATOM_INDEX_PROP = "cofkit_event_original_atom_idx"


@dataclass(frozen=True)
class LinkageEvent:
    """One indivisible, local linkage interpretation."""

    event_id: str
    family: str
    atoms: tuple[int, ...]
    bonds: tuple[int, ...]
    cut_bonds: tuple[int, ...]
    instance_ids: tuple[str | None, ...]
    confidence: str
    endpoint_roles: tuple[tuple[int, str], ...]
    site_id: str
    metadata: Mapping[str, object] = field(default_factory=dict)

    def to_dict(self) -> dict[str, object]:
        return {
            "event_id": self.event_id,
            "family": self.family,
            "atoms": list(self.atoms),
            "bonds": list(self.bonds),
            "cut_bonds": list(self.cut_bonds),
            "instance_ids": list(self.instance_ids),
            "confidence": self.confidence,
            "endpoint_roles": [
                {"atom": atom_idx, "role": role}
                for atom_idx, role in self.endpoint_roles
            ],
            "site_id": self.site_id,
            "metadata": dict(self.metadata),
        }


@dataclass(frozen=True)
class ReconstructedRole:
    """A reactive precursor role recovered on one cut fragment."""

    role: str
    fragment_id: int
    reactive_atoms: tuple[int, ...]
    connectivity: int
    validation_passed: bool
    metadata: Mapping[str, object] = field(default_factory=dict)

    def to_dict(self) -> dict[str, object]:
        return {
            "role": self.role,
            "fragment_id": self.fragment_id,
            "reactive_atoms": list(self.reactive_atoms),
            "connectivity": self.connectivity,
            "validation_passed": self.validation_passed,
            "metadata": dict(self.metadata),
        }


@dataclass(frozen=True)
class ReconstructionHypothesis:
    """One event set and the globally validated reconstruction it produces."""

    hypothesis_id: str
    events: tuple[LinkageEvent, ...]
    monomers: tuple[legacy.DecomposedMonomer, ...] = ()
    roles: tuple[ReconstructedRole, ...] = ()
    topology_graph: legacy.LinkageTopologyGraph | None = None
    topology: str | None = None
    cofid: str | None = None
    score: float = 0.0
    status: str = EVENT_STATUS_ENDPOINT
    validation_errors: tuple[str, ...] = ()
    metadata: Mapping[str, object] = field(default_factory=dict)
    fragment_by_atom: Mapping[int, int] = field(default_factory=dict, repr=False)

    @property
    def family(self) -> str:
        families = tuple(sorted({event.family for event in self.events}))
        return families[0] if len(families) == 1 else "mixed"

    @property
    def complete(self) -> bool:
        return self.status == EVENT_STATUS_COMPLETE and self.cofid is not None

    def to_dict(self) -> dict[str, object]:
        return {
            "hypothesis_id": self.hypothesis_id,
            "family": self.family,
            "event_ids": [event.event_id for event in self.events],
            "monomers": [
                {
                    "connectivity": monomer.connectivity,
                    "reactive_group": monomer.reactive_group,
                    "canonical_smiles": monomer.canonical_smiles,
                    "amount": monomer.amount,
                }
                for monomer in self.monomers
            ],
            "roles": [role.to_dict() for role in self.roles],
            "topology": self.topology,
            "cofid": self.cofid,
            "score": self.score,
            "status": self.status,
            "validation_errors": list(self.validation_errors),
            "metadata": dict(self.metadata),
        }


@dataclass(frozen=True)
class EventDetectionResult:
    events: tuple[LinkageEvent, ...]
    suppressed_events: tuple[Mapping[str, object], ...] = ()
    diagnostics: Mapping[str, object] = field(default_factory=dict)

    def to_dict(self) -> dict[str, object]:
        counts = Counter(event.family for event in self.events)
        return {
            "event_count": len(self.events),
            "event_counts_by_family": dict(sorted(counts.items())),
            "events": [event.to_dict() for event in self.events],
            "suppressed_events": [dict(record) for record in self.suppressed_events],
            "diagnostics": dict(self.diagnostics),
        }


def decompose_cif_to_cofid_event(
    cif_path: str | Path,
    *,
    topology: str | None = None,
    linkage: str = "auto",
    bond_mode: str = "auto",
) -> legacy.CifDecompositionResult:
    """Run the experimental event/hypothesis decomposition pipeline."""

    input_path = Path(cif_path)
    normalized_linkage = str(linkage).strip().lower()
    requested_family: str | None = None
    if normalized_linkage != "auto":
        spec = legacy._resolve_linkage_spec(normalized_linkage)
        if spec is None:
            return legacy.CifDecompositionResult(
                status="skipped",
                input_cif=str(input_path),
                topology=topology,
                linkage=str(linkage),
                reason=f"unsupported linkage {linkage!r} in event decomposition mode",
                metadata={
                    "decomposition_mode": "event",
                    "event_status": EVENT_STATUS_UNSUPPORTED,
                },
            )
        requested_family = spec.linkage_code

    normalized_topology: str | None = None
    try:
        if Chem is None:
            raise RuntimeError("RDKit is required for CIF decomposition.")
        if bond_mode not in {"auto", "distance"}:
            raise ValueError("bond_mode must be 'auto' or 'distance'")
        if topology is not None:
            normalized_topology = str(topology).strip().lower()
            if not normalized_topology:
                raise ValueError("topology must not be blank")
            legacy.default_topology_repository().get_index_entry(normalized_topology)
        atoms = read_periodic_cif_atoms(input_path)
        legacy._validate_supported_cif_periodicity(atoms)
        build_result = legacy._build_bonded_mol(atoms, bond_mode=bond_mode)
        detection = detect_linkage_events(build_result)
        hypotheses, generation_metadata = _generate_and_validate_hypotheses(
            input_path,
            atoms,
            build_result,
            detection.events,
            topology=normalized_topology,
            bond_mode=bond_mode,
        )
        hypotheses = _apply_triazine_motif_policy(hypotheses, detection.events)
        return _select_event_result(
            input_path,
            requested_family=requested_family,
            topology=normalized_topology,
            detection=detection,
            hypotheses=hypotheses,
            generation_metadata=generation_metadata,
        )
    except Exception as exc:
        return legacy.CifDecompositionResult(
            status="error",
            input_cif=str(input_path),
            topology=normalized_topology if normalized_topology is not None else (
                None if topology is None else str(topology)
            ),
            linkage=requested_family or normalized_linkage,
            reason=f"{type(exc).__name__}: {exc}",
            metadata={
                "decomposition_mode": "event",
                "event_status": EVENT_STATUS_CHEMICAL,
                "pipeline_stage": "graph_normalization_or_event_detection",
            },
        )


def detect_linkage_events(build_result: legacy.BondedMolBuildResult) -> EventDetectionResult:
    """Detect and locally resolve linkage events on one normalized CIF graph."""

    if Chem is None:
        raise RuntimeError("RDKit is required for event detection.")
    mol = Chem.Mol(build_result.mol)
    _clear_decomposition_roles(mol)
    events: list[LinkageEvent] = []
    detector_diagnostics: dict[str, object] = {}

    nitrogen_events, nitrogen_diagnostics = _detect_nitrogen_events(mol)
    events.extend(nitrogen_events)
    detector_diagnostics["nitrogen"] = nitrogen_diagnostics

    vinylene_events, vinylene_diagnostics = _detect_vinylene_events(mol)
    events.extend(vinylene_events)
    detector_diagnostics["vinylene"] = vinylene_diagnostics

    boest_events, boest_diagnostics = _detect_boronate_ester_events(mol)
    events.extend(boest_events)
    detector_diagnostics["boest"] = boest_diagnostics

    boroxine_events, boroxine_diagnostics = _detect_ring_events(
        mol,
        build_result.candidates,
        legacy._BOROXINE_SPEC,
    )
    events.extend(boroxine_events)
    detector_diagnostics["boroxine"] = boroxine_diagnostics

    triazine_events, triazine_diagnostics = _detect_ring_events(
        mol,
        build_result.candidates,
        legacy._TRIAZINE_SPEC,
    )
    events.extend(triazine_events)
    detector_diagnostics["triazine"] = triazine_diagnostics

    accepted, suppressed = _resolve_local_overlaps(tuple(events))
    return EventDetectionResult(
        events=tuple(sorted(accepted, key=lambda event: event.event_id)),
        suppressed_events=tuple(suppressed),
        diagnostics={
            "detectors": detector_diagnostics,
            "resolution_strategy": (
                "local event overlap only: azine > acylhydrazone > imine; "
                "beta-ketoenamine > overlapping vinylene"
            ),
        },
    )


def _clear_decomposition_roles(mol) -> None:
    for atom in mol.GetAtoms():
        if atom.HasProp("cofkit_decompose_role"):
            atom.ClearProp("cofkit_decompose_role")


def _event(
    mol,
    *,
    event_id: str,
    family: str,
    atoms: Iterable[int],
    bonds: Iterable[int],
    cut_bonds: Iterable[int],
    confidence: str,
    endpoint_roles: Iterable[tuple[int, str]],
    site_id: str,
    metadata: Mapping[str, object] | None = None,
) -> LinkageEvent:
    atom_indices = tuple(sorted(set(int(atom_idx) for atom_idx in atoms)))
    return LinkageEvent(
        event_id=event_id,
        family=family,
        atoms=atom_indices,
        bonds=tuple(sorted(set(int(bond_idx) for bond_idx in bonds))),
        cut_bonds=tuple(sorted(set(int(bond_idx) for bond_idx in cut_bonds))),
        instance_ids=tuple(_atom_instance_id(mol.GetAtomWithIdx(atom_idx)) for atom_idx in atom_indices),
        confidence=confidence,
        endpoint_roles=tuple(sorted(set((int(atom_idx), str(role)) for atom_idx, role in endpoint_roles))),
        site_id=site_id,
        metadata={} if metadata is None else dict(metadata),
    )


def _atom_instance_id(atom) -> str | None:
    if not atom.HasProp("instance_id"):
        return None
    value = atom.GetProp("instance_id").strip()
    return value or None


def _small_local_ring_bond_indices(mol, *, maximum_size: int = 6) -> set[int]:
    Chem.GetSymmSSSR(mol)
    return {
        int(bond_idx)
        for ring in mol.GetRingInfo().BondRings()
        if len(ring) <= maximum_size
        for bond_idx in ring
    }


def _external_heavy_neighbors(mol, atom_idx: int, *, excluding: set[int]) -> tuple[int, ...]:
    return tuple(
        sorted(
            neighbor.GetIdx()
            for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors()
            if neighbor.GetAtomicNum() > 1 and neighbor.GetIdx() not in excluding
        )
    )


def _is_aldehyde_derived_cn_carbon(mol, carbon_idx: int, nitrogen_idx: int) -> bool:
    carbon = mol.GetAtomWithIdx(carbon_idx)
    if carbon.GetAtomicNum() != 6 or carbon.GetFormalCharge() != 0:
        return False
    bond = mol.GetBondBetweenAtoms(carbon_idx, nitrogen_idx)
    if bond is None or abs(float(bond.GetBondTypeAsDouble()) - 2.0) > 1.0e-6:
        return False
    external = _external_heavy_neighbors(mol, carbon_idx, excluding={nitrogen_idx})
    if len(external) != 1 or mol.GetAtomWithIdx(external[0]).GetAtomicNum() != 6:
        return False
    heavy_valence = sum(
        float(candidate.GetBondTypeAsDouble())
        for candidate in carbon.GetBonds()
        if candidate.GetOtherAtom(carbon).GetAtomicNum() > 1
    )
    return heavy_valence <= 3.1


def _is_neutral_imine_nitrogen(
    mol,
    nitrogen_idx: int,
    carbon_idx: int,
    *,
    required_external_element: int | None,
) -> bool:
    nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
    if nitrogen.GetAtomicNum() != 7 or nitrogen.GetFormalCharge() != 0:
        return False
    external = _external_heavy_neighbors(mol, nitrogen_idx, excluding={carbon_idx})
    if len(external) != 1:
        return False
    return required_external_element is None or (
        mol.GetAtomWithIdx(external[0]).GetAtomicNum() == required_external_element
    )


def _cn_endpoints(mol, bond_idx: int) -> tuple[int, int]:
    bond = mol.GetBondWithIdx(int(bond_idx))
    first = bond.GetBeginAtom()
    second = bond.GetEndAtom()
    carbon = first if first.GetAtomicNum() == 6 else second
    nitrogen = first if first.GetAtomicNum() == 7 else second
    return int(carbon.GetIdx()), int(nitrogen.GetIdx())


def _detect_nitrogen_events(mol) -> tuple[tuple[LinkageEvent, ...], Mapping[str, object]]:
    eligible_double = legacy._eligible_carbon_nitrogen_double_bonds(mol)
    eligible_bken_single = legacy._eligible_beta_ketoenamine_single_bonds(mol)
    small_ring_bonds = _small_local_ring_bond_indices(mol)
    eligible_double = {
        bond_idx: endpoints
        for bond_idx, endpoints in eligible_double.items()
        if bond_idx not in small_ring_bonds
    }
    double_by_nitrogen: defaultdict[int, list[int]] = defaultdict(list)
    for bond_idx, (_carbon_idx, nitrogen_idx) in eligible_double.items():
        double_by_nitrogen[nitrogen_idx].append(bond_idx)

    events: list[LinkageEvent] = []
    claimed_azine: set[int] = set()
    claimed_hydrazone: set[int] = set()
    rejected_endpoint = 0

    for bond_idx, (carbon_idx, nitrogen_idx) in sorted(eligible_double.items()):
        if bond_idx in claimed_azine:
            continue
        if not _is_aldehyde_derived_cn_carbon(mol, carbon_idx, nitrogen_idx):
            continue
        nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
        partner_event: LinkageEvent | None = None
        for partner_nitrogen in nitrogen.GetNeighbors():
            partner_n_idx = int(partner_nitrogen.GetIdx())
            if partner_nitrogen.GetAtomicNum() != 7:
                continue
            nn_bond = mol.GetBondBetweenAtoms(nitrogen_idx, partner_n_idx)
            if nn_bond is None or abs(float(nn_bond.GetBondTypeAsDouble()) - 1.0) > 1.0e-6:
                continue
            partner_bonds = [
                candidate
                for candidate in double_by_nitrogen.get(partner_n_idx, ())
                if candidate != bond_idx
            ]
            if len(partner_bonds) == 1:
                partner_bond_idx = partner_bonds[0]
                partner_carbon_idx, _ = eligible_double[partner_bond_idx]
                if (
                    _is_aldehyde_derived_cn_carbon(mol, partner_carbon_idx, partner_n_idx)
                    and _is_neutral_imine_nitrogen(
                        mol,
                        nitrogen_idx,
                        carbon_idx,
                        required_external_element=7,
                    )
                    and _is_neutral_imine_nitrogen(
                        mol,
                        partner_n_idx,
                        partner_carbon_idx,
                        required_external_element=7,
                    )
                ):
                    pair = tuple(sorted((bond_idx, partner_bond_idx)))
                    partner_event = _event(
                        mol,
                        event_id=f"azine:{pair[0]}-{pair[1]}",
                        family="azine",
                        atoms=(carbon_idx, nitrogen_idx, partner_n_idx, partner_carbon_idx),
                        bonds=(bond_idx, int(nn_bond.GetIdx()), partner_bond_idx),
                        cut_bonds=pair,
                        confidence="high",
                        endpoint_roles=(
                            (carbon_idx, "aldehyde"),
                            (partner_carbon_idx, "aldehyde"),
                            (nitrogen_idx, "hydrazine"),
                            (partner_n_idx, "hydrazine"),
                        ),
                        site_id=f"azine:{pair[0]}-{pair[1]}",
                        metadata={"event_stoichiometry": "2 aldehyde sites + 1 hydrazine unit"},
                    )
                    claimed_azine.update(pair)
                    break

            if not partner_bonds and _is_acylhydrazone_partner(mol, partner_n_idx, nitrogen_idx):
                if _is_neutral_imine_nitrogen(
                    mol,
                    nitrogen_idx,
                    carbon_idx,
                    required_external_element=7,
                ):
                    carbonyl_idx = _hydrazide_carbonyl_neighbor(mol, partner_n_idx, nitrogen_idx)
                    atoms = [carbon_idx, nitrogen_idx, partner_n_idx]
                    bonds = [bond_idx, int(nn_bond.GetIdx())]
                    if carbonyl_idx is not None:
                        atoms.append(carbonyl_idx)
                        nc_bond = mol.GetBondBetweenAtoms(partner_n_idx, carbonyl_idx)
                        if nc_bond is not None:
                            bonds.append(int(nc_bond.GetIdx()))
                    partner_event = _event(
                        mol,
                        event_id=f"hydrazone:{bond_idx}",
                        family="hydrazone",
                        atoms=atoms,
                        bonds=bonds,
                        cut_bonds=(bond_idx,),
                        confidence="high",
                        endpoint_roles=((carbon_idx, "aldehyde"), (nitrogen_idx, "hydrazide")),
                        site_id=f"hydrazone:{bond_idx}",
                        metadata={"internal_family": "acylhydrazone"},
                    )
                    claimed_hydrazone.add(bond_idx)
                    break
        if partner_event is not None:
            if partner_event.event_id not in {event.event_id for event in events}:
                events.append(partner_event)

    for bond_idx, (carbon_idx, nitrogen_idx) in sorted(eligible_double.items()):
        if bond_idx in claimed_azine or bond_idx in claimed_hydrazone:
            continue
        if legacy._carbon_has_beta_ketoenamine_environment(
            mol,
            carbon_idx,
            exclude_atom_idx=nitrogen_idx,
            anchor_bond_order=1.0,
        ):
            events.append(_make_bken_event(mol, bond_idx, carbon_idx, nitrogen_idx, representation="enol_imine"))
            continue
        if not _is_aldehyde_derived_cn_carbon(mol, carbon_idx, nitrogen_idx) or not _is_neutral_imine_nitrogen(
            mol,
            nitrogen_idx,
            carbon_idx,
            required_external_element=6,
        ):
            rejected_endpoint += 1
            continue
        events.append(
            _event(
                mol,
                event_id=f"imine:{bond_idx}",
                family="imine",
                atoms=(carbon_idx, nitrogen_idx),
                bonds=(bond_idx,),
                cut_bonds=(bond_idx,),
                confidence="high",
                endpoint_roles=((carbon_idx, "aldehyde"), (nitrogen_idx, "amine")),
                site_id=f"imine:{bond_idx}",
                metadata={"event_stoichiometry": "1 aldehyde site + 1 amine site"},
            )
        )

    for bond_idx, (carbon_idx, nitrogen_idx) in sorted(eligible_bken_single.items()):
        events.append(_make_bken_event(mol, bond_idx, carbon_idx, nitrogen_idx, representation="keto_enamine"))

    counts = Counter(event.family for event in events)
    return tuple(events), {
        "eligible_cn_double_bond_count": len(eligible_double),
        "eligible_beta_ketoenamine_single_bond_count": len(eligible_bken_single),
        "small_ring_cn_rejected_count": len(
            set(legacy._eligible_carbon_nitrogen_double_bonds(mol)).intersection(small_ring_bonds)
        ),
        "strict_endpoint_rejected_count": rejected_endpoint,
        "event_counts": dict(sorted(counts.items())),
    }


def _hydrazide_carbonyl_neighbor(mol, nitrogen_idx: int, exclude_atom_idx: int) -> int | None:
    for neighbor in mol.GetAtomWithIdx(nitrogen_idx).GetNeighbors():
        if neighbor.GetIdx() == exclude_atom_idx or neighbor.GetAtomicNum() != 6:
            continue
        bond = mol.GetBondBetweenAtoms(nitrogen_idx, neighbor.GetIdx())
        if bond is None or abs(float(bond.GetBondTypeAsDouble()) - 1.0) > 1.0e-6:
            continue
        if legacy._atom_has_double_bonded_oxygen(
            mol,
            neighbor.GetIdx(),
            exclude_atom_idx=nitrogen_idx,
        ):
            return int(neighbor.GetIdx())
    return None


def _is_acylhydrazone_partner(mol, nitrogen_idx: int, imine_nitrogen_idx: int) -> bool:
    nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
    if nitrogen.GetFormalCharge() != 0:
        return False
    external = _external_heavy_neighbors(mol, nitrogen_idx, excluding={imine_nitrogen_idx})
    return len(external) == 1 and _hydrazide_carbonyl_neighbor(
        mol,
        nitrogen_idx,
        imine_nitrogen_idx,
    ) == external[0]


def _make_bken_event(
    mol,
    bond_idx: int,
    carbon_idx: int,
    nitrogen_idx: int,
    *,
    representation: str,
) -> LinkageEvent:
    context_atoms = {carbon_idx, nitrogen_idx}
    context_bonds = {bond_idx}
    for ring_anchor in mol.GetAtomWithIdx(carbon_idx).GetNeighbors():
        if ring_anchor.GetIdx() == nitrogen_idx or ring_anchor.GetAtomicNum() != 6:
            continue
        anchor_bond = mol.GetBondBetweenAtoms(carbon_idx, ring_anchor.GetIdx())
        if anchor_bond is None or not ring_anchor.IsInRing():
            continue
        if any(
            neighbor.GetAtomicNum() == 6
            and neighbor.GetIdx() != carbon_idx
            and legacy._atom_has_double_bonded_oxygen(
                mol,
                neighbor.GetIdx(),
                exclude_atom_idx=ring_anchor.GetIdx(),
            )
            for neighbor in ring_anchor.GetNeighbors()
        ):
            context_atoms.add(int(ring_anchor.GetIdx()))
            context_bonds.add(int(anchor_bond.GetIdx()))
            for neighbor in ring_anchor.GetNeighbors():
                if neighbor.GetAtomicNum() != 6 or neighbor.GetIdx() == carbon_idx:
                    continue
                if legacy._atom_has_double_bonded_oxygen(
                    mol,
                    neighbor.GetIdx(),
                    exclude_atom_idx=ring_anchor.GetIdx(),
                ):
                    context_atoms.add(int(neighbor.GetIdx()))
                    ring_bond = mol.GetBondBetweenAtoms(ring_anchor.GetIdx(), neighbor.GetIdx())
                    if ring_bond is not None:
                        context_bonds.add(int(ring_bond.GetIdx()))
            break
    return _event(
        mol,
        event_id=f"bken:{bond_idx}",
        family="bken",
        atoms=context_atoms,
        bonds=context_bonds,
        cut_bonds=(bond_idx,),
        confidence="high",
        endpoint_roles=((carbon_idx, "keto_aldehyde"), (nitrogen_idx, "amine")),
        site_id=f"bken:{bond_idx}",
        metadata={
            "representation": representation,
            "event_stoichiometry": "1 hydroxy-aldehyde site + 1 amine site",
        },
    )


def _detect_vinylene_events(mol) -> tuple[tuple[LinkageEvent, ...], Mapping[str, object]]:
    small_ring_bonds = _small_local_ring_bond_indices(mol)
    events: list[LinkageEvent] = []
    formal_double_count = 0
    same_instance_rejected = 0
    small_ring_rejected = 0
    endpoint_rejected = 0
    ambiguous_orientation_count = 0
    low_confidence_count = 0

    for bond in mol.GetBonds():
        if bond.GetIsAromatic() or abs(float(bond.GetBondTypeAsDouble()) - 2.0) > 1.0e-6:
            continue
        first = bond.GetBeginAtom()
        second = bond.GetEndAtom()
        if first.GetAtomicNum() != 6 or second.GetAtomicNum() != 6:
            continue
        formal_double_count += 1
        bond_idx = int(bond.GetIdx())
        first_idx = int(first.GetIdx())
        second_idx = int(second.GetIdx())
        first_instance = _atom_instance_id(first)
        second_instance = _atom_instance_id(second)
        if first_instance is not None and first_instance == second_instance:
            same_instance_rejected += 1
            continue
        if bond_idx in small_ring_bonds:
            small_ring_rejected += 1
            continue
        first_anchors = legacy._vinylene_carbon_anchor_indices(
            mol,
            first_idx,
            exclude_atom_idx=second_idx,
        )
        second_anchors = legacy._vinylene_carbon_anchor_indices(
            mol,
            second_idx,
            exclude_atom_idx=first_idx,
        )
        if not first_anchors or not second_anchors:
            endpoint_rejected += 1
            continue
        first_score = legacy._vinylene_activation_score(
            mol,
            first_anchors,
            exclude_atom_idx=first_idx,
        )
        second_score = legacy._vinylene_activation_score(
            mol,
            second_anchors,
            exclude_atom_idx=second_idx,
        )
        if first_score == second_score:
            orientations = (
                (second_idx, first_idx, second_score, first_score, "a"),
                (first_idx, second_idx, first_score, second_score, "b"),
            )
            ambiguous_orientation_count += 1
        elif first_score > second_score:
            orientations = ((second_idx, first_idx, first_score, second_score, "preferred"),)
        else:
            orientations = ((first_idx, second_idx, second_score, first_score, "preferred"),)

        accepted_for_bond = 0
        for aldehyde_idx, activated_idx, activated_score, other_score, variant in orientations:
            if not _is_aldehyde_derived_vinylene_endpoint(mol, aldehyde_idx, activated_idx):
                continue
            confidence = (
                "high"
                if activated_score > other_score and activated_score > 0
                else "medium"
                if activated_score > 0
                else "low"
            )
            if confidence == "low":
                low_confidence_count += 1
            atoms = {
                first_idx,
                second_idx,
                *legacy._vinylene_carbon_anchor_indices(
                    mol,
                    aldehyde_idx,
                    exclude_atom_idx=activated_idx,
                ),
                *legacy._vinylene_carbon_anchor_indices(
                    mol,
                    activated_idx,
                    exclude_atom_idx=aldehyde_idx,
                ),
            }
            context_bonds = {bond_idx}
            for atom_idx in (aldehyde_idx, activated_idx):
                partner_idx = activated_idx if atom_idx == aldehyde_idx else aldehyde_idx
                for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                    if neighbor.GetIdx() == partner_idx or neighbor.GetAtomicNum() != 6:
                        continue
                    context = mol.GetBondBetweenAtoms(atom_idx, neighbor.GetIdx())
                    if context is not None:
                        context_bonds.add(int(context.GetIdx()))
            site_id = f"vinylene:{bond_idx}"
            events.append(
                _event(
                    mol,
                    event_id=f"{site_id}:{variant}",
                    family="vinylene",
                    atoms=atoms,
                    bonds=context_bonds,
                    cut_bonds=(bond_idx,),
                    confidence=confidence,
                    endpoint_roles=(
                        (aldehyde_idx, "aldehyde"),
                        (activated_idx, "activated_methylene"),
                    ),
                    site_id=site_id,
                    metadata={
                        "orientation": variant,
                        "activated_score": activated_score,
                        "other_endpoint_score": other_score,
                        "event_stoichiometry": (
                            "1 aldehyde site + 1 activated-methylene site"
                        ),
                    },
                )
            )
            accepted_for_bond += 1
        if accepted_for_bond == 0:
            endpoint_rejected += 1

    return tuple(events), {
        "formal_cc_double_bond_count": formal_double_count,
        "same_instance_rejected_count": same_instance_rejected,
        "small_ring_rejected_count": small_ring_rejected,
        "endpoint_rejected_count": endpoint_rejected,
        "ambiguous_orientation_site_count": ambiguous_orientation_count,
        "low_confidence_event_count": low_confidence_count,
        "accepted_event_count": len(events),
        "orientation_policy": (
            "activation ranks hypotheses; tied orientations branch and zero-score "
            "structural candidates remain low confidence"
        ),
    }


def _is_aldehyde_derived_vinylene_endpoint(mol, carbon_idx: int, partner_idx: int) -> bool:
    carbon = mol.GetAtomWithIdx(carbon_idx)
    if carbon.GetAtomicNum() != 6 or carbon.GetFormalCharge() != 0:
        return False
    external = _external_heavy_neighbors(mol, carbon_idx, excluding={partner_idx})
    if len(external) != 1 or mol.GetAtomWithIdx(external[0]).GetAtomicNum() != 6:
        return False
    heavy_valence = sum(
        float(candidate.GetBondTypeAsDouble())
        for candidate in carbon.GetBonds()
        if candidate.GetOtherAtom(carbon).GetAtomicNum() > 1
    )
    return heavy_valence <= 3.1


def _detect_boronate_ester_events(mol) -> tuple[tuple[LinkageEvent, ...], Mapping[str, object]]:
    events: list[LinkageEvent] = []
    rejected_partial = 0
    seen_sites: set[tuple[int, int, int]] = set()
    for boron in mol.GetAtoms():
        if boron.GetAtomicNum() != 5:
            continue
        boron_idx = int(boron.GetIdx())
        qualifying_oxygens: list[tuple[int, int, int]] = []
        for oxygen in boron.GetNeighbors():
            if oxygen.GetAtomicNum() != 8:
                continue
            bo_bond = mol.GetBondBetweenAtoms(boron_idx, oxygen.GetIdx())
            if bo_bond is None or abs(float(bo_bond.GetBondTypeAsDouble()) - 1.0) > 1.0e-6:
                continue
            external = _external_heavy_neighbors(mol, oxygen.GetIdx(), excluding={boron_idx})
            if len(external) != 1 or mol.GetAtomWithIdx(external[0]).GetAtomicNum() != 6:
                continue
            if len(_external_heavy_neighbors(mol, oxygen.GetIdx(), excluding=set())) != 2:
                continue
            qualifying_oxygens.append((int(oxygen.GetIdx()), external[0], int(bo_bond.GetIdx())))
        found_for_boron = False
        for first_index in range(len(qualifying_oxygens)):
            for second_index in range(first_index + 1, len(qualifying_oxygens)):
                oxygen_1, carbon_1, bond_1 = qualifying_oxygens[first_index]
                oxygen_2, carbon_2, bond_2 = qualifying_oxygens[second_index]
                if carbon_1 == carbon_2:
                    continue
                carbon_bond = mol.GetBondBetweenAtoms(carbon_1, carbon_2)
                if carbon_bond is None:
                    continue
                external_boron = _external_heavy_neighbors(
                    mol,
                    boron_idx,
                    excluding={oxygen_1, oxygen_2},
                )
                if len(external_boron) != 1 or mol.GetAtomWithIdx(external_boron[0]).GetAtomicNum() != 6:
                    continue
                site_key = (boron_idx, min(oxygen_1, oxygen_2), max(oxygen_1, oxygen_2))
                if site_key in seen_sites:
                    continue
                seen_sites.add(site_key)
                found_for_boron = True
                oxygen_carbon_1 = mol.GetBondBetweenAtoms(oxygen_1, carbon_1)
                oxygen_carbon_2 = mol.GetBondBetweenAtoms(oxygen_2, carbon_2)
                ring_bonds = {
                    bond_1,
                    bond_2,
                    int(carbon_bond.GetIdx()),
                }
                if oxygen_carbon_1 is not None:
                    ring_bonds.add(int(oxygen_carbon_1.GetIdx()))
                if oxygen_carbon_2 is not None:
                    ring_bonds.add(int(oxygen_carbon_2.GetIdx()))
                site_id = f"boest:{boron_idx}:{site_key[1]}-{site_key[2]}"
                events.append(
                    _event(
                        mol,
                        event_id=site_id,
                        family="boest",
                        atoms=(boron_idx, oxygen_1, carbon_1, carbon_2, oxygen_2),
                        bonds=ring_bonds,
                        cut_bonds=(bond_1, bond_2),
                        confidence="high",
                        endpoint_roles=(
                            (boron_idx, "boronic_acid"),
                            (oxygen_1, "catechol"),
                            (oxygen_2, "catechol"),
                        ),
                        site_id=site_id,
                        metadata={
                            "ring_size": 5,
                            "event_stoichiometry": (
                                "2 B-O cuts + 1 boronic-acid site + 1 vicinal-diol site"
                            ),
                        },
                    )
                )
        if qualifying_oxygens and not found_for_boron:
            rejected_partial += 1
    return tuple(events), {
        "accepted_event_count": len(events),
        "partial_or_malformed_boron_site_count": rejected_partial,
        "atomic_event_policy": "paired B-O cuts in one five-membered B-O-C-C-O event",
    }


def _detect_ring_events(
    mol,
    candidates: tuple[legacy.BondCandidate, ...],
    spec: legacy.RingDecompositionSpec,
) -> tuple[tuple[LinkageEvent, ...], Mapping[str, object]]:
    Chem.GetSymmSSSR(mol)
    candidates_by_pair: defaultdict[frozenset[int], list[legacy.BondCandidate]] = defaultdict(list)
    for candidate in candidates:
        candidates_by_pair[frozenset((candidate.atom_idx_1, candidate.atom_idx_2))].append(candidate)
    events: list[LinkageEvent] = []
    seen_rings: set[frozenset[int]] = set()
    used_atoms: set[int] = set()
    periodic_rejected = 0
    overlap_rejected = 0
    for raw_ring in mol.GetRingInfo().AtomRings():
        ring = tuple(int(atom_idx) for atom_idx in raw_ring)
        ring_key = frozenset(ring)
        if ring_key in seen_rings or not legacy._matches_ring_decomposition_pattern(mol, ring, spec):
            continue
        seen_rings.add(ring_key)
        if used_atoms.intersection(ring):
            overlap_rejected += 1
            continue
        ring_bonds: list[int] = []
        for index, atom_idx in enumerate(ring):
            next_idx = ring[(index + 1) % len(ring)]
            bond = mol.GetBondBetweenAtoms(atom_idx, next_idx)
            if bond is None:
                ring_bonds = []
                break
            ring_bonds.append(int(bond.GetIdx()))
        if len(ring_bonds) != 6:
            continue
        atom_images = legacy._unwrap_ring_atom_images(ring, candidates_by_pair)
        if atom_images is None:
            periodic_rejected += 1
            continue
        anchors = tuple(
            atom_idx
            for atom_idx in ring
            if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == spec.anchor_atomic_num
        )
        if len(anchors) != 3:
            continue
        anchor_images = tuple(atom_images[atom_idx] for atom_idx in anchors)
        canonical_ring = "-".join(str(atom_idx) for atom_idx in sorted(ring))
        site_id = f"{spec.linkage_code}:ring:{canonical_ring}"
        endpoint_roles = tuple((atom_idx, spec.reactive_group) for atom_idx in anchors)
        base_metadata: dict[str, object] = {
            "ring_atoms": list(ring),
            "anchor_atom_indices": list(anchors),
            "anchor_images": [list(image) for image in anchor_images],
            "ring_size": 6,
        }
        if spec.linkage_code == "triazine":
            for parity, variant in ((0, "matching_a"), (1, "matching_b")):
                retained = tuple(ring_bonds[index] for index in range(6) if index % 2 == parity)
                cut = tuple(ring_bonds[index] for index in range(6) if index % 2 != parity)
                events.append(
                    _event(
                        mol,
                        event_id=f"{site_id}:{variant}",
                        family="triazine",
                        atoms=ring,
                        bonds=ring_bonds,
                        cut_bonds=cut,
                        confidence="high",
                        endpoint_roles=endpoint_roles,
                        site_id=site_id,
                        metadata={
                            **base_metadata,
                            "variant": variant,
                            "retained_bonds": list(retained),
                            "event_stoichiometry": "3 nitrile sites per triazine ring",
                        },
                    )
                )
        else:
            events.append(
                _event(
                    mol,
                    event_id=site_id,
                    family=spec.linkage_code,
                    atoms=ring,
                    bonds=ring_bonds,
                    cut_bonds=ring_bonds,
                    confidence="high",
                    endpoint_roles=endpoint_roles,
                    site_id=site_id,
                    metadata={
                        **base_metadata,
                        "event_stoichiometry": "3 boronic-acid sites per boroxine ring",
                    },
                )
            )
        used_atoms.update(ring)
    return tuple(events), {
        "structural_ring_count": len(seen_rings),
        "accepted_event_count": len(events),
        "periodic_closure_rejected_count": periodic_rejected,
        "overlapping_ring_rejected_count": overlap_rejected,
        "variant_policy": (
            "two alternating nitrile matchings"
            if spec.linkage_code == "triazine"
            else "complete ring is one atomic event"
        ),
    }


def _resolve_local_overlaps(
    events: tuple[LinkageEvent, ...],
) -> tuple[tuple[LinkageEvent, ...], tuple[Mapping[str, object], ...]]:
    bken_claimed_bonds = {
        bond_idx
        for event in events
        if event.family == "bken"
        for bond_idx in event.bonds
    }
    accepted: list[LinkageEvent] = []
    suppressed: list[Mapping[str, object]] = []
    for event in events:
        if event.family == "vinylene" and set(event.cut_bonds).intersection(bken_claimed_bonds):
            suppressed.append({
                "event": event.to_dict(),
                "status": EVENT_STATUS_SUPPRESSED,
                "reason": "overlaps canonical beta-ketoenamine event",
                "priority": "beta-ketoenamine > overlapping vinylene",
            })
            continue
        accepted.append(event)
    return tuple(accepted), tuple(suppressed)


def _generate_and_validate_hypotheses(
    input_path: Path,
    atoms: PeriodicCifAtoms,
    build_result: legacy.BondedMolBuildResult,
    events: tuple[LinkageEvent, ...],
    *,
    topology: str | None,
    bond_mode: str,
) -> tuple[tuple[ReconstructionHypothesis, ...], Mapping[str, object]]:
    events_by_family: defaultdict[str, list[LinkageEvent]] = defaultdict(list)
    for event in events:
        events_by_family[event.family].append(event)

    hypotheses: list[ReconstructionHypothesis] = []
    family_metadata: dict[str, object] = {}
    for family in _CANONICAL_FAMILIES:
        family_events = tuple(events_by_family.get(family, ()))
        if not family_events:
            family_metadata[family] = {
                "event_count": 0,
                "site_count": 0,
                "hypothesis_count": 0,
                "truncated": False,
            }
            continue
        event_sets, truncated, theoretical_count = _family_event_sets(family_events)
        family_hypotheses = [
            _reconstruct_and_validate_hypothesis(
                input_path,
                atoms,
                build_result,
                selected_events,
                topology=topology,
                bond_mode=bond_mode,
            )
            for selected_events in event_sets
        ]
        hypotheses.extend(family_hypotheses)
        family_metadata[family] = {
            "event_count": len(family_events),
            "site_count": len({event.site_id for event in family_events}),
            "hypothesis_count": len(family_hypotheses),
            "theoretical_hypothesis_count": theoretical_count,
            "truncated": truncated,
        }

    family_atom_sets = {
        family: set().union(*(set(event.atoms) for event in family_events))
        for family, family_events in events_by_family.items()
        if family_events
    }
    nonoverlapping_family_pairs = [
        [first, second]
        for first_index, first in enumerate(sorted(family_atom_sets))
        for second in sorted(family_atom_sets)[first_index + 1 :]
        if family_atom_sets[first].isdisjoint(family_atom_sets[second])
    ]
    return tuple(hypotheses), {
        "families": family_metadata,
        "maximum_hypotheses_per_family": _MAX_FAMILY_HYPOTHESES,
        "mixed_family_policy": (
            "non-overlapping family combinations are reported for future mixed-linkage "
            "serialization; current COFid output requires one linkage family"
        ),
        "nonoverlapping_family_pairs": nonoverlapping_family_pairs,
    }


def _family_event_sets(
    events: tuple[LinkageEvent, ...],
) -> tuple[tuple[tuple[LinkageEvent, ...], ...], bool, int]:
    groups: defaultdict[str, list[LinkageEvent]] = defaultdict(list)
    for event in events:
        groups[event.site_id].append(event)
    ordered_groups = tuple(
        tuple(
            sorted(
                variants,
                key=lambda event: (-_CONFIDENCE_SCORE.get(event.confidence, 0.0), event.event_id),
            )
        )
        for _site_id, variants in sorted(groups.items())
    )
    theoretical_count = 1
    for variants in ordered_groups:
        theoretical_count *= len(variants)
    combinations = product(*ordered_groups)
    selected = tuple(islice(combinations, _MAX_FAMILY_HYPOTHESES))
    return selected, theoretical_count > _MAX_FAMILY_HYPOTHESES, theoretical_count


def _reconstruct_and_validate_hypothesis(
    input_path: Path,
    atoms: PeriodicCifAtoms,
    build_result: legacy.BondedMolBuildResult,
    events: tuple[LinkageEvent, ...],
    *,
    topology: str | None,
    bond_mode: str,
) -> ReconstructionHypothesis:
    family = events[0].family
    hypothesis_id = f"{family}:" + "|".join(event.event_id for event in events)
    base_score = sum(_CONFIDENCE_SCORE.get(event.confidence, 0.0) for event in events)
    hypothesis = ReconstructionHypothesis(
        hypothesis_id=hypothesis_id,
        events=events,
        score=base_score,
        status=EVENT_STATUS_ENDPOINT,
    )
    try:
        spec = legacy._resolve_linkage_spec(family)
        if spec is None:
            return replace(
                hypothesis,
                status=EVENT_STATUS_UNSUPPORTED,
                validation_errors=(f"unsupported event family {family!r}",),
            )
        cut_result = _cut_and_reconstruct(build_result, events, spec)
        hypothesis = replace(
            hypothesis,
            monomers=cut_result.monomers,
            roles=cut_result.roles,
            topology_graph=cut_result.topology_graph,
            metadata=cut_result.metadata,
            fragment_by_atom=cut_result.fragment_by_atom,
            score=base_score + 25.0,
        )
        if cut_result.status != "ok":
            return replace(
                hypothesis,
                status=cut_result.status,
                validation_errors=cut_result.errors,
            )

        recovery_error, recovery_metadata = legacy._validate_recovered_precursors(
            cut_result.monomers,
            spec,
        )
        hypothesis_metadata = {
            **dict(hypothesis.metadata),
            "precursor_recovery_validation": recovery_metadata,
        }
        if recovery_error is not None:
            return replace(
                hypothesis,
                status=EVENT_STATUS_CHEMICAL,
                validation_errors=(recovery_error,),
                metadata=hypothesis_metadata,
                score=base_score + 50.0,
            )

        selected_topology = topology
        if selected_topology is not None:
            topology_error, topology_validation = legacy._validate_selected_topology(
                cut_result.topology_graph,
                selected_topology,
                cut_result.monomers,
                spec,
            )
            hypothesis_metadata["topology_validation"] = topology_validation
            if topology_error is not None:
                return replace(
                    hypothesis,
                    status=EVENT_STATUS_TOPOLOGY,
                    validation_errors=(topology_error,),
                    metadata=hypothesis_metadata,
                    score=base_score + 75.0,
                )
        else:
            details = legacy.LinkageDecompositionDetails(
                cut_result.monomers,
                dict(cut_result.metadata),
                cut_result.topology_graph,
            )
            detection = legacy.detect_cif_topology(
                input_path,
                linkage=family,
                bond_mode=bond_mode,
                atoms=atoms,
                details=details,
            )
            hypothesis_metadata["topology_detection"] = detection.to_dict()
            if not detection.ok:
                return replace(
                    hypothesis,
                    status=EVENT_STATUS_TOPOLOGY,
                    validation_errors=(detection.reason or "topology auto-detection failed",),
                    metadata=hypothesis_metadata,
                    score=base_score + 75.0,
                )
            selected_topology = str(detection.selected_topology)

        cofid_monomers = tuple(
            sorted(
                (monomer.to_cofid_monomer() for monomer in cut_result.monomers),
                key=lambda monomer: (-monomer.connectivity, monomer.canonical_smiles),
            )
        )
        cofid = serialize_cofid(
            monomers=cofid_monomers,
            topology=selected_topology,
            linkage=family,
        )
        parse_cofid(cofid)
        cofid_to_build_request(cofid)
        return replace(
            hypothesis,
            topology=selected_topology,
            cofid=cofid,
            status=EVENT_STATUS_COMPLETE,
            validation_errors=(),
            metadata=hypothesis_metadata,
            score=base_score + 1000.0,
        )
    except Exception as exc:
        return replace(
            hypothesis,
            status=EVENT_STATUS_CHEMICAL,
            validation_errors=(f"{type(exc).__name__}: {exc}",),
            score=base_score + 10.0,
        )


@dataclass(frozen=True)
class _CutReconstructionResult:
    status: str
    monomers: tuple[legacy.DecomposedMonomer, ...]
    roles: tuple[ReconstructedRole, ...]
    topology_graph: legacy.LinkageTopologyGraph | None
    fragment_by_atom: Mapping[int, int]
    errors: tuple[str, ...] = ()
    metadata: Mapping[str, object] = field(default_factory=dict)


def _cut_and_reconstruct(
    build_result: legacy.BondedMolBuildResult,
    events: tuple[LinkageEvent, ...],
    spec: legacy.DecompositionSpec,
) -> _CutReconstructionResult:
    mol = Chem.Mol(build_result.mol)
    _clear_decomposition_roles(mol)
    for atom in mol.GetAtoms():
        atom.SetIntProp(_ORIGINAL_ATOM_INDEX_PROP, int(atom.GetIdx()))

    endpoint_roles: dict[int, str] = {}
    for event in events:
        for atom_idx, role in event.endpoint_roles:
            previous = endpoint_roles.get(atom_idx)
            if previous is not None and previous != role:
                return _CutReconstructionResult(
                    status=EVENT_STATUS_ENDPOINT,
                    monomers=(),
                    roles=(),
                    topology_graph=None,
                    fragment_by_atom={},
                    errors=(
                        f"atom {atom_idx} is assigned incompatible endpoint roles "
                        f"{previous!r} and {role!r}",
                    ),
                )
            endpoint_roles[atom_idx] = role
            mol.GetAtomWithIdx(atom_idx).SetProp("cofkit_decompose_role", role)

    cut_bond_indices = tuple(sorted({bond_idx for event in events for bond_idx in event.cut_bonds}))
    if not cut_bond_indices:
        return _CutReconstructionResult(
            status=EVENT_STATUS_ENDPOINT,
            monomers=(),
            roles=(),
            topology_graph=None,
            fragment_by_atom={},
            errors=("hypothesis has no cut bonds",),
        )
    if any(bond_idx < 0 or bond_idx >= mol.GetNumBonds() for bond_idx in cut_bond_indices):
        return _CutReconstructionResult(
            status=EVENT_STATUS_ENDPOINT,
            monomers=(),
            roles=(),
            topology_graph=None,
            fragment_by_atom={},
            errors=("hypothesis references a bond outside the normalized graph",),
        )

    if events[0].family == "triazine":
        retained_bonds = {
            int(bond_idx)
            for event in events
            for bond_idx in event.metadata.get("retained_bonds", ())
        }
        ring_atoms = {atom_idx for event in events for atom_idx in event.atoms}
        for atom_idx in ring_atoms:
            mol.GetAtomWithIdx(atom_idx).SetIsAromatic(False)
        for event in events:
            for bond_idx in event.bonds:
                bond = mol.GetBondWithIdx(bond_idx)
                bond.SetIsAromatic(False)
                if bond_idx in retained_bonds:
                    bond.SetBondType(Chem.BondType.TRIPLE)

    cut_pairs = []
    for bond_idx in cut_bond_indices:
        bond = mol.GetBondWithIdx(bond_idx)
        cut_pairs.append((int(bond.GetBeginAtomIdx()), int(bond.GetEndAtomIdx())))

    framework_atoms = _event_framework_atoms(mol, events)
    allowed_residue_atoms = _allowed_linkage_residue_atoms(mol, events)
    editable = Chem.RWMol(mol)
    for bond_idx in sorted(cut_bond_indices, reverse=True):
        bond = editable.GetBondWithIdx(bond_idx)
        editable.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    cut_mol = editable.GetMol()
    fragment_members = tuple(
        tuple(int(atom_idx) for atom_idx in members)
        for members in Chem.GetMolFrags(cut_mol, asMols=False, sanitizeFrags=False)
    )
    fragment_by_atom = {
        atom_idx: fragment_idx
        for fragment_idx, members in enumerate(fragment_members)
        for atom_idx in members
    }
    fragments = Chem.GetMolFrags(cut_mol, asMols=True, sanitizeFrags=False)

    monomers_by_key: dict[
        tuple[str, str, int],
        tuple[legacy.DecomposedMonomer, int],
    ] = {}
    roles: list[ReconstructedRole] = []
    unexplained_fragments: list[dict[str, object]] = []
    ignored_guest_fragments = 0
    for fragment_idx, (fragment, members) in enumerate(zip(fragments, fragment_members)):
        fragment_role_atoms: defaultdict[str, list[int]] = defaultdict(list)
        for atom in fragment.GetAtoms():
            if atom.HasProp("cofkit_decompose_role"):
                original_idx = (
                    atom.GetIntProp(_ORIGINAL_ATOM_INDEX_PROP)
                    if atom.HasProp(_ORIGINAL_ATOM_INDEX_PROP)
                    else int(atom.GetIdx())
                )
                fragment_role_atoms[atom.GetProp("cofkit_decompose_role")].append(original_idx)
        if not fragment_role_atoms:
            member_set = set(members)
            if not member_set.intersection(framework_atoms):
                ignored_guest_fragments += 1
                continue
            if member_set.issubset(allowed_residue_atoms):
                continue
            unexplained_fragments.append({
                "fragment_id": fragment_idx,
                "atom_indices": list(members),
                "atomic_numbers": [mol.GetAtomWithIdx(atom_idx).GetAtomicNum() for atom_idx in members],
            })
            continue
        if len(fragment_role_atoms) != 1:
            return _CutReconstructionResult(
                status=EVENT_STATUS_ENDPOINT,
                monomers=(),
                roles=tuple(roles),
                topology_graph=None,
                fragment_by_atom=fragment_by_atom,
                errors=(
                    f"fragment {fragment_idx} owns multiple precursor roles: "
                    f"{sorted(fragment_role_atoms)!r}",
                ),
                metadata={"fragment_role_atoms": dict(fragment_role_atoms)},
            )
        role = next(iter(fragment_role_atoms))
        monomer = _repair_event_fragment(fragment, spec, events[0].family)
        if monomer is None:
            return _CutReconstructionResult(
                status=EVENT_STATUS_ENDPOINT,
                monomers=(),
                roles=tuple(roles),
                topology_graph=None,
                fragment_by_atom=fragment_by_atom,
                errors=(f"fragment {fragment_idx} could not be repaired as role {role!r}",),
            )
        roles.append(
            ReconstructedRole(
                role=role,
                fragment_id=fragment_idx,
                reactive_atoms=tuple(sorted(fragment_role_atoms[role])),
                connectivity=monomer.connectivity,
                validation_passed=True,
            )
        )
        key = (monomer.reactive_group, monomer.canonical_smiles, monomer.connectivity)
        previous = monomers_by_key.get(key)
        if previous is None:
            monomers_by_key[key] = (monomer, 1)
        else:
            monomers_by_key[key] = (previous[0], previous[1] + 1)

    monomers = tuple(
        sorted(
            (
                legacy.DecomposedMonomer(
                    connectivity=monomer.connectivity,
                    reactive_group=monomer.reactive_group,
                    canonical_smiles=monomer.canonical_smiles,
                    amount=amount,
                )
                for monomer, amount in monomers_by_key.values()
            ),
            key=lambda monomer: (-monomer.connectivity, monomer.reactive_group, monomer.canonical_smiles),
        )
    )
    topology_graph, topology_metadata = _event_topology_graph(
        mol,
        events,
        cut_pairs,
        fragment_by_atom,
        build_result.candidates,
    )
    metadata: dict[str, object] = {
        **dict(build_result.metadata),
        "decomposition_strategy": "event_hypothesis",
        "n_events": len(events),
        "n_cut_bonds": len(cut_bond_indices),
        "n_fragments_after_cut": len(fragments),
        "n_unique_monomers": len(monomers),
        "n_ignored_guest_fragments": ignored_guest_fragments,
        "n_unexplained_framework_fragments": len(unexplained_fragments),
        "event_topology": topology_metadata,
    }
    if topology_graph is not None:
        metadata["topology_graph"] = topology_graph.to_metadata()
    if unexplained_fragments:
        metadata["unexplained_framework_fragments"] = unexplained_fragments
        return _CutReconstructionResult(
            status=EVENT_STATUS_UNEXPLAINED,
            monomers=monomers,
            roles=tuple(roles),
            topology_graph=topology_graph,
            fragment_by_atom=fragment_by_atom,
            errors=(
                f"{len(unexplained_fragments)} framework fragment(s) were not assigned "
                "to reconstructed precursor roles",
            ),
            metadata=metadata,
        )
    return _CutReconstructionResult(
        status="ok",
        monomers=monomers,
        roles=tuple(roles),
        topology_graph=topology_graph,
        fragment_by_atom=fragment_by_atom,
        metadata=metadata,
    )


def _event_framework_atoms(mol, events: tuple[LinkageEvent, ...]) -> set[int]:
    event_atoms = {atom_idx for event in events for atom_idx in event.atoms}
    framework: set[int] = set()
    for members in Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False):
        member_set = {int(atom_idx) for atom_idx in members}
        if member_set.intersection(event_atoms):
            framework.update(member_set)
    return framework


def _allowed_linkage_residue_atoms(mol, events: tuple[LinkageEvent, ...]) -> set[int]:
    if not events or events[0].family != "boroxine":
        return set()
    return {
        atom_idx
        for event in events
        for atom_idx in event.atoms
        if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 8
    }


def _repair_event_fragment(
    fragment,
    spec: legacy.DecompositionSpec,
    family: str,
) -> legacy.DecomposedMonomer | None:
    if family == "triazine":
        endpoint_roles = {
            atom.GetProp("cofkit_decompose_role")
            for atom in fragment.GetAtoms()
            if atom.HasProp("cofkit_decompose_role")
        }
        if endpoint_roles != {"nitrile"}:
            return None
        return legacy._finalize_repaired_fragment(
            fragment,
            "nitrile",
            connectivity_counter=legacy._connectivity_count,
        )
    if isinstance(spec, legacy.RingDecompositionSpec):
        return legacy._repair_ring_fragment_to_monomer(fragment, spec)
    return legacy._repair_linkage_fragment_to_monomer(fragment, spec)


def _event_topology_graph(
    mol,
    events: tuple[LinkageEvent, ...],
    cut_pairs: list[tuple[int, int]],
    fragment_by_atom: Mapping[int, int],
    candidates: tuple[legacy.BondCandidate, ...],
) -> tuple[legacy.LinkageTopologyGraph | None, Mapping[str, object]]:
    family = events[0].family
    if family not in {"boroxine", "triazine"}:
        graph = legacy._linkage_topology_graph(
            mol,
            cut_pairs,
            fragment_by_atom,
            candidates,
        )
        return graph, {"strategy": "cut-edge linkage graph"}
    ring_events: list[legacy.RingDecompositionEvent] = []
    for event in events:
        anchors = tuple(int(value) for value in event.metadata.get("anchor_atom_indices", ()))
        raw_images = event.metadata.get("anchor_images", ())
        images = tuple(tuple(int(axis) for axis in image) for image in raw_images)
        if len(anchors) != 3 or len(images) != 3:
            return None, {
                "strategy": "virtual ring node",
                "error": f"event {event.event_id} lacks three anchor images",
            }
        ring_events.append(
            legacy.RingDecompositionEvent(
                anchor_atom_indices=anchors,  # type: ignore[arg-type]
                anchor_images=images,  # type: ignore[arg-type]
            )
        )
    atom_image_potentials = legacy._fragment_atom_image_potentials(fragment_by_atom, candidates)
    graph, component_count = legacy._ring_topology_graph(
        tuple(ring_events),
        fragment_by_atom,
        atom_image_potentials,
    )
    return graph, {
        "strategy": "virtual ring node",
        "n_topology_components": component_count,
    }


def _apply_triazine_motif_policy(
    hypotheses: tuple[ReconstructionHypothesis, ...],
    events: tuple[LinkageEvent, ...],
) -> tuple[ReconstructionHypothesis, ...]:
    triazine_sites: dict[str, tuple[int, ...]] = {}
    for event in events:
        if event.family != "triazine":
            continue
        triazine_sites[event.site_id] = tuple(
            int(atom_idx) for atom_idx in event.metadata.get("ring_atoms", event.atoms)
        )
    if not triazine_sites:
        return hypotheses

    complete_competitors = [
        hypothesis
        for hypothesis in hypotheses
        if hypothesis.complete and hypothesis.family != "triazine"
    ]
    containing_competitors: list[ReconstructionHypothesis] = []
    for hypothesis in complete_competitors:
        if all(
            _atoms_share_one_fragment(ring_atoms, hypothesis.fragment_by_atom)
            for ring_atoms in triazine_sites.values()
        ):
            containing_competitors.append(hypothesis)
    if not containing_competitors:
        return hypotheses

    competitor_families = tuple(sorted({hypothesis.family for hypothesis in containing_competitors}))
    updated: list[ReconstructionHypothesis] = []
    for hypothesis in hypotheses:
        if hypothesis.family != "triazine":
            updated.append(hypothesis)
            continue
        metadata = {
            **dict(hypothesis.metadata),
            "triazine_classification": {
                "classification": "triazine-containing motif",
                "competing_complete_linkages": list(competitor_families),
                "policy": (
                    "the intact triazine ring is wholly contained in a monomer from another "
                    "complete reconstruction"
                ),
                "candidate_cofid": hypothesis.cofid,
                "preclassification_status": hypothesis.status,
                "preclassification_errors": list(hypothesis.validation_errors),
            },
        }
        updated.append(
            replace(
                hypothesis,
                status=EVENT_STATUS_TRIAZINE_MOTIF,
                cofid=None,
                validation_errors=(
                    "triazine is an intact monomer motif under complete "
                    f"{', '.join(competitor_families)} reconstruction",
                ),
                metadata=metadata,
                score=hypothesis.score - 500.0,
            )
        )
    return tuple(updated)


def _atoms_share_one_fragment(
    atom_indices: Sequence[int],
    fragment_by_atom: Mapping[int, int],
) -> bool:
    fragments = {fragment_by_atom.get(int(atom_idx)) for atom_idx in atom_indices}
    return None not in fragments and len(fragments) == 1


def _select_event_result(
    input_path: Path,
    *,
    requested_family: str | None,
    topology: str | None,
    detection: EventDetectionResult,
    hypotheses: tuple[ReconstructionHypothesis, ...],
    generation_metadata: Mapping[str, object],
) -> legacy.CifDecompositionResult:
    considered = tuple(
        hypothesis
        for hypothesis in hypotheses
        if requested_family is None or hypothesis.family == requested_family
    )
    complete_by_key: dict[tuple[str, str], ReconstructionHypothesis] = {}
    for hypothesis in considered:
        if not hypothesis.complete or hypothesis.cofid is None:
            continue
        key = (hypothesis.family, hypothesis.cofid)
        previous = complete_by_key.get(key)
        if previous is None or hypothesis.score > previous.score:
            complete_by_key[key] = hypothesis
    complete = tuple(
        sorted(
            complete_by_key.values(),
            key=lambda hypothesis: (-hypothesis.score, hypothesis.family, hypothesis.cofid or ""),
        )
    )

    metadata: dict[str, object] = {
        "decomposition_mode": "event",
        "event_pipeline_version": 1,
        "event_detection": detection.to_dict(),
        "hypothesis_generation": dict(generation_metadata),
        "hypotheses": [hypothesis.to_dict() for hypothesis in hypotheses],
        "benchmark_contract": {
            "legacy_default_unchanged": True,
            "selection_unit": "globally validated reconstruction hypothesis",
            "external_labelled_benchmark_pending": True,
        },
    }
    requested_linkage = requested_family or "auto"
    if len(complete) == 1:
        selected = complete[0]
        metadata.update({
            "event_status": EVENT_STATUS_COMPLETE,
            "selected_hypothesis_id": selected.hypothesis_id,
            "successful_hypothesis_count": 1,
        })
        return legacy.CifDecompositionResult(
            status="success",
            input_cif=str(input_path),
            topology=selected.topology,
            linkage=selected.family,
            cofid=selected.cofid,
            monomers=selected.monomers,
            metadata=metadata,
        )
    if len(complete) > 1:
        labels = ", ".join(
            f"{hypothesis.family}:{hypothesis.hypothesis_id}"
            for hypothesis in complete[:8]
        )
        metadata.update({
            "event_status": EVENT_STATUS_AMBIGUOUS,
            "successful_hypothesis_count": len(complete),
            "successful_hypothesis_ids": [hypothesis.hypothesis_id for hypothesis in complete],
        })
        return legacy.CifDecompositionResult(
            status="ambiguous",
            input_cif=str(input_path),
            topology=topology,
            linkage=requested_linkage,
            reason=f"multiple complete event reconstructions remain: {labels}",
            metadata=metadata,
        )

    best_failure = _best_failed_hypothesis(considered)
    if best_failure is None:
        family_message = requested_family or "supported linkage"
        reason = f"no {family_message} events were detected"
        metadata.update({
            "event_status": EVENT_STATUS_UNSUPPORTED,
            "successful_hypothesis_count": 0,
        })
        return legacy.CifDecompositionResult(
            status="skipped",
            input_cif=str(input_path),
            topology=topology,
            linkage=requested_linkage,
            reason=reason,
            metadata=metadata,
        )

    reason = (
        best_failure.validation_errors[0]
        if best_failure.validation_errors
        else "no event reconstruction passed global validation"
    )
    metadata.update({
        "event_status": best_failure.status,
        "best_failed_hypothesis_id": best_failure.hypothesis_id,
        "successful_hypothesis_count": 0,
    })
    return legacy.CifDecompositionResult(
        status="skipped",
        input_cif=str(input_path),
        topology=best_failure.topology or topology,
        linkage=requested_linkage,
        monomers=best_failure.monomers,
        reason=reason,
        metadata=metadata,
    )


def _best_failed_hypothesis(
    hypotheses: tuple[ReconstructionHypothesis, ...],
) -> ReconstructionHypothesis | None:
    if not hypotheses:
        return None
    status_progress = {
        EVENT_STATUS_TRIAZINE_MOTIF: 6,
        EVENT_STATUS_TOPOLOGY: 5,
        EVENT_STATUS_CHEMICAL: 4,
        EVENT_STATUS_UNEXPLAINED: 3,
        EVENT_STATUS_ENDPOINT: 2,
        EVENT_STATUS_UNSUPPORTED: 1,
    }
    return max(
        hypotheses,
        key=lambda hypothesis: (
            status_progress.get(hypothesis.status, 0),
            hypothesis.score,
            hypothesis.hypothesis_id,
        ),
    )


__all__ = [
    "EventDetectionResult",
    "LinkageEvent",
    "ReconstructedRole",
    "ReconstructionHypothesis",
    "decompose_cif_to_cofid_event",
    "detect_linkage_events",
]
