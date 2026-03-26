from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Callable, Mapping

from .geometry import Vec3, add, matmul_vec, scale
from .model import Candidate, MonomerSpec, MotifRef


@dataclass(frozen=True)
class MonomerEligibilityAssessment:
    monomer_id: str
    role: str
    eligible: bool
    motif_count: int
    reason_codes: tuple[str, ...] = ()
    details: Mapping[str, object] = field(default_factory=dict)

    def to_metadata(self) -> dict[str, object]:
        return {
            "monomer_id": self.monomer_id,
            "role": self.role,
            "eligible": self.eligible,
            "motif_count": self.motif_count,
            "reason_codes": list(self.reason_codes),
            "details": dict(self.details),
        }


@dataclass(frozen=True)
class PostBuildConversionAssessment:
    profile_id: str
    requested: bool
    eligible: bool
    canonical_built_linkage: str | None = None
    required_external_conditions: tuple[str, ...] = ()
    reason_codes: tuple[str, ...] = ()
    changes_applied: bool = False
    monomer_assessments: Mapping[str, MonomerEligibilityAssessment] = field(default_factory=dict)
    details: Mapping[str, object] = field(default_factory=dict)

    def to_metadata(self) -> dict[str, object]:
        return {
            "requested": self.requested,
            "eligible": self.eligible,
            "canonical_built_linkage": self.canonical_built_linkage,
            "required_external_conditions": list(self.required_external_conditions),
            "reason_codes": list(self.reason_codes),
            "changes_applied": self.changes_applied,
            "monomer_assessments": {
                monomer_id: assessment.to_metadata()
                for monomer_id, assessment in sorted(self.monomer_assessments.items())
            },
            "details": dict(self.details),
        }


PostBuildConversionEvaluator = Callable[
    [Candidate, Mapping[str, MonomerSpec], bool, "PostBuildConversionProfile"],
    PostBuildConversionAssessment,
]


@dataclass(frozen=True)
class PostBuildConversionProfile:
    id: str
    description: str
    required_built_template_ids: tuple[str, ...] = ()
    required_external_conditions: tuple[str, ...] = ()
    conversion_product_family: str = ""
    evaluator: PostBuildConversionEvaluator | None = None


class PostBuildConversionRegistry:
    def __init__(self, profiles: Mapping[str, PostBuildConversionProfile] | None = None) -> None:
        self._profiles: dict[str, PostBuildConversionProfile] = dict(profiles or {})

    def register(self, profile: PostBuildConversionProfile) -> None:
        self._profiles[profile.id] = profile

    def get(self, profile_id: str) -> PostBuildConversionProfile:
        try:
            return self._profiles[profile_id]
        except KeyError as exc:
            raise KeyError(f"unknown post-build conversion profile {profile_id!r}") from exc

    def supported_profile_ids(self) -> tuple[str, ...]:
        return tuple(sorted(self._profiles))

    def assess_requested(
        self,
        candidate: Candidate,
        monomer_specs: Mapping[str, MonomerSpec],
        requested_ids: tuple[str, ...],
    ) -> dict[str, PostBuildConversionAssessment]:
        assessments: dict[str, PostBuildConversionAssessment] = {}
        for profile_id in dict.fromkeys(requested_ids):
            profile = self.get(profile_id)
            if profile.evaluator is None:
                raise ValueError(f"post-build conversion profile {profile_id!r} has no evaluator")
            assessments[profile_id] = profile.evaluator(candidate, monomer_specs, True, profile)
        return assessments

    @classmethod
    def builtin(cls) -> "PostBuildConversionRegistry":
        registry = cls()
        registry.register(
            PostBuildConversionProfile(
                id="sulfur_enabled_imine_conversion",
                description=(
                    "Builds the canonical imine graph, then enables a benzothiazole-like "
                    "external-sulfur annulation during atomistic realization for eligible events."
                ),
                required_built_template_ids=("imine_bridge",),
                required_external_conditions=("sulfur_enabled",),
                conversion_product_family="benzothiazole_linkage",
                evaluator=_assess_sulfur_enabled_imine_conversion,
            )
        )
        return registry


def builtin_post_build_conversion_registry() -> PostBuildConversionRegistry:
    return PostBuildConversionRegistry.builtin()


def annotate_post_build_conversions(
    candidate: Candidate,
    monomer_specs: Mapping[str, MonomerSpec],
    requested_ids: tuple[str, ...],
    *,
    registry: PostBuildConversionRegistry | None = None,
) -> Candidate:
    if not requested_ids:
        return candidate
    effective_registry = registry or builtin_post_build_conversion_registry()
    assessments = effective_registry.assess_requested(candidate, monomer_specs, requested_ids)
    metadata = dict(candidate.metadata)
    metadata["post_build_conversions"] = {
        profile_id: assessment.to_metadata()
        for profile_id, assessment in sorted(assessments.items())
    }
    return replace(candidate, metadata=metadata)


def is_imine_cof_amine_monomer(monomer: MonomerSpec) -> bool:
    return assess_imine_cof_amine_monomer(monomer).eligible


def is_imine_cof_aldehyde_monomer(monomer: MonomerSpec) -> bool:
    return assess_imine_cof_aldehyde_monomer(monomer).eligible


def assess_imine_cof_amine_monomer(monomer: MonomerSpec) -> MonomerEligibilityAssessment:
    reasons: list[str] = []
    motif_kinds = {motif.kind for motif in monomer.motifs}
    if not monomer.motifs:
        reasons.append("no_reactive_motifs")
    if len(monomer.motifs) < 2:
        reasons.append("insufficient_topicity")
    if motif_kinds and motif_kinds != {"amine"}:
        reasons.append("mixed_or_non_amine_motifs")
    if not monomer.atom_symbols or not monomer.bonds:
        reasons.append("non_atomistic_monomer")

    if not reasons:
        for motif in monomer.motifs:
            reasons.extend(_amine_site_reason_codes(monomer, motif))

    return MonomerEligibilityAssessment(
        monomer_id=monomer.id,
        role="amine",
        eligible=not reasons,
        motif_count=len(monomer.motifs),
        reason_codes=_dedupe(reasons),
        details={"motif_ids": tuple(motif.id for motif in monomer.motifs)},
    )


def assess_imine_cof_aldehyde_monomer(monomer: MonomerSpec) -> MonomerEligibilityAssessment:
    reasons: list[str] = []
    motif_kinds = {motif.kind for motif in monomer.motifs}
    if not monomer.motifs:
        reasons.append("no_reactive_motifs")
    if len(monomer.motifs) < 2:
        reasons.append("insufficient_topicity")
    if motif_kinds and motif_kinds != {"aldehyde"}:
        reasons.append("mixed_or_non_aldehyde_motifs")
    if not monomer.atom_symbols or not monomer.bonds:
        reasons.append("non_atomistic_monomer")

    if not reasons:
        for motif in monomer.motifs:
            reasons.extend(_aldehyde_site_reason_codes(monomer, motif))

    return MonomerEligibilityAssessment(
        monomer_id=monomer.id,
        role="aldehyde",
        eligible=not reasons,
        motif_count=len(monomer.motifs),
        reason_codes=_dedupe(reasons),
        details={"motif_ids": tuple(motif.id for motif in monomer.motifs)},
    )


def _assess_sulfur_enabled_imine_conversion(
    candidate: Candidate,
    monomer_specs: Mapping[str, MonomerSpec],
    requested: bool,
    profile: PostBuildConversionProfile,
) -> PostBuildConversionAssessment:
    reasons: list[str] = []
    built_template_ids = tuple(sorted({event.template_id for event in candidate.events}))
    if not candidate.events:
        reasons.append("no_built_linkages")
    if built_template_ids != profile.required_built_template_ids:
        reasons.append("not_built_from_pure_imine_linkage")

    participating_monomer_ids = tuple(
        sorted({participant.monomer_id for event in candidate.events for participant in event.participants})
    )
    monomer_assessments: dict[str, MonomerEligibilityAssessment] = {}
    amine_monomer_ids: list[str] = []
    aldehyde_monomer_ids: list[str] = []

    for monomer_id in participating_monomer_ids:
        monomer = monomer_specs[monomer_id]
        motif_kinds = {motif.kind for motif in monomer.motifs}
        if motif_kinds == {"amine"}:
            assessment = assess_imine_cof_amine_monomer(monomer)
            amine_monomer_ids.append(monomer_id)
        elif motif_kinds == {"aldehyde"}:
            assessment = assess_imine_cof_aldehyde_monomer(monomer)
            aldehyde_monomer_ids.append(monomer_id)
        else:
            assessment = MonomerEligibilityAssessment(
                monomer_id=monomer_id,
                role="unsupported",
                eligible=False,
                motif_count=len(monomer.motifs),
                reason_codes=("unsupported_motif_kinds",),
                details={"motif_kinds": tuple(sorted(motif_kinds))},
            )
        monomer_assessments[monomer_id] = assessment

    if not amine_monomer_ids:
        reasons.append("missing_amine_monomer")
    if not aldehyde_monomer_ids:
        reasons.append("missing_aldehyde_monomer")
    if any(not monomer_assessments[monomer_id].eligible for monomer_id in amine_monomer_ids):
        reasons.append("ineligible_amine_monomer")
    if any(not monomer_assessments[monomer_id].eligible for monomer_id in aldehyde_monomer_ids):
        reasons.append("ineligible_aldehyde_monomer")
    if any(assessment.role == "unsupported" for assessment in monomer_assessments.values()):
        reasons.append("unsupported_participating_monomer_role")

    used_ortho_sites: set[tuple[str, int]] = set()
    event_plans: dict[str, dict[str, object]] = {}
    event_failures: dict[str, tuple[str, ...]] = {}
    if requested and not reasons:
        for event in candidate.events:
            plan, failure_codes = _benzothiazole_event_plan(
                candidate,
                event,
                monomer_specs,
                used_ortho_sites=used_ortho_sites,
                monomer_assessments=monomer_assessments,
            )
            if plan is None:
                reasons.extend(failure_codes)
                event_failures[event.id] = failure_codes
                continue
            event_plans[event.id] = plan

    eligible = requested and not reasons and len(event_plans) == len(candidate.events)
    if eligible:
        reasons.extend(
            (
                "built_from_imine_bridge",
                "valid_polyfunctional_amine_monomer",
                "valid_polyaldehyde_monomer",
                "ortho_annulation_sites_available",
                "external_sulfur_required",
                "atomistic_benzothiazole_realization_enabled",
            )
        )

    return PostBuildConversionAssessment(
        profile_id=profile.id,
        requested=requested,
        eligible=eligible,
        canonical_built_linkage="imine_bridge" if "not_built_from_pure_imine_linkage" not in reasons else None,
        required_external_conditions=profile.required_external_conditions,
        reason_codes=_dedupe(reasons),
        changes_applied=eligible,
        monomer_assessments=monomer_assessments,
        details={
            "built_template_ids": list(built_template_ids),
            "participating_monomer_ids": list(participating_monomer_ids),
            "conversion_product_family": profile.conversion_product_family,
            "conversion_realization_id": "benzothiazole_external_sulfur",
            "n_event_plans": len(event_plans),
            "event_plans": event_plans,
            "event_failures": {
                event_id: list(reason_codes)
                for event_id, reason_codes in sorted(event_failures.items())
            },
        },
    )


def _benzothiazole_event_plan(
    candidate: Candidate,
    event,
    monomer_specs: Mapping[str, MonomerSpec],
    *,
    used_ortho_sites: set[tuple[str, int]],
    monomer_assessments: Mapping[str, MonomerEligibilityAssessment],
) -> tuple[dict[str, object] | None, tuple[str, ...]]:
    reasons: list[str] = []
    if event.template_id != "imine_bridge":
        return None, ("unsupported_event_template",)

    try:
        amine_ref, aldehyde_ref = _split_bridge_participants(event, monomer_specs, "amine", "aldehyde")
    except ValueError:
        return None, ("unsupported_event_participants",)

    if not monomer_assessments.get(amine_ref.monomer_id, MonomerEligibilityAssessment("", "", False, 0)).eligible:
        reasons.append("event_amine_not_imine_eligible")
    if not monomer_assessments.get(aldehyde_ref.monomer_id, MonomerEligibilityAssessment("", "", False, 0)).eligible:
        reasons.append("event_aldehyde_not_imine_eligible")
    if reasons:
        return None, tuple(reasons)

    amine_spec = monomer_specs[amine_ref.monomer_id]
    aldehyde_spec = monomer_specs[aldehyde_ref.monomer_id]
    amine_motif = amine_spec.motif_by_id(amine_ref.motif_id)
    aldehyde_motif = aldehyde_spec.motif_by_id(aldehyde_ref.motif_id)

    nitrogen_atom_id = _motif_atom_id(amine_motif, "reactive_atom_id", "amine motif missing reactive_atom_id")
    anchor_atom_id = _motif_atom_id(amine_motif, "anchor_atom_id", "amine motif missing anchor_atom_id")
    carbon_atom_id = _motif_atom_id(aldehyde_motif, "reactive_atom_id", "aldehyde motif missing reactive_atom_id")
    aldehydic_hydrogen_atom_ids = _hydrogen_neighbor_atom_ids(aldehyde_spec, carbon_atom_id)
    if len(aldehydic_hydrogen_atom_ids) != 1:
        return None, ("event_missing_aldehydic_hydrogen",)

    ortho_options = _ortho_annulation_options(
        amine_spec,
        anchor_atom_id=anchor_atom_id,
        excluded_atom_id=nitrogen_atom_id,
    )
    if not ortho_options:
        return None, ("event_missing_convertible_ortho_c_h",)

    carbon_world = _world_atom_position(candidate, aldehyde_ref, aldehyde_spec, carbon_atom_id)
    ranked_options = sorted(
        ortho_options,
        key=lambda option: _distance(
            _world_atom_position(candidate, amine_ref, amine_spec, option[0]),
            carbon_world,
        ),
    )
    selected_ortho_carbon_atom_id = None
    selected_ortho_hydrogen_atom_id = None
    selected_distance = None
    for ortho_carbon_atom_id, ortho_hydrogen_atom_id in ranked_options:
        site_key = (amine_ref.monomer_instance_id, ortho_carbon_atom_id)
        if site_key in used_ortho_sites:
            continue
        selected_ortho_carbon_atom_id = ortho_carbon_atom_id
        selected_ortho_hydrogen_atom_id = ortho_hydrogen_atom_id
        selected_distance = _distance(
            _world_atom_position(candidate, amine_ref, amine_spec, ortho_carbon_atom_id),
            carbon_world,
        )
        used_ortho_sites.add(site_key)
        break
    if selected_ortho_carbon_atom_id is None or selected_ortho_hydrogen_atom_id is None:
        return None, ("event_ortho_site_conflict",)

    return {
        "template_id": event.template_id,
        "amine_instance_id": amine_ref.monomer_instance_id,
        "amine_monomer_id": amine_ref.monomer_id,
        "amine_motif_id": amine_ref.motif_id,
        "aldehyde_instance_id": aldehyde_ref.monomer_instance_id,
        "aldehyde_monomer_id": aldehyde_ref.monomer_id,
        "aldehyde_motif_id": aldehyde_ref.motif_id,
        "nitrogen_atom_id": nitrogen_atom_id,
        "anchor_atom_id": anchor_atom_id,
        "imine_carbon_atom_id": carbon_atom_id,
        "aldehydic_hydrogen_atom_id": aldehydic_hydrogen_atom_ids[0],
        "ortho_carbon_atom_id": selected_ortho_carbon_atom_id,
        "ortho_hydrogen_atom_id": selected_ortho_hydrogen_atom_id,
        "ortho_candidate_atom_ids": [ortho_carbon_atom_id for ortho_carbon_atom_id, _ in ranked_options],
        "imine_to_ortho_distance": selected_distance,
    }, ()


def _split_bridge_participants(
    event,
    monomer_specs: Mapping[str, MonomerSpec],
    first_kind: str,
    second_kind: str,
) -> tuple[MotifRef, MotifRef]:
    if len(event.participants) != 2:
        raise ValueError(f"event {event.id!r} expected exactly two participants")
    first_ref, second_ref = event.participants
    first_motif = monomer_specs[first_ref.monomer_id].motif_by_id(first_ref.motif_id)
    second_motif = monomer_specs[second_ref.monomer_id].motif_by_id(second_ref.motif_id)
    if first_motif.kind == first_kind and second_motif.kind == second_kind:
        return first_ref, second_ref
    if first_motif.kind == second_kind and second_motif.kind == first_kind:
        return second_ref, first_ref
    raise ValueError(
        f"event {event.id!r} does not match expected participant kinds {(first_kind, second_kind)!r}"
    )


def _motif_atom_id(motif, key: str, error_message: str) -> int:
    atom_id = motif.metadata.get(key)
    if not isinstance(atom_id, int):
        raise ValueError(error_message)
    return atom_id


def _amine_site_reason_codes(monomer: MonomerSpec, motif) -> list[str]:
    reasons: list[str] = []
    if motif.kind != "amine":
        reasons.append("unexpected_motif_kind")
        return reasons
    reactive_atom_id = motif.metadata.get("reactive_atom_id")
    anchor_atom_id = motif.metadata.get("anchor_atom_id")
    if not isinstance(reactive_atom_id, int) or not isinstance(anchor_atom_id, int):
        reasons.append("missing_site_metadata")
        return reasons
    if not _valid_atom_id(monomer, reactive_atom_id) or not _valid_atom_id(monomer, anchor_atom_id):
        reasons.append("invalid_site_metadata")
        return reasons
    if monomer.atom_symbols[reactive_atom_id] != "N":
        reasons.append("reactive_atom_not_nitrogen")
        return reasons

    neighbors = _neighbors(monomer, reactive_atom_id)
    hydrogen_neighbors = [atom_id for atom_id, _order in neighbors if monomer.atom_symbols[atom_id] == "H"]
    heavy_neighbors = [atom_id for atom_id, _order in neighbors if monomer.atom_symbols[atom_id] != "H"]
    if len(hydrogen_neighbors) != 2 or len(heavy_neighbors) != 1:
        reasons.append("amine_not_primary")
        return reasons
    if heavy_neighbors[0] != anchor_atom_id:
        reasons.append("amine_anchor_mismatch")
        return reasons
    if monomer.atom_symbols[anchor_atom_id] != "C":
        reasons.append("amine_anchor_not_carbon")
        return reasons
    if not _is_aromatic_sp2_anchor(monomer, anchor_atom_id, excluded_atom_id=reactive_atom_id):
        reasons.append("amine_anchor_not_aromatic_sp2")
    return reasons


def _aldehyde_site_reason_codes(monomer: MonomerSpec, motif) -> list[str]:
    reasons: list[str] = []
    if motif.kind != "aldehyde":
        reasons.append("unexpected_motif_kind")
        return reasons
    reactive_atom_id = motif.metadata.get("reactive_atom_id")
    anchor_atom_id = motif.metadata.get("anchor_atom_id")
    if not isinstance(reactive_atom_id, int) or not isinstance(anchor_atom_id, int):
        reasons.append("missing_site_metadata")
        return reasons
    if not _valid_atom_id(monomer, reactive_atom_id) or not _valid_atom_id(monomer, anchor_atom_id):
        reasons.append("invalid_site_metadata")
        return reasons
    if monomer.atom_symbols[reactive_atom_id] != "C":
        reasons.append("reactive_atom_not_carbonyl_carbon")
        return reasons

    neighbors = _neighbors(monomer, reactive_atom_id)
    oxygen_neighbors = [
        atom_id
        for atom_id, order in neighbors
        if monomer.atom_symbols[atom_id] == "O" and order >= 1.5
    ]
    hydrogen_neighbors = [atom_id for atom_id, _order in neighbors if monomer.atom_symbols[atom_id] == "H"]
    scaffold_neighbors = [
        atom_id
        for atom_id, _order in neighbors
        if monomer.atom_symbols[atom_id] not in {"H", "O"}
    ]
    if len(oxygen_neighbors) != 1 or len(hydrogen_neighbors) != 1 or len(scaffold_neighbors) != 1:
        reasons.append("aldehyde_not_formyl")
        return reasons
    if scaffold_neighbors[0] != anchor_atom_id:
        reasons.append("aldehyde_anchor_mismatch")
    return reasons


def _ortho_annulation_options(
    monomer: MonomerSpec,
    *,
    anchor_atom_id: int,
    excluded_atom_id: int,
) -> tuple[tuple[int, int], ...]:
    options: list[tuple[int, int]] = []
    for neighbor_atom_id, bond_order in _neighbors(monomer, anchor_atom_id):
        if neighbor_atom_id == excluded_atom_id or bond_order < 1.5:
            continue
        if monomer.atom_symbols[neighbor_atom_id] != "C":
            continue
        if not _is_aromatic_sp2_ring_carbon(monomer, neighbor_atom_id):
            continue
        hydrogen_atom_ids = _hydrogen_neighbor_atom_ids(monomer, neighbor_atom_id)
        if len(hydrogen_atom_ids) != 1:
            continue
        options.append((neighbor_atom_id, hydrogen_atom_ids[0]))
    return tuple(options)


def _hydrogen_neighbor_atom_ids(monomer: MonomerSpec, atom_id: int) -> tuple[int, ...]:
    return tuple(
        neighbor_atom_id
        for neighbor_atom_id, _order in _neighbors(monomer, atom_id)
        if monomer.atom_symbols[neighbor_atom_id] == "H"
    )


def _valid_atom_id(monomer: MonomerSpec, atom_id: int) -> bool:
    return 0 <= atom_id < len(monomer.atom_symbols)


def _neighbors(monomer: MonomerSpec, atom_id: int) -> list[tuple[int, float]]:
    neighbors: list[tuple[int, float]] = []
    for atom_id_1, atom_id_2, bond_order in monomer.bonds:
        if atom_id_1 == atom_id:
            neighbors.append((atom_id_2, bond_order))
        elif atom_id_2 == atom_id:
            neighbors.append((atom_id_1, bond_order))
    return neighbors


def _is_aromatic_sp2_anchor(monomer: MonomerSpec, atom_id: int, *, excluded_atom_id: int) -> bool:
    heavy_neighbors = [
        (neighbor_atom_id, bond_order)
        for neighbor_atom_id, bond_order in _neighbors(monomer, atom_id)
        if neighbor_atom_id != excluded_atom_id and monomer.atom_symbols[neighbor_atom_id] != "H"
    ]
    if len(heavy_neighbors) < 2:
        return False
    return any(bond_order >= 1.5 for _neighbor_atom_id, bond_order in heavy_neighbors)


def _is_aromatic_sp2_ring_carbon(monomer: MonomerSpec, atom_id: int) -> bool:
    heavy_neighbors = [
        (neighbor_atom_id, bond_order)
        for neighbor_atom_id, bond_order in _neighbors(monomer, atom_id)
        if monomer.atom_symbols[neighbor_atom_id] != "H"
    ]
    if len(heavy_neighbors) < 2:
        return False
    return any(bond_order >= 1.5 for _neighbor_atom_id, bond_order in heavy_neighbors)


def _world_atom_position(
    candidate: Candidate,
    ref: MotifRef,
    monomer: MonomerSpec,
    atom_id: int,
) -> Vec3:
    pose = candidate.state.monomer_poses[ref.monomer_instance_id]
    image_shift = _periodic_offset(candidate.state.cell, ref.periodic_image)
    return add(add(pose.translation, matmul_vec(pose.rotation_matrix, monomer.atom_positions[atom_id])), image_shift)


def _periodic_offset(
    cell: tuple[Vec3, Vec3, Vec3],
    image: tuple[int, int, int],
) -> Vec3:
    return add(add(scale(cell[0], image[0]), scale(cell[1], image[1])), scale(cell[2], image[2]))


def _distance(left: Vec3, right: Vec3) -> float:
    dx = left[0] - right[0]
    dy = left[1] - right[1]
    dz = left[2] - right[2]
    return (dx * dx + dy * dy + dz * dz) ** 0.5


def _dedupe(reason_codes: list[str]) -> tuple[str, ...]:
    return tuple(dict.fromkeys(reason_codes))


__all__ = [
    "MonomerEligibilityAssessment",
    "PostBuildConversionAssessment",
    "PostBuildConversionProfile",
    "PostBuildConversionRegistry",
    "annotate_post_build_conversions",
    "assess_imine_cof_aldehyde_monomer",
    "assess_imine_cof_amine_monomer",
    "builtin_post_build_conversion_registry",
    "is_imine_cof_aldehyde_monomer",
    "is_imine_cof_amine_monomer",
]
