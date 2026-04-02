from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

from .model import MonomerSpec, ReactionTemplate


@dataclass(frozen=True)
class BinaryBridgeRole:
    role_id: str
    motif_kind: str
    library_prefix: str = ""


@dataclass(frozen=True)
class BinaryBridgePairOrder:
    roles: tuple[BinaryBridgeRole, BinaryBridgeRole]
    ordered_indices: tuple[int, int]


@dataclass(frozen=True)
class ReactionLinkageProfile:
    template_id: str
    bridge_target_distance: float
    binary_bridge_roles: tuple[BinaryBridgeRole, ...] = ()
    event_realizer: str | None = None
    workflow_family: str = "binary_bridge"
    topology_assignment_mode: str = "topology_guided"
    geometry_profile_id: str | None = None
    validation_profile_id: str | None = None
    library_layout: str = "role_count_files"

    @property
    def supports_binary_bridge_pair_generation(self) -> bool:
        return self.workflow_family == "binary_bridge" and len(self.binary_bridge_roles) == 2

    @property
    def supports_atomistic_realization(self) -> bool:
        return self.event_realizer is not None

    @property
    def supports_topology_guided_generation(self) -> bool:
        return self.topology_assignment_mode == "topology_guided"

    def resolve_pair_order(self, first_kind: str, second_kind: str) -> BinaryBridgePairOrder:
        if len(self.binary_bridge_roles) != 2:
            raise ValueError(f"template {self.template_id!r} is not configured for binary bridge pair generation")
        first_role, second_role = self.binary_bridge_roles
        if first_kind == first_role.motif_kind and second_kind == second_role.motif_kind:
            return BinaryBridgePairOrder(roles=(first_role, second_role), ordered_indices=(0, 1))
        if first_kind == second_role.motif_kind and second_kind == first_role.motif_kind:
            return BinaryBridgePairOrder(roles=(first_role, second_role), ordered_indices=(1, 0))
        raise ValueError(
            f"pair motif kinds {(first_kind, second_kind)!r} do not match template {self.template_id!r} "
            f"roles {(first_role.motif_kind, second_role.motif_kind)!r}"
        )


@dataclass
class ReactionLibrary:
    templates: dict[str, ReactionTemplate] = field(default_factory=dict)
    linkage_profiles: dict[str, ReactionLinkageProfile] = field(default_factory=dict)

    def register(self, template: ReactionTemplate) -> None:
        self.templates[template.id] = template

    def register_linkage_profile(self, profile: ReactionLinkageProfile) -> None:
        self.linkage_profiles[profile.template_id] = profile

    def get(self, template_id: str) -> ReactionTemplate:
        try:
            return self.templates[template_id]
        except KeyError as exc:
            raise KeyError(f"unknown reaction template {template_id!r}") from exc

    def linkage_profile(self, template: ReactionTemplate | str) -> ReactionLinkageProfile | None:
        return linkage_profile(template, profiles=self.linkage_profiles)

    def bridge_target_distance(self, template: ReactionTemplate | str, *, default_bridge_distance: float = 1.5) -> float:
        return bridge_target_distance(
            template,
            profiles=self.linkage_profiles,
            default_bridge_distance=default_bridge_distance,
        )

    def supports_binary_bridge_pair_generation(self, template: ReactionTemplate | str) -> bool:
        return supports_binary_bridge_pair_generation(template, profiles=self.linkage_profiles)

    def resolve_binary_bridge_pair_order(self, template: ReactionTemplate | str, first_kind: str, second_kind: str) -> BinaryBridgePairOrder:
        return resolve_binary_bridge_pair_order(
            template,
            first_kind,
            second_kind,
            profiles=self.linkage_profiles,
        )

    def resolve_binary_bridge_monomers(
        self,
        template: ReactionTemplate | str,
        first: MonomerSpec,
        second: MonomerSpec,
    ) -> tuple[MonomerSpec, MonomerSpec]:
        return resolve_binary_bridge_monomers(
            template,
            first,
            second,
            profiles=self.linkage_profiles,
        )

    def event_realizer(self, template: ReactionTemplate | str) -> str | None:
        return linkage_event_realizer(template, profiles=self.linkage_profiles)

    def workflow_family(self, template: ReactionTemplate | str) -> str:
        return workflow_family(template, profiles=self.linkage_profiles)

    def topology_assignment_mode(self, template: ReactionTemplate | str) -> str:
        return topology_assignment_mode(template, profiles=self.linkage_profiles)

    def selected(self, ids: Iterable[str], dimensionality: str) -> tuple[ReactionTemplate, ...]:
        selected: list[ReactionTemplate] = []
        for template_id in ids:
            template = self.get(template_id)
            if template.allows_dimensionality(dimensionality):
                selected.append(template)
        return tuple(selected)

    @classmethod
    def builtin(cls) -> "ReactionLibrary":
        lib = cls()
        for template in _builtin_templates():
            lib.register(template)
        for profile in _builtin_linkage_profiles():
            lib.register_linkage_profile(profile)
        return lib


def linkage_profile(
    template: ReactionTemplate | str,
    *,
    profiles: dict[str, ReactionLinkageProfile] | None = None,
) -> ReactionLinkageProfile | None:
    lookup = _linkage_profile_lookup(profiles)
    return lookup.get(_template_id(template))


def bridge_target_distance(
    template: ReactionTemplate | str,
    *,
    profiles: dict[str, ReactionLinkageProfile] | None = None,
    default_bridge_distance: float = 1.5,
) -> float:
    profile = linkage_profile(template, profiles=profiles)
    if profile is not None:
        return profile.bridge_target_distance
    if isinstance(template, ReactionTemplate) and template.topology_role == "ring":
        return 1.45
    return default_bridge_distance


def supports_binary_bridge_pair_generation(
    template: ReactionTemplate | str,
    *,
    profiles: dict[str, ReactionLinkageProfile] | None = None,
) -> bool:
    profile = linkage_profile(template, profiles=profiles)
    return profile.supports_binary_bridge_pair_generation if profile is not None else False


def resolve_binary_bridge_pair_order(
    template: ReactionTemplate | str,
    first_kind: str,
    second_kind: str,
    *,
    profiles: dict[str, ReactionLinkageProfile] | None = None,
) -> BinaryBridgePairOrder:
    profile = linkage_profile(template, profiles=profiles)
    if profile is None:
        raise ValueError(f"template {_template_id(template)!r} has no registered binary bridge profile")
    return profile.resolve_pair_order(first_kind, second_kind)


def resolve_binary_bridge_monomers(
    template: ReactionTemplate | str,
    first: MonomerSpec,
    second: MonomerSpec,
    *,
    profiles: dict[str, ReactionLinkageProfile] | None = None,
) -> tuple[MonomerSpec, MonomerSpec]:
    first_kind = _single_motif_kind(first)
    second_kind = _single_motif_kind(second)
    order = resolve_binary_bridge_pair_order(
        template,
        first_kind,
        second_kind,
        profiles=profiles,
    )
    monomers = (first, second)
    return tuple(monomers[index] for index in order.ordered_indices)  # type: ignore[return-value]


def linkage_event_realizer(
    template: ReactionTemplate | str,
    *,
    profiles: dict[str, ReactionLinkageProfile] | None = None,
) -> str | None:
    profile = linkage_profile(template, profiles=profiles)
    return None if profile is None else profile.event_realizer


def workflow_family(
    template: ReactionTemplate | str,
    *,
    profiles: dict[str, ReactionLinkageProfile] | None = None,
) -> str:
    profile = linkage_profile(template, profiles=profiles)
    if profile is not None:
        return profile.workflow_family
    if isinstance(template, ReactionTemplate) and template.topology_role == "ring":
        return "ring_forming"
    return "template_driven"


def topology_assignment_mode(
    template: ReactionTemplate | str,
    *,
    profiles: dict[str, ReactionLinkageProfile] | None = None,
) -> str:
    profile = linkage_profile(template, profiles=profiles)
    if profile is not None:
        return profile.topology_assignment_mode
    if isinstance(template, ReactionTemplate) and template.topology_role == "ring":
        return "topology_bypass"
    return "topology_guided"


def _builtin_templates() -> tuple[ReactionTemplate, ...]:
    return (
        ReactionTemplate(
            id="imine_bridge",
            arity=2,
            reactant_motif_kinds=("amine", "aldehyde"),
            product_name="imine",
            topology_role="bridge",
            planarity_prior="planar",
            torsion_prior="restricted",
            allowed_dimensionalities=("2D", "3D"),
        ),
        ReactionTemplate(
            id="hydrazone_bridge",
            arity=2,
            reactant_motif_kinds=("hydrazide", "aldehyde"),
            product_name="hydrazone",
            topology_role="bridge",
            planarity_prior="planar",
            torsion_prior="restricted",
            allowed_dimensionalities=("2D", "3D"),
        ),
        ReactionTemplate(
            id="azine_bridge",
            arity=2,
            reactant_motif_kinds=("hydrazine", "aldehyde"),
            product_name="azine",
            topology_role="bridge",
            planarity_prior="planar",
            torsion_prior="restricted",
            allowed_dimensionalities=("2D", "3D"),
        ),
        ReactionTemplate(
            id="boronate_ester_bridge",
            arity=2,
            reactant_motif_kinds=("boronic_acid", "catechol"),
            product_name="boronate_ester",
            topology_role="bridge",
            planarity_prior="semi-planar",
            torsion_prior="moderate",
            allowed_dimensionalities=("2D", "3D"),
        ),
        ReactionTemplate(
            id="keto_enamine_bridge",
            arity=2,
            reactant_motif_kinds=("amine", "keto_aldehyde"),
            product_name="beta_ketoenamine",
            topology_role="bridge",
            planarity_prior="planar",
            torsion_prior="restricted",
            allowed_dimensionalities=("2D", "3D"),
        ),
        ReactionTemplate(
            id="vinylene_bridge",
            arity=2,
            reactant_motif_kinds=("activated_methylene", "aldehyde"),
            product_name="vinylene",
            topology_role="bridge",
            planarity_prior="planar",
            torsion_prior="restricted",
            allowed_dimensionalities=("2D", "3D"),
        ),
        ReactionTemplate(
            id="boroxine_trimerization",
            arity=3,
            reactant_motif_kinds=("boronic_acid", "boronic_acid", "boronic_acid"),
            product_name="boroxine",
            topology_role="ring",
            planarity_prior="planar",
            torsion_prior="locked",
            allowed_dimensionalities=("2D",),
        ),
        ReactionTemplate(
            id="triazine_trimerization",
            arity=3,
            reactant_motif_kinds=("nitrile", "nitrile", "nitrile"),
            product_name="triazine",
            topology_role="ring",
            planarity_prior="planar",
            torsion_prior="locked",
            allowed_dimensionalities=("2D", "3D"),
        ),
    )


def _builtin_linkage_profiles() -> tuple[ReactionLinkageProfile, ...]:
    return (
        ReactionLinkageProfile(
            template_id="imine_bridge",
            bridge_target_distance=1.3,
            binary_bridge_roles=(
                BinaryBridgeRole(role_id="amine", motif_kind="amine", library_prefix="amines"),
                BinaryBridgeRole(role_id="aldehyde", motif_kind="aldehyde", library_prefix="aldehydes"),
            ),
            event_realizer="imine_bridge",
            geometry_profile_id="imine_bridge",
            validation_profile_id="imine_bridge",
        ),
        ReactionLinkageProfile(
            template_id="hydrazone_bridge",
            bridge_target_distance=1.3,
            binary_bridge_roles=(
                BinaryBridgeRole(role_id="hydrazide", motif_kind="hydrazide", library_prefix="hydrazides"),
                BinaryBridgeRole(role_id="aldehyde", motif_kind="aldehyde", library_prefix="aldehydes"),
            ),
            event_realizer="hydrazone_bridge",
            geometry_profile_id="hydrazone_bridge",
            validation_profile_id="hydrazone_bridge",
        ),
        ReactionLinkageProfile(
            template_id="azine_bridge",
            bridge_target_distance=1.3,
            binary_bridge_roles=(
                BinaryBridgeRole(role_id="hydrazine", motif_kind="hydrazine", library_prefix="hydrazines"),
                BinaryBridgeRole(role_id="aldehyde", motif_kind="aldehyde", library_prefix="aldehydes"),
            ),
            event_realizer="azine_bridge",
            geometry_profile_id="azine_bridge",
            validation_profile_id="azine_bridge",
        ),
        ReactionLinkageProfile(
            template_id="boronate_ester_bridge",
            bridge_target_distance=1.4,
            binary_bridge_roles=(
                BinaryBridgeRole(role_id="boronic_acid", motif_kind="boronic_acid", library_prefix="boronic_acids"),
                BinaryBridgeRole(role_id="catechol", motif_kind="catechol", library_prefix="catechols"),
            ),
            event_realizer="boronate_ester_bridge",
            geometry_profile_id="boronate_ester_bridge",
            validation_profile_id="boronate_ester_bridge",
        ),
        ReactionLinkageProfile(
            template_id="keto_enamine_bridge",
            bridge_target_distance=1.36,
            binary_bridge_roles=(
                BinaryBridgeRole(role_id="amine", motif_kind="amine", library_prefix="amines"),
                BinaryBridgeRole(role_id="keto_aldehyde", motif_kind="keto_aldehyde", library_prefix="keto_aldehydes"),
            ),
            event_realizer="keto_enamine_bridge",
            geometry_profile_id="keto_enamine_bridge",
            validation_profile_id="keto_enamine_bridge",
        ),
        ReactionLinkageProfile(
            template_id="vinylene_bridge",
            bridge_target_distance=1.34,
            binary_bridge_roles=(
                BinaryBridgeRole(
                    role_id="activated_methylene",
                    motif_kind="activated_methylene",
                    library_prefix="activated_methylenes",
                ),
                    BinaryBridgeRole(role_id="aldehyde", motif_kind="aldehyde", library_prefix="aldehydes"),
            ),
            event_realizer="vinylene_bridge",
            geometry_profile_id="vinylene_bridge",
            validation_profile_id="vinylene_bridge",
        ),
        ReactionLinkageProfile(
            template_id="boroxine_trimerization",
            bridge_target_distance=1.45,
            workflow_family="ring_forming",
            topology_assignment_mode="topology_bypass",
            geometry_profile_id="boroxine_trimerization",
            validation_profile_id="boroxine_trimerization",
            library_layout="role_count_files",
        ),
        ReactionLinkageProfile(
            template_id="triazine_trimerization",
            bridge_target_distance=1.45,
            workflow_family="ring_forming",
            topology_assignment_mode="topology_bypass",
            geometry_profile_id="triazine_trimerization",
            validation_profile_id="triazine_trimerization",
            library_layout="role_count_files",
        ),
    )


def _template_id(template: ReactionTemplate | str) -> str:
    return template if isinstance(template, str) else template.id


def _linkage_profile_lookup(
    profiles: dict[str, ReactionLinkageProfile] | None,
) -> dict[str, ReactionLinkageProfile]:
    if profiles is not None:
        return profiles
    return {profile.template_id: profile for profile in _builtin_linkage_profiles()}


def _single_motif_kind(monomer: MonomerSpec) -> str:
    kinds = {motif.kind for motif in monomer.motifs}
    if len(kinds) != 1:
        raise ValueError(
            f"monomer {monomer.id!r} must expose exactly one motif kind for binary bridge generation, got {tuple(sorted(kinds))!r}"
        )
    return next(iter(kinds))
