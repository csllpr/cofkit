from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field
from math import acos, cos, pi, sin
from typing import Callable, Mapping

from .geometry import Vec3, add, cross, dot, matmul_vec, norm, normalize, scale, sub, transpose
from .model import Candidate, MonomerSpec, MotifRef, Pose, ReactionEvent
from .reactions import bridge_target_distance, linkage_event_realizer


@dataclass(frozen=True)
class RealizedAtom:
    atom_id: int
    label: str
    symbol: str
    local_position: Vec3


@dataclass(frozen=True)
class RealizedBond:
    label_1: str
    label_2: str
    distance: float
    symmetry_1: str = "."
    symmetry_2: str = "."


@dataclass(frozen=True)
class EventRealization:
    removed_atom_ids: Mapping[str, tuple[int, ...]] = field(default_factory=dict)
    atom_position_overrides: Mapping[str, Mapping[int, Vec3]] = field(default_factory=dict)
    bonds: tuple[RealizedBond, ...] = ()
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class ReactionRealizationResult:
    atoms_by_instance: Mapping[str, tuple[RealizedAtom, ...]]
    removed_atom_ids_by_instance: Mapping[str, tuple[int, ...]] = field(default_factory=dict)
    bonds: tuple[RealizedBond, ...] = ()
    metadata: Mapping[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class PostBuildRealization:
    removed_atom_ids: Mapping[str, tuple[int, ...]] = field(default_factory=dict)
    atom_position_overrides: Mapping[str, Mapping[int, Vec3]] = field(default_factory=dict)
    added_atoms: Mapping[str, tuple[RealizedAtom, ...]] = field(default_factory=dict)
    bonds: tuple[RealizedBond, ...] = ()
    notes: tuple[str, ...] = ()
    metadata: Mapping[str, object] = field(default_factory=dict)


ReactionEventRealizationHandler = Callable[
    ["ReactionRealizer", ReactionEvent, Candidate, Mapping[str, MonomerSpec]],
    EventRealization | None,
]


class ReactionEventRealizationRegistry:
    def __init__(self, handlers: Mapping[str, ReactionEventRealizationHandler] | None = None) -> None:
        self._handlers: dict[str, ReactionEventRealizationHandler] = dict(handlers or {})

    def register(self, realizer_id: str, handler: ReactionEventRealizationHandler) -> None:
        self._handlers[realizer_id] = handler

    def get(self, realizer_id: str) -> ReactionEventRealizationHandler | None:
        return self._handlers.get(realizer_id)

    def supported_realizer_ids(self) -> tuple[str, ...]:
        return tuple(sorted(self._handlers))

    @classmethod
    def builtin(cls) -> "ReactionEventRealizationRegistry":
        registry = cls()
        registry.register("imine_bridge", _dispatch_imine_bridge_event)
        registry.register("hydrazone_bridge", _dispatch_hydrazone_bridge_event)
        registry.register("keto_enamine_bridge", _dispatch_keto_enamine_bridge_event)
        registry.register("boronate_ester_bridge", _dispatch_boronate_ester_bridge_event)
        registry.register("vinylene_bridge", _dispatch_vinylene_bridge_event)
        return registry


class ReactionRealizer:
    def __init__(self, registry: ReactionEventRealizationRegistry | None = None) -> None:
        self.registry = registry or ReactionEventRealizationRegistry.builtin()

    def realize(
        self,
        candidate: Candidate,
        monomer_specs: Mapping[str, MonomerSpec],
        instance_to_monomer: Mapping[str, str],
    ) -> ReactionRealizationResult | None:
        atomistic_instances = {
            instance_id
            for instance_id, monomer_id in instance_to_monomer.items()
            if monomer_specs[monomer_id].atom_symbols and monomer_specs[monomer_id].atom_positions
        }
        if not atomistic_instances:
            return None

        removed_atom_ids: dict[str, set[int]] = defaultdict(set)
        atom_position_overrides: dict[str, dict[int, Vec3]] = defaultdict(dict)
        added_atoms_by_instance: dict[str, list[RealizedAtom]] = defaultdict(list)
        bonds: list[RealizedBond] = []
        notes: list[str] = []
        applied_events = 0
        applied_templates: Counter[str] = Counter()

        for event in candidate.events:
            realization = self._realize_event(event, candidate, monomer_specs)
            if realization is None:
                continue
            applied_events += 1
            applied_templates[event.template_id] += 1
            for instance_id, atom_ids in realization.removed_atom_ids.items():
                removed_atom_ids[instance_id].update(atom_ids)
            for instance_id, overrides in realization.atom_position_overrides.items():
                atom_position_overrides[instance_id].update(overrides)
            bonds.extend(realization.bonds)
            notes.extend(realization.notes)

        post_build_metadata: dict[str, object] = {}
        post_build_realization = self._apply_post_build_conversions(
            candidate,
            monomer_specs,
            instance_to_monomer,
            removed_atom_ids,
            atom_position_overrides,
        )
        if post_build_realization is not None:
            for instance_id, atom_ids in post_build_realization.removed_atom_ids.items():
                removed_atom_ids[instance_id].update(atom_ids)
            for instance_id, overrides in post_build_realization.atom_position_overrides.items():
                atom_position_overrides[instance_id].update(overrides)
            for instance_id, atoms in post_build_realization.added_atoms.items():
                added_atoms_by_instance[instance_id].extend(atoms)
            bonds.extend(post_build_realization.bonds)
            notes.extend(post_build_realization.notes)
            post_build_metadata = dict(post_build_realization.metadata)

        hydrogen_cleanup_metadata = self._cleanup_reaction_hydrogens(
            candidate,
            monomer_specs,
            instance_to_monomer,
            removed_atom_ids,
            atom_position_overrides,
            extra_heavy_atoms=self._added_heavy_atom_world_positions(candidate, added_atoms_by_instance),
        )

        atoms_by_instance: dict[str, tuple[RealizedAtom, ...]] = {}
        removed_symbol_counts: Counter[str] = Counter()
        added_symbol_counts: Counter[str] = Counter()
        for instance_id, monomer_id in instance_to_monomer.items():
            monomer = monomer_specs[monomer_id]
            if not (monomer.atom_symbols and monomer.atom_positions):
                continue
            removed = removed_atom_ids.get(instance_id, set())
            realized_atoms: list[RealizedAtom] = []
            for atom_id, (symbol, position) in enumerate(zip(monomer.atom_symbols, monomer.atom_positions)):
                if atom_id in removed:
                    removed_symbol_counts[symbol] += 1
                    continue
                local_position = atom_position_overrides.get(instance_id, {}).get(atom_id, position)
                realized_atoms.append(
                    RealizedAtom(
                        atom_id=atom_id,
                        label=self.atom_label(instance_id, symbol, atom_id),
                        symbol=symbol,
                        local_position=local_position,
                    )
                )
            for atom in sorted(added_atoms_by_instance.get(instance_id, ()), key=lambda atom: atom.atom_id):
                added_symbol_counts[atom.symbol] += 1
                realized_atoms.append(atom)
            atoms_by_instance[instance_id] = tuple(realized_atoms)

        if applied_events == 0:
            return None

        metadata = {
            "applied_event_count": applied_events,
            "applied_templates": dict(applied_templates),
            "removed_atom_count": sum(removed_symbol_counts.values()),
            "removed_atom_symbols": dict(removed_symbol_counts),
            "added_atom_count": sum(added_symbol_counts.values()),
            "added_atom_symbols": dict(added_symbol_counts),
            "notes": tuple(dict.fromkeys(notes)),
        }
        if hydrogen_cleanup_metadata:
            metadata["hydrogen_cleanup"] = hydrogen_cleanup_metadata
        if post_build_metadata:
            metadata["post_build_conversions"] = post_build_metadata
        return ReactionRealizationResult(
            atoms_by_instance=atoms_by_instance,
            removed_atom_ids_by_instance={
                instance_id: tuple(sorted(atom_ids)) for instance_id, atom_ids in removed_atom_ids.items()
            },
            bonds=tuple(bonds),
            metadata=metadata,
        )

    def _realize_event(
        self,
        event: ReactionEvent,
        candidate: Candidate,
        monomer_specs: Mapping[str, MonomerSpec],
    ) -> EventRealization | None:
        realizer_id = linkage_event_realizer(event.template_id)
        if realizer_id is None:
            return None
        handler = self.registry.get(realizer_id)
        if handler is None:
            return None
        return handler(self, event, candidate, monomer_specs)

    def _apply_post_build_conversions(
        self,
        candidate: Candidate,
        monomer_specs: Mapping[str, MonomerSpec],
        instance_to_monomer: Mapping[str, str],
        removed_atom_ids: Mapping[str, set[int]],
        atom_position_overrides: Mapping[str, Mapping[int, Vec3]],
    ) -> PostBuildRealization | None:
        raw_profiles = candidate.metadata.get("post_build_conversions")
        if not isinstance(raw_profiles, Mapping):
            return None

        base_removed_atom_ids = removed_atom_ids
        base_atom_position_overrides = atom_position_overrides
        conversion_removed_atom_ids: dict[str, set[int]] = defaultdict(set)
        conversion_atom_position_overrides: dict[str, dict[int, Vec3]] = defaultdict(dict)
        added_atoms: dict[str, list[RealizedAtom]] = defaultdict(list)
        bonds: list[RealizedBond] = []
        notes: list[str] = []
        applied_profiles: dict[str, object] = {}

        for profile_id, profile_metadata in raw_profiles.items():
            if not isinstance(profile_metadata, Mapping):
                continue
            if not profile_metadata.get("requested") or not profile_metadata.get("eligible"):
                continue
            if str(profile_id) == "sulfur_enabled_imine_conversion":
                realization = self._realize_sulfur_enabled_imine_conversion(
                    candidate,
                    monomer_specs,
                    instance_to_monomer,
                    self._merge_removed_atom_ids(base_removed_atom_ids, conversion_removed_atom_ids),
                    self._merge_atom_position_overrides(base_atom_position_overrides, conversion_atom_position_overrides),
                    profile_metadata,
                )
                if realization is None:
                    continue
                for instance_id, atom_ids in realization.removed_atom_ids.items():
                    conversion_removed_atom_ids[instance_id].update(atom_ids)
                for instance_id, overrides in realization.atom_position_overrides.items():
                    conversion_atom_position_overrides[instance_id].update(overrides)
                for instance_id, atoms in realization.added_atoms.items():
                    added_atoms[instance_id].extend(atoms)
                bonds.extend(realization.bonds)
                notes.extend(realization.notes)
                applied_profiles[str(profile_id)] = dict(realization.metadata)

        if not bonds and not added_atoms and not conversion_removed_atom_ids and not conversion_atom_position_overrides:
            return None
        return PostBuildRealization(
            removed_atom_ids={
                instance_id: tuple(sorted(atom_ids))
                for instance_id, atom_ids in conversion_removed_atom_ids.items()
            },
            atom_position_overrides={
                instance_id: dict(overrides)
                for instance_id, overrides in conversion_atom_position_overrides.items()
            },
            added_atoms={
                instance_id: tuple(atoms)
                for instance_id, atoms in added_atoms.items()
            },
            bonds=tuple(bonds),
            notes=tuple(notes),
            metadata={"applied_profiles": applied_profiles},
        )

    def _realize_sulfur_enabled_imine_conversion(
        self,
        candidate: Candidate,
        monomer_specs: Mapping[str, MonomerSpec],
        instance_to_monomer: Mapping[str, str],
        removed_atom_ids: Mapping[str, set[int]],
        atom_position_overrides: Mapping[str, Mapping[int, Vec3]],
        profile_metadata: Mapping[str, object],
    ) -> PostBuildRealization | None:
        details = profile_metadata.get("details")
        if not isinstance(details, Mapping):
            return None
        raw_event_plans = details.get("event_plans")
        if not isinstance(raw_event_plans, Mapping) or not raw_event_plans:
            return None

        converted_removed_atom_ids: dict[str, set[int]] = defaultdict(set)
        fitted_atom_position_overrides: dict[str, dict[int, Vec3]] = defaultdict(dict)
        added_atoms: dict[str, list[RealizedAtom]] = defaultdict(list)
        bonds: list[RealizedBond] = []
        notes: list[str] = []
        applied_event_ids: list[str] = []
        used_ortho_sites: set[tuple[str, int]] = set()
        sulfur_geometry: dict[str, dict[str, float]] = {}

        for event in candidate.events:
            raw_plan = raw_event_plans.get(event.id)
            if not isinstance(raw_plan, Mapping):
                continue
            if event.template_id != "imine_bridge":
                continue

            amine_ref, aldehyde_ref = self._split_bridge_participants(
                event,
                monomer_specs,
                "amine",
                "aldehyde",
            )
            amine_instance_id = str(raw_plan.get("amine_instance_id", amine_ref.monomer_instance_id))
            aldehyde_instance_id = str(raw_plan.get("aldehyde_instance_id", aldehyde_ref.monomer_instance_id))
            if amine_instance_id != amine_ref.monomer_instance_id or aldehyde_instance_id != aldehyde_ref.monomer_instance_id:
                continue

            amine_spec = monomer_specs[amine_ref.monomer_id]
            aldehyde_spec = monomer_specs[aldehyde_ref.monomer_id]
            amine_motif = amine_spec.motif_by_id(amine_ref.motif_id)
            aldehyde_motif = aldehyde_spec.motif_by_id(aldehyde_ref.motif_id)

            nitrogen_atom_id = self._coerce_int(raw_plan.get("nitrogen_atom_id"))
            anchor_atom_id = self._coerce_int(raw_plan.get("anchor_atom_id"))
            carbon_atom_id = self._coerce_int(raw_plan.get("imine_carbon_atom_id"))
            aldehydic_hydrogen_atom_id = self._coerce_int(raw_plan.get("aldehydic_hydrogen_atom_id"))
            ortho_carbon_atom_id = self._coerce_int(raw_plan.get("ortho_carbon_atom_id"))
            ortho_hydrogen_atom_id = self._coerce_int(raw_plan.get("ortho_hydrogen_atom_id"))
            if None in {
                nitrogen_atom_id,
                anchor_atom_id,
                carbon_atom_id,
                aldehydic_hydrogen_atom_id,
                ortho_carbon_atom_id,
                ortho_hydrogen_atom_id,
            }:
                continue

            aldehyde_anchor_atom_id = self._motif_atom_id_from_metadata(
                aldehyde_motif,
                "anchor_atom_id",
                context=f"{event.id} benzothiazole aldehyde anchor",
            )
            oxygen_atom_id = self._unique_symbol_atom_id(
                aldehyde_spec,
                aldehyde_motif.atom_ids,
                "O",
                context=f"{event.id} benzothiazole aldehyde oxygen",
            )
            current_removed_atom_ids = self._merge_removed_atom_ids(
                removed_atom_ids,
                converted_removed_atom_ids,
            )
            current_atom_position_overrides = self._merge_atom_position_overrides(
                atom_position_overrides,
                fitted_atom_position_overrides,
            )
            blocker_atoms = (
                *self._heavy_atom_world_positions(
                    candidate,
                    monomer_specs,
                    instance_to_monomer,
                    current_removed_atom_ids,
                    current_atom_position_overrides,
                ),
                *self._added_heavy_atom_world_positions(candidate, added_atoms),
            )
            ortho_option_ids = self._ortho_option_ids(
                raw_plan,
                default_ortho_atom_id=ortho_carbon_atom_id,
            )
            sort_carbon_world = self._world_atom_position(
                candidate,
                aldehyde_ref,
                aldehyde_spec,
                carbon_atom_id,
                atom_position_overrides=current_atom_position_overrides,
            )
            ortho_option_ids = tuple(
                sorted(
                    ortho_option_ids,
                    key=lambda atom_id: self._distance(
                        sort_carbon_world,
                        self._world_atom_position(
                            candidate,
                            amine_ref,
                            amine_spec,
                            atom_id,
                            atom_position_overrides=current_atom_position_overrides,
                        ),
                    ),
                )
            )
            best_fit: tuple[int, int, Vec3, dict[str, float]] | None = None
            for ortho_option_atom_id in ortho_option_ids:
                if (amine_instance_id, ortho_option_atom_id) in used_ortho_sites:
                    continue
                ortho_hydrogen_candidates = self._directly_bonded_hydrogen_atom_ids(amine_spec, ortho_option_atom_id)
                if len(ortho_hydrogen_candidates) != 1:
                    continue
                carbon_world = self._world_atom_position(
                    candidate=candidate,
                    ref=aldehyde_ref,
                    monomer=aldehyde_spec,
                    atom_id=carbon_atom_id,
                    atom_position_overrides=current_atom_position_overrides,
                )
                nitrogen_world = self._world_atom_position(
                    candidate=candidate,
                    ref=amine_ref,
                    monomer=amine_spec,
                    atom_id=nitrogen_atom_id,
                    atom_position_overrides=current_atom_position_overrides,
                )
                anchor_world = self._world_atom_position(
                    candidate=candidate,
                    ref=amine_ref,
                    monomer=amine_spec,
                    atom_id=anchor_atom_id,
                    atom_position_overrides=current_atom_position_overrides,
                )
                ortho_world = self._world_atom_position(
                    candidate=candidate,
                    ref=amine_ref,
                    monomer=amine_spec,
                    atom_id=ortho_option_atom_id,
                    atom_position_overrides=current_atom_position_overrides,
                )
                sulfur_world, sulfur_metrics = self._benzothiazole_sulfur_world_position(
                    carbon_world=carbon_world,
                    ortho_world=ortho_world,
                    nitrogen_world=nitrogen_world,
                    anchor_world=anchor_world,
                    blockers=tuple(
                        (instance_id, atom_id, world_position)
                        for instance_id, atom_id, world_position in blocker_atoms
                        if (instance_id, atom_id)
                        not in {
                            (aldehyde_instance_id, carbon_atom_id),
                            (aldehyde_instance_id, aldehyde_anchor_atom_id),
                            (amine_instance_id, nitrogen_atom_id),
                            (amine_instance_id, anchor_atom_id),
                            (amine_instance_id, ortho_option_atom_id),
                        }
                    ),
                )
                geometry_metrics = self._benzothiazole_geometry_metrics(
                    candidate=candidate,
                    amine_ref=amine_ref,
                    amine_spec=amine_spec,
                    nitrogen_atom_id=nitrogen_atom_id,
                    anchor_atom_id=anchor_atom_id,
                    ortho_carbon_atom_id=ortho_option_atom_id,
                    aldehyde_ref=aldehyde_ref,
                    aldehyde_spec=aldehyde_spec,
                    carbon_atom_id=carbon_atom_id,
                    aldehyde_anchor_atom_id=aldehyde_anchor_atom_id,
                    atom_position_overrides=current_atom_position_overrides,
                    world_overrides={},
                    sulfur_world=sulfur_world,
                    blockers=tuple(
                        (instance_id, atom_id, world_position)
                        for instance_id, atom_id, world_position in blocker_atoms
                        if (instance_id, atom_id)
                        not in {
                            (aldehyde_instance_id, carbon_atom_id),
                            (aldehyde_instance_id, aldehyde_anchor_atom_id),
                            (amine_instance_id, nitrogen_atom_id),
                            (amine_instance_id, anchor_atom_id),
                            (amine_instance_id, ortho_option_atom_id),
                        }
                    ),
                    fit_objective=sulfur_metrics["fit_objective"],
                    local_relaxation_applied=False,
                )
                geometry_metrics["placement_offset"] = sulfur_metrics["placement_offset"]
                if best_fit is None or geometry_metrics["fit_objective"] < best_fit[3]["fit_objective"]:
                    best_fit = (
                        ortho_option_atom_id,
                        ortho_hydrogen_candidates[0],
                        sulfur_world,
                        geometry_metrics,
                    )

            if best_fit is None:
                carbon_world = self._world_atom_position(
                    candidate,
                    aldehyde_ref,
                    aldehyde_spec,
                    carbon_atom_id,
                    atom_position_overrides=current_atom_position_overrides,
                )
                nitrogen_world = self._world_atom_position(
                    candidate,
                    amine_ref,
                    amine_spec,
                    nitrogen_atom_id,
                    atom_position_overrides=current_atom_position_overrides,
                )
                anchor_world = self._world_atom_position(
                    candidate,
                    amine_ref,
                    amine_spec,
                    anchor_atom_id,
                    atom_position_overrides=current_atom_position_overrides,
                )
                ortho_world = self._world_atom_position(
                    candidate,
                    amine_ref,
                    amine_spec,
                    ortho_carbon_atom_id,
                    atom_position_overrides=current_atom_position_overrides,
                )
                sulfur_world, geometry_metrics = self._benzothiazole_sulfur_world_position(
                    carbon_world=carbon_world,
                    ortho_world=ortho_world,
                    nitrogen_world=nitrogen_world,
                    anchor_world=anchor_world,
                    blockers=tuple(
                        (instance_id, atom_id, world_position)
                        for instance_id, atom_id, world_position in blocker_atoms
                        if not (
                            (instance_id == aldehyde_instance_id and atom_id == carbon_atom_id)
                            or (instance_id == amine_instance_id and atom_id == ortho_carbon_atom_id)
                        )
                    ),
                )
                selected_ortho_carbon_atom_id = ortho_carbon_atom_id
                selected_ortho_hydrogen_atom_id = ortho_hydrogen_atom_id
                geometry_metrics = self._benzothiazole_geometry_metrics(
                    candidate=candidate,
                    amine_ref=amine_ref,
                    amine_spec=amine_spec,
                    nitrogen_atom_id=nitrogen_atom_id,
                    anchor_atom_id=anchor_atom_id,
                    ortho_carbon_atom_id=selected_ortho_carbon_atom_id,
                    aldehyde_ref=aldehyde_ref,
                    aldehyde_spec=aldehyde_spec,
                    carbon_atom_id=carbon_atom_id,
                    aldehyde_anchor_atom_id=aldehyde_anchor_atom_id,
                    atom_position_overrides=current_atom_position_overrides,
                    world_overrides={},
                    sulfur_world=sulfur_world,
                    blockers=tuple(
                        (instance_id, atom_id, world_position)
                        for instance_id, atom_id, world_position in blocker_atoms
                        if (instance_id, atom_id)
                        not in {
                            (aldehyde_instance_id, carbon_atom_id),
                            (aldehyde_instance_id, aldehyde_anchor_atom_id),
                            (amine_instance_id, nitrogen_atom_id),
                            (amine_instance_id, anchor_atom_id),
                            (amine_instance_id, selected_ortho_carbon_atom_id),
                        }
                    ),
                    fit_objective=float("inf"),
                    local_relaxation_applied=False,
                )
            else:
                (
                    selected_ortho_carbon_atom_id,
                    selected_ortho_hydrogen_atom_id,
                    sulfur_world,
                    geometry_metrics,
                ) = best_fit
                carbon_world = self._world_atom_position(
                    candidate,
                    aldehyde_ref,
                    aldehyde_spec,
                    carbon_atom_id,
                    atom_position_overrides=current_atom_position_overrides,
                )
                nitrogen_world = self._world_atom_position(
                    candidate,
                    amine_ref,
                    amine_spec,
                    nitrogen_atom_id,
                    atom_position_overrides=current_atom_position_overrides,
                )

            sulfur_atom_id = len(amine_spec.atom_symbols) + len(added_atoms[amine_instance_id])
            sulfur_label = self.atom_label(amine_instance_id, "S", sulfur_atom_id)
            sulfur_local = self._local_position_from_world(
                candidate.state.monomer_poses[amine_instance_id],
                sulfur_world,
                image_shift=self._periodic_offset(candidate.state.cell, amine_ref.periodic_image),
            )
            added_atoms[amine_instance_id].append(
                RealizedAtom(
                    atom_id=sulfur_atom_id,
                    label=sulfur_label,
                    symbol="S",
                    local_position=sulfur_local,
                )
            )

            carbon_label = self.atom_label(
                aldehyde_instance_id,
                aldehyde_spec.atom_symbols[carbon_atom_id],
                carbon_atom_id,
            )
            ortho_label = self.atom_label(
                amine_instance_id,
                amine_spec.atom_symbols[selected_ortho_carbon_atom_id],
                selected_ortho_carbon_atom_id,
            )
            bonds.extend(
                (
                    RealizedBond(
                        label_1=sulfur_label,
                        label_2=carbon_label,
                        distance=self._distance(sulfur_world, carbon_world),
                        symmetry_1=".",
                        symmetry_2=self._symmetry_code(aldehyde_ref.periodic_image),
                    ),
                    RealizedBond(
                        label_1=sulfur_label,
                        label_2=ortho_label,
                        distance=geometry_metrics["actual_ortho_s_distance"],
                    ),
                )
            )
            converted_removed_atom_ids[aldehyde_instance_id].add(aldehydic_hydrogen_atom_id)
            converted_removed_atom_ids[amine_instance_id].add(selected_ortho_hydrogen_atom_id)
            used_ortho_sites.add((amine_instance_id, selected_ortho_carbon_atom_id))
            applied_event_ids.append(event.id)
            sulfur_geometry[event.id] = geometry_metrics

        if not applied_event_ids:
            return None
        notes.extend(
            (
                "Requested sulfur-enabled imine conversion applies a benzothiazole-like annulation during atomistic realization only.",
                "The built candidate graph remains imine-linked; the exported atomistic structure adds one external sulfur atom per converted imine event.",
            )
        )
        return PostBuildRealization(
            removed_atom_ids={
                instance_id: tuple(sorted(atom_ids))
                for instance_id, atom_ids in converted_removed_atom_ids.items()
            },
            atom_position_overrides={
                instance_id: dict(overrides)
                for instance_id, overrides in fitted_atom_position_overrides.items()
            },
            added_atoms={
                instance_id: tuple(atoms)
                for instance_id, atoms in added_atoms.items()
            },
            bonds=tuple(bonds),
            notes=tuple(notes),
            metadata={
                "profile_id": "sulfur_enabled_imine_conversion",
                "applied_event_ids": tuple(applied_event_ids),
                "n_added_sulfur_atoms": sum(len(atoms) for atoms in added_atoms.values()),
                "sulfur_geometry": sulfur_geometry,
            },
        )

    def _added_heavy_atom_world_positions(
        self,
        candidate: Candidate,
        added_atoms_by_instance: Mapping[str, list[RealizedAtom]],
    ) -> tuple[tuple[str, int, Vec3], ...]:
        heavy_atoms: list[tuple[str, int, Vec3]] = []
        for instance_id, atoms in added_atoms_by_instance.items():
            pose = candidate.state.monomer_poses[instance_id]
            for atom in atoms:
                if atom.symbol == "H":
                    continue
                heavy_atoms.append((instance_id, atom.atom_id, self._world_position(pose, atom.local_position)))
        return tuple(heavy_atoms)

    def _benzothiazole_sulfur_world_position(
        self,
        *,
        carbon_world: Vec3,
        ortho_world: Vec3,
        nitrogen_world: Vec3,
        anchor_world: Vec3,
        blockers: tuple[tuple[str, int, Vec3], ...],
    ) -> tuple[Vec3, dict[str, float]]:
        carbon_s_distance = 1.90
        ortho_s_distance = 1.90
        nitrogen_s_distance = 2.35
        anchor_s_distance = 2.50
        clearance_target = 1.08
        baseline = sub(ortho_world, carbon_world)
        baseline_length = norm(baseline)
        if baseline_length < 1e-8:
            fallback = add(carbon_world, (carbon_s_distance, 0.0, 0.0))
            return fallback, {
                "carbon_s_distance": carbon_s_distance,
                "ortho_s_distance": self._distance(fallback, ortho_world),
                "min_nonbonded_clearance": self._minimum_distance(
                    fallback,
                    tuple(world for _instance_id, _atom_id, world in blockers),
                )
                or 0.0,
                "placement_offset": carbon_s_distance,
            }

        midpoint = scale(add(carbon_world, ortho_world), 0.5)
        baseline_unit = scale(baseline, 1.0 / baseline_length)
        bend_reference = sub(scale(add(nitrogen_world, anchor_world), 0.5), midpoint)
        side_vector = sub(bend_reference, scale(baseline_unit, dot(bend_reference, baseline_unit)))
        if norm(side_vector) < 1e-8:
            side_vector = cross(sub(anchor_world, carbon_world), baseline_unit)
        if norm(side_vector) < 1e-8:
            side_vector = cross(sub(nitrogen_world, carbon_world), baseline_unit)
        if norm(side_vector) < 1e-8:
            side_vector = self._any_orthogonal_unit(baseline_unit) or (0.0, 0.0, 1.0)
        opposite_bend_unit = scale(normalize(side_vector), -1.0)
        blocker_positions = tuple(world for _instance_id, _atom_id, world in blockers)

        best: tuple[float, Vec3, dict[str, float], float, float] | None = None
        best_x = 0.0
        best_y = max(
            0.40,
            (carbon_s_distance * carbon_s_distance - (0.5 * baseline_length) ** 2) ** 0.5
            if baseline_length < 2.0 * carbon_s_distance
            else 0.75,
        )
        search_passes = (
            (1.10, 0.15, 2.80, 0.12),
            (0.30, 0.05, 0.60, 0.03),
            (0.08, 0.05, 0.16, 0.01),
        )

        for span_x, min_y, span_y, step in search_passes:
            x_center = best_x
            y_center = max(min_y, best_y)
            x_start = x_center - span_x
            x_stop = x_center + span_x
            y_start = max(0.05, y_center - span_y)
            y_stop = y_center + span_y
            x_steps = max(1, int(round((x_stop - x_start) / step)))
            y_steps = max(1, int(round((y_stop - y_start) / step)))
            for x_index in range(x_steps + 1):
                x_offset = x_start + step * x_index
                for y_index in range(y_steps + 1):
                    y_offset = y_start + step * y_index
                    candidate_world = add(
                        midpoint,
                        add(scale(baseline_unit, x_offset), scale(opposite_bend_unit, y_offset)),
                    )
                    carbon_distance = self._distance(candidate_world, carbon_world)
                    ortho_distance = self._distance(candidate_world, ortho_world)
                    clearance = self._minimum_distance(candidate_world, blocker_positions) or 0.0
                    nitrogen_sulfur_distance = self._distance(candidate_world, nitrogen_world)
                    sulfur_anchor_distance = self._distance(candidate_world, anchor_world)
                    n_c_s_angle = self._angle(nitrogen_world, carbon_world, candidate_world)
                    c_s_ortho_angle = self._angle(carbon_world, candidate_world, ortho_world)
                    objective = 6.0 * (carbon_distance - carbon_s_distance) ** 2
                    objective += 6.0 * (ortho_distance - ortho_s_distance) ** 2
                    objective += 4.0 * (nitrogen_sulfur_distance - nitrogen_s_distance) ** 2
                    objective += 1.2 * (sulfur_anchor_distance - anchor_s_distance) ** 2
                    objective += 2.0 * (carbon_distance - ortho_distance) ** 2
                    objective += 0.04 * (n_c_s_angle - 100.0) ** 2
                    objective += 0.04 * (c_s_ortho_angle - 100.0) ** 2
                    if nitrogen_sulfur_distance < 2.10:
                        objective += 12.0 * (2.10 - nitrogen_sulfur_distance) ** 2
                    if sulfur_anchor_distance < 2.20:
                        objective += 6.0 * (2.20 - sulfur_anchor_distance) ** 2
                    if clearance < clearance_target:
                        objective += 10.0 * (clearance_target - clearance) ** 2
                    metrics = {
                        "carbon_s_distance": carbon_distance,
                        "ortho_s_distance": ortho_distance,
                        "min_nonbonded_clearance": clearance,
                        "placement_offset": self._distance(candidate_world, midpoint),
                        "n_c_s_angle": n_c_s_angle,
                        "c_s_ortho_angle": c_s_ortho_angle,
                        "nitrogen_sulfur_distance": nitrogen_sulfur_distance,
                        "sulfur_anchor_distance": sulfur_anchor_distance,
                        "fit_objective": objective,
                        "sulfur_opposite_bend": 1.0,
                    }
                    if best is None or objective < best[0]:
                        best = (objective, candidate_world, metrics, x_offset, y_offset)
                        best_x = x_offset
                        best_y = y_offset

        if best is not None:
            return best[1], best[2]

        fallback = add(midpoint, scale(opposite_bend_unit, carbon_s_distance))
        return fallback, {
            "carbon_s_distance": self._distance(fallback, carbon_world),
            "ortho_s_distance": self._distance(fallback, ortho_world),
            "min_nonbonded_clearance": self._minimum_distance(fallback, blocker_positions) or 0.0,
            "placement_offset": self._distance(fallback, midpoint),
            "n_c_s_angle": self._angle(nitrogen_world, carbon_world, fallback),
            "c_s_ortho_angle": self._angle(carbon_world, fallback, ortho_world),
            "nitrogen_sulfur_distance": self._distance(fallback, nitrogen_world),
            "sulfur_anchor_distance": self._distance(fallback, anchor_world),
            "fit_objective": float("inf"),
            "sulfur_opposite_bend": 1.0,
        }

    def _merge_removed_atom_ids(
        self,
        base_removed_atom_ids: Mapping[str, set[int]],
        extra_removed_atom_ids: Mapping[str, set[int]],
    ) -> dict[str, set[int]]:
        merged: dict[str, set[int]] = defaultdict(set)
        for source in (base_removed_atom_ids, extra_removed_atom_ids):
            for instance_id, atom_ids in source.items():
                merged[instance_id].update(atom_ids)
        return merged

    def _merge_atom_position_overrides(
        self,
        base_overrides: Mapping[str, Mapping[int, Vec3]],
        extra_overrides: Mapping[str, Mapping[int, Vec3]],
    ) -> dict[str, dict[int, Vec3]]:
        merged: dict[str, dict[int, Vec3]] = defaultdict(dict)
        for source in (base_overrides, extra_overrides):
            for instance_id, overrides in source.items():
                merged[instance_id].update(overrides)
        return merged

    def _ortho_option_ids(
        self,
        raw_plan: Mapping[str, object],
        *,
        default_ortho_atom_id: int,
    ) -> tuple[int, ...]:
        raw_candidate_ids = raw_plan.get("ortho_candidate_atom_ids")
        candidate_ids: list[int] = []
        if isinstance(raw_candidate_ids, (list, tuple)):
            candidate_ids.extend(
                atom_id
                for atom_id in (self._coerce_int(value) for value in raw_candidate_ids)
                if atom_id is not None
            )
        selected_atom_id = self._coerce_int(raw_plan.get("ortho_carbon_atom_id"))
        if selected_atom_id is not None:
            candidate_ids.append(selected_atom_id)
        candidate_ids.append(default_ortho_atom_id)
        deduped: list[int] = []
        seen: set[int] = set()
        for atom_id in candidate_ids:
            if atom_id in seen:
                continue
            seen.add(atom_id)
            deduped.append(atom_id)
        return tuple(deduped)

    def _fit_benzothiazole_annulation_positions(
        self,
        *,
        candidate: Candidate,
        amine_ref: MotifRef,
        amine_spec: MonomerSpec,
        amine_motif,
        nitrogen_atom_id: int,
        anchor_atom_id: int,
        ortho_carbon_atom_id: int,
        aldehyde_ref: MotifRef,
        aldehyde_spec: MonomerSpec,
        aldehyde_motif,
        carbon_atom_id: int,
        aldehyde_anchor_atom_id: int,
        oxygen_atom_id: int,
        atom_position_overrides: Mapping[str, Mapping[int, Vec3]],
        blockers: tuple[tuple[str, int, Vec3], ...],
    ) -> tuple[dict[str, dict[int, Vec3]], Vec3, dict[str, float]] | None:
        anchor_world = self._world_atom_position(
            candidate,
            amine_ref,
            amine_spec,
            anchor_atom_id,
            atom_position_overrides=atom_position_overrides,
        )
        ortho_world = self._world_atom_position(
            candidate,
            amine_ref,
            amine_spec,
            ortho_carbon_atom_id,
            atom_position_overrides=atom_position_overrides,
        )
        nitrogen_world = self._world_atom_position(
            candidate,
            amine_ref,
            amine_spec,
            nitrogen_atom_id,
            atom_position_overrides=atom_position_overrides,
        )
        carbon_world = self._world_atom_position(
            candidate,
            aldehyde_ref,
            aldehyde_spec,
            carbon_atom_id,
            atom_position_overrides=atom_position_overrides,
        )
        aldehyde_anchor_world = self._world_atom_position(
            candidate,
            aldehyde_ref,
            aldehyde_spec,
            aldehyde_anchor_atom_id,
            atom_position_overrides=atom_position_overrides,
        )
        oxygen_world = self._world_atom_position(candidate, aldehyde_ref, aldehyde_spec, oxygen_atom_id)

        plane_normal = add(
            cross(sub(ortho_world, anchor_world), sub(nitrogen_world, anchor_world)),
            cross(sub(carbon_world, aldehyde_anchor_world), sub(nitrogen_world, carbon_world)),
        )
        if norm(plane_normal) < 1e-8:
            plane_normal = cross(sub(aldehyde_anchor_world, anchor_world), sub(ortho_world, anchor_world))
        if norm(plane_normal) < 1e-8:
            plane_normal = cross(sub(nitrogen_world, carbon_world), sub(ortho_world, anchor_world))
        if norm(plane_normal) < 1e-8:
            return None
        plane_normal = normalize(plane_normal)

        basis_u = self._orthogonal_component(sub(ortho_world, anchor_world), plane_normal)
        if norm(basis_u) < 1e-8:
            basis_u = self._orthogonal_component(sub(aldehyde_anchor_world, anchor_world), plane_normal)
        if norm(basis_u) < 1e-8:
            basis_u = self._orthogonal_component(sub(nitrogen_world, anchor_world), plane_normal)
        if norm(basis_u) < 1e-8:
            return None
        basis_u = normalize(basis_u)
        basis_v = normalize(cross(plane_normal, basis_u))

        def to_plane(point: Vec3) -> tuple[float, float]:
            relative = sub(point, anchor_world)
            return (dot(relative, basis_u), dot(relative, basis_v))

        def to_world(point_2d: tuple[float, float]) -> Vec3:
            return add(anchor_world, add(scale(basis_u, point_2d[0]), scale(basis_v, point_2d[1])))

        anchor_2d = (0.0, 0.0)
        ortho_2d = to_plane(ortho_world)
        nitrogen_current_2d = to_plane(nitrogen_world)
        carbon_current_2d = to_plane(carbon_world)
        aldehyde_anchor_2d = to_plane(aldehyde_anchor_world)

        target_anchor_n = self._distance(
            amine_spec.atom_positions[anchor_atom_id],
            amine_spec.atom_positions[nitrogen_atom_id],
        )
        target_aldehyde_anchor_c = self._distance(
            aldehyde_spec.atom_positions[aldehyde_anchor_atom_id],
            aldehyde_spec.atom_positions[carbon_atom_id],
        )
        target_n_c = 1.33
        target_c_s = 1.76
        target_ortho_s = 1.74
        target_c_s_ortho_angle = 92.0
        target_c_ortho_chord = (
            target_c_s * target_c_s
            + target_ortho_s * target_ortho_s
            - 2.0 * target_c_s * target_ortho_s * cos(pi * target_c_s_ortho_angle / 180.0)
        ) ** 0.5
        target_anchor_angle = self._angle(
            amine_spec.atom_positions[ortho_carbon_atom_id],
            amine_spec.atom_positions[anchor_atom_id],
            amine_spec.atom_positions[nitrogen_atom_id],
        )
        target_carbon_angle = self._angle(
            aldehyde_spec.atom_positions[aldehyde_anchor_atom_id],
            aldehyde_spec.atom_positions[carbon_atom_id],
            aldehyde_spec.atom_positions[oxygen_atom_id],
        )

        sulfur_world_candidate, sulfur_metrics = self._benzothiazole_sulfur_world_position(
            carbon_world=carbon_world,
            ortho_world=ortho_world,
            nitrogen_world=nitrogen_world,
            anchor_world=anchor_world,
            blockers=blockers,
        )
        world_overrides: dict[str, dict[int, Vec3]] = {}
        metrics = self._benzothiazole_geometry_metrics(
            candidate=candidate,
            amine_ref=amine_ref,
            amine_spec=amine_spec,
            nitrogen_atom_id=nitrogen_atom_id,
            anchor_atom_id=anchor_atom_id,
            ortho_carbon_atom_id=ortho_carbon_atom_id,
            aldehyde_ref=aldehyde_ref,
            aldehyde_spec=aldehyde_spec,
            carbon_atom_id=carbon_atom_id,
            aldehyde_anchor_atom_id=aldehyde_anchor_atom_id,
            atom_position_overrides=atom_position_overrides,
            world_overrides=world_overrides,
            sulfur_world=sulfur_world_candidate,
            blockers=blockers,
            fit_objective=sulfur_metrics["fit_objective"],
            local_relaxation_applied=False,
        )
        metrics["placement_offset"] = sulfur_metrics["placement_offset"]
        return world_overrides, sulfur_world_candidate, metrics

    def _relax_benzothiazole_local_geometry(
        self,
        *,
        candidate: Candidate,
        amine_ref: MotifRef,
        amine_spec: MonomerSpec,
        nitrogen_atom_id: int,
        anchor_atom_id: int,
        ortho_carbon_atom_id: int,
        aldehyde_ref: MotifRef,
        aldehyde_spec: MonomerSpec,
        carbon_atom_id: int,
        aldehyde_anchor_atom_id: int,
        atom_position_overrides: Mapping[str, Mapping[int, Vec3]],
        initial_world_overrides: Mapping[str, Mapping[int, Vec3]],
        sulfur_world: Vec3,
        blockers: tuple[tuple[str, int, Vec3], ...],
    ) -> tuple[dict[str, dict[int, Vec3]], Vec3] | None:
        anchor_world = self._world_atom_position(
            candidate,
            amine_ref,
            amine_spec,
            anchor_atom_id,
            atom_position_overrides=atom_position_overrides,
        )
        ortho_world = initial_world_overrides.get(amine_ref.monomer_instance_id, {}).get(
            ortho_carbon_atom_id,
            self._world_atom_position(
                candidate,
                amine_ref,
                amine_spec,
                ortho_carbon_atom_id,
                atom_position_overrides=atom_position_overrides,
            ),
        )
        nitrogen_world = initial_world_overrides.get(amine_ref.monomer_instance_id, {}).get(
            nitrogen_atom_id,
            self._world_atom_position(
                candidate,
                amine_ref,
                amine_spec,
                nitrogen_atom_id,
                atom_position_overrides=atom_position_overrides,
            ),
        )
        carbon_world = initial_world_overrides.get(aldehyde_ref.monomer_instance_id, {}).get(
            carbon_atom_id,
            self._world_atom_position(
                candidate,
                aldehyde_ref,
                aldehyde_spec,
                carbon_atom_id,
                atom_position_overrides=atom_position_overrides,
            ),
        )
        aldehyde_anchor_world = initial_world_overrides.get(aldehyde_ref.monomer_instance_id, {}).get(
            aldehyde_anchor_atom_id,
            self._world_atom_position(
                candidate,
                aldehyde_ref,
                aldehyde_spec,
                aldehyde_anchor_atom_id,
                atom_position_overrides=atom_position_overrides,
            ),
        )
        plane_normal = add(
            cross(sub(ortho_world, anchor_world), sub(nitrogen_world, anchor_world)),
            cross(sub(carbon_world, aldehyde_anchor_world), sub(nitrogen_world, carbon_world)),
        )
        if norm(plane_normal) < 1e-8:
            plane_normal = cross(sub(aldehyde_anchor_world, anchor_world), sub(ortho_world, anchor_world))
        if norm(plane_normal) < 1e-8:
            return None
        plane_normal = normalize(plane_normal)
        plane_origin = anchor_world

        def world_position(
            ref: MotifRef,
            monomer: MonomerSpec,
            atom_id: int,
        ) -> Vec3:
            override = initial_world_overrides.get(ref.monomer_instance_id, {}).get(atom_id)
            if override is not None:
                return override
            return self._world_atom_position(
                candidate,
                ref,
                monomer,
                atom_id,
                atom_position_overrides=atom_position_overrides,
            )

        positions: dict[object, Vec3] = {}
        movable_keys: set[object] = set()
        fixed_keys: set[object] = set()
        sulfur_key: object = ("added", amine_ref.monomer_instance_id, "S")

        def monomer_key(instance_id: str, atom_id: int) -> tuple[str, int]:
            return (instance_id, atom_id)

        def add_movable(ref: MotifRef, monomer: MonomerSpec, atom_id: int) -> None:
            key = monomer_key(ref.monomer_instance_id, atom_id)
            if key not in positions:
                positions[key] = world_position(ref, monomer, atom_id)
            movable_keys.add(key)

        def add_fixed(ref: MotifRef, monomer: MonomerSpec, atom_id: int) -> None:
            key = monomer_key(ref.monomer_instance_id, atom_id)
            if key not in positions:
                positions[key] = world_position(ref, monomer, atom_id)
            fixed_keys.add(key)

        add_movable(amine_ref, amine_spec, nitrogen_atom_id)
        add_movable(amine_ref, amine_spec, ortho_carbon_atom_id)
        add_fixed(amine_ref, amine_spec, anchor_atom_id)
        for atom_id in self._bonded_heavy_neighbor_atom_ids(amine_spec, ortho_carbon_atom_id, set()):
            if atom_id == anchor_atom_id:
                continue
            add_movable(amine_ref, amine_spec, atom_id)
        add_movable(aldehyde_ref, aldehyde_spec, carbon_atom_id)
        add_movable(aldehyde_ref, aldehyde_spec, aldehyde_anchor_atom_id)
        for atom_id in self._bonded_heavy_neighbor_atom_ids(aldehyde_spec, aldehyde_anchor_atom_id, set()):
            if atom_id == carbon_atom_id:
                continue
            add_movable(aldehyde_ref, aldehyde_spec, atom_id)

        for atom_key in tuple(movable_keys):
            instance_id, atom_id = atom_key
            if not isinstance(atom_id, int):
                continue
            ref = amine_ref if instance_id == amine_ref.monomer_instance_id else aldehyde_ref
            monomer = amine_spec if instance_id == amine_ref.monomer_instance_id else aldehyde_spec
            for neighbor_atom_id in self._bonded_heavy_neighbor_atom_ids(monomer, atom_id, set()):
                neighbor_key = monomer_key(instance_id, neighbor_atom_id)
                if neighbor_key in movable_keys or neighbor_key in fixed_keys:
                    continue
                add_fixed(ref, monomer, neighbor_atom_id)

        positions[sulfur_key] = sulfur_world
        movable_keys.add(sulfur_key)
        initial_positions = dict(positions)

        springs: list[tuple[object, object, float, float]] = []
        for ref, monomer in ((amine_ref, amine_spec), (aldehyde_ref, aldehyde_spec)):
            instance_id = ref.monomer_instance_id
            for atom_id_1, atom_id_2, _bond_order in monomer.bonds:
                key_1 = monomer_key(instance_id, atom_id_1)
                key_2 = monomer_key(instance_id, atom_id_2)
                if key_1 not in positions or key_2 not in positions:
                    continue
                springs.append(
                    (
                        key_1,
                        key_2,
                        self._distance(monomer.atom_positions[atom_id_1], monomer.atom_positions[atom_id_2]),
                        1.0,
                    )
                )
        springs.extend(
            (
                (monomer_key(amine_ref.monomer_instance_id, nitrogen_atom_id), monomer_key(aldehyde_ref.monomer_instance_id, carbon_atom_id), 1.33, 1.6),
                (monomer_key(aldehyde_ref.monomer_instance_id, carbon_atom_id), sulfur_key, 1.76, 1.6),
                (sulfur_key, monomer_key(amine_ref.monomer_instance_id, ortho_carbon_atom_id), 1.74, 1.6),
                (monomer_key(amine_ref.monomer_instance_id, nitrogen_atom_id), sulfur_key, 2.45, 0.4),
                (monomer_key(aldehyde_ref.monomer_instance_id, carbon_atom_id), monomer_key(amine_ref.monomer_instance_id, ortho_carbon_atom_id), 2.75, 0.5),
            )
        )
        for key in tuple(movable_keys):
            rest_key = ("rest", key)
            positions[rest_key] = initial_positions[key]
            springs.append((key, rest_key, 0.0, 0.03))

        for _iteration in range(180):
            for key_1, key_2, target_distance, weight in springs:
                world_1 = positions[key_1]
                world_2 = positions[key_2]
                delta = sub(world_2, world_1)
                distance = norm(delta)
                if distance < 1e-8:
                    continue
                correction = scale(delta, 0.5 * weight * (distance - target_distance) / distance)
                key_1_movable = key_1 in movable_keys
                key_2_movable = key_2 in movable_keys
                if key_1_movable and key_2_movable:
                    positions[key_1] = self._project_onto_plane(add(world_1, correction), plane_origin, plane_normal)
                    positions[key_2] = self._project_onto_plane(sub(world_2, correction), plane_origin, plane_normal)
                elif key_1_movable:
                    positions[key_1] = self._project_onto_plane(
                        add(world_1, scale(delta, weight * (distance - target_distance) / distance)),
                        plane_origin,
                        plane_normal,
                    )
                elif key_2_movable:
                    positions[key_2] = self._project_onto_plane(
                        sub(world_2, scale(delta, weight * (distance - target_distance) / distance)),
                        plane_origin,
                        plane_normal,
                    )

        relaxed_world_overrides: dict[str, dict[int, Vec3]] = defaultdict(dict)
        for key in movable_keys:
            if key == sulfur_key:
                continue
            instance_id, atom_id = key
            assert isinstance(instance_id, str)
            assert isinstance(atom_id, int)
            world_position_value = positions[key]
            initial_world = initial_positions[key]
            if self._distance(world_position_value, initial_world) <= 1e-4:
                continue
            relaxed_world_overrides[instance_id][atom_id] = world_position_value
        relaxed_sulfur_world = positions[sulfur_key]
        local_blocker_keys = {
            key
            for key in positions
            if isinstance(key, tuple)
            and len(key) == 2
            and isinstance(key[0], str)
            and isinstance(key[1], int)
        }
        sulfur_clearance = self._minimum_distance(
            relaxed_sulfur_world,
            tuple(
                world_position
                for instance_id, atom_id, world_position in blockers
                if (instance_id, atom_id) not in local_blocker_keys
            ),
        ) or 0.0
        actual_anchor_c_distance = self._distance(
            positions[monomer_key(aldehyde_ref.monomer_instance_id, aldehyde_anchor_atom_id)],
            positions[monomer_key(aldehyde_ref.monomer_instance_id, carbon_atom_id)],
        )
        if sulfur_clearance < 0.95 or actual_anchor_c_distance > 1.65:
            return None
        return (
            {
                instance_id: dict(overrides)
                for instance_id, overrides in relaxed_world_overrides.items()
            },
            relaxed_sulfur_world,
        )

    def _benzothiazole_geometry_metrics(
        self,
        *,
        candidate: Candidate,
        amine_ref: MotifRef,
        amine_spec: MonomerSpec,
        nitrogen_atom_id: int,
        anchor_atom_id: int,
        ortho_carbon_atom_id: int,
        aldehyde_ref: MotifRef,
        aldehyde_spec: MonomerSpec,
        carbon_atom_id: int,
        aldehyde_anchor_atom_id: int,
        atom_position_overrides: Mapping[str, Mapping[int, Vec3]],
        world_overrides: Mapping[str, Mapping[int, Vec3]],
        sulfur_world: Vec3,
        blockers: tuple[tuple[str, int, Vec3], ...],
        fit_objective: float,
        local_relaxation_applied: bool,
    ) -> dict[str, float]:
        def world_position(ref: MotifRef, monomer: MonomerSpec, atom_id: int) -> Vec3:
            override = world_overrides.get(ref.monomer_instance_id, {}).get(atom_id)
            if override is not None:
                return override
            return self._world_atom_position(
                candidate,
                ref,
                monomer,
                atom_id,
                atom_position_overrides=atom_position_overrides,
            )

        anchor_world = world_position(amine_ref, amine_spec, anchor_atom_id)
        ortho_world = world_position(amine_ref, amine_spec, ortho_carbon_atom_id)
        nitrogen_world = world_position(amine_ref, amine_spec, nitrogen_atom_id)
        carbon_world = world_position(aldehyde_ref, aldehyde_spec, carbon_atom_id)
        aldehyde_anchor_world = world_position(aldehyde_ref, aldehyde_spec, aldehyde_anchor_atom_id)
        midpoint = scale(add(carbon_world, ortho_world), 0.5)
        bend_reference = sub(scale(add(nitrogen_world, anchor_world), 0.5), midpoint)
        sulfur_direction = sub(sulfur_world, midpoint)
        sulfur_same_side_as_bend = dot(bend_reference, sulfur_direction) > 0.0
        sulfur_clearance = self._minimum_distance(
            sulfur_world,
            tuple(world_position_value for _instance_id, _atom_id, world_position_value in blockers),
        ) or 0.0
        return {
            "fit_objective": fit_objective,
            "actual_anchor_n_distance": self._distance(anchor_world, nitrogen_world),
            "actual_n_c_distance": self._distance(nitrogen_world, carbon_world),
            "actual_aldehyde_anchor_c_distance": self._distance(aldehyde_anchor_world, carbon_world),
            "actual_carbon_s_distance": self._distance(carbon_world, sulfur_world),
            "actual_ortho_s_distance": self._distance(ortho_world, sulfur_world),
            "actual_c_ortho_distance": self._distance(carbon_world, ortho_world),
            "nitrogen_sulfur_distance": self._distance(nitrogen_world, sulfur_world),
            "sulfur_anchor_distance": self._distance(sulfur_world, anchor_world),
            "sulfur_clearance": sulfur_clearance,
            "anchor_angle": self._angle(ortho_world, anchor_world, nitrogen_world),
            "carbon_angle": self._angle(aldehyde_anchor_world, carbon_world, nitrogen_world),
            "n_c_s_angle": self._angle(nitrogen_world, carbon_world, sulfur_world),
            "c_s_ortho_angle": self._angle(carbon_world, sulfur_world, ortho_world),
            "selected_ortho_carbon_atom_id": float(ortho_carbon_atom_id),
            "sulfur_opposite_bend": 0.0 if sulfur_same_side_as_bend else 1.0,
            "local_relaxation_applied": 1.0 if local_relaxation_applied else 0.0,
        }

    def _project_onto_plane(self, point: Vec3, plane_origin: Vec3, plane_normal: Vec3) -> Vec3:
        offset = sub(point, plane_origin)
        return sub(point, scale(plane_normal, dot(offset, plane_normal)))

    def _distance_2d(self, first: tuple[float, float], second: tuple[float, float]) -> float:
        return self._squared_distance_2d(first, second) ** 0.5

    def _local_position_from_world(
        self,
        pose: Pose,
        world_position: Vec3,
        *,
        image_shift: Vec3 = (0.0, 0.0, 0.0),
    ) -> Vec3:
        return matmul_vec(transpose(pose.rotation_matrix), sub(sub(world_position, image_shift), pose.translation))

    def _fit_imine_bridge_positions(
        self,
        *,
        candidate: Candidate,
        amine_ref: MotifRef,
        amine_spec: MonomerSpec,
        amine_motif,
        nitrogen_world: Vec3,
        amine_anchor_world: Vec3,
        aldehyde_ref: MotifRef,
        aldehyde_spec: MonomerSpec,
        aldehyde_motif,
        carbon_world: Vec3,
        aldehyde_anchor_world: Vec3,
        oxygen_world: Vec3,
        target_cn_distance: float,
    ) -> tuple[Vec3, Vec3] | None:
        bridge_axis = sub(amine_anchor_world, aldehyde_anchor_world)
        if norm(bridge_axis) < 1e-8:
            return None
        axis = normalize(bridge_axis)

        amine_pose = candidate.state.monomer_poses[amine_ref.monomer_instance_id]
        aldehyde_pose = candidate.state.monomer_poses[aldehyde_ref.monomer_instance_id]
        amine_normal = matmul_vec(amine_pose.rotation_matrix, amine_motif.frame.normal)
        aldehyde_normal = matmul_vec(aldehyde_pose.rotation_matrix, aldehyde_motif.frame.normal)
        plane_normal = self._orthogonal_component(add(amine_normal, aldehyde_normal), axis)
        if norm(plane_normal) < 1e-8:
            plane_normal = self._orthogonal_component(cross(axis, sub(oxygen_world, carbon_world)), axis)
        if norm(plane_normal) < 1e-8:
            plane_normal = self._orthogonal_component(cross(axis, sub(nitrogen_world, carbon_world)), axis)
        if norm(plane_normal) < 1e-8:
            plane_normal = (0.0, 0.0, 1.0) if abs(axis[2]) < 0.9 else (1.0, 0.0, 0.0)
        plane_normal = normalize(plane_normal)
        lateral_axis = cross(plane_normal, axis)
        if norm(lateral_axis) < 1e-8:
            return None
        lateral_axis = normalize(lateral_axis)

        carbon_side = self._signed_value(dot(sub(oxygen_world, carbon_world), lateral_axis), default=1.0)
        nitrogen_side = -carbon_side
        anchor_distance = self._distance(aldehyde_anchor_world, amine_anchor_world)
        carbon_length = self._distance(aldehyde_anchor_world, carbon_world)
        nitrogen_length = self._distance(amine_anchor_world, nitrogen_world)
        if anchor_distance < 1e-8 or carbon_length < 1e-8 or nitrogen_length < 1e-8:
            return None

        current_carbon_2d = (
            dot(sub(carbon_world, aldehyde_anchor_world), axis),
            dot(sub(carbon_world, aldehyde_anchor_world), lateral_axis),
        )
        current_nitrogen_2d = (
            dot(sub(nitrogen_world, aldehyde_anchor_world), axis),
            dot(sub(nitrogen_world, aldehyde_anchor_world), lateral_axis),
        )
        target_carbon_angle = 127.2
        target_nitrogen_angle = 127.8
        best: tuple[float, tuple[float, float], tuple[float, float]] | None = None

        for step in range(721):
            alpha = pi * step / 720.0
            carbon_2d = (
                carbon_length * cos(alpha),
                carbon_side * carbon_length * sin(alpha),
            )
            intersections = self._circle_intersections_2d(
                carbon_2d,
                target_cn_distance,
                (anchor_distance, 0.0),
                nitrogen_length,
            )
            if not intersections:
                continue
            for nitrogen_2d in intersections:
                if nitrogen_side * nitrogen_2d[1] <= 1e-6:
                    continue
                carbon_angle = self._planar_angle_2d((0.0, 0.0), carbon_2d, nitrogen_2d)
                nitrogen_angle = self._planar_angle_2d(carbon_2d, nitrogen_2d, (anchor_distance, 0.0))
                objective = 0.12 * (carbon_angle - target_carbon_angle) ** 2
                objective += 0.12 * (nitrogen_angle - target_nitrogen_angle) ** 2
                objective += 0.5 * self._squared_distance_2d(carbon_2d, current_carbon_2d)
                objective += 0.5 * self._squared_distance_2d(nitrogen_2d, current_nitrogen_2d)
                if best is None or objective < best[0]:
                    best = (objective, carbon_2d, nitrogen_2d)

        if best is None:
            return None

        _, carbon_2d, nitrogen_2d = best
        return (
            add(
                aldehyde_anchor_world,
                add(scale(axis, carbon_2d[0]), scale(lateral_axis, carbon_2d[1])),
            ),
            add(
                aldehyde_anchor_world,
                add(scale(axis, nitrogen_2d[0]), scale(lateral_axis, nitrogen_2d[1])),
            ),
        )

    def _orthogonal_component(self, vector: Vec3, axis: Vec3) -> Vec3:
        return sub(vector, scale(axis, dot(vector, axis)))

    def _signed_value(self, value: float, *, default: float) -> float:
        if value > 1e-8:
            return 1.0
        if value < -1e-8:
            return -1.0
        return default

    def _circle_intersections_2d(
        self,
        center_a: tuple[float, float],
        radius_a: float,
        center_b: tuple[float, float],
        radius_b: float,
    ) -> tuple[tuple[float, float], ...]:
        delta_x = center_b[0] - center_a[0]
        delta_y = center_b[1] - center_a[1]
        distance = (delta_x * delta_x + delta_y * delta_y) ** 0.5
        if distance < 1e-8 or distance > radius_a + radius_b or distance < abs(radius_a - radius_b):
            return ()
        along = (radius_a * radius_a - radius_b * radius_b + distance * distance) / (2.0 * distance)
        height_sq = radius_a * radius_a - along * along
        if height_sq < -1e-8:
            return ()
        height = max(0.0, height_sq) ** 0.5
        unit_x = delta_x / distance
        unit_y = delta_y / distance
        base_x = center_a[0] + along * unit_x
        base_y = center_a[1] + along * unit_y
        perpendicular = (-unit_y, unit_x)
        first = (base_x + height * perpendicular[0], base_y + height * perpendicular[1])
        second = (base_x - height * perpendicular[0], base_y - height * perpendicular[1])
        if self._squared_distance_2d(first, second) <= 1e-12:
            return (first,)
        return (first, second)

    def _planar_angle_2d(
        self,
        point_a: tuple[float, float],
        point_b: tuple[float, float],
        point_c: tuple[float, float],
    ) -> float:
        vector_ba = (point_a[0] - point_b[0], point_a[1] - point_b[1])
        vector_bc = (point_c[0] - point_b[0], point_c[1] - point_b[1])
        norm_ba = (vector_ba[0] * vector_ba[0] + vector_ba[1] * vector_ba[1]) ** 0.5
        norm_bc = (vector_bc[0] * vector_bc[0] + vector_bc[1] * vector_bc[1]) ** 0.5
        if norm_ba < 1e-8 or norm_bc < 1e-8:
            return 180.0
        cosine = (vector_ba[0] * vector_bc[0] + vector_ba[1] * vector_bc[1]) / (norm_ba * norm_bc)
        return 180.0 * acos(max(-1.0, min(1.0, cosine))) / pi

    def _angle(self, point_a: Vec3, point_b: Vec3, point_c: Vec3) -> float:
        vector_ba = sub(point_a, point_b)
        vector_bc = sub(point_c, point_b)
        norm_ba = norm(vector_ba)
        norm_bc = norm(vector_bc)
        if norm_ba < 1e-8 or norm_bc < 1e-8:
            return 180.0
        cosine = dot(vector_ba, vector_bc) / (norm_ba * norm_bc)
        return 180.0 * acos(max(-1.0, min(1.0, cosine))) / pi

    def _squared_distance_2d(self, first: tuple[float, float], second: tuple[float, float]) -> float:
        delta_x = first[0] - second[0]
        delta_y = first[1] - second[1]
        return delta_x * delta_x + delta_y * delta_y

    def _coerce_int(self, value: object) -> int | None:
        if isinstance(value, int):
            return value
        if isinstance(value, float) and value.is_integer():
            return int(value)
        return None

    def _realize_imine_bridge(
        self,
        event: ReactionEvent,
        candidate: Candidate,
        monomer_specs: Mapping[str, MonomerSpec],
    ) -> EventRealization:
        amine_ref, aldehyde_ref = self._split_bridge_participants(event, monomer_specs, "amine", "aldehyde")
        amine_motif = monomer_specs[amine_ref.monomer_id].motif_by_id(amine_ref.motif_id)
        aldehyde_motif = monomer_specs[aldehyde_ref.monomer_id].motif_by_id(aldehyde_ref.motif_id)

        nitrogen_atom_id = self._reactive_atom_id(amine_motif)
        carbon_atom_id = self._reactive_atom_id(aldehyde_motif)
        oxygen_atom_id = self._unique_symbol_atom_id(
            monomer_specs[aldehyde_ref.monomer_id],
            aldehyde_motif.atom_ids,
            "O",
            context=f"{event.id} aldehyde oxygen",
        )
        hydrogen_atom_ids = self._select_hydrogens(
            candidate,
            amine_ref,
            monomer_specs[amine_ref.monomer_id],
            self._hydrogen_atom_ids_for_atom(monomer_specs[amine_ref.monomer_id], amine_motif, nitrogen_atom_id),
            target=self._world_atom_position(
                candidate,
                aldehyde_ref,
                monomer_specs[aldehyde_ref.monomer_id],
                carbon_atom_id,
            ),
            count=2,
            context="imine formation",
        )

        amine_spec = monomer_specs[amine_ref.monomer_id]
        aldehyde_spec = monomer_specs[aldehyde_ref.monomer_id]
        carbon_world = self._world_atom_position(candidate, aldehyde_ref, aldehyde_spec, carbon_atom_id)
        nitrogen_world = self._world_atom_position(candidate, amine_ref, amine_spec, nitrogen_atom_id)
        amine_anchor_atom_id = self._motif_atom_id_from_metadata(
            amine_motif,
            "anchor_atom_id",
            context=f"{event.id} imine amine anchor",
        )
        aldehyde_anchor_atom_id = self._motif_atom_id_from_metadata(
            aldehyde_motif,
            "anchor_atom_id",
            context=f"{event.id} imine aldehyde anchor",
        )
        amine_anchor_world = self._world_atom_position(candidate, amine_ref, amine_spec, amine_anchor_atom_id)
        aldehyde_anchor_world = self._world_atom_position(candidate, aldehyde_ref, aldehyde_spec, aldehyde_anchor_atom_id)
        oxygen_world = self._world_atom_position(candidate, aldehyde_ref, aldehyde_spec, oxygen_atom_id)
        target_cn_distance = bridge_target_distance("imine_bridge")
        closed_positions = self._fit_imine_bridge_positions(
            candidate=candidate,
            amine_ref=amine_ref,
            amine_spec=amine_spec,
            amine_motif=amine_motif,
            nitrogen_world=nitrogen_world,
            amine_anchor_world=amine_anchor_world,
            aldehyde_ref=aldehyde_ref,
            aldehyde_spec=aldehyde_spec,
            aldehyde_motif=aldehyde_motif,
            carbon_world=carbon_world,
            aldehyde_anchor_world=aldehyde_anchor_world,
            oxygen_world=oxygen_world,
            target_cn_distance=target_cn_distance,
        )
        atom_position_overrides: dict[str, dict[int, Vec3]] = {}
        if closed_positions is not None:
            carbon_world, nitrogen_world = closed_positions
            atom_position_overrides = {
                aldehyde_ref.monomer_instance_id: {
                    carbon_atom_id: self._local_position_from_world(
                        candidate.state.monomer_poses[aldehyde_ref.monomer_instance_id],
                        carbon_world,
                        image_shift=self._periodic_offset(candidate.state.cell, aldehyde_ref.periodic_image),
                    ),
                },
                amine_ref.monomer_instance_id: {
                    nitrogen_atom_id: self._local_position_from_world(
                        candidate.state.monomer_poses[amine_ref.monomer_instance_id],
                        nitrogen_world,
                        image_shift=self._periodic_offset(candidate.state.cell, amine_ref.periodic_image),
                    ),
                },
            }

        return EventRealization(
            removed_atom_ids={
                amine_ref.monomer_instance_id: tuple(sorted(hydrogen_atom_ids)),
                aldehyde_ref.monomer_instance_id: (oxygen_atom_id,),
            },
            atom_position_overrides=atom_position_overrides,
            bonds=(
                RealizedBond(
                    label_1=self.atom_label(aldehyde_ref.monomer_instance_id, aldehyde_spec.atom_symbols[carbon_atom_id], carbon_atom_id),
                    label_2=self.atom_label(amine_ref.monomer_instance_id, amine_spec.atom_symbols[nitrogen_atom_id], nitrogen_atom_id),
                    distance=self._distance(carbon_world, nitrogen_world),
                    symmetry_1=".",
                    symmetry_2=self._symmetry_code(amine_ref.periodic_image),
                ),
            ),
            notes=(
                "Imine realization removes the aldehyde oxygen and both amine hydrogens per reacting pair.",
                "The exported product applies a local imine chain-closure fit so the retained aryl-C/C=N/aryl-N segment is bent rather than left collinear.",
            ),
        )

    def _realize_hydrazone_bridge(
        self,
        event: ReactionEvent,
        candidate: Candidate,
        monomer_specs: Mapping[str, MonomerSpec],
    ) -> EventRealization:
        hydrazide_ref, aldehyde_ref = self._split_bridge_participants(event, monomer_specs, "hydrazide", "aldehyde")
        hydrazide_spec = monomer_specs[hydrazide_ref.monomer_id]
        aldehyde_spec = monomer_specs[aldehyde_ref.monomer_id]
        hydrazide_motif = hydrazide_spec.motif_by_id(hydrazide_ref.motif_id)
        aldehyde_motif = aldehyde_spec.motif_by_id(aldehyde_ref.motif_id)

        nitrogen_atom_id = self._reactive_atom_id(hydrazide_motif)
        carbon_atom_id = self._reactive_atom_id(aldehyde_motif)
        oxygen_atom_id = self._unique_symbol_atom_id(
            aldehyde_spec,
            aldehyde_motif.atom_ids,
            "O",
            context=f"{event.id} aldehyde oxygen",
        )
        carbon_world = self._world_atom_position(candidate, aldehyde_ref, aldehyde_spec, carbon_atom_id)
        hydrogen_atom_ids = self._select_hydrogens(
            candidate,
            hydrazide_ref,
            hydrazide_spec,
            self._hydrogen_atom_ids_for_atom(hydrazide_spec, hydrazide_motif, nitrogen_atom_id),
            target=carbon_world,
            count=2,
            context="hydrazone formation",
        )
        nitrogen_world = self._world_atom_position(candidate, hydrazide_ref, hydrazide_spec, nitrogen_atom_id)

        return EventRealization(
            removed_atom_ids={
                hydrazide_ref.monomer_instance_id: tuple(sorted(hydrogen_atom_ids)),
                aldehyde_ref.monomer_instance_id: (oxygen_atom_id,),
            },
            bonds=(
                RealizedBond(
                    label_1=self.atom_label(
                        aldehyde_ref.monomer_instance_id,
                        aldehyde_spec.atom_symbols[carbon_atom_id],
                        carbon_atom_id,
                    ),
                    label_2=self.atom_label(
                        hydrazide_ref.monomer_instance_id,
                        hydrazide_spec.atom_symbols[nitrogen_atom_id],
                        nitrogen_atom_id,
                    ),
                    distance=self._distance(carbon_world, nitrogen_world),
                    symmetry_1=".",
                    symmetry_2=self._symmetry_code(hydrazide_ref.periodic_image),
                ),
            ),
            notes=(
                "Hydrazone realization removes the aldehyde oxygen and both terminal hydrazide hydrogens per reacting pair.",
                "The hydrazide carbonyl is retained; the exported product realizes the inter-monomer C=N connectivity only.",
            ),
        )

    def _realize_keto_enamine_bridge(
        self,
        event: ReactionEvent,
        candidate: Candidate,
        monomer_specs: Mapping[str, MonomerSpec],
    ) -> EventRealization:
        amine_ref, keto_aldehyde_ref = self._split_bridge_participants(event, monomer_specs, "amine", "keto_aldehyde")
        amine_spec = monomer_specs[amine_ref.monomer_id]
        keto_aldehyde_spec = monomer_specs[keto_aldehyde_ref.monomer_id]
        amine_motif = amine_spec.motif_by_id(amine_ref.motif_id)
        keto_aldehyde_motif = keto_aldehyde_spec.motif_by_id(keto_aldehyde_ref.motif_id)

        nitrogen_atom_id = self._reactive_atom_id(amine_motif)
        carbon_atom_id = self._reactive_atom_id(keto_aldehyde_motif)
        oxygen_atom_id = self._motif_atom_id_from_metadata(
            keto_aldehyde_motif,
            "aldehyde_oxygen_atom_id",
            context=f"{event.id} keto-aldehyde oxygen",
        )
        tautomer_oxygen_atom_id = self._motif_atom_id_from_metadata(
            keto_aldehyde_motif,
            "ortho_hydroxyl_oxygen_atom_id",
            context=f"{event.id} keto-aldehyde ortho hydroxyl oxygen",
        )
        tautomer_hydrogen_atom_id = self._motif_atom_id_from_metadata(
            keto_aldehyde_motif,
            "ortho_hydroxyl_hydrogen_atom_id",
            context=f"{event.id} keto-aldehyde ortho hydroxyl hydrogen",
        )
        carbonyl_anchor_atom_id = self._keto_enamine_carbonyl_anchor_atom_id(
            keto_aldehyde_spec,
            keto_aldehyde_motif,
            tautomer_oxygen_atom_id,
        )
        carbon_world = self._world_atom_position(candidate, keto_aldehyde_ref, keto_aldehyde_spec, carbon_atom_id)
        hydrogen_atom_id = self._select_hydrogens(
            candidate,
            amine_ref,
            amine_spec,
            self._hydrogen_atom_ids_for_atom(amine_spec, amine_motif, nitrogen_atom_id),
            target=carbon_world,
            count=1,
            context="beta-ketoenamine formation",
        )[0]
        nitrogen_world = self._world_atom_position(candidate, amine_ref, amine_spec, nitrogen_atom_id)
        carbonyl_oxygen_local_position = self._shorten_bond_local_position(
            keto_aldehyde_spec.atom_positions[carbonyl_anchor_atom_id],
            keto_aldehyde_spec.atom_positions[tautomer_oxygen_atom_id],
            target_distance=1.24,
        )

        return EventRealization(
            removed_atom_ids={
                amine_ref.monomer_instance_id: (hydrogen_atom_id,),
                keto_aldehyde_ref.monomer_instance_id: (oxygen_atom_id, tautomer_hydrogen_atom_id),
            },
            atom_position_overrides={
                keto_aldehyde_ref.monomer_instance_id: {
                    tautomer_oxygen_atom_id: carbonyl_oxygen_local_position,
                },
            },
            bonds=(
                RealizedBond(
                    label_1=self.atom_label(
                        keto_aldehyde_ref.monomer_instance_id,
                        keto_aldehyde_spec.atom_symbols[carbon_atom_id],
                        carbon_atom_id,
                    ),
                    label_2=self.atom_label(
                        amine_ref.monomer_instance_id,
                        amine_spec.atom_symbols[nitrogen_atom_id],
                        nitrogen_atom_id,
                    ),
                    distance=self._distance(carbon_world, nitrogen_world),
                    symmetry_1=".",
                    symmetry_2=self._symmetry_code(amine_ref.periodic_image),
                ),
            ),
            notes=(
                "Beta-ketoenamine realization removes one amine hydrogen, the aldehydic oxygen, and the ortho-hydroxyl hydrogen per reacting pair.",
                "The exported product keeps the tautomerized ortho oxygen and shortens its local C=O distance to a carbonyl-like value.",
            ),
        )

    def _realize_boronate_ester_bridge(
        self,
        event: ReactionEvent,
        candidate: Candidate,
        monomer_specs: Mapping[str, MonomerSpec],
    ) -> EventRealization:
        boronic_ref, catechol_ref = self._split_bridge_participants(event, monomer_specs, "boronic_acid", "catechol")
        boronic_spec = monomer_specs[boronic_ref.monomer_id]
        catechol_spec = monomer_specs[catechol_ref.monomer_id]
        boronic_motif = boronic_spec.motif_by_id(boronic_ref.motif_id)
        catechol_motif = catechol_spec.motif_by_id(catechol_ref.motif_id)

        boron_atom_id = self._reactive_atom_id(boronic_motif)
        boronic_oxygen_atom_ids = self._motif_atom_ids_from_metadata(
            boronic_motif,
            "oxygen_atom_ids",
            fallback_symbol="O",
            monomer=boronic_spec,
        )
        catechol_oxygen_atom_ids = self._motif_atom_ids_from_metadata(
            catechol_motif,
            "reactive_atom_ids",
            fallback_symbol="O",
            monomer=catechol_spec,
        )
        if len(boronic_oxygen_atom_ids) != 2 or len(catechol_oxygen_atom_ids) != 2:
            raise ValueError(f"event {event.id!r} requires exactly two boronic and two catechol oxygen atoms")
        boronic_hydrogen_atom_ids = self._hydrogen_atom_ids_for_atoms(boronic_spec, boronic_motif, boronic_oxygen_atom_ids)
        catechol_hydrogen_atom_ids = self._hydrogen_atom_ids_for_atoms(catechol_spec, catechol_motif, catechol_oxygen_atom_ids)
        boron_world = self._world_atom_position(candidate, boronic_ref, boronic_spec, boron_atom_id)

        bonds = tuple(
            RealizedBond(
                label_1=self.atom_label(boronic_ref.monomer_instance_id, boronic_spec.atom_symbols[boron_atom_id], boron_atom_id),
                label_2=self.atom_label(catechol_ref.monomer_instance_id, catechol_spec.atom_symbols[oxygen_atom_id], oxygen_atom_id),
                distance=self._distance(
                    boron_world,
                    self._world_atom_position(candidate, catechol_ref, catechol_spec, oxygen_atom_id),
                ),
                symmetry_1=".",
                symmetry_2=self._symmetry_code(catechol_ref.periodic_image),
            )
            for oxygen_atom_id in catechol_oxygen_atom_ids
        )
        return EventRealization(
            removed_atom_ids={
                boronic_ref.monomer_instance_id: tuple(sorted((*boronic_oxygen_atom_ids, *boronic_hydrogen_atom_ids))),
                catechol_ref.monomer_instance_id: tuple(sorted(catechol_hydrogen_atom_ids)),
            },
            bonds=bonds,
            notes=(
                "Boronate ester realization removes both boronic-acid hydroxyl oxygens and all four condensation hydrogens per reacting pair.",
                "The exported product realizes the two inter-monomer B-O bonds to the catechol oxygens while keeping monomer-internal coordinates rigid.",
            ),
        )

    def _realize_vinylene_bridge(
        self,
        event: ReactionEvent,
        candidate: Candidate,
        monomer_specs: Mapping[str, MonomerSpec],
    ) -> EventRealization:
        activated_ref, aldehyde_ref = self._split_bridge_participants(
            event,
            monomer_specs,
            "activated_methylene",
            "aldehyde",
        )
        activated_spec = monomer_specs[activated_ref.monomer_id]
        aldehyde_spec = monomer_specs[aldehyde_ref.monomer_id]
        activated_motif = activated_spec.motif_by_id(activated_ref.motif_id)
        aldehyde_motif = aldehyde_spec.motif_by_id(aldehyde_ref.motif_id)

        activated_carbon_atom_id = self._reactive_atom_id(activated_motif)
        aldehyde_carbon_atom_id = self._reactive_atom_id(aldehyde_motif)
        aldehyde_oxygen_atom_id = self._unique_symbol_atom_id(
            aldehyde_spec,
            aldehyde_motif.atom_ids,
            "O",
            context=f"{event.id} aldehyde oxygen",
        )
        aldehyde_carbon_world = self._world_atom_position(candidate, aldehyde_ref, aldehyde_spec, aldehyde_carbon_atom_id)
        activated_hydrogen_atom_ids = self._select_hydrogens(
            candidate,
            activated_ref,
            activated_spec,
            self._hydrogen_atom_ids_for_atom(activated_spec, activated_motif, activated_carbon_atom_id),
            target=aldehyde_carbon_world,
            count=2,
            context="vinylene formation",
        )
        activated_carbon_world = self._world_atom_position(candidate, activated_ref, activated_spec, activated_carbon_atom_id)

        return EventRealization(
            removed_atom_ids={
                activated_ref.monomer_instance_id: activated_hydrogen_atom_ids,
                aldehyde_ref.monomer_instance_id: (aldehyde_oxygen_atom_id,),
            },
            bonds=(
                RealizedBond(
                    label_1=self.atom_label(
                        aldehyde_ref.monomer_instance_id,
                        aldehyde_spec.atom_symbols[aldehyde_carbon_atom_id],
                        aldehyde_carbon_atom_id,
                    ),
                    label_2=self.atom_label(
                        activated_ref.monomer_instance_id,
                        activated_spec.atom_symbols[activated_carbon_atom_id],
                        activated_carbon_atom_id,
                    ),
                    distance=self._distance(aldehyde_carbon_world, activated_carbon_world),
                    symmetry_1=".",
                    symmetry_2=self._symmetry_code(activated_ref.periodic_image),
                ),
            ),
            notes=(
                "Vinylene realization removes the aldehydic oxygen and two activated-carbon hydrogens per reacting pair.",
                "The exported product realizes the inter-monomer C=C connectivity only; internal bond-order redistribution within the activated precursor is not reoptimized.",
            ),
        )

    def _split_bridge_participants(
        self,
        event: ReactionEvent,
        monomer_specs: Mapping[str, MonomerSpec],
        first_kind: str,
        second_kind: str,
    ) -> tuple[MotifRef, MotifRef]:
        if len(event.participants) != 2:
            raise ValueError(f"event {event.id!r} must have exactly two participants")
        first, second = event.participants
        first_motif = monomer_specs[first.monomer_id].motif_by_id(first.motif_id)
        second_motif = monomer_specs[second.monomer_id].motif_by_id(second.motif_id)
        if first_motif.kind == first_kind and second_motif.kind == second_kind:
            return first, second
        if first_motif.kind == second_kind and second_motif.kind == first_kind:
            return second, first
        raise ValueError(
            f"event {event.id!r} does not match expected participant kinds {(first_kind, second_kind)!r}"
        )

    def _reactive_atom_id(self, motif) -> int:
        atom_id = motif.metadata.get("reactive_atom_id")
        if not isinstance(atom_id, int):
            raise ValueError(f"motif {motif.id!r} is missing integer reactive_atom_id metadata")
        return atom_id

    def _motif_atom_ids_from_metadata(
        self,
        motif,
        key: str,
        *,
        fallback_symbol: str,
        monomer: MonomerSpec,
    ) -> tuple[int, ...]:
        atom_ids = motif.metadata.get(key)
        if isinstance(atom_ids, tuple) and atom_ids and all(isinstance(atom_id, int) for atom_id in atom_ids):
            return tuple(atom_ids)
        return tuple(atom_id for atom_id in motif.atom_ids if monomer.atom_symbols[atom_id] == fallback_symbol)

    def _motif_atom_id_from_metadata(self, motif, key: str, *, context: str) -> int:
        atom_id = motif.metadata.get(key)
        if not isinstance(atom_id, int):
            raise ValueError(f"{context} is missing integer {key!r} metadata")
        return atom_id

    def _hydrogen_atom_ids_for_atom(self, monomer: MonomerSpec, motif, atom_id: int) -> tuple[int, ...]:
        hydrogen_atom_ids = motif.metadata.get("hydrogen_atom_ids")
        if isinstance(hydrogen_atom_ids, tuple) and hydrogen_atom_ids and all(isinstance(item, int) for item in hydrogen_atom_ids):
            attached = tuple(
                hydrogen_atom_id
                for hydrogen_atom_id in hydrogen_atom_ids
                if self._bond_exists(monomer, atom_id, hydrogen_atom_id)
            )
            if attached:
                return attached
            return hydrogen_atom_ids
        attached = tuple(
            other_atom_id
            for first_atom_id, second_atom_id, _order in monomer.bonds
            for other_atom_id in (
                (second_atom_id,) if first_atom_id == atom_id and monomer.atom_symbols[second_atom_id] == "H" else ()
            )
        )
        reverse = tuple(
            first_atom_id
            for first_atom_id, second_atom_id, _order in monomer.bonds
            if second_atom_id == atom_id and monomer.atom_symbols[first_atom_id] == "H"
        )
        bonded = tuple(sorted(set(attached + reverse)))
        if bonded:
            return bonded
        return tuple(
            candidate_atom_id
            for candidate_atom_id in motif.atom_ids
            if candidate_atom_id != atom_id and monomer.atom_symbols[candidate_atom_id] == "H"
        )

    def _hydrogen_atom_ids_for_atoms(self, monomer: MonomerSpec, motif, atom_ids: tuple[int, ...]) -> tuple[int, ...]:
        hydrogen_atom_ids = motif.metadata.get("hydrogen_atom_ids")
        if isinstance(hydrogen_atom_ids, tuple) and all(isinstance(item, int) for item in hydrogen_atom_ids):
            return tuple(hydrogen_atom_ids)
        collected: list[int] = []
        for atom_id in atom_ids:
            collected.extend(self._hydrogen_atom_ids_for_atom(monomer, motif, atom_id))
        return tuple(sorted(set(collected)))

    def _bond_exists(self, monomer: MonomerSpec, first_atom_id: int, second_atom_id: int) -> bool:
        return any(
            (left == first_atom_id and right == second_atom_id) or (left == second_atom_id and right == first_atom_id)
            for left, right, _order in monomer.bonds
        )

    def _select_hydrogens(
        self,
        candidate: Candidate,
        ref: MotifRef,
        monomer: MonomerSpec,
        hydrogen_atom_ids: tuple[int, ...],
        *,
        target: Vec3,
        count: int,
        context: str,
    ) -> tuple[int, ...]:
        if len(hydrogen_atom_ids) < count:
            raise ValueError(f"{context} requires {count} removable hydrogens, found {hydrogen_atom_ids!r}")
        ranked = sorted(
            hydrogen_atom_ids,
            key=lambda atom_id: self._distance(
                self._world_atom_position(candidate, ref, monomer, atom_id),
                target,
            ),
        )
        return tuple(ranked[:count])

    def _unique_symbol_atom_id(
        self,
        monomer: MonomerSpec,
        atom_ids: tuple[int, ...],
        symbol: str,
        *,
        context: str,
    ) -> int:
        matches = [atom_id for atom_id in atom_ids if monomer.atom_symbols[atom_id] == symbol]
        if len(matches) != 1:
            raise ValueError(f"{context} expected exactly one {symbol} atom, found {matches!r}")
        return matches[0]

    def _world_atom_position(
        self,
        candidate: Candidate,
        ref: MotifRef,
        monomer: MonomerSpec,
        atom_id: int,
        atom_position_overrides: Mapping[str, Mapping[int, Vec3]] | None = None,
    ) -> Vec3:
        pose = candidate.state.monomer_poses[ref.monomer_instance_id]
        image_shift = self._periodic_offset(candidate.state.cell, ref.periodic_image)
        if atom_position_overrides is not None:
            instance_overrides = atom_position_overrides.get(ref.monomer_instance_id, {})
            if atom_id in instance_overrides:
                return add(self._world_position(pose, instance_overrides[atom_id]), image_shift)
        return add(self._world_position(pose, monomer.atom_positions[atom_id]), image_shift)

    def _periodic_offset(
        self,
        cell: tuple[Vec3, Vec3, Vec3],
        image: tuple[int, int, int],
    ) -> Vec3:
        return add(add(scale(cell[0], image[0]), scale(cell[1], image[1])), scale(cell[2], image[2]))

    def _world_position(self, pose: Pose, local_position: Vec3) -> Vec3:
        return add(pose.translation, matmul_vec(pose.rotation_matrix, local_position))

    def _distance(self, left: Vec3, right: Vec3) -> float:
        dx = left[0] - right[0]
        dy = left[1] - right[1]
        dz = left[2] - right[2]
        return (dx * dx + dy * dy + dz * dz) ** 0.5

    def _shorten_bond_local_position(
        self,
        anchor_position: Vec3,
        atom_position: Vec3,
        *,
        target_distance: float,
    ) -> Vec3:
        vector = (
            atom_position[0] - anchor_position[0],
            atom_position[1] - anchor_position[1],
            atom_position[2] - anchor_position[2],
        )
        current_distance = self._distance(anchor_position, atom_position)
        if current_distance < 1e-8:
            return (anchor_position[0] + target_distance, anchor_position[1], anchor_position[2])
        scale_factor = target_distance / current_distance
        return (
            anchor_position[0] + vector[0] * scale_factor,
            anchor_position[1] + vector[1] * scale_factor,
            anchor_position[2] + vector[2] * scale_factor,
        )

    def _keto_enamine_carbonyl_anchor_atom_id(
        self,
        monomer: MonomerSpec,
        motif,
        oxygen_atom_id: int,
    ) -> int:
        metadata_atom_id = motif.metadata.get("ortho_hydroxyl_anchor_atom_id")
        if isinstance(metadata_atom_id, int):
            return metadata_atom_id
        heavy_neighbors = [
            candidate_atom_id
            for first_atom_id, second_atom_id, _order in monomer.bonds
            for candidate_atom_id in (
                (second_atom_id,) if first_atom_id == oxygen_atom_id and monomer.atom_symbols[second_atom_id] != "H" else ()
            )
        ]
        heavy_neighbors.extend(
            first_atom_id
            for first_atom_id, second_atom_id, _order in monomer.bonds
            if second_atom_id == oxygen_atom_id and monomer.atom_symbols[first_atom_id] != "H"
        )
        unique_neighbors = tuple(sorted(set(heavy_neighbors)))
        if len(unique_neighbors) != 1:
            raise ValueError(
                f"motif {motif.id!r} is missing ortho_hydroxyl_anchor_atom_id metadata and could not infer a unique carbonyl anchor from bonds"
            )
        return unique_neighbors[0]

    def _symmetry_code(self, image: tuple[int, int, int]) -> str:
        if image == (0, 0, 0):
            return "."
        return f"1_{image[0] + 5}{image[1] + 5}{image[2] + 5}"

    def _cleanup_reaction_hydrogens(
        self,
        candidate: Candidate,
        monomer_specs: Mapping[str, MonomerSpec],
        instance_to_monomer: Mapping[str, str],
        removed_atom_ids: Mapping[str, set[int]],
        atom_position_overrides: dict[str, dict[int, Vec3]],
        extra_heavy_atoms: tuple[tuple[str, int, Vec3], ...] = (),
    ) -> Mapping[str, object]:
        product_connections = self._product_connection_map(
            candidate,
            monomer_specs,
            atom_position_overrides=atom_position_overrides,
        )
        heavy_atoms = (
            *self._heavy_atom_world_positions(
                candidate,
                monomer_specs,
                instance_to_monomer,
                removed_atom_ids,
                atom_position_overrides,
            ),
            *extra_heavy_atoms,
        )
        moved_atom_labels: list[str] = []
        visited_parents: set[tuple[str, int]] = set()
        visited_hydrogen_atoms: set[tuple[str, int]] = set()

        for event in candidate.events:
            for ref in event.participants:
                monomer = monomer_specs[ref.monomer_id]
                if not (monomer.atom_symbols and monomer.atom_positions):
                    continue
                motif = monomer.motif_by_id(ref.motif_id)
                removed = removed_atom_ids.get(ref.monomer_instance_id, set())
                for parent_atom_id in self._hydrogen_cleanup_scope_atom_ids(monomer, motif, removed):
                    parent_key = (ref.monomer_instance_id, parent_atom_id)
                    if parent_key in visited_parents:
                        continue
                    visited_parents.add(parent_key)
                    hydrogen_atom_ids = self._attached_hydrogen_ids_for_parent(
                        monomer,
                        motif,
                        parent_atom_id,
                        removed,
                    )
                    if len(hydrogen_atom_ids) != 1:
                        continue
                    hydrogen_atom_id = hydrogen_atom_ids[0]
                    hydrogen_key = (ref.monomer_instance_id, hydrogen_atom_id)
                    if hydrogen_key in visited_hydrogen_atoms:
                        continue
                    new_local_position = self._cleaned_single_hydrogen_local_position(
                        candidate,
                        ref.monomer_instance_id,
                        monomer,
                        motif,
                        parent_atom_id,
                        hydrogen_atom_id,
                        removed,
                        atom_position_overrides.get(ref.monomer_instance_id, {}),
                        product_connections,
                        heavy_atoms,
                    )
                    if new_local_position is None:
                        continue
                    current_local_position = atom_position_overrides.get(ref.monomer_instance_id, {}).get(
                        hydrogen_atom_id,
                        monomer.atom_positions[hydrogen_atom_id],
                    )
                    if self._distance(new_local_position, current_local_position) <= 1e-4:
                        continue
                    atom_position_overrides[ref.monomer_instance_id][hydrogen_atom_id] = new_local_position
                    visited_hydrogen_atoms.add(hydrogen_key)
                    moved_atom_labels.append(
                        self.atom_label(
                            ref.monomer_instance_id,
                            monomer.atom_symbols[hydrogen_atom_id],
                            hydrogen_atom_id,
                        )
                    )

        if not moved_atom_labels:
            return {}
        return {
            "n_overrides": len(moved_atom_labels),
            "atom_labels": tuple(moved_atom_labels),
        }

    def _product_connection_map(
        self,
        candidate: Candidate,
        monomer_specs: Mapping[str, MonomerSpec],
        atom_position_overrides: Mapping[str, Mapping[int, Vec3]] | None = None,
    ) -> Mapping[tuple[str, int], tuple[tuple[str, int, Vec3], ...]]:
        connections: dict[tuple[str, int], list[tuple[str, int, Vec3]]] = defaultdict(list)
        for event in candidate.events:
            for first_ref, first_atom_id, second_ref, second_atom_id in self._event_connection_pairs(event, monomer_specs):
                first_monomer = monomer_specs[first_ref.monomer_id]
                second_monomer = monomer_specs[second_ref.monomer_id]
                first_world = self._world_atom_position(
                    candidate,
                    first_ref,
                    first_monomer,
                    first_atom_id,
                    atom_position_overrides=atom_position_overrides,
                )
                second_world = self._world_atom_position(
                    candidate,
                    second_ref,
                    second_monomer,
                    second_atom_id,
                    atom_position_overrides=atom_position_overrides,
                )
                first_pose = candidate.state.monomer_poses[first_ref.monomer_instance_id]
                second_pose = candidate.state.monomer_poses[second_ref.monomer_instance_id]
                first_local_vector = matmul_vec(
                    transpose(first_pose.rotation_matrix),
                    sub(second_world, first_world),
                )
                second_local_vector = matmul_vec(
                    transpose(second_pose.rotation_matrix),
                    sub(first_world, second_world),
                )
                connections[(first_ref.monomer_instance_id, first_atom_id)].append(
                    (second_ref.monomer_instance_id, second_atom_id, first_local_vector)
                )
                connections[(second_ref.monomer_instance_id, second_atom_id)].append(
                    (first_ref.monomer_instance_id, first_atom_id, second_local_vector)
                )
        return {key: tuple(values) for key, values in connections.items()}

    def _event_connection_pairs(
        self,
        event: ReactionEvent,
        monomer_specs: Mapping[str, MonomerSpec],
    ) -> tuple[tuple[MotifRef, int, MotifRef, int], ...]:
        realizer_id = linkage_event_realizer(event.template_id)
        if realizer_id in {"imine_bridge", "hydrazone_bridge", "keto_enamine_bridge", "vinylene_bridge"}:
            if realizer_id == "imine_bridge":
                first_kind, second_kind = "amine", "aldehyde"
            elif realizer_id == "hydrazone_bridge":
                first_kind, second_kind = "hydrazide", "aldehyde"
            elif realizer_id == "keto_enamine_bridge":
                first_kind, second_kind = "amine", "keto_aldehyde"
            else:
                first_kind, second_kind = "activated_methylene", "aldehyde"
            first_ref, second_ref = self._split_bridge_participants(event, monomer_specs, first_kind, second_kind)
            first_motif = monomer_specs[first_ref.monomer_id].motif_by_id(first_ref.motif_id)
            second_motif = monomer_specs[second_ref.monomer_id].motif_by_id(second_ref.motif_id)
            return (
                (
                    first_ref,
                    self._reactive_atom_id(first_motif),
                    second_ref,
                    self._reactive_atom_id(second_motif),
                ),
            )
        if realizer_id == "boronate_ester_bridge":
            boronic_ref, catechol_ref = self._split_bridge_participants(
                event,
                monomer_specs,
                "boronic_acid",
                "catechol",
            )
            boronic_spec = monomer_specs[boronic_ref.monomer_id]
            catechol_spec = monomer_specs[catechol_ref.monomer_id]
            boronic_motif = boronic_spec.motif_by_id(boronic_ref.motif_id)
            catechol_motif = catechol_spec.motif_by_id(catechol_ref.motif_id)
            boron_atom_id = self._reactive_atom_id(boronic_motif)
            catechol_oxygen_atom_ids = self._motif_atom_ids_from_metadata(
                catechol_motif,
                "reactive_atom_ids",
                fallback_symbol="O",
                monomer=catechol_spec,
            )
            return tuple(
                (boronic_ref, boron_atom_id, catechol_ref, oxygen_atom_id)
                for oxygen_atom_id in catechol_oxygen_atom_ids
            )
        return ()

    def _heavy_atom_world_positions(
        self,
        candidate: Candidate,
        monomer_specs: Mapping[str, MonomerSpec],
        instance_to_monomer: Mapping[str, str],
        removed_atom_ids: Mapping[str, set[int]],
        atom_position_overrides: Mapping[str, Mapping[int, Vec3]],
    ) -> tuple[tuple[str, int, Vec3], ...]:
        heavy_atoms: list[tuple[str, int, Vec3]] = []
        for instance_id, monomer_id in instance_to_monomer.items():
            monomer = monomer_specs[monomer_id]
            if not (monomer.atom_symbols and monomer.atom_positions):
                continue
            pose = candidate.state.monomer_poses[instance_id]
            removed = removed_atom_ids.get(instance_id, set())
            overrides = atom_position_overrides.get(instance_id, {})
            for atom_id, symbol in enumerate(monomer.atom_symbols):
                if symbol == "H" or atom_id in removed:
                    continue
                local_position = overrides.get(atom_id, monomer.atom_positions[atom_id])
                heavy_atoms.append((instance_id, atom_id, self._world_position(pose, local_position)))
        return tuple(heavy_atoms)

    def _hydrogen_cleanup_scope_atom_ids(
        self,
        monomer: MonomerSpec,
        motif,
        removed_atom_ids: set[int],
    ) -> tuple[int, ...]:
        seed_atom_ids = {
            atom_id
            for atom_id in motif.atom_ids
            if atom_id not in removed_atom_ids and monomer.atom_symbols[atom_id] != "H"
        }
        scope = set(seed_atom_ids)
        if monomer.bonds:
            for atom_id in tuple(seed_atom_ids):
                scope.update(self._bonded_heavy_neighbor_atom_ids(monomer, atom_id, removed_atom_ids))
        else:
            fallback_ids = []
            reactive_atom_id = motif.metadata.get("reactive_atom_id")
            if isinstance(reactive_atom_id, int):
                fallback_ids.append(reactive_atom_id)
            anchor_atom_id = motif.metadata.get("anchor_atom_id")
            if isinstance(anchor_atom_id, int):
                fallback_ids.append(anchor_atom_id)
            anchor_atom_ids = motif.metadata.get("anchor_atom_ids")
            if isinstance(anchor_atom_ids, tuple):
                fallback_ids.extend(
                    atom_id for atom_id in anchor_atom_ids if isinstance(atom_id, int)
                )
            scope.update(
                atom_id
                for atom_id in fallback_ids
                if atom_id not in removed_atom_ids and monomer.atom_symbols[atom_id] != "H"
            )
        return tuple(sorted(scope))

    def _attached_hydrogen_ids_for_parent(
        self,
        monomer: MonomerSpec,
        motif,
        parent_atom_id: int,
        removed_atom_ids: set[int],
    ) -> tuple[int, ...]:
        attached = tuple(
            atom_id
            for atom_id in self._directly_bonded_hydrogen_atom_ids(monomer, parent_atom_id)
            if atom_id not in removed_atom_ids
        )
        if attached:
            return attached
        fallback_parent_atom_ids = self._fallback_hydrogen_parent_atom_ids(motif)
        if parent_atom_id not in fallback_parent_atom_ids:
            return ()
        return tuple(
            atom_id
            for atom_id in self._hydrogen_atom_ids_for_atom(monomer, motif, parent_atom_id)
            if atom_id not in removed_atom_ids
        )

    def _cleaned_single_hydrogen_local_position(
        self,
        candidate: Candidate,
        instance_id: str,
        monomer: MonomerSpec,
        motif,
        parent_atom_id: int,
        hydrogen_atom_id: int,
        removed_atom_ids: set[int],
        atom_position_overrides: Mapping[int, Vec3],
        product_connections: Mapping[tuple[str, int], tuple[tuple[str, int, Vec3], ...]],
        heavy_atoms: tuple[tuple[str, int, Vec3], ...],
    ) -> Vec3 | None:
        pose = candidate.state.monomer_poses[instance_id]
        parent_local = atom_position_overrides.get(parent_atom_id, monomer.atom_positions[parent_atom_id])
        hydrogen_local = atom_position_overrides.get(hydrogen_atom_id, monomer.atom_positions[hydrogen_atom_id])
        bond_length = self._reference_hydrogen_bond_length(
            monomer,
            parent_atom_id,
            hydrogen_atom_id,
        )
        current_direction = self._normalize_or_none(sub(hydrogen_local, parent_local))
        if current_direction is None:
            return None

        internal_neighbor_atom_ids = self._internal_heavy_neighbor_atom_ids(
            monomer,
            motif,
            parent_atom_id,
            removed_atom_ids,
        )
        internal_neighbor_vectors = [
            sub(atom_position_overrides.get(atom_id, monomer.atom_positions[atom_id]), parent_local)
            for atom_id in internal_neighbor_atom_ids
        ]
        external_connections = product_connections.get((instance_id, parent_atom_id), ())
        external_neighbor_vectors = [vector for _other_instance, _other_atom_id, vector in external_connections]
        bonded_neighbor_vectors = tuple(
            vector
            for vector in (*internal_neighbor_vectors, *external_neighbor_vectors)
            if norm(vector) >= 1e-8
        )
        heavy_neighbor_directions = tuple(
            direction
            for direction in (
                self._normalize_or_none(vector)
                for vector in bonded_neighbor_vectors
            )
            if direction is not None
        )

        parent_world = self._world_position(pose, parent_local)
        bonded_atom_keys = {(instance_id, atom_id) for atom_id in internal_neighbor_atom_ids}
        bonded_atom_keys.update(
            (other_instance_id, other_atom_id) for other_instance_id, other_atom_id, _vector in external_connections
        )
        blockers = tuple(
            world_position
            for other_instance_id, other_atom_id, world_position in heavy_atoms
            if (other_instance_id, other_atom_id) != (instance_id, parent_atom_id)
            and (other_instance_id, other_atom_id) not in bonded_atom_keys
            and self._distance(parent_world, world_position) <= 4.0
        )
        current_world = add(parent_world, matmul_vec(pose.rotation_matrix, scale(current_direction, bond_length)))
        current_clearance = self._minimum_distance(current_world, blockers)
        if not external_connections and (current_clearance is None or current_clearance >= 1.9):
            return None

        repulsion_vector = self._repulsion_direction(parent_world, blockers)
        repulsion_local = None
        if repulsion_vector is not None:
            repulsion_local = self._normalize_or_none(matmul_vec(transpose(pose.rotation_matrix), repulsion_vector))

        candidate_directions, preferred_direction = self._single_hydrogen_candidate_directions(
            current_direction=current_direction,
            heavy_neighbor_directions=heavy_neighbor_directions,
            motif=motif,
            repulsion_local=repulsion_local,
        )
        if not candidate_directions:
            return None

        preferred = preferred_direction or current_direction
        best_direction = current_direction
        best_score = self._score_hydrogen_direction(
            candidate_direction=current_direction,
            preferred_direction=preferred,
            parent_world=parent_world,
            pose=pose,
            bond_length=bond_length,
            bonded_neighbor_vectors=bonded_neighbor_vectors,
            blockers=blockers,
        )
        for candidate_direction in candidate_directions:
            candidate_score = self._score_hydrogen_direction(
                candidate_direction=candidate_direction,
                preferred_direction=preferred,
                parent_world=parent_world,
                pose=pose,
                bond_length=bond_length,
                bonded_neighbor_vectors=bonded_neighbor_vectors,
                blockers=blockers,
            )
            if candidate_score > best_score:
                best_direction = candidate_direction
                best_score = candidate_score

        target_local_position = add(parent_local, scale(best_direction, bond_length))
        if self._distance(target_local_position, hydrogen_local) <= 1e-4:
            return None
        return target_local_position

    def _internal_heavy_neighbor_atom_ids(
        self,
        monomer: MonomerSpec,
        motif,
        parent_atom_id: int,
        removed_atom_ids: set[int],
    ) -> tuple[int, ...]:
        bonded_neighbors = self._bonded_heavy_neighbor_atom_ids(monomer, parent_atom_id, removed_atom_ids)
        if bonded_neighbors:
            return bonded_neighbors
        fallback_ids: list[int] = []
        reactive_atom_id = motif.metadata.get("reactive_atom_id")
        if isinstance(reactive_atom_id, int) and reactive_atom_id == parent_atom_id:
            anchor_atom_id = motif.metadata.get("anchor_atom_id")
            if isinstance(anchor_atom_id, int):
                fallback_ids.append(anchor_atom_id)
            anchor_atom_ids = motif.metadata.get("anchor_atom_ids")
            if isinstance(anchor_atom_ids, tuple):
                fallback_ids.extend(atom_id for atom_id in anchor_atom_ids if isinstance(atom_id, int))
        return tuple(
            sorted(
                atom_id
                for atom_id in set(fallback_ids)
                if atom_id != parent_atom_id
                and atom_id not in removed_atom_ids
                and monomer.atom_symbols[atom_id] != "H"
            )
        )

    def _bonded_heavy_neighbor_atom_ids(
        self,
        monomer: MonomerSpec,
        atom_id: int,
        removed_atom_ids: set[int],
    ) -> tuple[int, ...]:
        neighbors: set[int] = set()
        for first_atom_id, second_atom_id, _order in monomer.bonds:
            if first_atom_id == atom_id and second_atom_id not in removed_atom_ids and monomer.atom_symbols[second_atom_id] != "H":
                neighbors.add(second_atom_id)
            elif second_atom_id == atom_id and first_atom_id not in removed_atom_ids and monomer.atom_symbols[first_atom_id] != "H":
                neighbors.add(first_atom_id)
        return tuple(sorted(neighbors))

    def _directly_bonded_hydrogen_atom_ids(self, monomer: MonomerSpec, atom_id: int) -> tuple[int, ...]:
        neighbors: set[int] = set()
        for first_atom_id, second_atom_id, _order in monomer.bonds:
            if first_atom_id == atom_id and monomer.atom_symbols[second_atom_id] == "H":
                neighbors.add(second_atom_id)
            elif second_atom_id == atom_id and monomer.atom_symbols[first_atom_id] == "H":
                neighbors.add(first_atom_id)
        return tuple(sorted(neighbors))

    def _single_hydrogen_candidate_directions(
        self,
        *,
        current_direction: Vec3,
        heavy_neighbor_directions: tuple[Vec3, ...],
        motif,
        repulsion_local: Vec3 | None,
    ) -> tuple[tuple[Vec3, ...], Vec3 | None]:
        candidates: list[Vec3] = []
        preferred_direction: Vec3 | None = None

        def _add(direction: Vec3 | None) -> None:
            if direction is None:
                return
            normalized = self._normalize_or_none(direction)
            if normalized is None:
                return
            if any(self._distance(normalized, existing) <= 1e-4 for existing in candidates):
                return
            candidates.append(normalized)

        _add(current_direction)

        if len(heavy_neighbor_directions) == 1:
            axis = heavy_neighbor_directions[0]
            preferred_direction = scale(axis, -1.0)
            _add(preferred_direction)
            torsion_seed = None
            if repulsion_local is not None:
                torsion_seed = self._orthogonal_component(repulsion_local, axis)
            if torsion_seed is None or norm(torsion_seed) < 1e-8:
                torsion_seed = self._orthogonal_component(current_direction, axis)
            if torsion_seed is None or norm(torsion_seed) < 1e-8:
                torsion_seed = self._orthogonal_component(motif.frame.normal, axis)
            if torsion_seed is None or norm(torsion_seed) < 1e-8:
                torsion_seed = self._any_orthogonal_unit(axis)
            if torsion_seed is not None and norm(torsion_seed) >= 1e-8:
                axis_unit = normalize(axis)
                basis_u = normalize(torsion_seed)
                basis_v = normalize(cross(axis_unit, basis_u))
                parallel_component = dot(current_direction, axis_unit)
                radial_component = max(0.0, 1.0 - parallel_component * parallel_component) ** 0.5
                if radial_component > 1e-6:
                    for step in range(12):
                        angle = 2.0 * pi * step / 12.0
                        _add(
                            add(
                                scale(axis_unit, parallel_component),
                                add(
                                    scale(basis_u, radial_component * cos(angle)),
                                    scale(basis_v, radial_component * sin(angle)),
                                ),
                            )
                        )
        elif heavy_neighbor_directions:
            composite = (
                sum(direction[0] for direction in heavy_neighbor_directions),
                sum(direction[1] for direction in heavy_neighbor_directions),
                sum(direction[2] for direction in heavy_neighbor_directions),
            )
            if norm(composite) >= 1e-8:
                preferred_direction = scale(normalize(composite), -1.0)
                _add(preferred_direction)
                if repulsion_local is not None:
                    for weight in (0.35, 0.7, 1.05):
                        blended = add(preferred_direction, scale(repulsion_local, weight))
                        blended_direction = self._normalize_or_none(blended)
                        if blended_direction is None:
                            continue
                        if dot(blended_direction, preferred_direction) >= 0.55:
                            _add(blended_direction)
            else:
                axis = heavy_neighbor_directions[0]
                preferred_direction = None
                in_plane_seed = None
                if repulsion_local is not None:
                    in_plane_seed = self._orthogonal_component(repulsion_local, axis)
                if in_plane_seed is None or norm(in_plane_seed) < 1e-8:
                    in_plane_seed = self._orthogonal_component(current_direction, axis)
                if in_plane_seed is None or norm(in_plane_seed) < 1e-8:
                    in_plane_seed = cross(motif.frame.normal, axis)
                if in_plane_seed is None or norm(in_plane_seed) < 1e-8:
                    in_plane_seed = self._any_orthogonal_unit(axis)
                if in_plane_seed is not None:
                    perpendicular = self._normalize_or_none(in_plane_seed)
                    _add(perpendicular)
                    if perpendicular is not None:
                        _add(scale(perpendicular, -1.0))
        elif repulsion_local is not None:
            preferred_direction = repulsion_local
            _add(repulsion_local)

        return tuple(candidates), preferred_direction

    def _score_hydrogen_direction(
        self,
        *,
        candidate_direction: Vec3,
        preferred_direction: Vec3,
        parent_world: Vec3,
        pose: Pose,
        bond_length: float,
        bonded_neighbor_vectors: tuple[Vec3, ...],
        blockers: tuple[Vec3, ...],
    ) -> tuple[float, float, float]:
        candidate_local = scale(candidate_direction, bond_length)
        bonded_clearance = self._minimum_distance(candidate_local, bonded_neighbor_vectors)
        candidate_world = add(parent_world, matmul_vec(pose.rotation_matrix, scale(candidate_direction, bond_length)))
        clearance = self._minimum_distance(candidate_world, blockers)
        return (
            bonded_clearance if bonded_clearance is not None else float("inf"),
            clearance if clearance is not None else float("inf"),
            dot(candidate_direction, preferred_direction),
        )

    def _minimum_distance(self, point: Vec3, others: tuple[Vec3, ...]) -> float | None:
        if not others:
            return None
        return min(self._distance(point, other) for other in others)

    def _repulsion_direction(self, origin: Vec3, blockers: tuple[Vec3, ...]) -> Vec3 | None:
        if not blockers:
            return None
        repulsion = (0.0, 0.0, 0.0)
        for blocker in blockers:
            vector = sub(origin, blocker)
            distance = norm(vector)
            if distance < 1e-8:
                continue
            repulsion = add(repulsion, scale(vector, 1.0 / max(distance * distance, 1e-6)))
        return self._normalize_or_none(repulsion)

    def _orthogonal_component(self, vector: Vec3, axis: Vec3) -> Vec3:
        axis_unit = self._normalize_or_none(axis)
        if axis_unit is None:
            return vector
        return sub(vector, scale(axis_unit, dot(vector, axis_unit)))

    def _any_orthogonal_unit(self, axis: Vec3) -> Vec3 | None:
        axis_unit = self._normalize_or_none(axis)
        if axis_unit is None:
            return None
        seed = (0.0, 0.0, 1.0) if abs(axis_unit[2]) < 0.9 else (1.0, 0.0, 0.0)
        orthogonal = self._orthogonal_component(seed, axis_unit)
        return self._normalize_or_none(orthogonal)

    def _normalize_or_none(self, vector: Vec3) -> Vec3 | None:
        if norm(vector) < 1e-8:
            return None
        return normalize(vector)

    def _default_hydrogen_bond_length(self, parent_symbol: str) -> float:
        return {
            "B": 1.19,
            "C": 1.09,
            "N": 1.01,
            "O": 0.98,
            "S": 1.34,
        }.get(parent_symbol, 1.0)

    def _reference_hydrogen_bond_length(
        self,
        monomer: MonomerSpec,
        parent_atom_id: int,
        hydrogen_atom_id: int,
    ) -> float:
        original_bond_length = self._distance(
            monomer.atom_positions[parent_atom_id],
            monomer.atom_positions[hydrogen_atom_id],
        )
        default_bond_length = self._default_hydrogen_bond_length(monomer.atom_symbols[parent_atom_id])
        if original_bond_length < 1e-6:
            return default_bond_length
        if abs(original_bond_length - default_bond_length) > 0.35:
            return default_bond_length
        return original_bond_length

    def _fallback_hydrogen_parent_atom_ids(self, motif) -> set[int]:
        parent_atom_ids: set[int] = set()
        for key in ("reactive_atom_id", "ortho_hydroxyl_oxygen_atom_id"):
            atom_id = motif.metadata.get(key)
            if isinstance(atom_id, int):
                parent_atom_ids.add(atom_id)
        for key in ("reactive_atom_ids", "oxygen_atom_ids"):
            atom_ids = motif.metadata.get(key)
            if isinstance(atom_ids, tuple):
                parent_atom_ids.update(atom_id for atom_id in atom_ids if isinstance(atom_id, int))
        return parent_atom_ids

    @staticmethod
    def atom_label(instance_id: str, symbol: str, atom_id: int) -> str:
        return f"{instance_id}_{symbol}{atom_id + 1}"


def _dispatch_imine_bridge_event(
    realizer: ReactionRealizer,
    event: ReactionEvent,
    candidate: Candidate,
    monomer_specs: Mapping[str, MonomerSpec],
) -> EventRealization:
    return realizer._realize_imine_bridge(event, candidate, monomer_specs)


def _dispatch_hydrazone_bridge_event(
    realizer: ReactionRealizer,
    event: ReactionEvent,
    candidate: Candidate,
    monomer_specs: Mapping[str, MonomerSpec],
) -> EventRealization:
    return realizer._realize_hydrazone_bridge(event, candidate, monomer_specs)


def _dispatch_keto_enamine_bridge_event(
    realizer: ReactionRealizer,
    event: ReactionEvent,
    candidate: Candidate,
    monomer_specs: Mapping[str, MonomerSpec],
) -> EventRealization:
    return realizer._realize_keto_enamine_bridge(event, candidate, monomer_specs)


def _dispatch_boronate_ester_bridge_event(
    realizer: ReactionRealizer,
    event: ReactionEvent,
    candidate: Candidate,
    monomer_specs: Mapping[str, MonomerSpec],
) -> EventRealization:
    return realizer._realize_boronate_ester_bridge(event, candidate, monomer_specs)


def _dispatch_vinylene_bridge_event(
    realizer: ReactionRealizer,
    event: ReactionEvent,
    candidate: Candidate,
    monomer_specs: Mapping[str, MonomerSpec],
) -> EventRealization:
    return realizer._realize_vinylene_bridge(event, candidate, monomer_specs)
