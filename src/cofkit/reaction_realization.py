from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Callable, Mapping

from .geometry import Vec3, add, matmul_vec, scale
from .model import Candidate, MonomerSpec, MotifRef, Pose, ReactionEvent
from .reactions import linkage_event_realizer


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

        atoms_by_instance: dict[str, tuple[RealizedAtom, ...]] = {}
        removed_symbol_counts: Counter[str] = Counter()
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
            atoms_by_instance[instance_id] = tuple(realized_atoms)

        if applied_events == 0:
            return None

        metadata = {
            "applied_event_count": applied_events,
            "applied_templates": dict(applied_templates),
            "removed_atom_count": sum(removed_symbol_counts.values()),
            "removed_atom_symbols": dict(removed_symbol_counts),
            "notes": tuple(dict.fromkeys(notes)),
        }
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
        hydrogen_atom_id = self._select_amine_hydrogen(
            candidate,
            amine_ref,
            aldehyde_ref,
            monomer_specs,
            nitrogen_atom_id,
        )

        amine_spec = monomer_specs[amine_ref.monomer_id]
        aldehyde_spec = monomer_specs[aldehyde_ref.monomer_id]
        carbon_world = self._world_atom_position(candidate, aldehyde_ref, aldehyde_spec, carbon_atom_id)
        nitrogen_world = self._world_atom_position(candidate, amine_ref, amine_spec, nitrogen_atom_id)

        return EventRealization(
            removed_atom_ids={
                amine_ref.monomer_instance_id: (hydrogen_atom_id,),
                aldehyde_ref.monomer_instance_id: (oxygen_atom_id,),
            },
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
                "Imine realization removes the aldehyde oxygen and one amine hydrogen per reacting pair.",
                "The exported product keeps monomer-internal coordinates rigid; only atom deletion and inter-monomer C=N connectivity are realized.",
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

    def _select_amine_hydrogen(
        self,
        candidate: Candidate,
        amine_ref: MotifRef,
        aldehyde_ref: MotifRef,
        monomer_specs: Mapping[str, MonomerSpec],
        nitrogen_atom_id: int,
    ) -> int:
        amine_spec = monomer_specs[amine_ref.monomer_id]
        motif = amine_spec.motif_by_id(amine_ref.motif_id)
        aldehyde_carbon_id = self._reactive_atom_id(monomer_specs[aldehyde_ref.monomer_id].motif_by_id(aldehyde_ref.motif_id))
        target = self._world_atom_position(candidate, aldehyde_ref, monomer_specs[aldehyde_ref.monomer_id], aldehyde_carbon_id)
        return self._select_hydrogens(
            candidate,
            amine_ref,
            amine_spec,
            self._hydrogen_atom_ids_for_atom(amine_spec, motif, nitrogen_atom_id),
            target=target,
            count=1,
            context="imine formation",
        )[0]

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
    ) -> Vec3:
        pose = candidate.state.monomer_poses[ref.monomer_instance_id]
        local = monomer.atom_positions[atom_id]
        image_shift = self._periodic_offset(candidate.state.cell, ref.periodic_image)
        return add(self._world_position(pose, local), image_shift)

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
