from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, permutations
from typing import Mapping

from .model import MonomerInstance, MonomerSpec, MotifRef, ReactionEvent, ReactionTemplate
from .planner import AssignmentPlan, NetPlan
from .single_node_topologies import resolve_single_node_topology_layout
from .single_node_topologies_3d import resolve_three_d_single_node_topology_layout


@dataclass(frozen=True)
class AssignmentOutcome:
    assignment_plan: AssignmentPlan
    monomer_instances: tuple[MonomerInstance, ...]
    events: tuple[ReactionEvent, ...]
    unreacted_motifs: tuple[MotifRef, ...]
    consumed_count: int


def _resolve_supported_single_node_layout(topology_id: str):
    try:
        return resolve_single_node_topology_layout(topology_id)
    except (KeyError, ValueError):
        return resolve_three_d_single_node_topology_layout(topology_id)


class AssignmentSolver:
    """Builds concrete monomer assignments and reaction events."""

    def build_assignment_plans(
        self,
        net_plan: NetPlan,
        monomer_specs: Mapping[str, MonomerSpec],
    ) -> tuple[AssignmentPlan, ...]:
        monomer_ids = net_plan.monomer_ids
        if net_plan.topology is None or not net_plan.topology.node_coordination:
            return (
                AssignmentPlan(
                    net_plan=net_plan,
                    slot_to_monomer={f"unit{i+1}": mid for i, mid in enumerate(monomer_ids)},
                    slot_to_conformer={
                        f"unit{i+1}": monomer_specs[mid].conformer_ids[0]
                        for i, mid in enumerate(monomer_ids)
                        if monomer_specs[mid].conformer_ids
                    },
                    metadata={"assignment_mode": "identity"},
                ),
            )

        if len(net_plan.topology.node_coordination) == 1 and len(monomer_ids) > 1:
            slot_count = len(monomer_ids)
        else:
            slot_count = len(net_plan.topology.node_coordination)
        slot_ids = tuple(f"slot{i+1}" for i in range(slot_count))
        unique_plans: list[AssignmentPlan] = []
        seen: set[tuple[tuple[str, str], ...]] = set()

        for ordering in permutations(monomer_ids, len(slot_ids)):
            if not self._ordering_matches_topology(ordering, slot_ids, net_plan, monomer_specs):
                continue
            mapping = {slot: monomer_id for slot, monomer_id in zip(slot_ids, ordering)}
            key = tuple(sorted(mapping.items()))
            if key in seen:
                continue
            seen.add(key)
            unique_plans.append(
                AssignmentPlan(
                    net_plan=net_plan,
                    slot_to_monomer=mapping,
                    slot_to_conformer={
                        slot: monomer_specs[mid].conformer_ids[0]
                        for slot, mid in mapping.items()
                        if monomer_specs[mid].conformer_ids
                    },
                    metadata={"assignment_mode": "coordination-match"},
                )
            )

        if unique_plans:
            return tuple(unique_plans)

        return (
            AssignmentPlan(
                net_plan=net_plan,
                slot_to_monomer={f"unit{i+1}": mid for i, mid in enumerate(monomer_ids)},
                metadata={"assignment_mode": "fallback-identity"},
            ),
        )

    def instantiate_monomers(self, assignment_plan: AssignmentPlan) -> tuple[MonomerInstance, ...]:
        instances: list[MonomerInstance] = []
        for idx, (slot_id, monomer_id) in enumerate(assignment_plan.slot_to_monomer.items(), start=1):
            conformer = assignment_plan.slot_to_conformer.get(slot_id)
            instances.append(
                MonomerInstance(
                    id=f"m{idx}",
                    monomer_id=monomer_id,
                    periodic_image=(0, 0, 0),
                    conformer_id=conformer,
                    metadata={"slot_id": slot_id},
                )
            )
        return tuple(instances)

    def solve_events(
        self,
        assignment_plan: AssignmentPlan,
        monomer_instances: tuple[MonomerInstance, ...],
        monomer_specs: Mapping[str, MonomerSpec],
        templates: tuple[ReactionTemplate, ...],
    ) -> AssignmentOutcome:
        motif_pool = tuple(self._motif_pool(monomer_instances, monomer_specs))
        best_events: tuple[ReactionEvent, ...] = ()
        best_consumed = -1

        def backtrack(used_keys: frozenset[tuple[str, str]], built_events: tuple[ReactionEvent, ...]) -> None:
            nonlocal best_events, best_consumed

            if len(used_keys) > best_consumed or (
                len(used_keys) == best_consumed and len(built_events) > len(best_events)
            ):
                best_events = built_events
                best_consumed = len(used_keys)

            remaining = [ref for ref in motif_pool if (ref.monomer_instance_id, ref.motif_id) not in used_keys]
            if len(used_keys) + len(remaining) < best_consumed:
                return

            next_event_id = len(built_events) + 1
            for template in templates:
                for combo in combinations(remaining, template.arity):
                    if not self._combo_matches_template(combo, template, monomer_specs):
                        continue
                    combo_keys = frozenset((ref.monomer_instance_id, ref.motif_id) for ref in combo)
                    if combo_keys & used_keys:
                        continue
                    event = ReactionEvent(
                        id=f"rxn{next_event_id}",
                        template_id=template.id,
                        participants=tuple(combo),
                    )
                    backtrack(used_keys | combo_keys, built_events + (event,))

        backtrack(frozenset(), tuple())

        used = {(ref.monomer_instance_id, ref.motif_id) for event in best_events for ref in event.participants}
        unreacted = tuple(
            ref for ref in motif_pool if (ref.monomer_instance_id, ref.motif_id) not in used
        )
        return AssignmentOutcome(
            assignment_plan=assignment_plan,
            monomer_instances=monomer_instances,
            events=self._decorate_single_node_bipartite_events(best_events, assignment_plan, monomer_instances, monomer_specs),
            unreacted_motifs=unreacted,
            consumed_count=best_consumed if best_consumed >= 0 else 0,
        )

    def _ordering_matches_topology(
        self,
        ordering: tuple[str, ...],
        slot_ids: tuple[str, ...],
        net_plan: NetPlan,
        monomer_specs: Mapping[str, MonomerSpec],
    ) -> bool:
        assert net_plan.topology is not None
        if len(net_plan.topology.node_coordination) == 1 and len(slot_ids) > 1:
            coordinations = (net_plan.topology.node_coordination[0],) * len(slot_ids)
        else:
            coordinations = net_plan.topology.node_coordination
        for slot_id, monomer_id, coordination in zip(slot_ids, ordering, coordinations):
            _ = slot_id
            if len(monomer_specs[monomer_id].motifs) != coordination:
                return False
        return True

    def _motif_pool(
        self,
        monomer_instances: tuple[MonomerInstance, ...],
        monomer_specs: Mapping[str, MonomerSpec],
    ) -> list[MotifRef]:
        pool: list[MotifRef] = []
        for instance in monomer_instances:
            monomer = monomer_specs[instance.monomer_id]
            for motif in monomer.motifs:
                pool.append(
                    MotifRef(
                        monomer_instance_id=instance.id,
                        monomer_id=monomer.id,
                        motif_id=motif.id,
                        periodic_image=instance.periodic_image,
                    )
                )
        return pool

    def _combo_matches_template(
        self,
        combo: tuple[MotifRef, ...],
        template: ReactionTemplate,
        monomer_specs: Mapping[str, MonomerSpec],
    ) -> bool:
        kinds = []
        for ref in combo:
            motif = monomer_specs[ref.monomer_id].motif_by_id(ref.motif_id)
            if not motif.accepts(template.id):
                return False
            kinds.append(motif.kind)
        return template.matches(tuple(kinds))

    def _decorate_single_node_bipartite_events(
        self,
        events: tuple[ReactionEvent, ...],
        assignment_plan: AssignmentPlan,
        monomer_instances: tuple[MonomerInstance, ...],
        monomer_specs: Mapping[str, MonomerSpec],
    ) -> tuple[ReactionEvent, ...]:
        topology = assignment_plan.net_plan.topology
        if topology is None or len(events) != 3 or len(monomer_instances) != 2:
            return events
        if int(topology.metadata.get("n_node_definitions", len(topology.node_coordination) or 0)) != 1:
            return events
        if len({instance.monomer_id for instance in monomer_instances}) != 2:
            return events
        if any(len(monomer_specs[instance.monomer_id].motifs) != 3 for instance in monomer_instances):
            return events
        layout = _resolve_supported_single_node_layout(topology.id)
        if layout.placement_model != "p1-two-node":
            raise ValueError(
                f"topology {topology.id!r} requires an expanded single-node builder; direct engine support is "
                "currently limited to two-node P1 single-node nets"
            )

        image_sequence = ((0, 0, 0), (-1, 0, 0), (0, -1, 0))
        decorated: list[ReactionEvent] = []
        second_instance_id = monomer_instances[1].id
        for image, event in zip(image_sequence, events):
            participants = []
            for participant in event.participants:
                periodic_image = image if participant.monomer_instance_id == second_instance_id else (0, 0, 0)
                participants.append(
                    MotifRef(
                        monomer_instance_id=participant.monomer_instance_id,
                        monomer_id=participant.monomer_id,
                        motif_id=participant.motif_id,
                        periodic_image=periodic_image,
                    )
                )
            decorated.append(
                ReactionEvent(
                    id=event.id,
                    template_id=event.template_id,
                    participants=tuple(participants),
                    product_state=event.product_state,
                    metadata=event.metadata,
                )
            )
        return tuple(decorated)
