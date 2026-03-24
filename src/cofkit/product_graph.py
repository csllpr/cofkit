from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Iterable, Mapping

from .model import MonomerInstance, MonomerSpec, ReactionEvent, ReactionTemplate


@dataclass
class PeriodicProductGraph:
    monomer_instances: dict[str, MonomerInstance] = field(default_factory=dict)
    reaction_events: list[ReactionEvent] = field(default_factory=list)
    _motif_usage: Counter[tuple[str, str]] = field(default_factory=Counter)

    def add_monomer(self, instance: MonomerInstance) -> None:
        if instance.id in self.monomer_instances:
            raise ValueError(f"duplicate monomer instance {instance.id!r}")
        self.monomer_instances[instance.id] = instance

    def add_reaction_event(
        self,
        event: ReactionEvent,
        monomer_specs: Mapping[str, MonomerSpec],
        templates: Mapping[str, ReactionTemplate],
    ) -> None:
        if event.template_id not in templates:
            raise KeyError(f"unknown reaction template {event.template_id!r}")

        template = templates[event.template_id]
        if len(event.participants) != template.arity:
            raise ValueError(
                f"event {event.id!r} has arity {len(event.participants)} but template {template.id!r} requires {template.arity}"
            )

        motif_kinds: list[str] = []
        for ref in event.participants:
            instance = self.monomer_instances.get(ref.monomer_instance_id)
            if instance is None:
                raise KeyError(f"unknown monomer instance {ref.monomer_instance_id!r}")
            if instance.monomer_id != ref.monomer_id:
                raise ValueError(
                    f"participant {ref.monomer_instance_id!r} claims monomer {ref.monomer_id!r} but instance stores {instance.monomer_id!r}"
                )

            monomer = monomer_specs[ref.monomer_id]
            motif = monomer.motif_by_id(ref.motif_id)
            if not motif.accepts(template.id):
                raise ValueError(
                    f"motif {ref.motif_id!r} in monomer {monomer.id!r} does not accept template {template.id!r}"
                )

            key = (ref.monomer_instance_id, ref.motif_id)
            if self._motif_usage[key] > 0:
                raise ValueError(
                    f"motif {ref.motif_id!r} on monomer instance {ref.monomer_instance_id!r} is already consumed"
                )

            motif_kinds.append(motif.kind)

        if not template.matches(tuple(motif_kinds)):
            raise ValueError(
                f"event {event.id!r} motif kinds {tuple(motif_kinds)!r} do not match template {template.id!r}"
            )

        self.reaction_events.append(event)
        for ref in event.participants:
            self._motif_usage[(ref.monomer_instance_id, ref.motif_id)] += 1

    def motif_usage_count(self, monomer_instance_id: str, motif_id: str) -> int:
        return self._motif_usage[(monomer_instance_id, motif_id)]

    def neighbors(self, monomer_instance_id: str) -> set[str]:
        neighbors: set[str] = set()
        for event in self.reaction_events:
            participants = {ref.monomer_instance_id for ref in event.participants}
            if monomer_instance_id in participants:
                neighbors |= participants
        neighbors.discard(monomer_instance_id)
        return neighbors

    def hyperedges(self) -> tuple[tuple[str, ...], ...]:
        return tuple(
            tuple(ref.monomer_instance_id for ref in event.participants)
            for event in self.reaction_events
        )

    def summary(self) -> dict[str, object]:
        by_template = Counter(event.template_id for event in self.reaction_events)
        return {
            "n_monomer_instances": len(self.monomer_instances),
            "n_reaction_events": len(self.reaction_events),
            "reaction_templates": dict(by_template),
        }

    def validate_all_motifs_known(self, monomer_specs: Mapping[str, MonomerSpec]) -> None:
        for instance in self.monomer_instances.values():
            if instance.monomer_id not in monomer_specs:
                raise KeyError(f"missing monomer spec for {instance.monomer_id!r}")

    def unreacted_motifs(self, monomer_specs: Mapping[str, MonomerSpec]) -> tuple[tuple[str, str], ...]:
        remaining: list[tuple[str, str]] = []
        for instance in self.monomer_instances.values():
            monomer = monomer_specs[instance.monomer_id]
            for motif in monomer.motifs:
                key = (instance.id, motif.id)
                if self._motif_usage[key] == 0:
                    remaining.append(key)
        return tuple(remaining)
