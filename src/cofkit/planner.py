from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping

from .model import MonomerSpec, ReactionTemplate


@dataclass(frozen=True)
class TopologyHint:
    id: str
    dimensionality: str
    node_coordination: tuple[int, ...] = ()
    metadata: Mapping[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class NetPlan:
    topology: TopologyHint | None
    monomer_ids: tuple[str, ...]
    reaction_ids: tuple[str, ...]
    metadata: Mapping[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class AssignmentPlan:
    net_plan: NetPlan
    slot_to_monomer: Mapping[str, str]
    slot_to_conformer: Mapping[str, str] = field(default_factory=dict)
    metadata: Mapping[str, object] = field(default_factory=dict)


class NetPlanner:
    """Proposes topology/net plans from monomer valence and allowed reactions."""

    def __init__(self, topology_repository: object | None = None):
        self._topology_repository = topology_repository

    def propose(
        self,
        monomers: tuple[MonomerSpec, ...],
        templates: tuple[ReactionTemplate, ...],
        target_dimensionality: str,
        target_topologies: tuple[str, ...] = (),
    ) -> tuple[NetPlan, ...]:
        monomer_ids = tuple(m.id for m in monomers)
        reaction_ids = tuple(t.id for t in templates)

        if self._contains_ring_forming_template(templates):
            return (
                NetPlan(
                    topology=None,
                    monomer_ids=monomer_ids,
                    reaction_ids=reaction_ids,
                    metadata={
                        "planning_mode": "ring-forming",
                        "reason": "ring-forming reactions currently bypass explicit topology assignment",
                    },
                ),
            )

        if target_topologies:
            hints = tuple(self._resolve_requested_topologies(target_topologies))
        else:
            hints = tuple(self._infer_repository_topologies(monomers, target_dimensionality))

        compatible_hints = tuple(
            hint for hint in hints if self._is_topology_compatible(hint, monomers, target_dimensionality)
        )
        if compatible_hints:
            return tuple(
                NetPlan(
                    topology=hint,
                    monomer_ids=monomer_ids,
                    reaction_ids=reaction_ids,
                    metadata={
                        "planning_mode": "topology-guided",
                        "connectivities": tuple(len(m.motifs) for m in monomers),
                        "topology_metadata": dict(hint.metadata),
                    },
                )
                for hint in compatible_hints
            )

        return (
            NetPlan(
                topology=None,
                monomer_ids=monomer_ids,
                reaction_ids=reaction_ids,
                metadata={
                    "planning_mode": "topology-free",
                    "reason": "no compatible topology identified",
                    "connectivities": tuple(len(m.motifs) for m in monomers),
                },
            ),
        )

    def _contains_ring_forming_template(self, templates: tuple[ReactionTemplate, ...]) -> bool:
        return any(t.topology_role == "ring" for t in templates)

    def _resolve_requested_topologies(self, target_topologies: tuple[str, ...]) -> tuple[TopologyHint, ...]:
        repository = self._repository()

        return tuple(repository.get_hint(tid) for tid in target_topologies)

    def _infer_repository_topologies(
        self,
        monomers: tuple[MonomerSpec, ...],
        target_dimensionality: str,
    ) -> tuple[TopologyHint, ...]:
        repository = self._repository()
        connectivities = tuple(len(m.motifs) for m in monomers)
        return tuple(
            hint
            for hint in repository.list_hints(
                dimensionality=target_dimensionality,
                node_connectivities=connectivities,
            )
            if self._is_topology_compatible(hint, monomers, target_dimensionality)
        )

    def _is_topology_compatible(
        self,
        hint: TopologyHint,
        monomers: tuple[MonomerSpec, ...],
        target_dimensionality: str,
    ) -> bool:
        if hint.dimensionality != target_dimensionality:
            return False
        if not hint.node_coordination:
            return True

        monomer_connectivities = sorted(len(m.motifs) for m in monomers)
        topology_connectivities = sorted(hint.node_coordination)
        if (
            len(topology_connectivities) == 1
            and len(monomer_connectivities) > 1
            and len(set(monomer_connectivities)) == 1
            and topology_connectivities[0] == monomer_connectivities[0]
        ):
            return True
        return monomer_connectivities == topology_connectivities

    def _repository(self):
        if self._topology_repository is not None:
            return self._topology_repository

        from .topologies import default_topology_repository

        self._topology_repository = default_topology_repository()
        return self._topology_repository
