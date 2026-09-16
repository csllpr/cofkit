from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping

from .model import MonomerSpec, ReactionTemplate
from .node_shape import classify_monomer_node_shape, classify_topology_node_shape, shapes_compatible


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

    def __init__(
        self,
        topology_repository: object | None = None,
        *,
        shape_aware_topology_filter: bool = True,
    ):
        self._topology_repository = topology_repository
        self.shape_aware_topology_filter = shape_aware_topology_filter

    def propose(
        self,
        monomers: tuple[MonomerSpec, ...],
        templates: tuple[ReactionTemplate, ...],
        target_dimensionality: str,
        target_topologies: tuple[str, ...] = (),
    ) -> tuple[NetPlan, ...]:
        monomer_ids = tuple(m.id for m in monomers)
        reaction_ids = tuple(t.id for t in templates)

        if any(template.topology_role == "ring" for template in templates):
            topology_ids = target_topologies or ("hcb",)
            hints = tuple(self._resolve_requested_topologies(topology_ids))
            compatible = tuple(
                hint
                for hint in hints
                if hint.dimensionality == target_dimensionality and 3 in hint.node_coordination
            )
            if not compatible:
                raise ValueError(
                    "ring-forming planning requires a target topology with 3-connected virtual product nodes"
                )
            return tuple(
                NetPlan(
                    topology=hint,
                    monomer_ids=monomer_ids,
                    reaction_ids=reaction_ids,
                    metadata={
                        "planning_mode": "virtual-node-topology",
                        "precursor_connectivities": tuple(len(monomer.motifs) for monomer in monomers),
                        "product_node_connectivity": 3,
                    },
                )
                for hint in compatible
            )

        if target_topologies:
            hints = tuple(self._resolve_requested_topologies(target_topologies))
        else:
            hints = tuple(self._infer_repository_topologies(monomers, target_dimensionality))

        # Explicitly requested topologies are planned even when the node-shape
        # classifier disagrees; the conflict is reported as a warning instead.
        explicit_request = bool(target_topologies)
        compatible_hints = tuple(
            hint
            for hint in hints
            if self._is_topology_compatible(
                hint,
                monomers,
                target_dimensionality,
                shape_aware=not explicit_request,
            )
        )
        shape_warnings = (
            self._requested_shape_warnings(compatible_hints, monomers) if explicit_request else ()
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
                        **({"shape_warnings": shape_warnings} if shape_warnings else {}),
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
        *,
        shape_aware: bool = True,
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
            connectivity_match = True
        else:
            connectivity_match = monomer_connectivities == topology_connectivities
        if not connectivity_match:
            return False

        if not shape_aware or not self.shape_aware_topology_filter:
            return True

        topology_shape = classify_topology_node_shape(hint.id)
        for monomer in monomers:
            if len(monomer.motifs) != 4:
                continue
            if shapes_compatible(classify_monomer_node_shape(monomer), topology_shape) is False:
                return False
        return True

    def _requested_shape_warnings(
        self,
        hints: tuple[TopologyHint, ...],
        monomers: tuple[MonomerSpec, ...],
    ) -> tuple[str, ...]:
        """Non-blocking node-shape conflicts for explicitly requested topologies."""
        if not self.shape_aware_topology_filter:
            return ()
        warnings: list[str] = []
        for hint in hints:
            if not hint.node_coordination:
                continue
            topology_shape = classify_topology_node_shape(hint.id)
            for monomer in monomers:
                if len(monomer.motifs) != 4:
                    continue
                monomer_shape = classify_monomer_node_shape(monomer)
                if shapes_compatible(monomer_shape, topology_shape) is False:
                    warnings.append(
                        f"requested topology {hint.id!r} conflicts with the classified tetratopic "
                        f"node shape and was kept anyway: monomer {monomer.id!r} node shape "
                        f"{monomer_shape.label!r} is incompatible with topology {hint.id!r} "
                        f"node shape {topology_shape.label!r}"
                    )
        return tuple(warnings)

    def _repository(self):
        if self._topology_repository is not None:
            return self._topology_repository

        from .topologies import default_topology_repository

        self._topology_repository = default_topology_repository()
        return self._topology_repository
