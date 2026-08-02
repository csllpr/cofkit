from __future__ import annotations

from dataclasses import dataclass, replace
from math import atan2, cos, sin
from typing import Mapping

from ..geometry import Frame, Vec3, add, centroid, dot, matmul_vec, norm, rotation_from_frame_to_axes, scale, sub
from ..indexed_topology_layouts import expand_indexed_topology
from ..model import (
    AssemblyState,
    Candidate,
    MonomerInstance,
    MonomerSpec,
    MotifRef,
    Pose,
    ReactionEvent,
)
from ..planner import AssignmentPlan, NetPlan, TopologyHint
from ..product_graph import PeriodicProductGraph
from ..reactions import ReactionLibrary
from ..ring_geometry import (
    RingGeometryOptimizer,
    fractional_to_cartesian,
    ring_geometry_profile,
    ring_geometry_report,
    validate_ring_geometry,
)
from ..search import AssignmentOutcome
from ..single_node_topologies import expand_single_node_topology
from ..stacking import enumerate_candidate_stackings

WORKFLOW_ID = "ring_forming"
IMPLEMENTATION_STATUS = "available"


@dataclass(frozen=True)
class RingFormationConfig:
    topology_id: str = "hcb"
    layer_spacing: float = 3.4
    optimize_geometry: bool = True
    stacking_ids: tuple[str, ...] = ()


@dataclass(frozen=True)
class RingBuildResult:
    candidate: Candidate
    graph: PeriodicProductGraph
    outcome: AssignmentOutcome


class RingFormingStructureGenerator:
    """Builds topology-creating cyclotrimerizations around virtual 3-connected nodes."""

    def __init__(
        self,
        config: RingFormationConfig | None = None,
        reaction_library: ReactionLibrary | None = None,
        optimizer: RingGeometryOptimizer | None = None,
    ) -> None:
        self.config = config or RingFormationConfig()
        self.reaction_library = reaction_library or ReactionLibrary.builtin()
        self.optimizer = optimizer or RingGeometryOptimizer()

    def generate(
        self,
        monomer: MonomerSpec,
        template_id: str,
        *,
        candidate_id: str = "ring-candidate-1",
        topology_id: str | None = None,
    ) -> Candidate:
        candidates = self.generate_candidates(
            monomer,
            template_id,
            candidate_id=candidate_id,
            topology_id=topology_id,
        )
        if len(candidates) != 1:
            raise ValueError(
                "multiple stacking registries were requested; use generate_candidates() to receive all variants"
            )
        return candidates[0]

    def generate_candidates(
        self,
        monomer: MonomerSpec,
        template_id: str,
        *,
        candidate_id: str = "ring-candidate-1",
        topology_id: str | None = None,
        stacking_ids: tuple[str, ...] | None = None,
    ) -> tuple[Candidate, ...]:
        result = self.build(
            monomer,
            template_id,
            candidate_id=candidate_id,
            topology_id=topology_id,
        )
        return enumerate_candidate_stackings(
            result.candidate,
            registry_ids=self.config.stacking_ids if stacking_ids is None else stacking_ids,
        )

    def build(
        self,
        monomer: MonomerSpec,
        template_id: str,
        *,
        candidate_id: str = "ring-candidate-1",
        topology_id: str | None = None,
    ) -> RingBuildResult:
        template = self.reaction_library.get(template_id)
        profile = self.reaction_library.linkage_profile(template)
        if profile is None or not profile.supports_ring_forming_generation:
            raise ValueError(f"template {template_id!r} is not configured for ring-forming generation")
        if template.arity != 3 or template.topology_role != "ring":
            raise ValueError(f"template {template_id!r} is not a three-participant ring-forming template")
        motif_kind = profile.ring_participant_motif_kind
        if any(motif.kind != motif_kind for motif in monomer.motifs):
            raise ValueError(
                f"monomer {monomer.id!r} must expose only {motif_kind!r} motifs for {template_id!r}"
            )
        if not monomer.atom_positions:
            raise ValueError("ring-forming atomistic generation requires monomer atom positions")

        if len(monomer.motifs) != 2:
            return self._build_indexed_node_precursor(
                monomer,
                template_id,
                candidate_id=candidate_id,
                topology_id=topology_id or self.config.topology_id,
            )

        selected_topology_id = topology_id or self.config.topology_id
        expanded = expand_single_node_topology(selected_topology_id)
        if expanded.connectivity != 3:
            raise ValueError(
                f"ring product nodes require connectivity 3, but topology {selected_topology_id!r} "
                f"has connectivity {expanded.connectivity}"
            )

        ring_profile = ring_geometry_profile(template_id)
        first_motif, second_motif = monomer.motifs
        precursor_span = norm(sub(second_motif.frame.origin, first_motif.frame.origin))
        if precursor_span < 1e-6:
            raise ValueError(f"monomer {monomer.id!r} has coincident reactive motif origins")
        base_edge_lengths = tuple(norm(edge.base_vector) for edge in expanded.edge_sites)
        base_edge_length = sum(base_edge_lengths) / len(base_edge_lengths)
        target_center_separation = precursor_span + 2.0 * ring_profile.participant_radius
        lateral_scale = target_center_separation / base_edge_length
        cell = (
            scale(expanded.cell[0], lateral_scale),
            scale(expanded.cell[1], lateral_scale),
            (0.0, 0.0, self.config.layer_spacing),
        )

        source_frame = Frame(
            origin=first_motif.frame.origin,
            primary=sub(second_motif.frame.origin, first_motif.frame.origin),
            normal=_monomer_normal(monomer),
        )
        node_by_id = {node.id: node for node in expanded.node_sites}
        instances: list[MonomerInstance] = []
        poses: dict[str, Pose] = {}
        endpoint_refs: dict[str, MotifRef] = {}
        for edge_index, edge in enumerate(expanded.edge_sites, start=1):
            instance_id = f"p{edge_index}"
            start_center = fractional_to_cartesian(node_by_id[edge.start_node_id].fractional_position, cell)
            end_center = add(
                fractional_to_cartesian(node_by_id[edge.end_node_id].fractional_position, cell),
                fractional_to_cartesian(tuple(float(v) for v in edge.end_image), cell),
            )
            direction = scale(sub(end_center, start_center), 1.0 / norm(sub(end_center, start_center)))
            target_first = add(start_center, scale(direction, ring_profile.participant_radius))
            target_second = sub(end_center, scale(direction, ring_profile.participant_radius))
            rotation = rotation_from_frame_to_axes(source_frame, direction, (0.0, 0.0, 1.0))
            translation = sub(target_first, matmul_vec(rotation, first_motif.frame.origin))
            placed_second = add(matmul_vec(rotation, second_motif.frame.origin), translation)
            if norm(sub(placed_second, target_second)) > 2e-4:
                raise ValueError("failed to place both precursor reactive sites on the virtual-node edge")
            instances.append(
                MonomerInstance(
                    id=instance_id,
                    monomer_id=monomer.id,
                    conformer_id=monomer.conformer_ids[0] if monomer.conformer_ids else None,
                    metadata={"topology_edge_id": edge.id},
                )
            )
            poses[instance_id] = Pose(translation=translation, rotation_matrix=rotation)
            endpoint_refs[f"{edge.id}:start"] = MotifRef(instance_id, monomer.id, first_motif.id)
            endpoint_refs[f"{edge.id}:end"] = MotifRef(
                instance_id,
                monomer.id,
                second_motif.id,
                periodic_image=tuple(-value for value in edge.end_image),
            )

        events: list[ReactionEvent] = []
        for event_index, node in enumerate(expanded.node_sites, start=1):
            participants = tuple(endpoint_refs[endpoint_id] for endpoint_id in node.edge_ids)
            if len(participants) != 3:
                raise ValueError(f"virtual ring node {node.id!r} has {len(participants)} incidences instead of 3")
            if len({ref.monomer_instance_id for ref in participants}) != 3:
                raise ValueError(f"virtual ring node {node.id!r} does not use three distinct precursor instances")
            events.append(
                ReactionEvent(
                    id=f"r{event_index}",
                    template_id=template_id,
                    participants=participants,
                    product_state=template.product_name,
                    metadata={
                        "topology_node_id": node.id,
                        "ring_center_fractional": node.fractional_position,
                        "ring_normal": (0.0, 0.0, 1.0),
                        "ring_atom_bond_length": ring_profile.ring_atom_bond_length,
                        "participant_order": "counterclockwise",
                    },
                )
            )

        topology_hint = TopologyHint(
            id=selected_topology_id,
            dimensionality="2D",
            node_coordination=(3,),
            metadata={"virtual_node_product": template.product_name},
        )
        net_plan = NetPlan(
            topology=topology_hint,
            monomer_ids=(monomer.id,),
            reaction_ids=(template_id,),
            metadata={
                "planning_mode": "virtual-node-topology",
                "precursor_connectivity": 2,
                "product_node_connectivity": 3,
            },
        )
        assignment_plan = AssignmentPlan(
            net_plan=net_plan,
            slot_to_monomer={"edge_precursor": monomer.id},
            slot_to_conformer={
                "edge_precursor": monomer.conformer_ids[0]
            } if monomer.conformer_ids else {},
            metadata={"assignment_mode": "virtual-ring-nodes/precursor-edges"},
        )
        outcome = AssignmentOutcome(
            assignment_plan=assignment_plan,
            monomer_instances=tuple(instances),
            events=tuple(events),
            unreacted_motifs=(),
            consumed_count=sum(len(event.participants) for event in events),
        )
        graph = PeriodicProductGraph()
        for instance in instances:
            graph.add_monomer(instance)
        for event in events:
            graph.add_reaction_event(event, {monomer.id: monomer}, {template.id: template})
        if graph.unreacted_motifs({monomer.id: monomer}):
            raise ValueError("virtual-node assembly left unreacted precursor motifs")

        initial_state = AssemblyState(cell=cell, monomer_poses=poses, stacking_state="disabled")
        optimization = (
            self.optimizer.optimize(tuple(events), initial_state, {monomer.id: monomer})
            if self.config.optimize_geometry
            else None
        )
        state = optimization.state if optimization is not None else initial_state
        geometry = ring_geometry_report(tuple(events), state, {monomer.id: monomer})
        validation = validate_ring_geometry(tuple(events), state, {monomer.id: monomer})
        flags = ["contains_ring_forming_event", "stacking_disabled"]
        if validation.classification != "accepted":
            flags.append("ring_geometry_rejected")
        candidate = Candidate(
            id=candidate_id,
            score=20.0 - geometry.total_residual,
            state=state,
            events=tuple(events),
            flags=tuple(flags),
            metadata={
                "workflow": WORKFLOW_ID,
                "graph_summary": graph.summary(),
                "net_plan": {"topology": selected_topology_id, "metadata": dict(net_plan.metadata)},
                "assignment": dict(assignment_plan.slot_to_monomer),
                "instance_to_monomer": {instance.id: monomer.id for instance in instances},
                "instance_metadata": {instance.id: dict(instance.metadata) for instance in instances},
                "score_metadata": {
                    "ring_geometry": geometry.as_dict(),
                    "n_unreacted_motifs": 0,
                },
                "optimization": dict(optimization.metrics) if optimization is not None else {"enabled": False},
                "ring_validation": {
                    "classification": validation.classification,
                    "reasons": validation.reasons,
                    "metrics": dict(validation.metrics),
                },
                "topology_assignment_mode": "virtual_node_topology",
            },
        )
        candidate = self._annotate_ring_embedding(
            candidate,
            monomer,
            cell_kind=expanded.metric_family,
            placement_mode="virtual-ring-nodes/precursor-edges",
        )
        return RingBuildResult(candidate=candidate, graph=graph, outcome=outcome)

    def _build_indexed_node_precursor(
        self,
        monomer: MonomerSpec,
        template_id: str,
        *,
        candidate_id: str,
        topology_id: str,
    ) -> RingBuildResult:
        if len(monomer.motifs) < 3:
            raise ValueError("indexed ring-forming node placement requires a precursor with at least three motifs")
        if len(monomer.motifs) == 3:
            try:
                expanded = expand_single_node_topology(topology_id)
                expanded_dimensionality = "2D"
            except (KeyError, ValueError):
                expanded = expand_indexed_topology(topology_id)
                expanded_dimensionality = expanded.dimensionality
        else:
            expanded = expand_indexed_topology(topology_id)
            expanded_dimensionality = expanded.dimensionality
        if expanded_dimensionality != "2D":
            raise ValueError("indexed ring-forming generation currently supports 2D topologies")
        template = self.reaction_library.get(template_id)
        ring_profile = ring_geometry_profile(template_id)
        precursor_connectivity = len(monomer.motifs)

        if precursor_connectivity == 3:
            sublattices = {node.sublattice for node in expanded.node_sites}
            if not expanded.is_bipartite or sublattices != {0, 1}:
                raise ValueError(
                    f"3-connected precursor topology {topology_id!r} must expose two bipartite sublattices "
                    "so precursor and virtual-ring roles are unambiguous"
                )
            real_nodes = tuple(node for node in expanded.node_sites if node.sublattice == 0)
            ring_nodes = tuple(node for node in expanded.node_sites if node.sublattice == 1)
        else:
            unsupported = sorted(
                {node.connectivity for node in expanded.node_sites} - {3, precursor_connectivity}
            )
            if unsupported:
                raise ValueError(
                    f"topology {topology_id!r} contains unsupported node connectivities {tuple(unsupported)!r}; "
                    f"expected only 3 and {precursor_connectivity}"
                )
            real_nodes = tuple(node for node in expanded.node_sites if node.connectivity == precursor_connectivity)
            ring_nodes = tuple(node for node in expanded.node_sites if node.connectivity == 3)
        if not real_nodes or not ring_nodes:
            raise ValueError(
                f"topology {topology_id!r} does not contain both precursor nodes and virtual 3-connected ring nodes"
            )
        real_node_ids = {node.id for node in real_nodes}
        ring_node_ids = {node.id for node in ring_nodes}
        if any(
            {edge.start_node_id, edge.end_node_id}.issubset(real_node_ids)
            or {edge.start_node_id, edge.end_node_id}.issubset(ring_node_ids)
            for edge in expanded.edge_sites
        ):
            raise ValueError(f"topology {topology_id!r} has edges that do not alternate precursor and virtual-ring nodes")

        placements = {
            node.id: self._planar_node_placement(monomer, node.directions, node.edge_ids)
            for node in real_nodes
        }
        node_by_id = {node.id: node for node in expanded.node_sites}
        target_edge_vectors: dict[str, Vec3] = {}
        for edge in expanded.edge_sites:
            if edge.start_node_id in real_node_ids:
                real_node_id = edge.start_node_id
                endpoint_key = f"{edge.id}:start"
            else:
                real_node_id = edge.end_node_id
                endpoint_key = f"{edge.id}:end"
            offset = placements[real_node_id][2][endpoint_key]
            if offset <= 0.1:
                raise ValueError(
                    f"precursor {monomer.id!r} does not point motif {placements[real_node_id][1][endpoint_key]!r} "
                    f"toward topology edge {edge.id!r}"
                )
            direction = scale(edge.base_vector, 1.0 / norm(edge.base_vector))
            target_edge_vectors[edge.id] = scale(direction, offset + ring_profile.participant_radius)
        cell = self._fit_planar_cell(expanded, target_edge_vectors)

        instances: list[MonomerInstance] = []
        poses: dict[str, Pose] = {}
        instance_by_node: dict[str, str] = {}
        for index, node in enumerate(real_nodes, start=1):
            instance_id = f"p{index}"
            rotation, _motif_by_edge, _offsets, local_center = placements[node.id]
            node_center = fractional_to_cartesian(node.fractional_position, cell)
            translation = sub(node_center, matmul_vec(rotation, local_center))
            instance_by_node[node.id] = instance_id
            instances.append(
                MonomerInstance(
                    id=instance_id,
                    monomer_id=monomer.id,
                    conformer_id=monomer.conformer_ids[0] if monomer.conformer_ids else None,
                    metadata={"topology_node_id": node.id},
                )
            )
            poses[instance_id] = Pose(translation=translation, rotation_matrix=rotation)

        ring_endpoint_refs: dict[str, MotifRef] = {}
        for edge in expanded.edge_sites:
            if edge.start_node_id in real_node_ids:
                real_node_id = edge.start_node_id
                ring_endpoint_key = f"{edge.id}:end"
                real_endpoint_key = f"{edge.id}:start"
                image = tuple(-value for value in edge.end_image)
            else:
                real_node_id = edge.end_node_id
                ring_endpoint_key = f"{edge.id}:start"
                real_endpoint_key = f"{edge.id}:end"
                image = edge.end_image
            ring_endpoint_refs[ring_endpoint_key] = MotifRef(
                instance_by_node[real_node_id],
                monomer.id,
                placements[real_node_id][1][real_endpoint_key],
                periodic_image=image,
            )

        events: list[ReactionEvent] = []
        for event_index, node in enumerate(ring_nodes, start=1):
            participants = tuple(ring_endpoint_refs[endpoint_key] for endpoint_key in node.edge_ids)
            if len(participants) != 3:
                raise ValueError(f"virtual ring node {node.id!r} has {len(participants)} incidences instead of 3")
            physical_copies = {(ref.monomer_instance_id, ref.periodic_image) for ref in participants}
            if len(physical_copies) != 3:
                raise ValueError(f"virtual ring node {node.id!r} does not use three distinct periodic precursor copies")
            events.append(
                ReactionEvent(
                    id=f"r{event_index}",
                    template_id=template_id,
                    participants=participants,
                    product_state=template.product_name,
                    metadata={
                        "topology_node_id": node.id,
                        "ring_center_fractional": node.fractional_position,
                        "ring_normal": (0.0, 0.0, 1.0),
                        "ring_atom_bond_length": ring_profile.ring_atom_bond_length,
                        "participant_order": "counterclockwise",
                    },
                )
            )

        topology_hint = TopologyHint(
            id=topology_id,
            dimensionality="2D",
            node_coordination=(3, precursor_connectivity),
            metadata={"virtual_node_product": template.product_name},
        )
        net_plan = NetPlan(
            topology=topology_hint,
            monomer_ids=(monomer.id,),
            reaction_ids=(template_id,),
            metadata={
                "planning_mode": "virtual-node-indexed-topology",
                "precursor_connectivity": precursor_connectivity,
                "product_node_connectivity": 3,
            },
        )
        assignment_plan = AssignmentPlan(
            net_plan=net_plan,
            slot_to_monomer={"precursor_node": monomer.id},
            metadata={"assignment_mode": "precursor-nodes/virtual-ring-nodes"},
        )
        outcome = AssignmentOutcome(
            assignment_plan=assignment_plan,
            monomer_instances=tuple(instances),
            events=tuple(events),
            unreacted_motifs=(),
            consumed_count=sum(len(event.participants) for event in events),
        )
        graph = PeriodicProductGraph()
        for instance in instances:
            graph.add_monomer(instance)
        for event in events:
            graph.add_reaction_event(event, {monomer.id: monomer}, {template.id: template})
        unreacted = graph.unreacted_motifs({monomer.id: monomer})
        if unreacted:
            raise ValueError(f"indexed virtual-node assembly left unreacted precursor motifs: {unreacted!r}")

        state = AssemblyState(cell=cell, monomer_poses=poses, stacking_state="disabled")
        optimization = self.optimizer.optimize(tuple(events), state, {monomer.id: monomer}) if self.config.optimize_geometry else None
        state = optimization.state if optimization is not None else state
        geometry = ring_geometry_report(tuple(events), state, {monomer.id: monomer})
        validation = validate_ring_geometry(tuple(events), state, {monomer.id: monomer})
        flags = ["contains_ring_forming_event", "stacking_disabled", "indexed_virtual_node_topology"]
        if validation.classification != "accepted":
            flags.append("ring_geometry_rejected")
        candidate = Candidate(
            id=candidate_id,
            score=20.0 - geometry.total_residual,
            state=state,
            events=tuple(events),
            flags=tuple(flags),
            metadata={
                "workflow": WORKFLOW_ID,
                "graph_summary": graph.summary(),
                "net_plan": {"topology": topology_id, "metadata": dict(net_plan.metadata)},
                "assignment": dict(assignment_plan.slot_to_monomer),
                "instance_to_monomer": {instance.id: monomer.id for instance in instances},
                "score_metadata": {"ring_geometry": geometry.as_dict(), "n_unreacted_motifs": 0},
                "optimization": dict(optimization.metrics) if optimization is not None else {"enabled": False},
                "ring_validation": {
                    "classification": validation.classification,
                    "reasons": validation.reasons,
                    "metrics": dict(validation.metrics),
                },
                "topology_assignment_mode": "virtual_node_topology",
            },
        )
        candidate = self._annotate_ring_embedding(
            candidate,
            monomer,
            cell_kind=expanded.metric_family,
            placement_mode="precursor-nodes/virtual-ring-nodes",
        )
        return RingBuildResult(candidate=candidate, graph=graph, outcome=outcome)

    def _annotate_ring_embedding(
        self,
        candidate: Candidate,
        monomer: MonomerSpec,
        *,
        cell_kind: str,
        placement_mode: str,
    ) -> Candidate:
        from ..reaction_realization import ReactionRealizer

        try:
            realization = ReactionRealizer().realize(
                candidate,
                {monomer.id: monomer},
                candidate.metadata["instance_to_monomer"],
            )
        except (KeyError, ValueError):
            # Coarse/synthetic MonomerSpec inputs can still be embedded even when
            # they do not carry enough atom metadata for product realization.
            realization = None
        z_values: list[float] = []
        for instance_id, pose in candidate.state.monomer_poses.items():
            realized_atoms = None if realization is None else realization.atoms_by_instance.get(instance_id)
            local_positions = (
                tuple(atom.local_position for atom in realized_atoms)
                if realized_atoms is not None
                else monomer.atom_positions
            )
            z_values.extend(add(matmul_vec(pose.rotation_matrix, position), pose.translation)[2] for position in local_positions)
        layer_z_span = max(z_values) - min(z_values) if z_values else 0.0
        pose_details = {
            instance_id: {
                "translation": pose.translation,
                "rotation_matrix": pose.rotation_matrix,
            }
            for instance_id, pose in candidate.state.monomer_poses.items()
        }
        embedding = {
            "mode": "virtual-node-topology",
            "topology": candidate.metadata["net_plan"]["topology"],
            "cell_kind": cell_kind,
            "placement_mode": placement_mode,
            "layer_z_span": layer_z_span,
            "layer_z_span_mode": "atomistic_product" if realization is not None else "precursor_coordinates",
            "stacking_enabled": False,
            "poses": pose_details,
        }
        return replace(
            candidate,
            metadata={
                **dict(candidate.metadata),
                "embedding": embedding,
                "stacking_mode": "disabled",
            },
        )

    def _planar_node_placement(
        self,
        monomer: MonomerSpec,
        target_directions: tuple[Vec3, ...],
        endpoint_keys: tuple[str, ...],
    ) -> tuple[tuple[Vec3, Vec3, Vec3], dict[str, str], dict[str, float], Vec3]:
        local_center = centroid(motif.frame.origin for motif in monomer.motifs)
        motifs = tuple(
            sorted(
                monomer.motifs,
                key=lambda motif: atan2(
                    motif.frame.origin[1] - local_center[1],
                    motif.frame.origin[0] - local_center[0],
                ),
            )
        )
        if len(motifs) != len(target_directions):
            raise ValueError("precursor motif count does not match topology-node incidence")
        best = None
        for reverse in (False, True):
            ordered = tuple(reversed(motifs)) if reverse else motifs
            for shift_index in range(len(ordered)):
                shifted = ordered[shift_index:] + ordered[:shift_index]
                local_angles = tuple(
                    atan2(motif.frame.origin[1] - local_center[1], motif.frame.origin[0] - local_center[0])
                    for motif in shifted
                )
                target_angles = tuple(atan2(direction[1], direction[0]) for direction in target_directions)
                differences = tuple(target - local for target, local in zip(target_angles, local_angles))
                angle = atan2(sum(sin(value) for value in differences), sum(cos(value) for value in differences))
                score = sum(cos(target - (local + angle)) for target, local in zip(target_angles, local_angles))
                if best is None or score > best[0]:
                    best = (score, angle, shifted)
        assert best is not None
        angle = best[1]
        rotation = (
            (cos(angle), -sin(angle), 0.0),
            (sin(angle), cos(angle), 0.0),
            (0.0, 0.0, 1.0),
        )
        motif_by_edge = {endpoint: motif.id for endpoint, motif in zip(endpoint_keys, best[2])}
        offsets = {
            endpoint: dot(
                matmul_vec(rotation, sub(motif.frame.origin, local_center)),
                direction,
            )
            for endpoint, motif, direction in zip(endpoint_keys, best[2], target_directions)
        }
        return rotation, motif_by_edge, offsets, local_center

    def _fit_planar_cell(self, expanded, target_edge_vectors: Mapping[str, Vec3]):
        node_positions = {node.id: node.fractional_position for node in expanded.node_sites}
        s_xx = s_xy = s_yy = 0.0
        rhs_x_0 = rhs_x_1 = rhs_y_0 = rhs_y_1 = 0.0
        for edge in expanded.edge_sites:
            delta_x = node_positions[edge.end_node_id][0] + edge.end_image[0] - node_positions[edge.start_node_id][0]
            delta_y = node_positions[edge.end_node_id][1] + edge.end_image[1] - node_positions[edge.start_node_id][1]
            target = target_edge_vectors[edge.id]
            s_xx += delta_x * delta_x
            s_xy += delta_x * delta_y
            s_yy += delta_y * delta_y
            rhs_x_0 += delta_x * target[0]
            rhs_x_1 += delta_y * target[0]
            rhs_y_0 += delta_x * target[1]
            rhs_y_1 += delta_y * target[1]
        determinant = s_xx * s_yy - s_xy * s_xy
        if abs(determinant) < 1e-10:
            raise ValueError(f"topology {expanded.topology_id!r} does not determine a non-singular 2D cell")
        return (
            (
                (rhs_x_0 * s_yy - rhs_x_1 * s_xy) / determinant,
                (rhs_y_0 * s_yy - rhs_y_1 * s_xy) / determinant,
                0.0,
            ),
            (
                (s_xx * rhs_x_1 - s_xy * rhs_x_0) / determinant,
                (s_xx * rhs_y_1 - s_xy * rhs_y_0) / determinant,
                0.0,
            ),
            (0.0, 0.0, self.config.layer_spacing),
        )


def _monomer_normal(monomer: MonomerSpec) -> Vec3:
    heavy_atom_positions = tuple(
        position
        for symbol, position in zip(monomer.atom_symbols, monomer.atom_positions)
        if symbol != "H"
    )
    fitted_normal = _smallest_covariance_axis(heavy_atom_positions)
    if fitted_normal is not None:
        return fitted_normal
    normals = tuple(motif.frame.normal for motif in monomer.motifs if norm(motif.frame.normal) > 1e-8)
    if not normals:
        return (0.0, 0.0, 1.0)
    summed = (
        sum(vector[0] for vector in normals),
        sum(vector[1] for vector in normals),
        sum(vector[2] for vector in normals),
    )
    return summed if norm(summed) > 1e-8 else normals[0]


def _smallest_covariance_axis(points: tuple[Vec3, ...]) -> Vec3 | None:
    """Returns the least-variance axis of a small point cloud via Jacobi rotations."""
    if len(points) < 3:
        return None
    center = centroid(points)
    covariance = [[0.0, 0.0, 0.0] for _ in range(3)]
    for point in points:
        offset = sub(point, center)
        for row in range(3):
            for column in range(3):
                covariance[row][column] += offset[row] * offset[column]
    eigenvectors = [[1.0 if row == column else 0.0 for column in range(3)] for row in range(3)]
    for _ in range(24):
        row, column = max(((0, 1), (0, 2), (1, 2)), key=lambda pair: abs(covariance[pair[0]][pair[1]]))
        if abs(covariance[row][column]) < 1e-12:
            break
        angle = 0.5 * atan2(
            2.0 * covariance[row][column],
            covariance[column][column] - covariance[row][row],
        )
        cosine = cos(angle)
        sine = sin(angle)
        rotation = [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
        rotation[row][row] = cosine
        rotation[column][column] = cosine
        rotation[row][column] = sine
        rotation[column][row] = -sine
        covariance = _matrix_multiply(_matrix_transpose(rotation), _matrix_multiply(covariance, rotation))
        eigenvectors = _matrix_multiply(eigenvectors, rotation)
    axis_index = min(range(3), key=lambda index: covariance[index][index])
    axis = tuple(eigenvectors[row][axis_index] for row in range(3))
    if norm(axis) < 1e-8:
        return None
    normalized = scale(axis, 1.0 / norm(axis))
    return scale(normalized, -1.0) if normalized[2] < 0.0 else normalized


def _matrix_multiply(left, right):
    return [
        [sum(left[row][inner] * right[inner][column] for inner in range(3)) for column in range(3)]
        for row in range(3)
    ]


def _matrix_transpose(matrix):
    return [[matrix[column][row] for column in range(3)] for row in range(3)]


__all__ = [
    "IMPLEMENTATION_STATUS",
    "RingBuildResult",
    "RingFormationConfig",
    "RingFormingStructureGenerator",
    "WORKFLOW_ID",
]
