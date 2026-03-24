from __future__ import annotations

from dataclasses import dataclass, field

from .embedding import EmbeddingConfig, PeriodicEmbedder
from .model import Candidate, CandidateEnsemble, MonomerSpec, ReactionTemplate
from .optimizer import ContinuousOptimizer, OptimizerConfig
from .planner import NetPlanner
from .product_graph import PeriodicProductGraph
from .reactions import ReactionLibrary
from .scoring import CandidateScorer
from .search import AssignmentOutcome, AssignmentSolver
from .single_node_topologies_3d import (
    list_supported_3d_single_node_topology_ids,
)
from .topologies import get_topology_hint


@dataclass(frozen=True)
class COFProject:
    monomers: tuple[MonomerSpec, ...]
    allowed_reactions: tuple[str, ...]
    target_dimensionality: str = "2D"
    target_topologies: tuple[str, ...] = ()
    stacking_mode: str = "disabled"
    metadata: dict[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class COFEngineConfig:
    default_layer_spacing: float = 8.0
    default_lateral_span: float = 30.0
    max_candidates: int = 16
    optimization_max_iterations: int = 8


class COFEngine:
    def __init__(
        self,
        reaction_library: ReactionLibrary | None = None,
        config: COFEngineConfig | None = None,
        net_planner: NetPlanner | None = None,
        assignment_solver: AssignmentSolver | None = None,
        embedder: PeriodicEmbedder | None = None,
        scorer: CandidateScorer | None = None,
        optimizer: ContinuousOptimizer | None = None,
    ):
        self.reaction_library = reaction_library or ReactionLibrary.builtin()
        self.config = config or COFEngineConfig()
        self.net_planner = net_planner or NetPlanner()
        self.assignment_solver = assignment_solver or AssignmentSolver()
        self.embedder = embedder or PeriodicEmbedder(
            EmbeddingConfig(
                default_layer_spacing=self.config.default_layer_spacing,
                default_lateral_span=self.config.default_lateral_span,
            )
        )
        self.scorer = scorer or CandidateScorer()
        self.optimizer = optimizer or ContinuousOptimizer(
            OptimizerConfig(max_iterations=self.config.optimization_max_iterations),
            scorer=self.scorer,
        )

    def run(self, project: COFProject) -> CandidateEnsemble:
        if project.stacking_mode != "disabled":
            raise ValueError("stacking exploration is out of scope; set stacking_mode='disabled'")

        templates = self.reaction_library.selected(project.allowed_reactions, project.target_dimensionality)
        if not templates:
            raise ValueError("project selected no reaction templates")

        pair_ensemble = self._run_supported_single_node_pair_project(project, templates)
        if pair_ensemble is not None:
            return pair_ensemble

        monomer_specs = {m.id: m for m in project.monomers}
        net_plans = self.net_planner.propose(
            monomers=project.monomers,
            templates=templates,
            target_dimensionality=project.target_dimensionality,
            target_topologies=project.target_topologies,
        )

        ensemble = CandidateEnsemble()
        candidate_index = 0
        template_map = {t.id: t for t in templates}
        for net_plan in net_plans:
            assignment_plans = self.assignment_solver.build_assignment_plans(net_plan, monomer_specs)
            for assignment_plan in assignment_plans:
                monomer_instances = self.assignment_solver.instantiate_monomers(assignment_plan)
                outcome = self.assignment_solver.solve_events(
                    assignment_plan=assignment_plan,
                    monomer_instances=monomer_instances,
                    monomer_specs=monomer_specs,
                    templates=templates,
                )
                graph = PeriodicProductGraph()
                for instance in monomer_instances:
                    graph.add_monomer(instance)
                for event in outcome.events:
                    graph.add_reaction_event(event, monomer_specs, template_map)

                candidate_index += 1
                candidate = self._make_candidate(
                    project=project,
                    templates=templates,
                    graph=graph,
                    outcome=outcome,
                    candidate_id=f"candidate-{candidate_index}",
                )
                ensemble.add(candidate)
                if len(ensemble.candidates) >= self.config.max_candidates:
                    return ensemble
        return ensemble

    def _run_supported_single_node_pair_project(
        self,
        project: COFProject,
        templates: tuple[ReactionTemplate, ...],
    ) -> CandidateEnsemble | None:
        topology_ids = self._supported_single_node_pair_topology_ids(project, templates)
        if not topology_ids:
            return None

        from .batch import BatchGenerationConfig, BatchStructureGenerator

        generator = BatchStructureGenerator(
            config=BatchGenerationConfig(
                allowed_reactions=project.allowed_reactions,
                target_dimensionality=project.target_dimensionality,
                topology_ids=topology_ids,
                enumerate_all_topologies=True,
                write_cif=False,
                embedding_config=self.embedder.config,
                engine_config=self.config,
                optimizer_config=self.optimizer.config,
            ),
            reaction_library=self.reaction_library,
        )
        summaries, candidates, _attempted_structures = generator.generate_monomer_pair_candidates(
            project.monomers[0],
            project.monomers[1],
            write_cif=False,
        )
        if not candidates:
            if summaries:
                error = summaries[0].metadata.get("error")
                if error is None:
                    error = summaries[0].metadata.get("failed_topologies")
                raise ValueError(str(error or "single-pair topology generation failed"))
            raise ValueError("single-pair topology generation failed")

        ensemble = CandidateEnsemble()
        for candidate in sorted(candidates, key=lambda candidate: candidate.score, reverse=True)[: self.config.max_candidates]:
            ensemble.add(candidate)
        return ensemble

    def _supported_single_node_pair_topology_ids(
        self,
        project: COFProject,
        templates: tuple[ReactionTemplate, ...],
    ) -> tuple[str, ...]:
        if len(project.monomers) != 2:
            return ()
        if len(templates) != 1 or not self.reaction_library.supports_binary_bridge_pair_generation(templates[0]):
            return ()

        connectivities = tuple(len(monomer.motifs) for monomer in project.monomers)
        if connectivities[0] >= 3 and connectivities[1] >= 3:
            pair_mode = "node-node"
        elif min(connectivities) == 2 and max(connectivities) >= 3:
            connectivity = max(connectivities)
            pair_mode = "node-linker"
        else:
            return ()

        if project.target_topologies:
            validated: list[str] = []
            mode_token = (
                f"{max(connectivities)}+2"
                if pair_mode == "node-linker"
                else f"{min(connectivities)}+{max(connectivities)}"
            )
            mode_key = (
                "two_monomer_node_linker_modes"
                if pair_mode == "node-linker"
                else "two_monomer_node_node_modes"
            )
            for topology_id in project.target_topologies:
                hint = get_topology_hint(topology_id)
                if hint.dimensionality != project.target_dimensionality:
                    return ()
                topology_modes = tuple(str(value) for value in hint.metadata.get(mode_key, ()))
                if mode_token not in topology_modes:
                    return ()
                validated.append(topology_id)
            return tuple(validated)

        if project.target_dimensionality == "3D":
            if pair_mode == "node-node" and connectivities[0] != connectivities[1]:
                return ()
            return list_supported_3d_single_node_topology_ids(connectivity, pair_mode=pair_mode)
        return ()

    def _make_candidate(
        self,
        project: COFProject,
        templates: tuple[ReactionTemplate, ...],
        graph: PeriodicProductGraph,
        outcome: AssignmentOutcome,
        candidate_id: str,
    ) -> Candidate:
        monomer_specs = {m.id: m for m in project.monomers}
        template_map = {t.id: t for t in templates}
        embedding = self.embedder.embed(outcome, monomer_specs, template_map)
        optimization = self.optimizer.optimize(outcome, embedding.state, monomer_specs, template_map)
        scoring = self.scorer.score(outcome, optimization.state, monomer_specs, template_map)

        flags: list[str] = []
        if outcome.unreacted_motifs:
            flags.append(f"unreacted_motifs:{len(outcome.unreacted_motifs)}")
        if outcome.assignment_plan.net_plan.topology is None:
            flags.append("no_topology_hint")
        if any(t.topology_role == "ring" for t in templates):
            flags.append("contains_ring_forming_event")
        if project.stacking_mode == "disabled":
            flags.append("stacking_disabled")

        metadata = {
            "graph_summary": graph.summary(),
            "target_topologies": project.target_topologies,
            "allowed_reactions": project.allowed_reactions,
            "stacking_mode": project.stacking_mode,
            "net_plan": {
                "topology": outcome.assignment_plan.net_plan.topology.id
                if outcome.assignment_plan.net_plan.topology is not None
                else None,
                "metadata": dict(outcome.assignment_plan.net_plan.metadata),
            },
            "assignment": dict(outcome.assignment_plan.slot_to_monomer),
            "instance_to_monomer": {
                instance.id: instance.monomer_id for instance in outcome.monomer_instances
            },
            "instance_to_slot": {
                instance.id: instance.metadata.get("slot_id") for instance in outcome.monomer_instances
            },
            "embedding": dict(embedding.metadata),
            "optimization": dict(optimization.metrics),
            "score_breakdown": dict(scoring.breakdown),
            "score_metadata": dict(scoring.metadata),
        }
        return Candidate(
            id=candidate_id,
            score=scoring.total,
            state=optimization.state,
            events=outcome.events,
            flags=tuple(flags),
            metadata=metadata,
        )
