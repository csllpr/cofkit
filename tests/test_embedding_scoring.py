import math
import sys
import unittest
from math import cos, pi, sin
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit.geometry import matmul, matmul_vec, normalize
from cofkit import (
    AssignmentPlan,
    AssignmentOutcome,
    AssignmentSolver,
    AssemblyState,
    CandidateScorer,
    COFEngine,
    COFProject,
    ContinuousOptimizer,
    EmbeddingConfig,
    Frame,
    MonomerInstance,
    MonomerSpec,
    MotifRef,
    NetPlan,
    NetPlanner,
    PeriodicEmbedder,
    Pose,
    ReactiveMotif,
    ReactionEvent,
    ReactionLibrary,
)
from cofkit.single_node_topologies import resolve_single_node_topology_layout
from cofkit.single_node_topologies_3d import resolve_three_d_single_node_topology_layout


def trigonal_motifs(prefix: str, kind: str, radius: float) -> tuple[ReactiveMotif, ...]:
    motifs = []
    for idx, angle in enumerate((0.0, 2.0 * pi / 3.0, -2.0 * pi / 3.0), start=1):
        origin = (radius * cos(angle), radius * sin(angle), 0.0)
        primary = (cos(angle), sin(angle), 0.0)
        motifs.append(
            ReactiveMotif(
                id=f"{prefix}{idx}",
                kind=kind,
                atom_ids=(idx,),
                frame=Frame(origin=origin, primary=primary, normal=(0.0, 0.0, 1.0)),
            )
        )
    return tuple(motifs)


def irregular_trigonal_motifs(prefix: str, kind: str, radii: tuple[float, float, float]) -> tuple[ReactiveMotif, ...]:
    motifs = []
    for idx, (angle, radius) in enumerate(zip((0.0, 2.0 * pi / 3.0, -2.0 * pi / 3.0), radii), start=1):
        origin = (radius * cos(angle), radius * sin(angle), 0.0)
        primary = (cos(angle), sin(angle), 0.0)
        motifs.append(
            ReactiveMotif(
                id=f"{prefix}{idx}",
                kind=kind,
                atom_ids=(idx,),
                frame=Frame(origin=origin, primary=primary, normal=(0.0, 0.0, 1.0)),
            )
        )
    return tuple(motifs)


def cyclic_motifs(prefix: str, kind: str, count: int, radius: float) -> tuple[ReactiveMotif, ...]:
    motifs = []
    for idx in range(count):
        angle = 2.0 * pi * idx / count
        origin = (radius * cos(angle), radius * sin(angle), 0.0)
        primary = (cos(angle), sin(angle), 0.0)
        motifs.append(
            ReactiveMotif(
                id=f"{prefix}{idx + 1}",
                kind=kind,
                atom_ids=(idx + 1,),
                frame=Frame(origin=origin, primary=primary, normal=(0.0, 0.0, 1.0)),
            )
        )
    return tuple(motifs)


def build_imine_case():
    tri_amine = MonomerSpec(
        id="tapb",
        name="TAPB-like triamine",
        motifs=(
            ReactiveMotif(id="n1", kind="amine", atom_ids=(1,), frame=Frame.xy()),
            ReactiveMotif(id="n2", kind="amine", atom_ids=(2,), frame=Frame.yz()),
            ReactiveMotif(id="n3", kind="amine", atom_ids=(3,), frame=Frame.zx()),
        ),
    )
    tri_aldehyde = MonomerSpec(
        id="tfp",
        name="TFP-like trialdehyde",
        motifs=(
            ReactiveMotif(id="c1", kind="aldehyde", atom_ids=(4,), frame=Frame.xy()),
            ReactiveMotif(id="c2", kind="aldehyde", atom_ids=(5,), frame=Frame.yz()),
            ReactiveMotif(id="c3", kind="aldehyde", atom_ids=(6,), frame=Frame.zx()),
        ),
    )
    specs = {tri_amine.id: tri_amine, tri_aldehyde.id: tri_aldehyde}
    templates = {"imine_bridge": ReactionLibrary.builtin().get("imine_bridge")}
    planner = NetPlanner()
    solver = AssignmentSolver()
    net_plan = planner.propose((tri_amine, tri_aldehyde), tuple(templates.values()), "2D")[0]
    assignment = solver.build_assignment_plans(net_plan, specs)[0]
    instances = solver.instantiate_monomers(assignment)
    outcome = solver.solve_events(
        assignment_plan=assignment,
        monomer_instances=instances,
        monomer_specs=specs,
        templates=tuple(templates.values()),
    )
    return specs, templates, outcome


def build_hcb_case():
    tri_amine = MonomerSpec(
        id="tapb",
        name="TAPB-like triamine",
        motifs=trigonal_motifs("n", "amine", radius=4.5),
    )
    tri_aldehyde = MonomerSpec(
        id="tfb",
        name="TFB-like trialdehyde",
        motifs=trigonal_motifs("c", "aldehyde", radius=2.4),
    )
    specs = {tri_amine.id: tri_amine, tri_aldehyde.id: tri_aldehyde}
    templates = {"imine_bridge": ReactionLibrary.builtin().get("imine_bridge")}
    planner = NetPlanner()
    solver = AssignmentSolver()
    net_plan = planner.propose(
        (tri_amine, tri_aldehyde),
        tuple(templates.values()),
        "2D",
        target_topologies=("hcb",),
    )[0]
    assignment = solver.build_assignment_plans(net_plan, specs)[0]
    instances = solver.instantiate_monomers(assignment)
    outcome = solver.solve_events(
        assignment_plan=assignment,
        monomer_instances=instances,
        monomer_specs=specs,
        templates=tuple(templates.values()),
    )
    return specs, templates, outcome


def build_asymmetric_hcb_case():
    tri_amine = MonomerSpec(
        id="asym_amine",
        name="asymmetric triamine",
        motifs=irregular_trigonal_motifs("n", "amine", radii=(7.0, 4.2, 5.1)),
    )
    tri_aldehyde = MonomerSpec(
        id="asym_aldehyde",
        name="asymmetric trialdehyde",
        motifs=irregular_trigonal_motifs("c", "aldehyde", radii=(2.4, 6.6, 3.3)),
    )
    specs = {tri_amine.id: tri_amine, tri_aldehyde.id: tri_aldehyde}
    templates = {"imine_bridge": ReactionLibrary.builtin().get("imine_bridge")}
    planner = NetPlanner()
    solver = AssignmentSolver()
    net_plan = planner.propose(
        (tri_amine, tri_aldehyde),
        tuple(templates.values()),
        "2D",
        target_topologies=("hcb",),
    )[0]
    assignment = solver.build_assignment_plans(net_plan, specs)[0]
    instances = solver.instantiate_monomers(assignment)
    outcome = solver.solve_events(
        assignment_plan=assignment,
        monomer_instances=instances,
        monomer_specs=specs,
        templates=tuple(templates.values()),
    )
    return specs, templates, outcome


def build_single_node_topology_case(topology_id: str):
    tri_amine = MonomerSpec(
        id="tapb",
        name="TAPB-like triamine",
        motifs=trigonal_motifs("n", "amine", radius=4.5),
    )
    tri_aldehyde = MonomerSpec(
        id="tfb",
        name="TFB-like trialdehyde",
        motifs=trigonal_motifs("c", "aldehyde", radius=2.4),
    )
    specs = {tri_amine.id: tri_amine, tri_aldehyde.id: tri_aldehyde}
    templates = {"imine_bridge": ReactionLibrary.builtin().get("imine_bridge")}
    planner = NetPlanner()
    solver = AssignmentSolver()
    net_plan = planner.propose(
        (tri_amine, tri_aldehyde),
        tuple(templates.values()),
        "2D",
        target_topologies=(topology_id,),
    )[0]
    assignment = solver.build_assignment_plans(net_plan, specs)[0]
    instances = solver.instantiate_monomers(assignment)
    outcome = solver.solve_events(
        assignment_plan=assignment,
        monomer_instances=instances,
        monomer_specs=specs,
        templates=tuple(templates.values()),
    )
    return specs, templates, outcome


def build_single_imine_bridge_case():
    amine = MonomerSpec(
        id="amine",
        name="single amine",
        motifs=(ReactiveMotif(id="n1", kind="amine", atom_ids=(1,), frame=Frame.xy()),),
    )
    aldehyde = MonomerSpec(
        id="aldehyde",
        name="single aldehyde",
        motifs=(ReactiveMotif(id="c1", kind="aldehyde", atom_ids=(2,), frame=Frame.xy()),),
    )
    template = ReactionLibrary.builtin().get("imine_bridge")
    assignment_plan = AssignmentPlan(
        net_plan=NetPlan(topology=None, monomer_ids=(amine.id, aldehyde.id), reaction_ids=(template.id,)),
        slot_to_monomer={"slot1": amine.id, "slot2": aldehyde.id},
    )
    outcome = AssignmentOutcome(
        assignment_plan=assignment_plan,
        monomer_instances=(
            MonomerInstance(id="m1", monomer_id=amine.id),
            MonomerInstance(id="m2", monomer_id=aldehyde.id),
        ),
        events=(
            ReactionEvent(
                id="rxn1",
                template_id=template.id,
                participants=(
                    MotifRef(monomer_instance_id="m1", monomer_id=amine.id, motif_id="n1"),
                    MotifRef(monomer_instance_id="m2", monomer_id=aldehyde.id, motif_id="c1"),
                ),
            ),
        ),
        unreacted_motifs=(),
        consumed_count=2,
    )
    specs = {amine.id: amine, aldehyde.id: aldehyde}
    templates = {template.id: template}
    base_state = AssemblyState(
        cell=((10.0, 0.0, 0.0), (0.0, 10.0, 0.0), (0.0, 0.0, 8.0)),
        monomer_poses={
            "m1": Pose(),
            "m2": Pose(
                translation=(1.3, 0.0, 0.0),
                rotation_matrix=(
                    (-1.0, 0.0, 0.0),
                    (0.0, -1.0, 0.0),
                    (0.0, 0.0, 1.0),
                ),
            ),
        },
        stacking_state="disabled",
    )
    return specs, templates, outcome, base_state


def rotation_about_axis(axis: tuple[float, float, float], angle: float) -> tuple[tuple[float, float, float], ...]:
    x, y, z = normalize(axis)
    c = math.cos(angle)
    s = math.sin(angle)
    one_minus_c = 1.0 - c
    return (
        (
            c + x * x * one_minus_c,
            x * y * one_minus_c - z * s,
            x * z * one_minus_c + y * s,
        ),
        (
            y * x * one_minus_c + z * s,
            c + y * y * one_minus_c,
            y * z * one_minus_c - x * s,
        ),
        (
            z * x * one_minus_c - y * s,
            z * y * one_minus_c + x * s,
            c + z * z * one_minus_c,
        ),
    )


class EmbeddingTests(unittest.TestCase):
    def test_embedder_uses_repository_selected_cell_for_workspace_topology_case(self):
        specs, templates, outcome = build_imine_case()
        embedder = PeriodicEmbedder(EmbeddingConfig(default_lateral_span=24.0, default_layer_spacing=7.5))

        embedding = embedder.embed(outcome, specs, templates)

        self.assertEqual(embedding.metadata["mode"], "topology-guided")
        self.assertEqual(embedding.metadata["topology"], "car")
        self.assertEqual(embedding.metadata["cell_kind"], "orthogonal")
        self.assertEqual(embedding.state.stacking_state, "disabled")
        self.assertEqual(len(embedding.state.monomer_poses), 2)
        self.assertAlmostEqual(embedding.state.cell[1][0], 0.0)
        self.assertGreater(embedding.state.cell[1][1], 0.0)
        self.assertNotEqual(
            embedding.state.monomer_poses["m1"].translation,
            embedding.state.monomer_poses["m2"].translation,
        )

    def test_hcb_embedding_uses_alternating_nodes_and_motif_radial_offsets(self):
        specs, templates, outcome = build_hcb_case()
        embedder = PeriodicEmbedder()

        embedding = embedder.embed(outcome, specs, templates)
        pose_a = embedding.state.monomer_poses["m1"]
        pose_b = embedding.state.monomer_poses["m2"]
        report = CandidateScorer().bridge_geometry_report(outcome, embedding.state, specs, templates)
        center_distance = math.dist(pose_a.translation, pose_b.translation)

        self.assertEqual(embedding.metadata["topology"], "hcb")
        self.assertEqual(embedding.metadata["placement_mode"], "single-node-bipartite")
        self.assertGreater(center_distance, embedding.metadata["target_distance"])
        self.assertAlmostEqual(center_distance, embedding.metadata["reactive_site_distance"], places=6)
        self.assertTrue(all(metric.actual_distance >= 1.29 for metric in report.event_metrics))
        self.assertTrue(all(metric.actual_distance <= 1.31 for metric in report.event_metrics))
        self.assertEqual(embedding.metadata["poses"]["m1"]["sublattice"], "A")
        self.assertEqual(embedding.metadata["poses"]["m2"]["sublattice"], "B")

    def test_hcb_embedding_uses_oblique_cell_for_asymmetric_trigonal_nodes(self):
        specs, templates, outcome = build_asymmetric_hcb_case()
        embedder = PeriodicEmbedder()

        embedding = embedder.embed(outcome, specs, templates)
        report = CandidateScorer().bridge_geometry_report(outcome, embedding.state, specs, templates)

        self.assertEqual(embedding.metadata["topology"], "hcb")
        self.assertEqual(embedding.metadata["placement_mode"], "single-node-bipartite")
        self.assertEqual(embedding.metadata["cell_kind"], "oblique")
        self.assertTrue(all(abs(metric.actual_distance - metric.target_distance) < 1e-6 for metric in report.event_metrics))

    def test_fes_layout_uses_expanded_single_node_topology(self):
        layout = resolve_single_node_topology_layout("fes")

        self.assertEqual(layout.symmetry_orbit_size, 4)
        self.assertTrue(layout.supports_current_builder)
        self.assertEqual(layout.placement_model, "p1-expanded")
        self.assertTrue(layout.supports_node_node)
        self.assertTrue(layout.supports_node_linker)

    def test_sql_layout_uses_expanded_single_node_topology(self):
        layout = resolve_single_node_topology_layout("sql")

        self.assertEqual(layout.connectivity, 4)
        self.assertEqual(layout.symmetry_orbit_size, 1)
        self.assertTrue(layout.supports_current_builder)
        self.assertEqual(layout.placement_model, "p1-expanded")
        self.assertTrue(layout.supports_node_node)
        self.assertTrue(layout.supports_node_linker)

    def test_hxl_layout_allows_node_linker_only(self):
        layout = resolve_single_node_topology_layout("hxl")

        self.assertEqual(layout.connectivity, 6)
        self.assertEqual(layout.symmetry_orbit_size, 1)
        self.assertTrue(layout.supports_current_builder)
        self.assertEqual(layout.placement_model, "p1-expanded")
        self.assertFalse(layout.supports_node_node)
        self.assertTrue(layout.supports_node_linker)

    def test_dia_layout_uses_three_d_single_node_topology(self):
        layout = resolve_three_d_single_node_topology_layout("dia")

        self.assertEqual(layout.connectivity, 4)
        self.assertEqual(layout.symmetry_orbit_size, 2)
        self.assertTrue(layout.supports_current_builder)
        self.assertEqual(layout.placement_model, "p1-two-node-3d")
        self.assertTrue(layout.supports_node_node)
        self.assertTrue(layout.supports_node_linker)

    def test_pcu_layout_allows_node_linker_only_in_three_d(self):
        layout = resolve_three_d_single_node_topology_layout("pcu")

        self.assertEqual(layout.connectivity, 6)
        self.assertEqual(layout.symmetry_orbit_size, 1)
        self.assertTrue(layout.supports_current_builder)
        self.assertEqual(layout.placement_model, "p1-self-edge-3d")
        self.assertFalse(layout.supports_node_node)
        self.assertTrue(layout.supports_node_linker)

    def test_engine_supports_fes_direct_single_pair_generation(self):
        tri_amine = MonomerSpec(
            id="tapb",
            name="TAPB-like triamine",
            motifs=trigonal_motifs("n", "amine", radius=4.5),
        )
        tri_aldehyde = MonomerSpec(
            id="tfb",
            name="TFB-like trialdehyde",
            motifs=trigonal_motifs("c", "aldehyde", radius=2.4),
        )
        project = COFProject(
            monomers=(tri_amine, tri_aldehyde),
            allowed_reactions=("imine_bridge",),
            target_dimensionality="2D",
            target_topologies=("fes",),
        )

        candidate = COFEngine().run(project).top(1)[0]

        self.assertEqual(candidate.metadata["net_plan"]["topology"], "fes")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "single-node-expanded-node-node")
        self.assertEqual(candidate.metadata["embedding"]["topology_family"], "single-node-2d")

    def test_engine_supports_sql_direct_single_pair_generation(self):
        tetra_amine = MonomerSpec(
            id="tetra_amine",
            name="square tetramine",
            motifs=cyclic_motifs("n", "amine", count=4, radius=4.5),
        )
        tetra_aldehyde = MonomerSpec(
            id="tetra_aldehyde",
            name="square tetraaldehyde",
            motifs=cyclic_motifs("c", "aldehyde", count=4, radius=2.4),
        )
        project = COFProject(
            monomers=(tetra_amine, tetra_aldehyde),
            allowed_reactions=("imine_bridge",),
            target_dimensionality="2D",
            target_topologies=("sql",),
        )

        candidate = COFEngine().run(project).top(1)[0]

        self.assertEqual(candidate.metadata["net_plan"]["topology"], "sql")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "single-node-expanded-node-node")
        self.assertEqual(candidate.metadata["embedding"]["topology_family"], "single-node-2d")

    def test_engine_supports_dia_direct_single_pair_generation(self):
        tetra_amine = MonomerSpec(
            id="tetra_amine",
            name="tetramine",
            motifs=cyclic_motifs("n", "amine", count=4, radius=4.5),
        )
        tetra_aldehyde = MonomerSpec(
            id="tetra_aldehyde",
            name="tetraaldehyde",
            motifs=cyclic_motifs("c", "aldehyde", count=4, radius=2.4),
        )
        project = COFProject(
            monomers=(tetra_amine, tetra_aldehyde),
            allowed_reactions=("imine_bridge",),
            target_dimensionality="3D",
            target_topologies=("dia",),
        )

        candidate = COFEngine().run(project).top(1)[0]

        self.assertEqual(candidate.metadata["net_plan"]["topology"], "dia")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "single-node-3d-node-node")
        self.assertEqual(candidate.metadata["embedding"]["topology_family"], "single-node-3d")

    def test_engine_enumerates_requested_two_d_single_node_topologies_for_single_pair(self):
        tri_amine = MonomerSpec(
            id="tapb",
            name="TAPB-like triamine",
            motifs=trigonal_motifs("n", "amine", radius=4.5),
        )
        di_aldehyde = MonomerSpec(
            id="dialdehyde",
            name="dialdehyde",
            motifs=cyclic_motifs("c", "aldehyde", count=2, radius=2.4),
        )
        project = COFProject(
            monomers=(tri_amine, di_aldehyde),
            allowed_reactions=("imine_bridge",),
            target_dimensionality="2D",
            target_topologies=("hcb", "hca", "fes", "fxt"),
        )

        ensemble = COFEngine().run(project)

        self.assertEqual(
            {candidate.metadata["net_plan"]["topology"] for candidate in ensemble.candidates},
            {"hcb", "hca", "fes", "fxt"},
        )


class ScoringTests(unittest.TestCase):
    def test_topology_guided_assignment_scores_above_topology_free_assignment(self):
        specs, templates, outcome = build_imine_case()
        embedder = PeriodicEmbedder()
        scorer = CandidateScorer()

        topology_state = embedder.embed(outcome, specs, templates).state
        topology_score = scorer.score(outcome, topology_state, specs, templates)

        topology_free_assignment = AssignmentPlan(
            net_plan=NetPlan(
                topology=None,
                monomer_ids=outcome.assignment_plan.net_plan.monomer_ids,
                reaction_ids=outcome.assignment_plan.net_plan.reaction_ids,
                metadata={"planning_mode": "topology-free"},
            ),
            slot_to_monomer=outcome.assignment_plan.slot_to_monomer,
            metadata={"assignment_mode": "identity"},
        )
        topology_free_outcome = type(outcome)(
            assignment_plan=topology_free_assignment,
            monomer_instances=outcome.monomer_instances,
            events=outcome.events,
            unreacted_motifs=outcome.unreacted_motifs,
            consumed_count=outcome.consumed_count,
        )
        topology_free_state = embedder.embed(topology_free_outcome, specs, templates).state
        topology_free_score = scorer.score(topology_free_outcome, topology_free_state, specs, templates)

        self.assertGreater(topology_score.total, topology_free_score.total)
        self.assertGreater(topology_score.breakdown["bridge_geometry"], 0.0)
        self.assertEqual(topology_score.breakdown["stacking_penalty"], 0.0)

    def test_optimizer_runs_and_does_not_worsen_imine_bridge_geometry(self):
        specs, templates, outcome = build_imine_case()
        embedder = PeriodicEmbedder()
        scorer = CandidateScorer()
        optimizer = ContinuousOptimizer(scorer=scorer)

        initial_state = embedder.embed(outcome, specs, templates).state
        initial_report = scorer.bridge_geometry_report(outcome, initial_state, specs, templates)

        optimized = optimizer.optimize(outcome, initial_state, specs, templates)
        final_report = scorer.bridge_geometry_report(outcome, optimized.state, specs, templates)

        self.assertTrue(optimized.metrics["enabled"])
        self.assertLessEqual(final_report.total_residual, initial_report.total_residual + 1e-9)
        self.assertGreaterEqual(final_report.score + 1e-9, initial_report.score)

    def test_normal_misalignment_residual_guides_twisted_imine_refinement(self):
        specs, templates, outcome, base_state = build_single_imine_bridge_case()
        scorer = CandidateScorer()
        optimizer = ContinuousOptimizer(scorer=scorer)

        event = outcome.events[0]
        participant = event.participants[1]
        pose = base_state.monomer_poses[participant.monomer_instance_id]
        motif = specs[participant.monomer_id].motif_by_id(participant.motif_id)
        bridge_axis = matmul_vec(pose.rotation_matrix, motif.frame.primary)
        twisted_pose = Pose(
            translation=pose.translation,
            rotation_matrix=matmul(rotation_about_axis(bridge_axis, pi / 2.0), pose.rotation_matrix),
        )
        twisted_poses = dict(base_state.monomer_poses)
        twisted_poses[participant.monomer_instance_id] = twisted_pose
        twisted_state = AssemblyState(
            cell=base_state.cell,
            monomer_poses=twisted_poses,
            torsions=base_state.torsions,
            layer_offsets=base_state.layer_offsets,
            stacking_state=base_state.stacking_state,
        )

        base_report = scorer.bridge_geometry_report(outcome, base_state, specs, templates)
        twisted_report = scorer.bridge_geometry_report(outcome, twisted_state, specs, templates)
        optimized = optimizer.optimize(outcome, twisted_state, specs, templates)
        final_report = scorer.bridge_geometry_report(outcome, optimized.state, specs, templates)

        self.assertAlmostEqual(
            twisted_report.event_metrics[0].alignment_residual,
            base_report.event_metrics[0].alignment_residual,
            places=6,
        )
        self.assertGreater(twisted_report.event_metrics[0].normal_misalignment_residual, 0.45)
        self.assertGreater(twisted_report.total_residual, base_report.total_residual + 0.45)
        self.assertLess(final_report.total_residual, twisted_report.total_residual)
        self.assertLess(
            final_report.event_metrics[0].normal_misalignment_residual,
            twisted_report.event_metrics[0].normal_misalignment_residual,
        )


class EngineIntegrationTests(unittest.TestCase):
    def test_engine_exposes_embedding_and_score_breakdown_metadata(self):
        specs, _, _ = build_imine_case()
        project = COFProject(
            monomers=tuple(specs.values()),
            allowed_reactions=("imine_bridge",),
            target_dimensionality="2D",
        )
        best = COFEngine().run(project).top(1)[0]

        self.assertIn("stacking_disabled", best.flags)
        self.assertEqual(best.metadata["embedding"]["mode"], "topology-guided")
        self.assertIn("optimization", best.metadata)
        self.assertIn("final_residual", best.metadata["optimization"])
        self.assertIn("bridge_geometry", best.metadata["score_breakdown"])
        self.assertIn("bridge_geometry_residual", best.metadata["score_metadata"])
        self.assertIn("normal_misalignment_residual", best.metadata["score_metadata"]["bridge_event_metrics"][0])
        self.assertFalse(best.metadata["score_metadata"]["stacking_considered"])

    def test_engine_rejects_non_disabled_stacking_modes(self):
        specs, _, _ = build_imine_case()
        project = COFProject(
            monomers=tuple(specs.values()),
            allowed_reactions=("imine_bridge",),
            stacking_mode="explore",
        )

        with self.assertRaises(ValueError):
            COFEngine().run(project)


if __name__ == "__main__":
    unittest.main()
