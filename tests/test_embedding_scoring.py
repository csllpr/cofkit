import math
import sys
import unittest
from math import cos, pi, sin, sqrt
from pathlib import Path


from cofkit.geometry import matmul, matmul_vec, normalize, rotation_from_frame_to_axes
from cofkit.linkage_geometry import effective_motif_origin
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
    OptimizerConfig,
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


def tetrahedral_motifs(prefix: str, kind: str, radius: float) -> tuple[ReactiveMotif, ...]:
    corners = ((1.0, 1.0, 1.0), (1.0, -1.0, -1.0), (-1.0, 1.0, -1.0), (-1.0, -1.0, 1.0))
    motifs = []
    for idx, corner in enumerate(corners):
        length = sqrt(3.0)
        origin = tuple(radius * component / length for component in corner)
        primary = tuple(component / length for component in corner)
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


def build_offcenter_imine_bridge_case():
    amine = MonomerSpec(
        id="amine",
        name="single amine, off-center motif",
        motifs=(
            ReactiveMotif(
                id="n1",
                kind="amine",
                atom_ids=(1,),
                frame=Frame(origin=(0.4, 0.2, 0.0), primary=(1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)),
            ),
        ),
    )
    aldehyde = MonomerSpec(
        id="aldehyde",
        name="single aldehyde, off-center motif",
        motifs=(
            ReactiveMotif(
                id="c1",
                kind="aldehyde",
                atom_ids=(2,),
                frame=Frame(origin=(-0.5, 0.3, 0.0), primary=(-1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)),
            ),
        ),
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
                translation=(1.8, 0.1, 0.0),
                rotation_matrix=rotation_about_axis((0.0, 0.0, 1.0), 0.4),
            ),
        },
        stacking_state="disabled",
    )
    return specs, templates, outcome, base_state


def build_two_event_linker_case():
    linker = MonomerSpec(
        id="linker",
        name="ditopic dialdehyde, shared motif origin",
        motifs=(
            ReactiveMotif(
                id="c1",
                kind="aldehyde",
                atom_ids=(1,),
                frame=Frame(origin=(0.0, 0.0, 0.0), primary=(-1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)),
            ),
            ReactiveMotif(
                id="c2",
                kind="aldehyde",
                atom_ids=(2,),
                frame=Frame(origin=(0.0, 0.0, 0.0), primary=(1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)),
            ),
        ),
    )
    amine_b = MonomerSpec(
        id="amine_b",
        name="amine at the c1 side",
        motifs=(
            ReactiveMotif(
                id="n1",
                kind="amine",
                atom_ids=(1,),
                frame=Frame(origin=(0.0, 0.0, 0.0), primary=(1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)),
            ),
        ),
    )
    amine_c = MonomerSpec(
        id="amine_c",
        name="amine at the c2 side",
        motifs=(
            ReactiveMotif(
                id="n1",
                kind="amine",
                atom_ids=(1,),
                frame=Frame(origin=(0.0, 0.0, 0.0), primary=(-1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)),
            ),
        ),
    )
    template = ReactionLibrary.builtin().get("imine_bridge")
    assignment_plan = AssignmentPlan(
        net_plan=NetPlan(
            topology=None,
            monomer_ids=(linker.id, amine_b.id, amine_c.id),
            reaction_ids=(template.id,),
        ),
        slot_to_monomer={"slot1": linker.id, "slot2": amine_b.id, "slot3": amine_c.id},
    )
    outcome = AssignmentOutcome(
        assignment_plan=assignment_plan,
        monomer_instances=(
            MonomerInstance(id="mL", monomer_id=linker.id),
            MonomerInstance(id="mB", monomer_id=amine_b.id),
            MonomerInstance(id="mC", monomer_id=amine_c.id),
        ),
        events=(
            ReactionEvent(
                id="rxn_b",
                template_id=template.id,
                participants=(
                    MotifRef(monomer_instance_id="mL", monomer_id=linker.id, motif_id="c1"),
                    MotifRef(monomer_instance_id="mB", monomer_id=amine_b.id, motif_id="n1"),
                ),
            ),
            ReactionEvent(
                id="rxn_c",
                template_id=template.id,
                participants=(
                    MotifRef(monomer_instance_id="mL", monomer_id=linker.id, motif_id="c2"),
                    MotifRef(monomer_instance_id="mC", monomer_id=amine_c.id, motif_id="n1"),
                ),
            ),
        ),
        unreacted_motifs=(),
        consumed_count=4,
    )
    specs = {linker.id: linker, amine_b.id: amine_b, amine_c.id: amine_c}
    templates = {template.id: template}
    # Both amines sit at exactly the 1.3 target distance along asymmetric
    # off-axis directions, so the initial residual is purely orientational and
    # worst-event-only steering of the linker over-rotates it away from the
    # lower-residual event. The twisted linker pose (+5.7 deg) makes the two
    # events disagree about the correction.
    direction_b = (cos(math.radians(171.0)), sin(math.radians(171.0)), 0.0)
    direction_c = (cos(math.radians(5.0)), sin(math.radians(5.0)), 0.0)
    base_state = AssemblyState(
        cell=((12.0, 0.0, 0.0), (0.0, 12.0, 0.0), (0.0, 0.0, 8.0)),
        monomer_poses={
            "mL": Pose(rotation_matrix=rotation_about_axis((0.0, 0.0, 1.0), 0.1)),
            "mB": Pose(translation=tuple(1.3 * component for component in direction_b)),
            "mC": Pose(translation=tuple(1.3 * component for component in direction_c)),
        },
        stacking_state="disabled",
    )
    return specs, templates, outcome, base_state


class EmbeddingTests(unittest.TestCase):
    def test_embedder_uses_repository_selected_cell_for_workspace_topology_case(self):
        specs, templates, outcome = build_imine_case()
        embedder = PeriodicEmbedder(EmbeddingConfig(default_lateral_span=24.0, default_layer_spacing=7.5))

        embedding = embedder.embed(outcome, specs, templates)

        self.assertEqual(embedding.metadata["mode"], "topology-guided")
        self.assertEqual(embedding.metadata["topology"], "car")
        # The built cell has equal-length, 90° in-plane vectors, so the
        # shared classifier honestly reports "square" (previously the
        # topology-id-derived answer claimed "orthogonal").
        self.assertEqual(embedding.metadata["cell_kind"], "square")
        self.assertEqual(embedding.state.stacking_state, "disabled")
        self.assertEqual(len(embedding.state.monomer_poses), 2)
        self.assertAlmostEqual(embedding.state.cell[1][0], 0.0)
        self.assertGreater(embedding.state.cell[1][1], 0.0)
        self.assertNotEqual(
            embedding.state.monomer_poses["m1"].translation,
            embedding.state.monomer_poses["m2"].translation,
        )

    def test_2d_embedding_labels_c_axis_as_vacuum_slab(self):
        """W1.4: single-layer embedding exports pad c with vacuum and say so."""
        specs, templates, outcome = build_hcb_case()
        embedding = PeriodicEmbedder().embed(outcome, specs, templates)

        self.assertEqual(embedding.metadata["c_axis_semantics"], "vacuum_slab")
        self.assertAlmostEqual(
            embedding.state.cell[2][2],
            EmbeddingConfig().default_layer_spacing,
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
            motifs=tetrahedral_motifs("n", "amine", radius=4.5),
        )
        tetra_aldehyde = MonomerSpec(
            id="tetra_aldehyde",
            name="tetraaldehyde",
            motifs=tetrahedral_motifs("c", "aldehyde", radius=2.4),
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

    def test_mixed_linkage_embedding_applies_per_template_origin_retraction(self):
        from cofkit.planner import TopologyHint

        def atomistic_radial_motifs(prefix: str, kinds: tuple[str, ...], radius: float):
            motifs = []
            atom_symbols = []
            atom_positions = []
            for idx, (kind, angle) in enumerate(zip(kinds, (0.0, 2.0 * pi / 3.0, -2.0 * pi / 3.0))):
                direction = (cos(angle), sin(angle), 0.0)
                reactive_atom_id = 2 * idx
                anchor_atom_id = 2 * idx + 1
                atom_symbols.extend(("C", "C"))
                atom_positions.append((radius * direction[0], radius * direction[1], 0.0))
                atom_positions.append(((radius - 1.0) * direction[0], (radius - 1.0) * direction[1], 0.0))
                motifs.append(
                    ReactiveMotif(
                        id=f"{prefix}{idx + 1}",
                        kind=kind,
                        atom_ids=(reactive_atom_id, anchor_atom_id),
                        frame=Frame(
                            origin=(radius * direction[0], radius * direction[1], 0.0),
                            primary=direction,
                            normal=(0.0, 0.0, 1.0),
                        ),
                        metadata={"reactive_atom_id": reactive_atom_id, "anchor_atom_id": anchor_atom_id},
                    )
                )
            return tuple(motifs), tuple(atom_symbols), tuple(atom_positions)

        def build_outcome(node_kinds, event_template_ids):
            node_motifs, node_symbols, node_positions = atomistic_radial_motifs("n", node_kinds, radius=4.5)
            linker_motifs, linker_symbols, linker_positions = atomistic_radial_motifs(
                "c", ("aldehyde", "aldehyde", "aldehyde"), radius=2.4
            )
            node = MonomerSpec(
                id="node",
                name="mixed node",
                motifs=node_motifs,
                atom_symbols=node_symbols,
                atom_positions=node_positions,
            )
            linker = MonomerSpec(
                id="linker",
                name="trialdehyde linker",
                motifs=linker_motifs,
                atom_symbols=linker_symbols,
                atom_positions=linker_positions,
            )
            topology = TopologyHint(
                id="hcb",
                dimensionality="2D",
                node_coordination=(3,),
                metadata={"n_node_definitions": 1},
            )
            assignment_plan = AssignmentPlan(
                net_plan=NetPlan(topology=topology, monomer_ids=("node", "linker"), reaction_ids=tuple(event_template_ids)),
                slot_to_monomer={"slot1": "node", "slot2": "linker"},
            )
            images = ((0, 0, 0), (-1, 0, 0), (0, -1, 0))
            events = tuple(
                ReactionEvent(
                    id=f"rxn{idx + 1}",
                    template_id=template_id,
                    participants=(
                        MotifRef(monomer_instance_id="m1", monomer_id="node", motif_id=f"n{idx + 1}"),
                        MotifRef(monomer_instance_id="m2", monomer_id="linker", motif_id=f"c{idx + 1}", periodic_image=images[idx]),
                    ),
                )
                for idx, template_id in enumerate(event_template_ids)
            )
            outcome = AssignmentOutcome(
                assignment_plan=assignment_plan,
                monomer_instances=(
                    MonomerInstance(id="m1", monomer_id="node"),
                    MonomerInstance(id="m2", monomer_id="linker"),
                ),
                events=events,
                unreacted_motifs=(),
                consumed_count=6,
            )
            templates = {template_id: ReactionLibrary.builtin().get(template_id) for template_id in set(event_template_ids)}
            return {"node": node, "linker": linker}, templates, outcome

        mixed_specs, mixed_templates, mixed_outcome = build_outcome(
            ("amine", "amine", "hydrazine"),
            ("imine_bridge", "imine_bridge", "azine_bridge"),
        )
        mixed_embedding = PeriodicEmbedder().embed(mixed_outcome, mixed_specs, mixed_templates)

        self.assertEqual(mixed_embedding.metadata["placement_mode"], "single-node-bipartite")
        node_offsets = mixed_embedding.metadata["poses"]["m1"]["radial_offsets"]
        linker_offsets = mixed_embedding.metadata["poses"]["m2"]["radial_offsets"]
        # Each bridge gets its own template's retraction: 0.11 for the two
        # imine motifs, 0.08 for the azine motif — before the per-template fix
        # the mixed build retracted nothing (all offsets stayed at 4.5/2.4).
        self.assertAlmostEqual(node_offsets[0], 4.5 - 0.11, places=6)
        self.assertAlmostEqual(node_offsets[1], 4.5 - 0.11, places=6)
        self.assertAlmostEqual(node_offsets[2], 4.5 - 0.08, places=6)
        self.assertAlmostEqual(linker_offsets[0], 2.4 - 0.11, places=6)
        self.assertAlmostEqual(linker_offsets[1], 2.4 - 0.11, places=6)
        self.assertAlmostEqual(linker_offsets[2], 2.4 - 0.08, places=6)

        # A pure imine build over the same geometry is bit-identical to the
        # shared-template behavior: its offsets equal a direct shared-template
        # rotation computation, and the mixed build's imine motifs match it.
        pure_specs, pure_templates, pure_outcome = build_outcome(
            ("amine", "amine", "amine"),
            ("imine_bridge", "imine_bridge", "imine_bridge"),
        )
        pure_embedding = PeriodicEmbedder().embed(pure_outcome, pure_specs, pure_templates)
        reference_rotation, reference_offsets = PeriodicEmbedder()._rotation_for_planar_motifs(
            pure_specs["node"],
            tuple((cos(pi / 6.0 + offset), sin(pi / 6.0 + offset), 0.0) for offset in (0.0, 2.0 * pi / 3.0, -2.0 * pi / 3.0)),
            template_id="imine_bridge",
        )
        self.assertEqual(
            pure_embedding.metadata["poses"]["m1"]["radial_offsets"],
            tuple(round(offset, 6) for offset in reference_offsets),
        )
        self.assertAlmostEqual(node_offsets[0], pure_embedding.metadata["poses"]["m1"]["radial_offsets"][0], places=12)


class ScoringTests(unittest.TestCase):
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

    def test_optimizer_translation_step_moves_long_bridge_toward_target(self):
        specs, templates, outcome, base_state = build_single_imine_bridge_case()
        stretched_poses = dict(base_state.monomer_poses)
        stretched_poses["m2"] = Pose(
            translation=(4.0, 0.0, 0.0),
            rotation_matrix=base_state.monomer_poses["m2"].rotation_matrix,
        )
        stretched_state = AssemblyState(
            cell=base_state.cell,
            monomer_poses=stretched_poses,
            torsions=base_state.torsions,
            layer_offsets=base_state.layer_offsets,
            stacking_state=base_state.stacking_state,
        )
        scorer = CandidateScorer()
        optimizer = ContinuousOptimizer(
            config=OptimizerConfig(max_iterations=1),
            scorer=scorer,
        )

        initial_report = scorer.bridge_geometry_report(outcome, stretched_state, specs, templates)
        optimized = optimizer.optimize(outcome, stretched_state, specs, templates)
        final_report = scorer.bridge_geometry_report(outcome, optimized.state, specs, templates)

        self.assertGreater(initial_report.event_metrics[0].actual_distance, 3.9)
        self.assertLess(final_report.event_metrics[0].actual_distance, initial_report.event_metrics[0].actual_distance)
        self.assertLess(final_report.total_residual, initial_report.total_residual)

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

    def test_orientation_proposal_rotates_about_attachment_point(self):
        specs, templates, outcome, base_state = build_offcenter_imine_bridge_case()
        optimizer = ContinuousOptimizer(scorer=CandidateScorer())

        proposal = optimizer._refine_orientations(base_state, outcome, specs, templates)

        event = outcome.events[0]
        rotations_changed = False
        for participant in event.participants:
            instance_id = participant.monomer_instance_id
            motif = specs[participant.monomer_id].motif_by_id(participant.motif_id)
            local_origin = effective_motif_origin(event.template_id, specs[participant.monomer_id], motif)
            before = optimizer._world_motif_origin(
                base_state.cell,
                base_state.monomer_poses[instance_id],
                local_origin,
                participant.periodic_image,
            )
            after = optimizer._world_motif_origin(
                base_state.cell,
                proposal.monomer_poses[instance_id],
                local_origin,
                participant.periodic_image,
            )
            for before_component, after_component in zip(before, after):
                self.assertAlmostEqual(before_component, after_component, places=9)
            if base_state.monomer_poses[instance_id].rotation_matrix != proposal.monomer_poses[instance_id].rotation_matrix:
                rotations_changed = True
        # The fixture is built so both proposals are genuine rotations; with
        # off-center motif origins the pre-fix code (rotation applied at a
        # fixed translation) would have displaced these attachment points.
        self.assertTrue(rotations_changed)

    def test_orientation_proposal_weights_all_incident_events(self):
        specs, templates, outcome, base_state = build_two_event_linker_case()
        scorer = CandidateScorer()
        optimizer = ContinuousOptimizer(scorer=scorer)

        initial_report = scorer.bridge_geometry_report(outcome, base_state, specs, templates)
        proposal = optimizer._refine_orientations(base_state, outcome, specs, templates)
        proposal_report = scorer.bridge_geometry_report(outcome, proposal, specs, templates)

        self.assertLess(proposal_report.total_residual, initial_report.total_residual)

        # Reconstruct the pre-fix behavior for the linker — steer by the single
        # worst-residual event with the translation held fixed — and confirm
        # that the blended, attachment-preserving proposal scores better.
        residuals = {metrics.event_id: metrics.total_residual for metrics in initial_report.event_metrics}
        worst_event = max(outcome.events, key=lambda event: residuals[event.id])
        linker_participant = next(p for p in worst_event.participants if p.monomer_instance_id == "mL")
        other_participant = next(p for p in worst_event.participants if p.monomer_instance_id != "mL")
        linker_pose = base_state.monomer_poses["mL"]
        other_pose = base_state.monomer_poses[other_participant.monomer_instance_id]
        motif = specs[linker_participant.monomer_id].motif_by_id(linker_participant.motif_id)
        origin = optimizer._world_motif_origin(
            base_state.cell,
            linker_pose,
            effective_motif_origin(worst_event.template_id, specs[linker_participant.monomer_id], motif),
            linker_participant.periodic_image,
        )
        other_origin = optimizer._world_motif_origin(
            base_state.cell,
            other_pose,
            effective_motif_origin(
                worst_event.template_id,
                specs[other_participant.monomer_id],
                specs[other_participant.monomer_id].motif_by_id(other_participant.motif_id),
            ),
            other_participant.periodic_image,
        )
        target_primary = normalize(tuple(o2 - o1 for o1, o2 in zip(origin, other_origin)))
        worst_only_rotation = rotation_from_frame_to_axes(motif.frame, target_primary, (0.0, 0.0, 1.0))
        worst_only_poses = dict(proposal.monomer_poses)
        worst_only_poses["mL"] = Pose(
            translation=linker_pose.translation,
            rotation_matrix=worst_only_rotation,
        )
        worst_only_state = AssemblyState(
            cell=base_state.cell,
            monomer_poses=worst_only_poses,
            torsions=base_state.torsions,
            layer_offsets=base_state.layer_offsets,
            stacking_state=base_state.stacking_state,
        )
        worst_only_report = scorer.bridge_geometry_report(outcome, worst_only_state, specs, templates)

        self.assertGreater(worst_only_report.total_residual, proposal_report.total_residual)


class DegenerateVectorPolicyTests(unittest.TestCase):
    """Audit A10 follow-up: degenerate vectors raise instead of fabricating +x.

    Probe evidence (2026-09-26): the old silent (1,0,0) fallback never fired
    on any live path in optimizer/embedding/batch across real builds and the
    full test suite, so degenerate geometry is now a visible ValueError.
    """

    def test_optimizer_safe_normalize_raises_on_degenerate_vector(self):
        optimizer = ContinuousOptimizer(scorer=CandidateScorer())
        with self.assertRaises(ValueError) as ctx:
            optimizer._safe_normalize((0.0, 0.0, 0.0))
        self.assertIn("ContinuousOptimizer", str(ctx.exception))
        self.assertIn("degenerate", str(ctx.exception))

    def test_embedder_safe_normalize_raises_on_degenerate_vector(self):
        embedder = PeriodicEmbedder()
        with self.assertRaises(ValueError) as ctx:
            embedder._safe_normalize((0.0, 0.0, 0.0))
        self.assertIn("PeriodicEmbedder", str(ctx.exception))
        self.assertIn("degenerate", str(ctx.exception))

    def test_refine_translations_uses_normal1_for_antiparallel_normals(self):
        # A flat bridge has antiparallel motif normals whose sum is
        # (near-)zero; the planarity correction must steer along normal1
        # rather than any fabricated direction, and must not raise.
        specs, templates, outcome, base_state = build_single_imine_bridge_case()
        flip_about_x = (
            (1.0, 0.0, 0.0),
            (0.0, -1.0, 0.0),
            (0.0, 0.0, -1.0),
        )
        poses = dict(base_state.monomer_poses)
        poses["m2"] = Pose(translation=(0.0, 0.0, 1.3), rotation_matrix=flip_about_x)
        state = AssemblyState(
            cell=base_state.cell,
            monomer_poses=poses,
            torsions=base_state.torsions,
            layer_offsets=base_state.layer_offsets,
            stacking_state=base_state.stacking_state,
        )
        optimizer = ContinuousOptimizer(scorer=CandidateScorer())

        refined = optimizer._refine_translations(state, outcome, specs, templates)

        # The fixture monomers carry no atom positions, so both motif origins
        # are their frame origins: delta = (0, 0, 1.3), normal1 = +z and
        # normal2 = -z sum to zero, and the guard picks plane_normal =
        # normal1 = +z. A fabricated in-plane plane_normal would zero the
        # planarity term (dot(delta, plane_normal) == 0), so the expected z
        # update discriminates the guard from any fabricated direction.
        step = optimizer.config.translation_step
        target = optimizer.scorer._target_distance(templates["imine_bridge"])
        expected_z = (0.5 * (1.3 - target) + 0.35 * 1.3) * step
        self.assertAlmostEqual(refined.monomer_poses["m1"].translation[2], expected_z)
        self.assertAlmostEqual(refined.monomer_poses["m2"].translation[2], 1.3 - expected_z)


class EngineIntegrationTests(unittest.TestCase):
    def test_engine_exposes_embedding_and_score_metadata(self):
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
        self.assertIsNone(best.score)
        self.assertIn("bridge_geometry_residual", best.metadata["score_metadata"])
        self.assertIn("normal_misalignment_residual", best.metadata["score_metadata"]["bridge_event_metrics"][0])


if __name__ == "__main__":
    unittest.main()
