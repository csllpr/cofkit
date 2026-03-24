import sys
import unittest
from math import cos, pi, sin
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit import (
    AssignmentSolver,
    BatchGenerationConfig,
    BatchStructureGenerator,
    COFEngine,
    COFProject,
    Frame,
    MonomerInstance,
    MonomerSpec,
    MotifRef,
    NetPlanner,
    PeriodicProductGraph,
    ReactiveMotif,
    ReactionEvent,
    ReactionLibrary,
)


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


class ReactionLibraryTests(unittest.TestCase):
    def test_builtin_imine_template_matches_order_independently(self):
        lib = ReactionLibrary.builtin()
        imine = lib.get("imine_bridge")
        self.assertTrue(imine.matches(("amine", "aldehyde")))
        self.assertTrue(imine.matches(("aldehyde", "amine")))
        self.assertFalse(imine.matches(("amine", "amine")))

    def test_binary_bridge_profiles_expose_pair_order_and_target_distance(self):
        lib = ReactionLibrary.builtin()

        order = lib.resolve_binary_bridge_pair_order("hydrazone_bridge", "aldehyde", "hydrazide")

        self.assertEqual(order.ordered_indices, (1, 0))
        self.assertEqual(tuple(role.role_id for role in order.roles), ("hydrazide", "aldehyde"))
        self.assertAlmostEqual(lib.bridge_target_distance("hydrazone_bridge"), 1.3, places=6)
        self.assertTrue(lib.supports_binary_bridge_pair_generation("hydrazone_bridge"))
        self.assertAlmostEqual(lib.bridge_target_distance("vinylene_bridge"), 1.34, places=6)


class ProductGraphTests(unittest.TestCase):
    def setUp(self):
        self.lib = ReactionLibrary.builtin()
        self.templates = {"imine_bridge": self.lib.get("imine_bridge")}
        self.m1 = MonomerSpec(
            id="amine",
            name="amine monomer",
            motifs=(ReactiveMotif(id="n1", kind="amine", atom_ids=(1,), frame=Frame.xy()),),
        )
        self.m2 = MonomerSpec(
            id="aldehyde",
            name="aldehyde monomer",
            motifs=(ReactiveMotif(id="c1", kind="aldehyde", atom_ids=(2,), frame=Frame.xy()),),
        )
        self.specs = {self.m1.id: self.m1, self.m2.id: self.m2}

    def test_rejects_motif_reuse(self):
        graph = PeriodicProductGraph()
        graph.add_monomer(MonomerInstance(id="m1", monomer_id="amine"))
        graph.add_monomer(MonomerInstance(id="m2", monomer_id="aldehyde"))

        first = ReactionEvent(
            id="rxn1",
            template_id="imine_bridge",
            participants=(
                MotifRef("m1", "amine", "n1"),
                MotifRef("m2", "aldehyde", "c1"),
            ),
        )
        graph.add_reaction_event(first, self.specs, self.templates)

        second = ReactionEvent(
            id="rxn2",
            template_id="imine_bridge",
            participants=(
                MotifRef("m1", "amine", "n1"),
                MotifRef("m2", "aldehyde", "c1"),
            ),
        )
        with self.assertRaises(ValueError):
            graph.add_reaction_event(second, self.specs, self.templates)


class PlannerAndSolverTests(unittest.TestCase):
    def test_net_planner_prefers_workspace_repository_match_for_two_tritopic_monomers(self):
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
        planner = NetPlanner()
        templates = (ReactionLibrary.builtin().get("imine_bridge"),)
        plans = planner.propose((tri_amine, tri_aldehyde), templates, "2D")
        self.assertEqual(plans[0].topology.id, "car")
        self.assertEqual(plans[0].topology.metadata["n_edge_definitions"], 4)

    def test_assignment_solver_consumes_all_bridgeable_motifs(self):
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
        planner = NetPlanner()
        solver = AssignmentSolver()
        templates = (ReactionLibrary.builtin().get("imine_bridge"),)
        plan = planner.propose((tri_amine, tri_aldehyde), templates, "2D")[0]
        assignment = solver.build_assignment_plans(plan, specs)[0]
        outcome = solver.solve_events(
            assignment_plan=assignment,
            monomer_instances=solver.instantiate_monomers(assignment),
            monomer_specs=specs,
            templates=templates,
        )
        self.assertEqual(len(outcome.events), 3)
        self.assertEqual(len(outcome.unreacted_motifs), 0)

    def test_explicit_hcb_request_accepts_single_node_bipartite_three_plus_three_case(self):
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
        planner = NetPlanner()
        templates = (ReactionLibrary.builtin().get("imine_bridge"),)

        plans = planner.propose((tri_amine, tri_aldehyde), templates, "2D", target_topologies=("hcb",))

        self.assertEqual(len(plans), 1)
        self.assertEqual(plans[0].topology.id, "hcb")
        self.assertEqual(plans[0].topology.metadata["n_node_definitions"], 1)


class EngineTests(unittest.TestCase):
    def test_imine_project_gets_topology_guided_candidate(self):
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
        project = COFProject(
            monomers=(tri_amine, tri_aldehyde),
            allowed_reactions=("imine_bridge",),
            target_dimensionality="2D",
        )
        engine = COFEngine()
        best = engine.run(project).top(1)[0]

        self.assertEqual(len(best.events), 3)
        self.assertGreater(best.score, 30.0)
        self.assertEqual(best.metadata["graph_summary"]["n_reaction_events"], 3)
        self.assertEqual(best.metadata["net_plan"]["topology"], "car")
        self.assertNotIn("no_topology_hint", best.flags)

    def test_ring_forming_project_bypasses_topology(self):
        tri_boronic = MonomerSpec(
            id="boronic",
            name="tri-boronic monomer",
            motifs=(
                ReactiveMotif(id="b1", kind="boronic_acid", atom_ids=(1,), frame=Frame.xy()),
                ReactiveMotif(id="b2", kind="boronic_acid", atom_ids=(2,), frame=Frame.yz()),
                ReactiveMotif(id="b3", kind="boronic_acid", atom_ids=(3,), frame=Frame.zx()),
            ),
        )
        project = COFProject(
            monomers=(tri_boronic,),
            allowed_reactions=("boroxine_trimerization",),
            target_dimensionality="2D",
        )
        engine = COFEngine()
        best = engine.run(project).top(1)[0]

        self.assertEqual(len(best.events), 1)
        self.assertIn("contains_ring_forming_event", best.flags)
        self.assertIn("no_topology_hint", best.flags)

    def test_engine_supports_explicit_indexed_mixed_connectivity_topology(self):
        tetra_amine = MonomerSpec(
            id="tetra_amine",
            name="tetramine",
            motifs=cyclic_motifs("n", "amine", count=4, radius=4.5),
        )
        hexa_aldehyde = MonomerSpec(
            id="hexa_aldehyde",
            name="hexaaldehyde",
            motifs=cyclic_motifs("c", "aldehyde", count=6, radius=2.4),
        )
        project = COFProject(
            monomers=(tetra_amine, hexa_aldehyde),
            allowed_reactions=("imine_bridge",),
            target_dimensionality="3D",
            target_topologies=("cor",),
        )

        best = COFEngine().run(project).top(1)[0]

        self.assertEqual(best.metadata["net_plan"]["topology"], "cor")
        self.assertEqual(best.metadata["embedding"]["placement_mode"], "indexed-topology-node-node-3d")
        self.assertEqual(best.metadata["graph_summary"]["n_monomer_instances"], 30)
        self.assertEqual(best.metadata["graph_summary"]["n_reaction_events"], 72)


class GenericBinaryBridgeBatchTests(unittest.TestCase):
    def test_batch_generator_handles_non_imine_binary_bridge_template(self):
        tri_hydrazide = MonomerSpec(
            id="tri_hydrazide",
            name="tri hydrazide",
            motifs=trigonal_motifs("h", "hydrazide", radius=4.5),
        )
        tri_aldehyde = MonomerSpec(
            id="tri_aldehyde",
            name="tri aldehyde",
            motifs=trigonal_motifs("c", "aldehyde", radius=2.4),
        )
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                allowed_reactions=("hydrazone_bridge",),
                single_node_topology_ids=("hcb",),
                write_cif=False,
            )
        )

        summary, candidate = generator.generate_monomer_pair_candidate(tri_aldehyde, tri_hydrazide)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "3+3-node-node")
        self.assertEqual(summary.reactant_record_ids, {"hydrazide": "tri_hydrazide", "aldehyde": "tri_aldehyde"})
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "hcb")
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 3)
        self.assertEqual(candidate.metadata["topology_selection"]["selected_topology"], "hcb")

    def test_batch_generator_handles_keto_enamine_binary_bridge_template(self):
        tri_amine = MonomerSpec(
            id="tri_amine",
            name="tri amine",
            motifs=trigonal_motifs("n", "amine", radius=4.5),
        )
        tri_keto_aldehyde = MonomerSpec(
            id="tri_keto_aldehyde",
            name="tri keto aldehyde",
            motifs=trigonal_motifs("c", "keto_aldehyde", radius=2.6),
        )
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                allowed_reactions=("keto_enamine_bridge",),
                single_node_topology_ids=("hcb",),
                write_cif=False,
            )
        )

        summary, candidate = generator.generate_monomer_pair_candidate(tri_amine, tri_keto_aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "3+3-node-node")
        self.assertEqual(summary.reactant_record_ids, {"amine": "tri_amine", "keto_aldehyde": "tri_keto_aldehyde"})
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "hcb")
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 3)

    def test_batch_generator_handles_boronate_ester_binary_bridge_template(self):
        tri_boronic_acid = MonomerSpec(
            id="tri_boronic_acid",
            name="tri boronic acid",
            motifs=trigonal_motifs("b", "boronic_acid", radius=4.5),
        )
        di_catechol = MonomerSpec(
            id="di_catechol",
            name="di catechol",
            motifs=cyclic_motifs("o", "catechol", count=2, radius=2.4),
        )
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                allowed_reactions=("boronate_ester_bridge",),
                single_node_topology_ids=("hcb",),
                write_cif=False,
            )
        )

        summary, candidate = generator.generate_monomer_pair_candidate(tri_boronic_acid, di_catechol)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "3+2-node-linker")
        self.assertEqual(summary.reactant_record_ids, {"boronic_acid": "tri_boronic_acid", "catechol": "di_catechol"})
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "hcb")
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 6)

    def test_batch_generator_handles_vinylene_binary_bridge_template(self):
        tri_activated_methylene = MonomerSpec(
            id="tri_activated_methylene",
            name="tri activated methylene",
            motifs=trigonal_motifs("v", "activated_methylene", radius=4.5),
        )
        di_aldehyde = MonomerSpec(
            id="di_aldehyde",
            name="di aldehyde",
            motifs=cyclic_motifs("c", "aldehyde", count=2, radius=2.4),
        )
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                allowed_reactions=("vinylene_bridge",),
                single_node_topology_ids=("hcb",),
                write_cif=False,
            )
        )

        summary, candidate = generator.generate_monomer_pair_candidate(tri_activated_methylene, di_aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "3+2-node-linker")
        self.assertEqual(
            summary.reactant_record_ids,
            {"activated_methylene": "tri_activated_methylene", "aldehyde": "di_aldehyde"},
        )
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "hcb")
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 6)

    def test_batch_generator_supports_explicit_indexed_node_linker_topology(self):
        tetra_amine = MonomerSpec(
            id="tetra_amine",
            name="tetramine",
            motifs=cyclic_motifs("n", "amine", count=4, radius=4.5),
        )
        di_aldehyde = MonomerSpec(
            id="di_aldehyde",
            name="dialdehyde",
            motifs=cyclic_motifs("c", "aldehyde", count=2, radius=2.4),
        )
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                topology_ids=("bew",),
                write_cif=False,
            )
        )

        summary, candidate = generator.generate_monomer_pair_candidate(tetra_amine, di_aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "4+2-node-linker")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "bew")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "node-linker-indexed-topology")
        self.assertEqual(candidate.metadata["graph_summary"]["n_monomer_instances"], 9)
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 12)

    def test_batch_generator_supports_explicit_indexed_mixed_connectivity_node_node_topology(self):
        tetra_amine = MonomerSpec(
            id="tetra_amine",
            name="tetramine",
            motifs=cyclic_motifs("n", "amine", count=4, radius=4.5),
        )
        hexa_aldehyde = MonomerSpec(
            id="hexa_aldehyde",
            name="hexaaldehyde",
            motifs=cyclic_motifs("c", "aldehyde", count=6, radius=2.4),
        )
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                topology_ids=("cor",),
                write_cif=False,
            )
        )

        summary, candidate = generator.generate_monomer_pair_candidate(tetra_amine, hexa_aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "4+6-node-node")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "cor")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "indexed-topology-node-node-3d")
        self.assertEqual(candidate.metadata["graph_summary"]["n_monomer_instances"], 30)
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 72)


if __name__ == "__main__":
    unittest.main()
