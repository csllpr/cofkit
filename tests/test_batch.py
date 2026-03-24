import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit import BatchGenerationConfig, BatchMonomerRecord, BatchStructureGenerator, Frame, MonomerSpec, ReactiveMotif

try:
    from rdkit import Chem  # noqa: F401
except ImportError:  # pragma: no cover - environment-dependent
    Chem = None


TAPB = "C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N"
TFB = "C1=C(C=C(C=C1C=O)C=O)C=O"
TEREPHTHALALDEHYDE = "O=Cc1ccc(C=O)cc1"
PPD = "Nc1ccc(N)cc1"
ASYMMETRIC_TRIAMINE = (
    "COC(=O)[C@@H]1C[C@H](OCc2c(-c3ccc(N)cc3)cc(-c3ccc(N)cc3)cc2-c2ccc(N)cc2)CN1C(=O)OC(C)(C)C"
)
LONG_DIALDEHYDE = "O=Cc1ccc(C#Cc2c3ccccc3c(C#Cc3ccc(C=O)cc3)c3ccccc23)cc1"
ASYMMETRIC_TRIALDEHYDE = "O=Cc1ccccc1C#Cc1cc(C#Cc2ccccc2C=O)cc(C#Cc2ccccc2C=O)c1"
TETRA_AMINE = "Nc1ccc(C(c2ccc(N)cc2)(c2ccc(N)cc2)c2ccc(N)cc2)cc1"
TETRA_ALDEHYDE = "O=Cc1ccc(C(c2ccc(C=O)cc2)(c2ccc(C=O)cc2)c2ccc(C=O)cc2)cc1"
HEXA_AMINE = "Nc1cc2c3cc(N)c(N)cc3c3cc(N)c(N)cc3c2cc1N"
HEXA_ALDEHYDE = "O=Cc1ccc(-c2cc3c(cc2-c2ccc(C=O)cc2)C2c4cc(-c5ccc(C=O)cc5)c(-c5ccc(C=O)cc5)cc4C3c3cc(-c4ccc(C=O)cc4)c(-c4ccc(C=O)cc4)cc32)cc1"
TP = "O=Cc1c(O)c(C=O)c(O)c(C=O)c1O"
COF42_HYDRAZIDE = "CCOc1cc(C(=O)NN)cc(C(=O)NN)c1OCC"


@unittest.skipIf(Chem is None, "RDKit is not available")
class BatchStructureGeneratorTests(unittest.TestCase):
    def setUp(self):
        self.generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
                single_node_topology_ids=("hcb",),
            )
        )

    def test_load_smiles_library_parses_fixture_format(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "amines_count_2.txt"
            path.write_text("smiles\nNc1ccc(N)cc1\nNc1ccnc(N)c1\n", encoding="utf-8")

            records = self.generator.load_smiles_library(
                path,
                motif_kind="amine",
                expected_connectivity=2,
            )

        self.assertEqual([record.id for record in records], ["amines_count_2_0001", "amines_count_2_0002"])
        self.assertEqual(records[0].smiles, "Nc1ccc(N)cc1")
        self.assertEqual(records[0].expected_connectivity, 2)

    def test_infer_monomer_record_detects_amine(self):
        record = self.generator.infer_monomer_record(
            PPD,
            record_id="generic_0001",
        )

        self.assertEqual(record.motif_kind, "amine")
        self.assertEqual(record.expected_connectivity, 2)
        self.assertTrue(record.metadata["auto_detected"])
        self.assertEqual(record.metadata["detected_motif_kinds"], ("amine",))

    def test_infer_monomer_record_prefers_keto_aldehyde_over_generic_aldehyde(self):
        record = self.generator.infer_monomer_record(
            TP,
            record_id="generic_0001",
        )

        self.assertEqual(record.motif_kind, "keto_aldehyde")
        self.assertEqual(record.expected_connectivity, 3)
        self.assertEqual(record.metadata["detected_motif_kinds"], ("aldehyde", "keto_aldehyde"))
        self.assertEqual(record.metadata["detected_connectivities"]["keto_aldehyde"], 3)

    def test_load_binary_bridge_test_set_auto_groups_generic_files_by_detected_role_and_connectivity(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                allowed_reactions=("keto_enamine_bridge",),
                rdkit_num_conformers=2,
                retain_top_results=5,
                single_node_topology_ids=("hcb",),
            )
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            (temp_path / "set_a.txt").write_text(f"smiles\n{PPD}\n", encoding="utf-8")
            (temp_path / "set_b.txt").write_text(f"smiles\n{TP}\n", encoding="utf-8")

            libraries = generator.load_binary_bridge_test_set_auto(
                temp_path,
                template_id="keto_enamine_bridge",
            )

        self.assertEqual(tuple(sorted(libraries)), ("amines_count_2", "keto_aldehydes_count_3"))
        self.assertEqual(libraries["amines_count_2"][0].motif_kind, "amine")
        self.assertEqual(libraries["keto_aldehydes_count_3"][0].motif_kind, "keto_aldehyde")

    def test_available_binary_bridge_template_ids_discovers_supported_templates(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            (temp_path / "amines_count_3.txt").write_text(f"smiles\n{TAPB}\n", encoding="utf-8")
            (temp_path / "aldehydes_count_3.txt").write_text(f"smiles\n{TFB}\n", encoding="utf-8")
            (temp_path / "hydrazides_count_2.txt").write_text(f"smiles\n{COF42_HYDRAZIDE}\n", encoding="utf-8")
            (temp_path / "keto_aldehydes_count_3.txt").write_text(f"smiles\n{TP}\n", encoding="utf-8")

            template_ids = self.generator.available_binary_bridge_template_ids(temp_path)

        self.assertEqual(template_ids, ("hydrazone_bridge", "imine_bridge", "keto_enamine_bridge"))

    def test_custom_smiles_monomer_builder_can_be_injected(self):
        calls = []

        def fake_builder(monomer_id, name, smiles, motif_kind, *, num_conformers, random_seed):
            del name, smiles
            calls.append((monomer_id, motif_kind, num_conformers, random_seed))
            return MonomerSpec(
                id=monomer_id,
                name=monomer_id,
                motifs=(
                    ReactiveMotif(id=f"{motif_kind}1", kind=motif_kind, atom_ids=(1,), frame=Frame.xy()),
                    ReactiveMotif(id=f"{motif_kind}2", kind=motif_kind, atom_ids=(2,), frame=Frame.yz()),
                ),
            )

        generator = BatchStructureGenerator(
            BatchGenerationConfig(rdkit_num_conformers=3, rdkit_random_seed=123),
            smiles_monomer_builder=fake_builder,
        )
        record = BatchMonomerRecord(
            id="fake_amine",
            name="fake_amine",
            smiles="N",
            motif_kind="amine",
            expected_connectivity=2,
        )

        built = generator.build_monomer(record)

        self.assertTrue(built.ok)
        self.assertEqual(calls, [("fake_amine", "amine", 3, 123)])

    def test_three_plus_three_pair_builds_hcb_candidate(self):
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tfb",
            name="tfb",
            smiles=TFB,
            motif_kind="aldehyde",
            expected_connectivity=3,
        )

        summary, candidate = self.generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "3+3-node-node")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "hcb")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "single-node-bipartite")
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 3)

    def test_generate_monomer_pair_candidate_accepts_direct_monomer_specs(self):
        amine_record = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde_record = BatchMonomerRecord(
            id="tfb",
            name="tfb",
            smiles=TFB,
            motif_kind="aldehyde",
            expected_connectivity=3,
        )

        built_amine = self.generator.build_monomer(amine_record)
        built_aldehyde = self.generator.build_monomer(aldehyde_record)

        assert built_amine.monomer is not None
        assert built_aldehyde.monomer is not None
        summary, candidate = self.generator.generate_monomer_pair_candidate(
            built_aldehyde.monomer,
            built_amine.monomer,
        )

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_id, "tapb__tfb")
        self.assertEqual(summary.amine_record_id, "tapb")
        self.assertEqual(summary.aldehyde_record_id, "tfb")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "hcb")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "single-node-bipartite")

    def test_three_plus_two_node_linker_pair_builds_hcb_candidate(self):
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        summary, candidate = self.generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "3+2-node-linker")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "hcb")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "node-linker-single-node")
        self.assertEqual(candidate.metadata["graph_summary"]["n_monomer_instances"], 5)
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 6)
        self.assertNotIn("unreacted_motifs:1", candidate.flags)
        self.assertLess(self._max_bridge_distance_error(candidate), 0.02)
        self.assertLess(candidate.metadata["score_metadata"]["bridge_geometry_residual"], 0.2)

    def test_reverse_polarity_two_plus_three_pair_builds_hcb_candidate(self):
        amine = BatchMonomerRecord(
            id="ppd",
            name="ppd",
            smiles=PPD,
            motif_kind="amine",
            expected_connectivity=2,
        )
        aldehyde = BatchMonomerRecord(
            id="tfb",
            name="tfb",
            smiles=TFB,
            motif_kind="aldehyde",
            expected_connectivity=3,
        )

        summary, candidate = self.generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "3+2-node-linker")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "hcb")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "node-linker-single-node")
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 6)
        self.assertLess(self._max_bridge_distance_error(candidate), 0.02)
        self.assertLess(candidate.metadata["score_metadata"]["bridge_geometry_residual"], 0.2)

    def test_three_plus_two_pair_handles_asymmetric_tritopic_node(self):
        amine = BatchMonomerRecord(
            id="asymmetric_triamine",
            name="asymmetric_triamine",
            smiles=ASYMMETRIC_TRIAMINE,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="long_dialdehyde",
            name="long_dialdehyde",
            smiles=LONG_DIALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        summary, candidate = self.generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "3+2-node-linker")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["embedding"]["cell_kind"], "oblique")
        self.assertLess(self._max_bridge_distance_error(candidate), 0.02)
        self.assertLess(self._max_linker_translation_mismatch(candidate), 0.02)

    def test_three_plus_three_pair_handles_asymmetric_tritopic_nodes(self):
        amine = BatchMonomerRecord(
            id="asymmetric_triamine",
            name="asymmetric_triamine",
            smiles=ASYMMETRIC_TRIAMINE,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="asymmetric_trialdehyde",
            name="asymmetric_trialdehyde",
            smiles=ASYMMETRIC_TRIALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=3,
        )

        summary, candidate = self.generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "3+3-node-node")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["embedding"]["cell_kind"], "oblique")
        self.assertLess(self._max_bridge_distance_error(candidate), 0.02)

    def test_three_plus_two_pair_can_export_cif(self):
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            summary, candidate = self.generator.generate_pair_candidate(
                amine,
                aldehyde,
                out_dir=temp_dir,
                write_cif=True,
            )

            self.assertEqual(summary.status, "ok")
            self.assertIsNotNone(candidate)
            self.assertIsNotNone(summary.cif_path)
            self.assertIn(f"{Path(temp_dir) / 'valid'}", summary.cif_path)
            cif_text = Path(summary.cif_path).read_text(encoding="utf-8")

        self.assertIn("_geom_bond_atom_site_label_1", cif_text)
        self.assertIn("_atom_site_label", cif_text)

    def test_three_plus_two_pair_exports_cif_by_default_when_output_dir_is_provided(self):
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            summary, candidate = self.generator.generate_pair_candidate(
                amine,
                aldehyde,
                out_dir=temp_dir,
            )

            self.assertEqual(summary.status, "ok")
            self.assertIsNotNone(candidate)
            self.assertIsNotNone(summary.cif_path)
            self.assertTrue(Path(summary.cif_path).is_file())
            self.assertIn(f"{Path(temp_dir) / 'valid'}", summary.cif_path)

    def test_hard_hard_invalid_bridge_blocks_cif_export(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
                single_node_topology_ids=("hcb",),
                hard_hard_max_bridge_distance=1.0,
            )
        )
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            summary, candidate = generator.generate_pair_candidate(
                amine,
                aldehyde,
                out_dir=temp_dir,
                write_cif=True,
            )

        self.assertEqual(summary.status, "ok")
        self.assertIsNotNone(candidate)
        self.assertIsNone(summary.cif_path)
        self.assertEqual(
            tuple(summary.metadata["hard_hard_invalid_reasons"]),
            ("bridge_distance_exceeds_cif_export_limit",),
        )
        self.assertTrue(summary.metadata["cif_export_blocked"])

    def test_three_plus_two_pair_can_target_fes_single_node_topology(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
                single_node_topology_ids=("fes",),
            )
        )
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "fes")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "node-linker-single-node-expanded")
        self.assertEqual(candidate.metadata["graph_summary"]["n_monomer_instances"], 10)
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 12)

    def test_three_plus_two_pair_can_target_hca_single_node_topology(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
                single_node_topology_ids=("hca",),
            )
        )
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "hca")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "node-linker-single-node-expanded")
        self.assertEqual(candidate.metadata["graph_summary"]["n_monomer_instances"], 15)
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 18)

    def test_three_plus_two_pair_can_target_fxt_single_node_topology(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
                single_node_topology_ids=("fxt",),
            )
        )
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "fxt")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "node-linker-single-node-expanded")
        self.assertEqual(candidate.metadata["graph_summary"]["n_monomer_instances"], 30)
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 36)

    def test_three_plus_three_pair_can_target_fes_single_node_topology(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
                single_node_topology_ids=("fes",),
            )
        )
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tfb",
            name="tfb",
            smiles=TFB,
            motif_kind="aldehyde",
            expected_connectivity=3,
        )

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "fes")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "single-node-expanded-node-node")
        self.assertEqual(candidate.metadata["graph_summary"]["n_monomer_instances"], 16)
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 24)

    def test_three_plus_three_pair_can_target_fxt_single_node_topology(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
                single_node_topology_ids=("fxt",),
            )
        )
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tfb",
            name="tfb",
            smiles=TFB,
            motif_kind="aldehyde",
            expected_connectivity=3,
        )

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "fxt")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "single-node-expanded-node-node")
        self.assertEqual(candidate.metadata["graph_summary"]["n_monomer_instances"], 12)
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 18)

    def test_three_plus_three_pair_rejects_non_bipartite_single_node_topologies(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
                single_node_topology_ids=("hca",),
            )
        )
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tfb",
            name="tfb",
            smiles=TFB,
            motif_kind="aldehyde",
            expected_connectivity=3,
        )

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "generation-failed")
        self.assertIsNone(candidate)
        self.assertIn("not bipartite", summary.metadata["failed_topologies"]["hca"])

    def test_default_topology_selection_respects_single_node_category_rules(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
            )
        )
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )

        three_plus_two_summary, _ = generator.generate_pair_candidate(
            amine,
            BatchMonomerRecord(
                id="tpal",
                name="tpal",
                smiles=TEREPHTHALALDEHYDE,
                motif_kind="aldehyde",
                expected_connectivity=2,
            ),
        )
        three_plus_three_summary, _ = generator.generate_pair_candidate(
            amine,
            BatchMonomerRecord(
                id="tfb",
                name="tfb",
                smiles=TFB,
                motif_kind="aldehyde",
                expected_connectivity=3,
            ),
        )

        self.assertEqual(
            three_plus_two_summary.metadata["topology_selection"]["available_topologies"],
            ("hcb", "hca", "fes", "fxt", "srs"),
        )
        self.assertEqual(
            three_plus_three_summary.metadata["topology_selection"]["available_topologies"],
            ("hcb", "fes", "fxt", "srs"),
        )

    def test_generate_pair_candidates_returns_all_supported_topologies(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
            )
        )
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        summaries, candidates, attempted_structures = generator.generate_pair_candidates(amine, aldehyde)

        self.assertEqual(attempted_structures, 5)
        self.assertEqual(len(summaries), 5)
        self.assertEqual(len(candidates), 5)
        self.assertEqual(
            {summary.topology_id for summary in summaries},
            {"hcb", "hca", "fes", "fxt", "srs"},
        )
        instance_counts = {
            candidate.metadata["net_plan"]["topology"]: candidate.metadata["graph_summary"]["n_monomer_instances"]
            for candidate in candidates
        }
        self.assertEqual(instance_counts, {"hcb": 5, "hca": 15, "fes": 10, "fxt": 30, "srs": 20})
        self.assertTrue(all(summary.cif_path is None for summary in summaries))

    def test_four_plus_two_pair_uses_expanded_default_topology_pool(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=5,
            )
        )
        amine = BatchMonomerRecord(
            id="tetra_amine",
            name="tetra_amine",
            smiles=TETRA_AMINE,
            motif_kind="amine",
            expected_connectivity=4,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "4+2-node-linker")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(
            set(summary.metadata["topology_selection"]["available_topologies"]),
            {"dia", "sql", "kgm", "pts", "lon", "qtz"},
        )
        self.assertIn(candidate.metadata["net_plan"]["topology"], {"dia", "sql", "kgm", "pts", "lon", "qtz"})
        self.assertGreater(candidate.metadata["graph_summary"]["n_reaction_events"], 0)

    def test_four_plus_two_pair_returns_explicit_2d_topologies_when_requested(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=5,
                single_node_topology_ids=("sql", "kgm", "htb"),
            )
        )
        amine = BatchMonomerRecord(
            id="tetra_amine",
            name="tetra_amine",
            smiles=TETRA_AMINE,
            motif_kind="amine",
            expected_connectivity=4,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        summaries, candidates, attempted_structures = generator.generate_pair_candidates(amine, aldehyde)

        self.assertEqual(attempted_structures, 3)
        self.assertEqual({summary.topology_id for summary in summaries}, {"sql", "kgm", "htb"})
        self.assertTrue(all(summary.pair_mode == "4+2-node-linker" for summary in summaries))
        instance_counts = {
            candidate.metadata["net_plan"]["topology"]: candidate.metadata["graph_summary"]["n_monomer_instances"]
            for candidate in candidates
        }
        self.assertEqual(instance_counts, {"sql": 3, "kgm": 9, "htb": 18})

    def test_reverse_polarity_two_plus_four_pair_returns_explicit_2d_topologies_when_requested(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=5,
                single_node_topology_ids=("sql", "kgm", "htb"),
            )
        )
        amine = BatchMonomerRecord(
            id="ppd",
            name="ppd",
            smiles=PPD,
            motif_kind="amine",
            expected_connectivity=2,
        )
        aldehyde = BatchMonomerRecord(
            id="tetra_aldehyde",
            name="tetra_aldehyde",
            smiles=TETRA_ALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=4,
        )

        summaries, candidates, attempted_structures = generator.generate_pair_candidates(amine, aldehyde)

        self.assertEqual(attempted_structures, 3)
        self.assertEqual({summary.topology_id for summary in summaries}, {"sql", "kgm", "htb"})
        self.assertTrue(all(summary.pair_mode == "4+2-node-linker" for summary in summaries))
        self.assertEqual({candidate.metadata["net_plan"]["topology"] for candidate in candidates}, {"sql", "kgm", "htb"})

    def test_four_plus_four_pair_uses_expanded_default_topology_pool(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=5,
            )
        )
        amine = BatchMonomerRecord(
            id="tetra_amine",
            name="tetra_amine",
            smiles=TETRA_AMINE,
            motif_kind="amine",
            expected_connectivity=4,
        )
        aldehyde = BatchMonomerRecord(
            id="tetra_aldehyde",
            name="tetra_aldehyde",
            smiles=TETRA_ALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=4,
        )

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "4+4-node-node")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(
            set(summary.metadata["topology_selection"]["available_topologies"]),
            {"dia", "sql", "pts", "lon", "qtz"},
        )
        self.assertIn(candidate.metadata["net_plan"]["topology"], {"dia", "sql", "pts", "lon", "qtz"})
        self.assertGreater(candidate.metadata["graph_summary"]["n_reaction_events"], 0)

    def test_four_plus_four_pair_can_target_sql_bipartite_topology(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=5,
                single_node_topology_ids=("sql",),
            )
        )
        amine = BatchMonomerRecord(
            id="tetra_amine",
            name="tetra_amine",
            smiles=TETRA_AMINE,
            motif_kind="amine",
            expected_connectivity=4,
        )
        aldehyde = BatchMonomerRecord(
            id="tetra_aldehyde",
            name="tetra_aldehyde",
            smiles=TETRA_ALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=4,
        )

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "4+4-node-node")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "sql")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "single-node-expanded-node-node")
        self.assertEqual(candidate.metadata["graph_summary"]["n_monomer_instances"], 4)
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 8)

    def test_six_plus_two_pair_uses_expanded_default_topology_pool(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=5,
            )
        )
        amine = BatchMonomerRecord(
            id="hexa_amine",
            name="hexa_amine",
            smiles=HEXA_AMINE,
            motif_kind="amine",
            expected_connectivity=6,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "6+2-node-linker")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(
            set(summary.metadata["topology_selection"]["available_topologies"]),
            {"pcu", "hxl", "acs"},
        )
        self.assertIn(candidate.metadata["net_plan"]["topology"], {"pcu", "hxl", "acs"})
        self.assertGreater(candidate.metadata["graph_summary"]["n_reaction_events"], 0)

    def test_six_plus_two_pair_can_target_hxl_when_requested(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=5,
                single_node_topology_ids=("hxl",),
            )
        )
        amine = BatchMonomerRecord(
            id="hexa_amine",
            name="hexa_amine",
            smiles=HEXA_AMINE,
            motif_kind="amine",
            expected_connectivity=6,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "6+2-node-linker")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "hxl")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "node-linker-single-node-expanded")
        self.assertEqual(candidate.metadata["graph_summary"]["n_monomer_instances"], 4)
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 6)

    def test_extended_topology_selection_respects_single_node_category_rules(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=5,
            )
        )
        amine_4 = BatchMonomerRecord(
            id="tetra_amine",
            name="tetra_amine",
            smiles=TETRA_AMINE,
            motif_kind="amine",
            expected_connectivity=4,
        )
        aldehyde_4 = BatchMonomerRecord(
            id="tetra_aldehyde",
            name="tetra_aldehyde",
            smiles=TETRA_ALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=4,
        )
        amine_6 = BatchMonomerRecord(
            id="hexa_amine",
            name="hexa_amine",
            smiles=HEXA_AMINE,
            motif_kind="amine",
            expected_connectivity=6,
        )

        four_plus_two_summary, _ = generator.generate_pair_candidate(
            amine_4,
            BatchMonomerRecord(
                id="tpal",
                name="tpal",
                smiles=TEREPHTHALALDEHYDE,
                motif_kind="aldehyde",
                expected_connectivity=2,
            ),
        )
        four_plus_four_summary, _ = generator.generate_pair_candidate(amine_4, aldehyde_4)
        six_plus_two_summary, _ = generator.generate_pair_candidate(
            amine_6,
            BatchMonomerRecord(
                id="tpal",
                name="tpal",
                smiles=TEREPHTHALALDEHYDE,
                motif_kind="aldehyde",
                expected_connectivity=2,
            ),
        )

        self.assertEqual(
            set(four_plus_two_summary.metadata["topology_selection"]["available_topologies"]),
            {"dia", "sql", "kgm", "pts", "lon", "qtz"},
        )
        self.assertEqual(
            set(four_plus_four_summary.metadata["topology_selection"]["available_topologies"]),
            {"dia", "sql", "pts", "lon", "qtz"},
        )
        self.assertEqual(
            set(six_plus_two_summary.metadata["topology_selection"]["available_topologies"]),
            {"pcu", "hxl", "acs"},
        )

    def test_default_selector_includes_curated_compatible_topologies(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=5,
            )
        )

        three_plus_two = generator._topology_ids_for_pair(connectivities=(3, 2), pair_mode="node-linker")
        three_plus_three = generator._topology_ids_for_pair(connectivities=(3, 3), pair_mode="node-node")
        three_plus_four = generator._topology_ids_for_pair(connectivities=(3, 4), pair_mode="node-node")
        three_plus_six = generator._topology_ids_for_pair(connectivities=(3, 6), pair_mode="node-node")
        four_plus_two = generator._topology_ids_for_pair(connectivities=(4, 2), pair_mode="node-linker")
        four_plus_four = generator._topology_ids_for_pair(connectivities=(4, 4), pair_mode="node-node")
        six_plus_two = generator._topology_ids_for_pair(connectivities=(6, 2), pair_mode="node-linker")
        six_plus_six = generator._topology_ids_for_pair(connectivities=(6, 6), pair_mode="node-node")

        self.assertIn("srs", three_plus_two)
        self.assertIn("srs", three_plus_three)
        self.assertIn("ctn", three_plus_four)
        self.assertIn("bor", three_plus_four)
        self.assertIn("tbo", three_plus_four)
        self.assertIn("kgd", three_plus_six)
        self.assertIn("sql", four_plus_two)
        self.assertIn("kgm", four_plus_two)
        self.assertIn("pts", four_plus_two)
        self.assertIn("lon", four_plus_two)
        self.assertIn("qtz", four_plus_two)
        self.assertIn("sql", four_plus_four)
        self.assertIn("pts", four_plus_four)
        self.assertIn("lon", four_plus_four)
        self.assertIn("qtz", four_plus_four)
        self.assertIn("hxl", six_plus_two)
        self.assertIn("acs", six_plus_two)
        self.assertIn("acs", six_plus_six)

    def test_file_driven_batch_runner_writes_manifest(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_root = Path(temp_dir)
            input_dir = temp_root / "input"
            output_dir = temp_root / "out"
            input_dir.mkdir()
            (input_dir / "amines_count_2.txt").write_text("smiles\nNc1ccc(N)cc1\n", encoding="utf-8")
            (input_dir / "amines_count_3.txt").write_text(f"smiles\n{TAPB}\n", encoding="utf-8")
            (input_dir / "aldehydes_count_2.txt").write_text(f"smiles\n{TEREPHTHALALDEHYDE}\n", encoding="utf-8")
            (input_dir / "aldehydes_count_3.txt").write_text(f"smiles\n{TFB}\n", encoding="utf-8")

            summary = self.generator.run_imine_batch(
                input_dir,
                output_dir,
                max_pairs=3,
                write_cif=False,
            )

            self.assertEqual(summary.attempted_pairs, 3)
            self.assertEqual(summary.successful_pairs, 3)
            self.assertTrue((output_dir / "manifest.jsonl").exists())
            self.assertTrue((output_dir / "summary.md").exists())
            self.assertEqual(sum(1 for _ in (output_dir / "manifest.jsonl").open()), 3)

    def test_file_driven_batch_runner_parallelizes_pair_generation(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
                max_workers=2,
            )
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_root = Path(temp_dir)
            input_dir = temp_root / "input"
            output_dir = temp_root / "out"
            input_dir.mkdir()
            (input_dir / "amines_count_2.txt").write_text("smiles\nNc1ccc(N)cc1\n", encoding="utf-8")
            (input_dir / "amines_count_3.txt").write_text(f"smiles\n{TAPB}\n", encoding="utf-8")
            (input_dir / "aldehydes_count_2.txt").write_text(f"smiles\n{TEREPHTHALALDEHYDE}\n", encoding="utf-8")
            (input_dir / "aldehydes_count_3.txt").write_text(f"smiles\n{TFB}\n", encoding="utf-8")

            summary = generator.run_imine_batch(
                input_dir,
                output_dir,
                max_pairs=4,
                write_cif=False,
            )

            self.assertEqual(summary.attempted_pairs, 3)
            self.assertEqual(summary.successful_pairs, 3)
            self.assertTrue((output_dir / "manifest.jsonl").exists())
            self.assertEqual(
                len(
                    [
                        line
                        for line in (output_dir / "manifest.jsonl").read_text(encoding="utf-8").splitlines()
                        if line
                    ]
                ),
                summary.successful_structures,
            )

    def test_file_driven_batch_runner_enumerates_all_topologies_by_default(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
            )
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_root = Path(temp_dir)
            input_dir = temp_root / "input"
            output_dir = temp_root / "out"
            input_dir.mkdir()
            (input_dir / "amines_count_2.txt").write_text("smiles\n", encoding="utf-8")
            (input_dir / "amines_count_3.txt").write_text(f"smiles\n{TAPB}\n", encoding="utf-8")
            (input_dir / "aldehydes_count_2.txt").write_text(f"smiles\n{TEREPHTHALALDEHYDE}\n", encoding="utf-8")
            (input_dir / "aldehydes_count_3.txt").write_text("smiles\n", encoding="utf-8")

            summary = generator.run_imine_batch(
                input_dir,
                output_dir,
            )

            manifest_rows = [line for line in (output_dir / "manifest.jsonl").read_text(encoding="utf-8").splitlines() if line]
            cif_count = sum(1 for _ in (output_dir / "cifs").rglob("*.cif"))

            self.assertEqual(summary.attempted_pairs, 1)
            self.assertEqual(summary.successful_pairs, 1)
            self.assertEqual(summary.attempted_structures, 5)
            self.assertEqual(summary.successful_structures, 5)
            self.assertEqual(len(manifest_rows), 5)
            self.assertEqual(cif_count, 5)
            self.assertTrue((output_dir / "cifs" / "valid").is_dir())

    def test_file_driven_batch_runner_discovers_four_and_six_connectivity_libraries(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=5,
            )
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_root = Path(temp_dir)
            input_dir = temp_root / "input"
            output_dir = temp_root / "out"
            input_dir.mkdir()
            (input_dir / "amines_count_2.txt").write_text("smiles\nNc1ccc(N)cc1\n", encoding="utf-8")
            (input_dir / "amines_count_4.txt").write_text(f"smiles\n{TETRA_AMINE}\n", encoding="utf-8")
            (input_dir / "amines_count_6.txt").write_text(f"smiles\n{HEXA_AMINE}\n", encoding="utf-8")
            (input_dir / "aldehydes_count_2.txt").write_text(f"smiles\n{TEREPHTHALALDEHYDE}\n", encoding="utf-8")
            (input_dir / "aldehydes_count_4.txt").write_text(f"smiles\n{TETRA_ALDEHYDE}\n", encoding="utf-8")
            (input_dir / "aldehydes_count_6.txt").write_text(f"smiles\n{HEXA_ALDEHYDE}\n", encoding="utf-8")

            summary = generator.run_imine_batch(
                input_dir,
                output_dir,
                write_cif=False,
            )

            self.assertEqual(summary.attempted_pairs, 8)
            self.assertEqual(summary.successful_pairs, 8)
            self.assertGreaterEqual(summary.attempted_structures, 8)
            self.assertGreater(summary.attempted_structures, summary.successful_structures)
            self.assertGreater(summary.successful_structures, 8)
            self.assertIn("4+4-node-node", summary.mode_counts)
            self.assertIn("4+6-node-node", summary.mode_counts)
            self.assertIn("4+2-node-linker", summary.mode_counts)
            self.assertIn("6+2-node-linker", summary.mode_counts)
            self.assertIn("dia", summary.topology_counts)
            self.assertIn("pcu", summary.topology_counts)
            self.assertIn("cor", summary.topology_counts)

    def _max_bridge_distance_error(self, candidate) -> float:
        return max(
            abs(float(metric["actual_distance"]) - float(metric["target_distance"]))
            for metric in candidate.metadata["score_metadata"]["bridge_event_metrics"]
        )

    def _max_linker_translation_mismatch(self, candidate) -> float:
        return max(
            float(details["translation_mismatch"])
            for details in candidate.metadata["embedding"]["poses"].values()
            if details.get("role") == "linker"
        )


if __name__ == "__main__":
    unittest.main()
