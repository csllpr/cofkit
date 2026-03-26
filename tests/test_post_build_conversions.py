import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit import (
    BatchGenerationConfig,
    BatchStructureGenerator,
    COFEngine,
    COFProject,
    CIFWriter,
    assess_imine_cof_aldehyde_monomer,
    assess_imine_cof_amine_monomer,
    build_rdkit_monomer,
)

try:
    from rdkit import Chem  # noqa: F401
except ImportError:  # pragma: no cover - environment-dependent
    Chem = None


TAPB = "C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N"
TFB = "C1=C(C=C(C=C1C=O)C=O)C=O"
TAPT = "C1=CC(=CC=C1C2=NC(=NC(=N2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N"
BTT = "C1=C(SC2=C1C3=C(C=C(S3)C=O)C4=C2C=C(S4)C=O)C=O"
PPD = "Nc1ccc(N)cc1"
TEREPHTHALALDEHYDE = "O=Cc1ccc(C=O)cc1"
ETHYLENEDIAMINE = "NCCN"


@unittest.skipIf(Chem is None, "RDKit is not available")
class PostBuildConversionTests(unittest.TestCase):
    def test_primary_aryl_diamine_is_imine_eligible(self):
        monomer = build_rdkit_monomer(
            "ppd",
            "ppd",
            PPD,
            "amine",
            num_conformers=2,
        )

        assessment = assess_imine_cof_amine_monomer(monomer)

        self.assertTrue(assessment.eligible)
        self.assertEqual(assessment.motif_count, 2)
        self.assertEqual(assessment.reason_codes, ())

    def test_aliphatic_diamine_is_not_imine_cof_eligible(self):
        monomer = build_rdkit_monomer(
            "eda",
            "eda",
            ETHYLENEDIAMINE,
            "amine",
            num_conformers=2,
        )

        assessment = assess_imine_cof_amine_monomer(monomer)

        self.assertFalse(assessment.eligible)
        self.assertIn("amine_anchor_not_aromatic_sp2", assessment.reason_codes)

    def test_dialdehyde_is_imine_eligible(self):
        monomer = build_rdkit_monomer(
            "tpal",
            "tpal",
            TEREPHTHALALDEHYDE,
            "aldehyde",
            num_conformers=2,
        )

        assessment = assess_imine_cof_aldehyde_monomer(monomer)

        self.assertTrue(assessment.eligible)
        self.assertEqual(assessment.motif_count, 2)
        self.assertEqual(assessment.reason_codes, ())

    def test_engine_can_annotate_sulfur_enabled_imine_conversion_candidates(self):
        candidate = self._build_engine_candidate(
            post_build_conversions=("sulfur_enabled_imine_conversion",),
        )

        self.assertEqual(candidate.metadata["graph_summary"]["reaction_templates"], {"imine_bridge": 3})
        self.assertIn("post_build_conversions", candidate.metadata)
        annotation = candidate.metadata["post_build_conversions"]["sulfur_enabled_imine_conversion"]
        self.assertTrue(annotation["eligible"])
        self.assertEqual(annotation["canonical_built_linkage"], "imine_bridge")
        self.assertTrue(annotation["changes_applied"])
        self.assertEqual(annotation["required_external_conditions"], ["sulfur_enabled"])
        self.assertEqual(annotation["details"]["conversion_product_family"], "benzothiazole_linkage")
        self.assertEqual(annotation["details"]["n_event_plans"], 3)

    def test_engine_does_not_annotate_when_not_requested(self):
        candidate = self._build_engine_candidate()

        self.assertNotIn("post_build_conversions", candidate.metadata)
        self.assertEqual(candidate.metadata["graph_summary"]["reaction_templates"], {"imine_bridge": 3})

    def test_batch_generator_surfaces_post_build_conversion_annotation_in_summary(self):
        tri_amine = build_rdkit_monomer(
            "tapb",
            "tapb",
            TAPB,
            "amine",
            num_conformers=2,
        )
        tri_aldehyde = build_rdkit_monomer(
            "tfb",
            "tfb",
            TFB,
            "aldehyde",
            num_conformers=2,
        )
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                allowed_reactions=("imine_bridge",),
                single_node_topology_ids=("hcb",),
                rdkit_num_conformers=2,
                post_build_conversions=("sulfur_enabled_imine_conversion",),
                write_cif=False,
            )
        )

        summary, candidate = generator.generate_monomer_pair_candidate(
            tri_aldehyde,
            tri_amine,
            write_cif=False,
        )

        self.assertEqual(summary.status, "ok")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertIn("post_build_conversions", candidate.metadata)
        self.assertIn("post_build_conversions", summary.metadata)
        annotation = summary.metadata["post_build_conversions"]["sulfur_enabled_imine_conversion"]
        self.assertTrue(annotation["eligible"])
        self.assertEqual(annotation["canonical_built_linkage"], "imine_bridge")
        self.assertTrue(annotation["changes_applied"])

    def test_cif_export_applies_atomistic_sulfur_conversion_for_requested_candidate(self):
        candidate = self._build_engine_candidate(
            post_build_conversions=("sulfur_enabled_imine_conversion",),
        )
        tri_amine = build_rdkit_monomer(
            "tapb",
            "tapb",
            TAPB,
            "amine",
            num_conformers=2,
        )
        tri_aldehyde = build_rdkit_monomer(
            "tfb",
            "tfb",
            TFB,
            "aldehyde",
            num_conformers=2,
        )

        export = CIFWriter().export_candidate(
            candidate,
            {tri_amine.id: tri_amine, tri_aldehyde.id: tri_aldehyde},
        )

        realization = export.metadata["reaction_realization"]
        self.assertEqual(realization["applied_templates"], {"imine_bridge": 3})
        self.assertEqual(realization["added_atom_symbols"], {"S": 3})
        self.assertEqual(realization["post_build_conversions"]["applied_profiles"]["sulfur_enabled_imine_conversion"]["n_added_sulfur_atoms"], 3)
        self.assertIn("benzothiazole-like annulation", " ".join(realization["notes"]))
        self.assertIn("_geom_bond_atom_site_label_1", export.text)
        self.assertIn(" S ", export.text)

    def test_benzothiazole_conversion_keeps_periodic_edge_closures_consistent(self):
        tapt = build_rdkit_monomer(
            "tapt",
            "tapt",
            TAPT,
            "amine",
            num_conformers=2,
        )
        btt = build_rdkit_monomer(
            "btt",
            "btt",
            BTT,
            "aldehyde",
            num_conformers=2,
        )
        candidate = COFEngine().run(
            COFProject(
                monomers=(tapt, btt),
                allowed_reactions=("imine_bridge",),
                target_dimensionality="2D",
                target_topologies=("hcb",),
                post_build_conversions=("sulfur_enabled_imine_conversion",),
            )
        ).top(1)[0]

        monomer_specs = {tapt.id: tapt, btt.id: btt}
        instance_to_monomer = {
            participant.monomer_instance_id: participant.monomer_id
            for event in candidate.events
            for participant in event.participants
        }
        realizer = CIFWriter().reaction_realizer
        result = realizer.realize(candidate, monomer_specs, instance_to_monomer)

        self.assertIsNotNone(result)
        assert result is not None
        realized_local_positions = {
            instance_id: {atom.atom_id: atom.local_position for atom in atoms}
            for instance_id, atoms in result.atoms_by_instance.items()
        }
        sulfur_atom_ids = sorted(
            atom.atom_id
            for atom in result.atoms_by_instance["m1"]
            if atom.symbol == "S" and atom.atom_id >= len(tapt.atom_symbols)
        )
        sulfur_geometry = result.metadata["post_build_conversions"]["applied_profiles"]["sulfur_enabled_imine_conversion"]["sulfur_geometry"]

        for sulfur_atom_id, event in zip(sulfur_atom_ids, candidate.events, strict=True):
            amine_ref, aldehyde_ref = event.participants
            if amine_ref.monomer_id != "tapt":
                amine_ref, aldehyde_ref = aldehyde_ref, amine_ref
            event_geometry = sulfur_geometry[event.id]
            sulfur_world = realizer._world_position(
                candidate.state.monomer_poses[amine_ref.monomer_instance_id],
                realized_local_positions[amine_ref.monomer_instance_id][sulfur_atom_id],
            )
            self.assertIsNotNone(sulfur_world)
            self.assertEqual(event_geometry["sulfur_opposite_bend"], 1.0)
            self.assertEqual(event_geometry["local_relaxation_applied"], 0.0)
            self.assertLess(event_geometry["actual_aldehyde_anchor_c_distance"], 1.6)
            self.assertLess(event_geometry["actual_carbon_s_distance"], 2.0)
            self.assertLess(event_geometry["actual_ortho_s_distance"], 2.5)
            self.assertGreater(event_geometry["nitrogen_sulfur_distance"], 2.2)
            self.assertGreater(event_geometry["n_c_s_angle"], 90.0)
            self.assertGreater(event_geometry["c_s_ortho_angle"], 90.0)

    def _build_engine_candidate(self, *, post_build_conversions: tuple[str, ...] = ()):
        tri_amine = build_rdkit_monomer(
            "tapb",
            "tapb",
            TAPB,
            "amine",
            num_conformers=2,
        )
        tri_aldehyde = build_rdkit_monomer(
            "tfb",
            "tfb",
            TFB,
            "aldehyde",
            num_conformers=2,
        )
        project = COFProject(
            monomers=(tri_amine, tri_aldehyde),
            allowed_reactions=("imine_bridge",),
            target_dimensionality="2D",
            target_topologies=("hcb",),
            post_build_conversions=post_build_conversions,
        )

        return COFEngine().run(project).top(1)[0]


if __name__ == "__main__":
    unittest.main()
