import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit import CIFWriter, COFEngine, COFProject, build_rdkit_monomer
from cofkit.reaction_realization import ReactionRealizer

try:
    from rdkit import Chem  # noqa: F401
except ImportError:  # pragma: no cover - environment-dependent
    Chem = None


@unittest.skipIf(Chem is None, "RDKit is not available")
class RDKitMonomerTests(unittest.TestCase):
    def test_build_rdkit_monomer_detects_tapb_contextual_amines(self):
        monomer = build_rdkit_monomer(
            "tapb",
            "1,3,5-tris(4-aminophenyl)benzene",
            "C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N",
            "amine",
        )

        self.assertEqual(len(monomer.motifs), 3)
        self.assertEqual(monomer.metadata["geometry_mode"], "rdkit-etkdg")
        self.assertTrue(all(motif.kind == "amine" for motif in monomer.motifs))
        self.assertTrue(all("imine_bridge" in motif.allowed_reaction_templates for motif in monomer.motifs))
        self.assertEqual(len(monomer.atom_symbols), len(monomer.atom_positions))

    def test_build_rdkit_monomer_detects_tfb_contextual_aldehydes(self):
        monomer = build_rdkit_monomer(
            "tfb",
            "1,3,5-benzenetricarbaldehyde",
            "C1=C(C=C(C=C1C=O)C=O)C=O",
            "aldehyde",
        )

        self.assertEqual(len(monomer.motifs), 3)
        self.assertTrue(all(motif.kind == "aldehyde" for motif in monomer.motifs))
        self.assertTrue(all("imine_bridge" in motif.allowed_reaction_templates for motif in monomer.motifs))
        self.assertEqual(len(monomer.atom_symbols), len(monomer.atom_positions))

    def test_build_rdkit_monomer_detects_hydrazide(self):
        monomer = build_rdkit_monomer(
            "benzohydrazide",
            "benzohydrazide",
            "NNC(=O)c1ccccc1",
            "hydrazide",
        )

        self.assertEqual(len(monomer.motifs), 1)
        motif = monomer.motifs[0]
        self.assertEqual(motif.kind, "hydrazide")
        self.assertIn("hydrazone_bridge", motif.allowed_reaction_templates)
        self.assertIn("carbonyl_carbon_atom_id", motif.metadata)

    def test_build_rdkit_monomer_detects_boronic_acid_and_catechol(self):
        boronic = build_rdkit_monomer(
            "phenyl_boronic_acid",
            "phenyl boronic acid",
            "OB(O)c1ccccc1",
            "boronic_acid",
        )
        catechol = build_rdkit_monomer(
            "catechol",
            "catechol",
            "Oc1ccccc1O",
            "catechol",
        )

        self.assertEqual(len(boronic.motifs), 1)
        self.assertEqual(boronic.motifs[0].kind, "boronic_acid")
        self.assertEqual(boronic.motifs[0].metadata["oxygen_atom_ids"].__len__(), 2)
        self.assertEqual(len(catechol.motifs), 1)
        self.assertEqual(catechol.motifs[0].kind, "catechol")
        self.assertEqual(catechol.motifs[0].metadata["reactive_atom_ids"], tuple(sorted(catechol.motifs[0].metadata["reactive_atom_ids"])))

    def test_build_rdkit_monomer_detects_corrected_hhtp_as_tritopic_catechol(self):
        monomer = build_rdkit_monomer(
            "HHTP",
            "hexahydroxytriphenylene",
            "OC1=C(O)C=C2C(=C1)C1=CC(O)=C(O)C=C1C1=CC(O)=C(O)C=C21",
            "catechol",
        )

        self.assertEqual(len(monomer.motifs), 3)
        self.assertTrue(all(motif.kind == "catechol" for motif in monomer.motifs))

    def test_build_rdkit_monomer_detects_keto_aldehyde_and_activated_methylene(self):
        keto_aldehyde = build_rdkit_monomer(
            "salicylaldehyde",
            "salicylaldehyde",
            "O=Cc1ccccc1O",
            "keto_aldehyde",
        )
        activated_methylene = build_rdkit_monomer(
            "malononitrile",
            "malononitrile",
            "N#CCC#N",
            "activated_methylene",
        )

        self.assertEqual(len(keto_aldehyde.motifs), 1)
        self.assertEqual(keto_aldehyde.motifs[0].kind, "keto_aldehyde")
        self.assertIn("ortho_hydroxyl_oxygen_atom_id", keto_aldehyde.motifs[0].metadata)
        self.assertIn("ortho_hydroxyl_anchor_atom_id", keto_aldehyde.motifs[0].metadata)
        self.assertEqual(len(activated_methylene.motifs), 1)
        self.assertEqual(activated_methylene.motifs[0].kind, "activated_methylene")
        self.assertGreaterEqual(len(activated_methylene.motifs[0].metadata["hydrogen_atom_ids"]), 2)

    def test_build_rdkit_monomer_assigns_distinct_keto_aldehyde_hydroxyls_for_tp(self):
        monomer = build_rdkit_monomer(
            "tp",
            "1,3,5-triformylphloroglucinol",
            "O=Cc1c(O)c(C=O)c(O)c(C=O)c1O",
            "keto_aldehyde",
        )

        self.assertEqual(len(monomer.motifs), 3)
        hydroxyl_oxygen_ids = {
            motif.metadata["ortho_hydroxyl_oxygen_atom_id"]
            for motif in monomer.motifs
        }
        self.assertEqual(len(hydroxyl_oxygen_ids), 3)

    def test_build_rdkit_monomer_detects_tmt_as_tritopic_activated_methylene(self):
        monomer = build_rdkit_monomer(
            "TMT",
            "2,4,6-trimethyl-s-triazine",
            "Cc1nc(C)nc(C)n1",
            "activated_methylene",
        )

        self.assertEqual(len(monomer.motifs), 3)
        self.assertTrue(all(motif.kind == "activated_methylene" for motif in monomer.motifs))

    def test_rdkit_hcb_case_keeps_bipartite_periodic_images_and_atomistic_export(self):
        tapb = build_rdkit_monomer(
            "tapb",
            "1,3,5-tris(4-aminophenyl)benzene",
            "C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N",
            "amine",
        )
        tfb = build_rdkit_monomer(
            "tfb",
            "1,3,5-benzenetricarbaldehyde",
            "C1=C(C=C(C=C1C=O)C=O)C=O",
            "aldehyde",
        )
        project = COFProject(
            monomers=(tapb, tfb),
            allowed_reactions=("imine_bridge",),
            target_dimensionality="2D",
            target_topologies=("hcb",),
        )

        candidate = COFEngine().run(project).top(1)[0]
        result = CIFWriter().export_candidate(candidate, project.monomers, data_name="tapb_tfb_hcb_rdkit_test")

        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "single-node-bipartite")
        self.assertEqual(
            [event.participants[1].periodic_image for event in candidate.events],
            [(0, 0, 0), (-1, 0, 0), (0, -1, 0)],
        )
        self.assertEqual(result.mode, "atomistic")
        self.assertGreater(result.n_sites, 40)

    def test_rdkit_hcb_case_exports_reacted_imine_product(self):
        tapb = build_rdkit_monomer(
            "tapb",
            "1,3,5-tris(4-aminophenyl)benzene",
            "C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N",
            "amine",
        )
        tfb = build_rdkit_monomer(
            "tfb",
            "1,3,5-benzenetricarbaldehyde",
            "C1=C(C=C(C=C1C=O)C=O)C=O",
            "aldehyde",
        )
        project = COFProject(
            monomers=(tapb, tfb),
            allowed_reactions=("imine_bridge",),
            target_dimensionality="2D",
            target_topologies=("hcb",),
        )

        candidate = COFEngine().run(project).top(1)[0]
        result = CIFWriter().export_candidate(candidate, project.monomers, data_name="tapb_tfb_hcb_rdkit_imine_product_test")

        realization = result.metadata["reaction_realization"]
        self.assertEqual(realization["applied_event_count"], 3)
        self.assertEqual(realization["applied_templates"], {"imine_bridge": 3})
        self.assertEqual(realization["removed_atom_symbols"], {"H": 3, "O": 3})
        self.assertEqual(result.n_sites, len(tapb.atom_symbols) + len(tfb.atom_symbols) - 6)
        self.assertGreater(len(tapb.bonds), 0)
        self.assertGreater(len(tfb.bonds), 0)

        atom_lines = [
            line for line in result.text.splitlines()
            if line.startswith("m") and len(line.split()) == 6
        ]
        self.assertTrue(atom_lines)
        self.assertFalse(any(line.split()[1] == "O" for line in atom_lines))
        self.assertIn("_geom_bond_atom_site_label_1", result.text)
        first_event = candidate.events[0]
        first_ref, second_ref = first_event.participants
        if first_ref.monomer_id == "tapb":
            amine_ref, aldehyde_ref = first_ref, second_ref
            amine_spec, aldehyde_spec = tapb, tfb
        else:
            amine_ref, aldehyde_ref = second_ref, first_ref
            amine_spec, aldehyde_spec = tapb, tfb
        amine_motif = amine_spec.motif_by_id(amine_ref.motif_id)
        aldehyde_motif = aldehyde_spec.motif_by_id(aldehyde_ref.motif_id)
        expected_prefix = (
            f"{ReactionRealizer.atom_label(aldehyde_ref.monomer_instance_id, 'C', aldehyde_motif.metadata['reactive_atom_id'])} "
            f"{ReactionRealizer.atom_label(amine_ref.monomer_instance_id, 'N', amine_motif.metadata['reactive_atom_id'])}"
        )
        bond_lines = [line for line in result.text.splitlines() if line.startswith("m") and len(line.split()) == 5]
        matching_bonds = [line for line in bond_lines if line.startswith(expected_prefix)]
        self.assertEqual(len(matching_bonds), 1)
        self.assertAlmostEqual(float(matching_bonds[0].split()[-1]), 1.3, places=2)
        self.assertGreater(len(bond_lines), len(candidate.events))
        self.assertTrue(any(line.startswith("m1_C1 m1_C2 ") for line in bond_lines))


if __name__ == "__main__":
    unittest.main()
