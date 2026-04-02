import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit.chem import Molecule, MotifDetector, default_motif_kind_registry


class MotifDetectorTests(unittest.TestCase):
    def test_builtin_registry_exposes_supported_motif_kinds(self):
        registry = default_motif_kind_registry()

        self.assertIn("amine", registry.supported_kinds())
        self.assertIn("aldehyde", registry.supported_kinds())
        self.assertIn("hydrazine", registry.supported_kinds())
        self.assertIn("hydrazide", registry.supported_kinds())
        self.assertIn("activated_methylene", registry.supported_kinds())
        self.assertEqual(registry.get("boronic_acid").cif_symbol, "B")
        self.assertEqual(registry.get("amine").allowed_reaction_templates, ("imine_bridge", "keto_enamine_bridge"))
        self.assertEqual(registry.get("hydrazine").allowed_reaction_templates, ("azine_bridge",))
        self.assertEqual(registry.get("hydrazide").allowed_reaction_templates, ("hydrazone_bridge",))

    def test_detect_amine(self):
        xyz = """14
Aniline
C  0.0  0.0  0.0
C  1.4  0.0  0.0
C  2.1  1.2  0.0
C  1.4  2.4  0.0
C  0.0  2.4  0.0
C -0.7  1.2  0.0
N  2.8 -1.2  0.0
H -0.5 -0.9  0.0
H  3.2  1.2  0.0
H  1.9  3.3  0.0
H -0.5  3.3  0.0
H -1.8  1.2  0.0
H  3.8 -1.1 -0.1
H  2.5 -2.1  0.1
"""
        mol = Molecule.from_xyz(xyz)
        spec = MotifDetector.detect_amine(mol, "aniline")

        self.assertEqual(len(spec.motifs), 1)
        motif = spec.motifs[0]
        self.assertEqual(motif.kind, "amine")
        self.assertEqual(motif.atom_ids[0], 6)
        self.assertEqual(motif.allowed_reaction_templates, ("imine_bridge", "keto_enamine_bridge"))
        self.assertEqual(motif.metadata["reactive_atom_id"], 6)
        self.assertEqual(len(spec.atom_symbols), len(spec.atom_positions))
        self.assertEqual(spec.atom_symbols[6], "N")
        self.assertEqual(spec.bonds, ())

    def test_detect_aldehyde(self):
        xyz = """14
Benzaldehyde
C  0.0  0.0  0.0
C  1.4  0.0  0.0
C  2.1  1.2  0.0
C  1.4  2.4  0.0
C  0.0  2.4  0.0
C -0.7  1.2  0.0
C  2.1 -1.4  0.0
O  3.3 -1.5  0.0
H  1.4 -2.3  0.0
H -0.5 -0.9  0.0
H  3.2  1.2  0.0
H  1.9  3.3  0.0
H -0.5  3.3  0.0
H -1.8  1.2  0.0
"""
        mol = Molecule.from_xyz(xyz)
        spec = MotifDetector.detect_aldehyde(mol, "benzaldehyde")

        self.assertEqual(len(spec.motifs), 1)
        motif = spec.motifs[0]
        self.assertEqual(motif.kind, "aldehyde")
        self.assertEqual(motif.atom_ids[0], 6)
        self.assertEqual(motif.atom_ids[1], 7)
        self.assertEqual(motif.atom_ids[2], 8)
        self.assertEqual(
            motif.allowed_reaction_templates,
            ("imine_bridge", "hydrazone_bridge", "azine_bridge", "vinylene_bridge"),
        )
        self.assertEqual(motif.metadata["reactive_atom_id"], 6)
        self.assertEqual(len(spec.atom_symbols), len(spec.atom_positions))
        self.assertEqual(spec.atom_symbols[7], "O")
        self.assertEqual(spec.bonds, ())

    def test_generic_detect_api_uses_registered_handler(self):
        xyz = """14
Aniline
C  0.0  0.0  0.0
C  1.4  0.0  0.0
C  2.1  1.2  0.0
C  1.4  2.4  0.0
C  0.0  2.4  0.0
C -0.7  1.2  0.0
N  2.8 -1.2  0.0
H -0.5 -0.9  0.0
H  3.2  1.2  0.0
H  1.9  3.3  0.0
H -0.5  3.3  0.0
H -1.8  1.2  0.0
H  3.8 -1.1 -0.1
H  2.5 -2.1  0.1
"""
        mol = Molecule.from_xyz(xyz)
        spec = MotifDetector.builtin().detect(mol, "aniline", "amine")

        self.assertEqual(len(spec.motifs), 1)
        self.assertEqual(spec.motifs[0].kind, "amine")


if __name__ == "__main__":
    unittest.main()
