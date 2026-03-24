import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import cofkit.topologies as topologies_module

from cofkit import RCSRArchiveImporter, default_topology_repository
from cofkit.topology_symmetry import has_gemmi


_EXPLICIT_2D_TOPOLOGIES = """\
CRYSTAL
NAME ladder
GROUP P1
CELL 1.0 1.0 90.0
NODE a 0.0 0.0
NODE b 0.5 0.0
EDGE 0.0 0.0 0.5 0.0
EDGE 0.0 0.0 0.5 1.0
END

CRYSTAL
NAME triangle
GROUP P1
CELL 1.0 1.0 90.0
NODE a 0.0 0.0
NODE b 0.333333 0.0
NODE c 0.666667 0.0
EDGE 0.0 0.0 0.333333 0.0
EDGE 0.333333 0.0 0.666667 0.0
EDGE 0.666667 0.0 0.0 0.0
END
"""

_EXPLICIT_3D_SELF_EDGE = """\
CRYSTAL
NAME cubic_self
GROUP P1
CELL 1.0 1.0 1.0 90.0 90.0 90.0
NODE a 0.0 0.0 0.0
EDGE 0.0 0.0 0.0 1.0 0.0 0.0
END
"""


class ZeroLinkerTopologyAnalysisTests(unittest.TestCase):
    def test_workspace_supported_topologies_expose_zero_linker_metadata(self):
        topologies_module._DEFAULT_REPOSITORY = None
        try:
            repository = default_topology_repository()
            hcb = repository.get_index_entry("hcb")
            hca = repository.get_index_entry("hca")
            dia = repository.get_index_entry("dia")
            pcu = repository.get_index_entry("pcu")
            hcb_definition = repository.load("hcb")
        finally:
            topologies_module._DEFAULT_REPOSITORY = None

        self.assertTrue(hcb.metadata["zero_linker_compatible"])
        self.assertTrue(hcb.metadata["zero_linker_bipartite"])
        self.assertEqual(hcb.metadata["zero_linker_role_count_lower_bound"], 2)
        self.assertEqual(hcb.metadata["zero_linker_scan_method"], "expanded-p1-periodic-bipartite-scan")
        self.assertTrue(hcb.metadata["zero_linker_builder_supported"])
        self.assertTrue(hcb.metadata["two_monomer_compatible"])
        self.assertEqual(hcb.metadata["two_monomer_node_node_modes"], ["3+3"])
        self.assertEqual(hcb.metadata["two_monomer_node_linker_modes"], ["3+2"])

        self.assertFalse(hca.metadata["zero_linker_compatible"])
        self.assertFalse(hca.metadata["zero_linker_bipartite"])
        self.assertEqual(hca.metadata["zero_linker_role_count_lower_bound"], 3)
        self.assertFalse(hca.metadata["zero_linker_builder_supported"])
        self.assertTrue(hca.metadata["two_monomer_compatible"])
        self.assertEqual(hca.metadata["two_monomer_node_node_modes"], [])
        self.assertEqual(hca.metadata["two_monomer_node_linker_modes"], ["3+2"])

        self.assertTrue(dia.metadata["zero_linker_compatible"])
        self.assertTrue(dia.metadata["zero_linker_builder_supported"])
        self.assertEqual(dia.metadata["two_monomer_node_node_modes"], ["4+4"])
        self.assertEqual(dia.metadata["two_monomer_node_linker_modes"], ["4+2"])
        self.assertTrue(pcu.metadata["zero_linker_compatible"])
        self.assertFalse(pcu.metadata["zero_linker_builder_supported"])
        self.assertEqual(pcu.metadata["two_monomer_node_node_modes"], ["6+6"])
        self.assertEqual(pcu.metadata["two_monomer_node_linker_modes"], ["6+2"])

        self.assertTrue(hcb_definition.metadata["zero_linker_compatible"])
        self.assertEqual(hcb_definition.metadata["zero_linker_scan_method"], "expanded-p1-periodic-bipartite-scan")

    def test_importer_scans_explicit_quotient_graph_topologies(self):
        importer = RCSRArchiveImporter()

        definitions = {
            definition.id: definition
            for definition in importer.iter_text_definitions(
                text=_EXPLICIT_2D_TOPOLOGIES,
                dimensionality="2D",
                source_archive="synthetic://2d",
                source_member="synthetic_2d.cgd",
            )
        }

        ladder = definitions["ladder"]
        triangle = definitions["triangle"]

        self.assertTrue(ladder.metadata["zero_linker_compatible"])
        self.assertTrue(ladder.metadata["zero_linker_bipartite"])
        self.assertEqual(ladder.metadata["zero_linker_role_count_lower_bound"], 2)
        self.assertEqual(ladder.metadata["zero_linker_scan_method"], "quotient-graph-parity-scan")
        self.assertTrue(ladder.metadata["zero_linker_builder_supported"])

        self.assertFalse(triangle.metadata["zero_linker_compatible"])
        self.assertFalse(triangle.metadata["zero_linker_bipartite"])
        self.assertEqual(triangle.metadata["zero_linker_role_count_lower_bound"], 3)
        self.assertEqual(triangle.metadata["zero_linker_scan_method"], "quotient-graph-parity-scan")

        ladder_entry = importer.definition_to_entry(ladder)
        self.assertTrue(ladder_entry.metadata["zero_linker_compatible"])
        self.assertEqual(ladder_entry.metadata["zero_linker_scan_method"], "quotient-graph-parity-scan")

    def test_explicit_self_edge_shift_is_detected_as_bipartite(self):
        importer = RCSRArchiveImporter()

        definition = next(
            importer.iter_text_definitions(
                text=_EXPLICIT_3D_SELF_EDGE,
                dimensionality="3D",
                source_archive="synthetic://3d",
                source_member="synthetic_3d.cgd",
            )
        )

        self.assertTrue(definition.metadata["zero_linker_compatible"])
        self.assertTrue(definition.metadata["zero_linker_bipartite"])
        self.assertEqual(definition.metadata["zero_linker_scan_method"], "quotient-graph-parity-scan")
        self.assertFalse(definition.metadata["zero_linker_builder_supported"])

    @unittest.skipUnless(has_gemmi(), "gemmi is required for generic 3D symmetry expansion coverage")
    def test_workspace_three_d_topologies_can_scan_via_generic_symmetry_expansion(self):
        topologies_module._DEFAULT_REPOSITORY = None
        try:
            repository = default_topology_repository()
            for topology_id in ("acs", "afi", "ana", "ast", "ant"):
                definition = repository.load(topology_id)
                self.assertEqual(definition.metadata["zero_linker_scan_method"], "gemmi-space-group-expanded-p1-scan")
                self.assertIsNotNone(definition.metadata["zero_linker_compatible"])
        finally:
            topologies_module._DEFAULT_REPOSITORY = None

    @unittest.skipUnless(has_gemmi(), "gemmi is required for generic 2D symmetry expansion coverage")
    def test_workspace_two_d_topologies_can_scan_via_generic_symmetry_expansion(self):
        topologies_module._DEFAULT_REPOSITORY = None
        try:
            repository = default_topology_repository()
            for topology_id in ("car", "bew", "cpa"):
                definition = repository.load(topology_id)
                self.assertEqual(definition.metadata["zero_linker_scan_method"], "gemmi-space-group-expanded-p1-scan")
                self.assertIsNotNone(definition.metadata["zero_linker_compatible"])
            self.assertEqual(repository.load("car").metadata["two_monomer_node_node_modes"], [])
            self.assertEqual(repository.load("car").metadata["two_monomer_node_linker_modes"], ["3+2"])
            self.assertEqual(repository.load("bew").metadata["two_monomer_node_node_modes"], [])
            self.assertEqual(repository.load("bew").metadata["two_monomer_node_linker_modes"], ["4+2"])
            self.assertFalse(repository.load("cpa").metadata["two_monomer_compatible"])
            cor = repository.load("cor")
            self.assertTrue(cor.metadata["two_monomer_compatible"])
            self.assertEqual(cor.metadata["two_monomer_node_node_modes"], ["4+6"])
            self.assertEqual(cor.metadata["two_monomer_node_linker_modes"], [])
        finally:
            topologies_module._DEFAULT_REPOSITORY = None


if __name__ == "__main__":
    unittest.main()
