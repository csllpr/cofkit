import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import cofkit.topologies as topologies_module
from cofkit.topology_data import write_workspace_topology_index

from cofkit import (
    Frame,
    MonomerSpec,
    NetPlanner,
    RCSRArchiveImporter,
    ReactiveMotif,
    TopologyRepository,
    default_topology_repository,
    imported_topology_data_metadata,
    load_topology,
)


_BUNDLE_2D = """\
CRYSTAL
NAME hcb
GROUP p6mm
CELL 1.0 1.0 120.0
NODE tri_a 0.0 0.0
NODE tri_b 0.333333 0.666667
EDGE 0.0 0.0 0.333333 0.666667
EDGE 0.0 0.0 0.333333 -0.333333
EDGE 1.0 0.0 0.333333 0.666667
END

CRYSTAL
NAME kgm
GROUP p6mm
CELL 1.0 1.0 120.0
NODE tri 3 0.0 0.0
NODE hex 6 0.5 0.5
EDGE tri hex
EDGE tri hex
EDGE tri hex
EDGE tri hex
EDGE tri hex
EDGE tri hex
END
"""

_BUNDLE_3D = """\
CRYSTAL
NAME dia
GROUP Fd-3m
CELL 1.0 1.0 1.0 90.0 90.0 90.0
NODE tetra_a 0.0 0.0 0.0
NODE tetra_b 0.25 0.25 0.25
EDGE tetra_a tetra_b
EDGE tetra_a tetra_b
EDGE tetra_a tetra_b
EDGE tetra_a tetra_b
END
"""

_WORKSPACE_HCB = """\
CRYSTAL
NAME hcb
GROUP P6/mmm
CELL 1.73205 1.73205 1.00000 90.0000 90.0000 120.0000
NODE 1 3 0.33333 0.66667 0.00000
EDGE 0.33333 0.66667 0.00000 0.66667 0.33333 0.00000
END
"""

_WORKSPACE_PCU = """\
CRYSTAL
NAME pcu
GROUP Pm-3m
CELL 1.00000 1.00000 1.00000 90.0000 90.0000 90.0000
NODE 1 6 0.00000 0.00000 0.00000
EDGE 0.00000 0.00000 0.00000 0.00000 0.00000 1.00000
END
"""


class TopologyRepositoryTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.root = Path(self.temp_dir.name)
        self.bundle_2d = self.root / "topo_2d.zip"
        self.bundle_3d = self.root / "topologies_3d_bundle.zip"
        with zipfile.ZipFile(self.bundle_2d, "w") as archive:
            archive.writestr("nets2d.cgd", _BUNDLE_2D)
        with zipfile.ZipFile(self.bundle_3d, "w") as archive:
            archive.writestr("nets3d.cgd", _BUNDLE_3D)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_importer_builds_index_with_connectivities(self):
        importer = RCSRArchiveImporter()

        entries = importer.build_index(self.bundle_2d, dimensionality="2D")

        self.assertEqual([entry.id for entry in entries], ["hcb", "kgm"])
        hcb = entries[0]
        self.assertEqual(hcb.n_node_definitions, 2)
        self.assertEqual(hcb.n_edge_definitions, 3)
        self.assertEqual(hcb.node_connectivities, (3, 3))
        self.assertEqual(hcb.space_group, "p6mm")

    def test_importer_handles_real_2d_cgd_coordinate_triplets(self):
        importer = RCSRArchiveImporter()

        with zipfile.ZipFile(self.bundle_2d, "w") as archive:
            archive.writestr(
                "real_2d.cgd",
                """\
CRYSTAL
NAME hcb
GROUP P6/mmm
CELL 1.73205 1.73205 1.00000 90.0000 90.0000 120.0000
NODE 1 3 0.33333 0.66667 0.00000
EDGE 0.33333 0.66667 0.00000 0.66667 0.33333 0.00000
END
""",
            )

        definitions = list(importer.iter_definitions(self.bundle_2d, dimensionality="2D"))
        self.assertEqual(definitions[0].node_definitions[0].position, (0.33333, 0.66667))
        self.assertEqual(definitions[0].edge_definitions[0].start, (0.33333, 0.66667))
        self.assertEqual(definitions[0].edge_definitions[0].end, (0.66667, 0.33333))

    def test_repository_lists_and_loads_real_topologies(self):
        repository = TopologyRepository.from_rcsr_archives(
            archive_paths={"2D": self.bundle_2d, "3D": self.bundle_3d}
        )

        index_entries = repository.list_index(dimensionality="2D", node_connectivities=(3, 6))
        self.assertEqual([entry.id for entry in index_entries], ["kgm"])
        self.assertTrue(index_entries[0].source_archive.endswith(".zip"))
        self.assertEqual(index_entries[0].n_edge_definitions, 6)

        kgm = repository.load("kgm")
        self.assertEqual(kgm.name, "kgm")
        self.assertEqual([node.connectivity for node in kgm.node_definitions], [3, 6])
        self.assertEqual(len(kgm.edge_definitions), 6)

        dia = repository.get_index_entry("dia")
        self.assertEqual(dia.dimensionality, "3D")
        self.assertEqual(dia.node_connectivities, (4, 4))

    def test_planner_prefers_repository_metadata_and_falls_back_when_empty(self):
        repository = TopologyRepository.from_rcsr_archives(archive_paths={"2D": self.bundle_2d})
        planner = NetPlanner(topology_repository=repository)

        tri = MonomerSpec(
            id="tri",
            name="tritopic",
            motifs=tuple(
                ReactiveMotif(id=f"m{i}", kind="amine", atom_ids=(i,), frame=Frame.xy()) for i in range(1, 4)
            ),
        )
        hexa = MonomerSpec(
            id="hexa",
            name="hexatopic",
            motifs=tuple(
                ReactiveMotif(id=f"h{i}", kind="aldehyde", atom_ids=(i,), frame=Frame.xy()) for i in range(1, 7)
            ),
        )

        plans = planner.propose((tri, hexa), tuple(), "2D")
        self.assertEqual(plans[0].topology.id, "kgm")
        self.assertEqual(plans[0].topology.metadata["n_edge_definitions"], 6)

        fallback_planner = NetPlanner(topology_repository=TopologyRepository())
        fallback = fallback_planner._infer_repository_topologies((tri, tri), "2D")
        self.assertEqual([hint.id for hint in fallback], ["hcb"])

    def test_default_repository_prefers_workspace_imported_data(self):
        topologies_module._DEFAULT_REPOSITORY = None
        try:
            repository = default_topology_repository()
            hcb = repository.get_index_entry("hcb")
            definition = load_topology("hcb")
            metadata = imported_topology_data_metadata()
        finally:
            topologies_module._DEFAULT_REPOSITORY = None

        self.assertEqual(hcb.dimensionality, "2D")
        self.assertEqual(hcb.node_connectivities, (3,))
        self.assertEqual(str(definition.metadata["source_archive"]), "bundled-topology://2d/hcb.cgd")
        self.assertEqual(metadata["import_summary"]["2D"]["indexed_topologies"], 194)
        self.assertGreater(metadata["import_summary"]["3D"]["indexed_topologies"], 0)
        self.assertTrue(str(metadata["archives_found"]["3D"]).endswith(".zip"))
        self.assertFalse(str(metadata["archives_found"]["3D"]).startswith("/"))

    def test_bundled_index_precomputes_two_d_symmetry_scans_for_bulk_listing(self):
        topologies_module._DEFAULT_REPOSITORY = None
        try:
            repository = default_topology_repository()
            entries = {entry.id: entry for entry in repository.list_index(dimensionality="2D")}
        finally:
            topologies_module._DEFAULT_REPOSITORY = None

        self.assertEqual(entries["car"].metadata["zero_linker_scan_method"], "gemmi-space-group-expanded-p1-scan")
        self.assertIsNotNone(entries["car"].metadata["zero_linker_compatible"])
        self.assertEqual(entries["car"].metadata["two_monomer_node_node_modes"], [])
        self.assertEqual(entries["car"].metadata["two_monomer_node_linker_modes"], ["3+2"])

    def test_discover_rcsr_archives_picks_nearest_matching_3d_bundle(self):
        from cofkit.topology_importers import discover_rcsr_archives

        topologies_dir = self.root / "downloads"
        topologies_dir.mkdir()
        (topologies_dir / "topo_2d.zip").write_bytes(b"")
        (topologies_dir / "rcsr_topology_3d_bundle.zip").write_bytes(b"")
        (topologies_dir / "some_other_3d.zip").write_bytes(b"")

        discovered = discover_rcsr_archives(topologies_dir)

        self.assertEqual(Path(discovered["2D"]).name, "topo_2d.zip")
        self.assertEqual(Path(discovered["3D"]).name, "rcsr_topology_3d_bundle.zip")

    def test_workspace_index_writer_precomputes_zero_linker_metadata_for_bulk_listing(self):
        workspace_root = self.root / "workspace_data"
        (workspace_root / "2d").mkdir(parents=True)
        (workspace_root / "3d").mkdir(parents=True)
        (workspace_root / "2d" / "hcb.cgd").write_text(_WORKSPACE_HCB, encoding="utf-8")
        (workspace_root / "3d" / "pcu.cgd").write_text(_WORKSPACE_PCU, encoding="utf-8")

        write_workspace_topology_index(data_dir=workspace_root)
        repository = TopologyRepository.from_workspace_data(data_dir=workspace_root)
        entries = {entry.id: entry for entry in repository.list_index()}

        self.assertTrue(entries["hcb"].metadata["zero_linker_compatible"])
        self.assertTrue(entries["hcb"].metadata["zero_linker_bipartite"])
        self.assertTrue(entries["hcb"].metadata["zero_linker_builder_supported"])
        self.assertEqual(entries["hcb"].metadata["two_monomer_node_node_modes"], ["3+3"])
        self.assertEqual(entries["hcb"].metadata["two_monomer_node_linker_modes"], ["3+2"])
        self.assertTrue(entries["pcu"].metadata["zero_linker_compatible"])
        self.assertTrue(entries["pcu"].metadata["zero_linker_bipartite"])
        self.assertFalse(entries["pcu"].metadata["zero_linker_builder_supported"])
        self.assertEqual(entries["pcu"].metadata["two_monomer_node_node_modes"], ["6+6"])
        self.assertEqual(entries["pcu"].metadata["two_monomer_node_linker_modes"], ["6+2"])


if __name__ == "__main__":
    unittest.main()
