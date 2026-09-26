import contextlib
import io
import unittest
from math import atan2, pi
from unittest.mock import patch


from cofkit import single_node_topologies as snt
from cofkit.topologies import load_topology


def _quotient_edge(
    edge_id: str,
    start: str,
    end: str,
    image: tuple[int, int, int],
) -> snt.ExpandedSingleNodeEdge:
    return snt.ExpandedSingleNodeEdge(
        id=edge_id,
        start_node_id=start,
        end_node_id=end,
        end_image=image,
        base_vector=(1.0, 0.0, 0.0),
        center_fractional=(0.0, 0.0, 0.0),
    )


class ExactBipartiteColoringTests(unittest.TestCase):
    def test_curated_topology_verdicts_unchanged(self):
        expected = {
            "hcb": True,
            "hca": False,
            "fes": True,
            "fxt": True,
            "sql": True,
            "kgm": False,
            "htb": False,
            "hxl": False,
        }
        for topology_id, verdict in expected.items():
            with self.subTest(topology_id=topology_id):
                self.assertIs(snt.expand_single_node_topology(topology_id).is_bipartite, verdict)

    def test_odd_cycle_closing_through_periodic_image_is_non_bipartite(self):
        # Quotient triangle A-B-C-A whose closing edge carries an even image
        # shift (2, 0). Every finite patch of the lifted net is consistently
        # 2-colorable (the retired finite-patch coloring returned bipartite
        # here), but no periodic 2-coloring exists: the zero-shift path forces
        # color(A) = color(C) while the closing edge forces them apart.
        edges = (
            _quotient_edge("e1", "A", "B", (0, 0, 0)),
            _quotient_edge("e2", "B", "C", (0, 0, 0)),
            _quotient_edge("e3", "C", "A", (2, 0, 0)),
        )

        is_bipartite, sublattices = snt._bipartite_coloring(edges, ("A", "B", "C"))

        self.assertFalse(is_bipartite)
        self.assertEqual(sublattices, {})

    def test_odd_parity_self_loop_is_a_bipartite_chain(self):
        # One node with a single edge to its own (1, 0) image lifts to a
        # bi-infinite chain, which is genuinely bipartite (the two endpoints
        # of an odd-parity self-loop live in different image-parity states).
        # sql/hxl-style curated nets rely on exactly this semantics.
        edges = (_quotient_edge("e1", "n1", "n1", (1, 0, 0)),)

        is_bipartite, sublattices = snt._bipartite_coloring(edges, ("n1",))

        self.assertTrue(is_bipartite)
        self.assertEqual(sublattices, {"n1": 0})

    def test_even_parity_self_loop_is_an_immediate_conflict(self):
        edges = (_quotient_edge("e1", "n1", "n1", (2, 0, 0)),)

        is_bipartite, sublattices = snt._bipartite_coloring(edges, ("n1",))

        self.assertFalse(is_bipartite)
        self.assertEqual(sublattices, {})

    def test_bipartite_multi_node_quotient_with_image_shifts(self):
        edges = (
            _quotient_edge("e1", "A", "B", (0, 0, 0)),
            _quotient_edge("e2", "A", "B", (1, 0, 0)),
            _quotient_edge("e3", "A", "B", (0, 1, 0)),
            _quotient_edge("e4", "A", "B", (1, 1, 0)),
        )

        is_bipartite, sublattices = snt._bipartite_coloring(edges, ("A", "B"))

        self.assertTrue(is_bipartite)
        self.assertEqual(sublattices, {"A": 0, "B": 1})


class MeasuredDirectionStarTests(unittest.TestCase):
    """hcb/hca/fes layout directions are measured from the space-group-expanded
    site (``_layout_directions_from_expanded_topology``), not fabricated by the
    retired rotate/mirror-trigonal heuristics. Equivalence audit 2026-09-26:
    the expanded path picks the asymmetric-node site itself (it sorts first in
    each curated expansion), so no monomer-orientation flip occurs; hcb and fes
    agree with the retired heuristics to the cgd coordinate precision, while
    hca's mirror heuristic was genuinely wrong — it fabricated the third arm at
    150° (mirroring across a metric axis that is not a site symmetry) where the
    measured neighbor direction is 60°."""

    def assertStarAlmostEqual(self, directions, expected_angles_degrees, places=6):
        angles = sorted(atan2(d[1], d[0]) * 180.0 / pi for d in directions)
        for actual, expected in zip(angles, sorted(expected_angles_degrees)):
            self.assertAlmostEqual(actual, expected, places=places)

    def test_hcb_star_measured_from_expanded_site(self):
        layout = snt.resolve_single_node_topology_layout("hcb")
        self.assertEqual(layout.inference_mode, "expanded")
        # Ideal star is (-150, -30, 90); the cgd node coordinates are rounded
        # to 1e-5, which shifts the measured angles by ~5e-4 degrees.
        self.assertStarAlmostEqual(layout.directions, (-150.0, -30.0, 90.0), places=3)

    def test_hca_star_measured_corrects_mirror_heuristic(self):
        layout = snt.resolve_single_node_topology_layout("hca")
        self.assertEqual(layout.inference_mode, "expanded")
        self.assertStarAlmostEqual(layout.directions, (-150.0, 0.0, 60.0))

    def test_fes_star_measured_from_expanded_site(self):
        layout = snt.resolve_single_node_topology_layout("fes")
        self.assertEqual(layout.inference_mode, "expanded")
        self.assertStarAlmostEqual(layout.directions, (-135.0, -45.0, 90.0))


class GeneralPlaneGroupExpansionTests(unittest.TestCase):
    """_point_operations_2d derives the layer point operations for arbitrary
    plane groups from gemmi's 3D space-group tables (operations preserving the
    z layer, projected to their in-plane fractional matrix and translation)."""

    def test_curated_topologies_expand_with_expected_counts(self):
        expected = {
            "hcb": (2, 3, True),
            "hca": (6, 9, False),
            "fes": (4, 6, True),
            "fxt": (12, 18, True),
            "sql": (1, 2, True),
            "kgm": (3, 6, False),
            "htb": (6, 12, False),
            "hxl": (1, 3, False),
        }
        for topology_id, (n_nodes, n_edges, is_bipartite) in expected.items():
            with self.subTest(topology_id=topology_id):
                expanded = snt.expand_single_node_topology(topology_id)
                self.assertEqual(len(expanded.node_sites), n_nodes)
                self.assertEqual(len(expanded.edge_sites), n_edges)
                self.assertIs(expanded.is_bipartite, is_bipartite)
                for node in expanded.node_sites:
                    self.assertEqual(len(node.directions), expanded.connectivity)

    def test_cmmm_topology_expands_with_centering(self):
        # cem is a real single-node RCSR 2D net in Cmmm; the C centering
        # doubles the mirror orbit of the asymmetric node (0, 0.13397).
        expanded = snt.expand_single_node_topology_definition(load_topology("cem"))

        positions = {
            (round(node.fractional_position[0], 5), round(node.fractional_position[1], 5))
            for node in expanded.node_sites
        }
        self.assertEqual(
            positions,
            {(0.0, 0.13397), (0.0, 0.86603), (0.5, 0.36603), (0.5, 0.63397)},
        )
        self.assertEqual(len(expanded.edge_sites), 10)
        for node in expanded.node_sites:
            self.assertEqual(len(node.directions), 5)

    def test_p4mbm_topology_expands(self):
        # tts is a real single-node RCSR 2D net in P4/mbm.
        expanded = snt.expand_single_node_topology_definition(load_topology("tts"))

        self.assertEqual(len(expanded.node_sites), 4)
        self.assertEqual(len(expanded.edge_sites), 10)
        for node in expanded.node_sites:
            self.assertEqual(len(node.directions), 5)

    def test_new_plane_group_topologies_list_without_omission_warning(self):
        snt._OMITTED_TOPOLOGY_WARNING_SEEN.clear()
        self.addCleanup(snt._OMITTED_TOPOLOGY_WARNING_SEEN.clear)
        entries = {
            "cem": {"inference_mode": "expanded"},
            "tts": {"inference_mode": "expanded"},
        }
        with patch.dict(snt._SUPPORTED_SINGLE_NODE_TOPOLOGIES, entries):
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                ids = snt.list_supported_single_node_topology_ids(5, pair_mode="node-linker")

        self.assertIn("cem", ids)
        self.assertIn("tts", ids)
        self.assertEqual(stderr.getvalue(), "")


class SupportedSingleNodeTopologyListingTests(unittest.TestCase):
    def setUp(self):
        snt._OMITTED_TOPOLOGY_WARNING_SEEN.clear()
        self.addCleanup(snt._OMITTED_TOPOLOGY_WARNING_SEEN.clear)

    def test_curated_topologies_list_without_warnings(self):
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            ids = snt.list_supported_single_node_topology_ids(3, pair_mode="node-node")

        self.assertIn("hcb", ids)
        self.assertEqual(stderr.getvalue(), "")

    def test_unsupported_space_group_omission_warns_once_per_topology(self):
        real_point_operations = snt._point_operations_2d

        def drop_p6mmm(space_group):
            if space_group == "P6/mmm":
                return ()
            return real_point_operations(space_group)

        snt.resolve_single_node_topology_layout.cache_clear()
        snt.expand_single_node_topology.cache_clear()
        self.addCleanup(snt.resolve_single_node_topology_layout.cache_clear)
        self.addCleanup(snt.expand_single_node_topology.cache_clear)
        with patch.object(snt, "_point_operations_2d", side_effect=drop_p6mmm):
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                ids = snt.list_supported_single_node_topology_ids(3, pair_mode="node-node")

            self.assertNotIn("hcb", ids)
            self.assertIn(
                "warning: single-node topology 'hcb' omitted from the supported list",
                stderr.getvalue(),
            )
            self.assertIn("space group", stderr.getvalue())

            repeat = io.StringIO()
            with contextlib.redirect_stderr(repeat):
                snt.list_supported_single_node_topology_ids(3, pair_mode="node-node")
                snt.list_supported_single_node_topology_ids(3, pair_mode="node-linker")

            self.assertEqual(repeat.getvalue(), "")

    def test_non_2d_registry_entry_warns_with_reason(self):
        with patch.dict(snt._SUPPORTED_SINGLE_NODE_TOPOLOGIES, {"dia": {"inference_mode": "expanded"}}):
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                ids = snt.list_supported_single_node_topology_ids(4, pair_mode="node-linker")

            self.assertNotIn("dia", ids)
            self.assertIn("warning: single-node topology 'dia' omitted from the supported list", stderr.getvalue())
            self.assertIn("must be 2D", stderr.getvalue())


if __name__ == "__main__":
    unittest.main()
