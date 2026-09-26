import contextlib
import io
import unittest
from unittest.mock import patch


from cofkit import single_node_topologies as snt


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
