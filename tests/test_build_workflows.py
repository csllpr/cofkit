import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit.build_workflows import builtin_build_workflow_registry
from cofkit.reactions import ReactionLibrary


class BuildWorkflowRegistryTests(unittest.TestCase):
    def setUp(self):
        self.registry = builtin_build_workflow_registry()
        self.reaction_library = ReactionLibrary.builtin()

    def test_builtin_registry_lists_isolated_workflow_families(self):
        workflows = {workflow.id: workflow for workflow in self.registry.list()}

        self.assertEqual(tuple(workflows), ("binary_bridge", "ring_forming", "composite_topology_unit"))
        self.assertEqual(workflows["binary_bridge"].implementation_status, "available")
        self.assertEqual(workflows["ring_forming"].implementation_status, "planned")
        self.assertTrue(workflows["composite_topology_unit"].explicit_selection_required)

    def test_binary_bridge_templates_resolve_to_binary_bridge_workflow(self):
        workflow = self.registry.resolve_template_workflow(
            self.reaction_library,
            ("imine_bridge",),
        )

        self.assertEqual(workflow.id, "binary_bridge")
        self.assertEqual(workflow.template_workflow_family, "binary_bridge")

    def test_ring_forming_templates_resolve_to_ring_forming_workflow(self):
        workflow = self.registry.resolve_template_workflow(
            self.reaction_library,
            ("boroxine_trimerization",),
        )

        self.assertEqual(workflow.id, "ring_forming")
        self.assertEqual(workflow.template_workflow_family, "ring_forming")

    def test_mixed_template_families_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "mixed families"):
            self.registry.resolve_template_workflow(
                self.reaction_library,
                ("imine_bridge", "boroxine_trimerization"),
            )


if __name__ == "__main__":
    unittest.main()
