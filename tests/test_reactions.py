"""Contract tests for ReactionLinkageProfile field semantics (action A22).

Pins the documented enforced-vs-informational split: informational fields
(``geometry_profile_id``, ``validation_profile_id``, ``library_layout``,
``require_distinct_participant_copies``) exist as annotations but must not
gate any ``supports_*`` outcome.
"""

from __future__ import annotations

import dataclasses
import unittest

from cofkit.reactions import (
    DEFAULT_BRIDGE_TARGET_DISTANCE,
    ReactionLibrary,
    ReactionTemplate,
    bridge_target_distance,
    supports_binary_bridge_pair_generation,
    topology_assignment_mode,
    workflow_family,
)


class LinkageProfileFieldContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.library = ReactionLibrary.builtin()

    def test_informational_fields_present_on_builtin_profiles(self):
        for template_id, profile in self.library.linkage_profiles.items():
            self.assertEqual(profile.geometry_profile_id, template_id)
            self.assertEqual(profile.validation_profile_id, template_id)
            self.assertEqual(profile.library_layout, "role_count_files")

    def test_distinct_participant_copies_is_annotation_only(self):
        # Ring profiles carry the annotation; binary-bridge profiles do not.
        for template_id in ("boroxine_trimerization", "triazine_trimerization"):
            self.assertTrue(self.library.linkage_profiles[template_id].require_distinct_participant_copies)
        for template_id, profile in self.library.linkage_profiles.items():
            if profile.workflow_family == "binary_bridge":
                self.assertFalse(profile.require_distinct_participant_copies, template_id)

    def test_informational_fields_do_not_gate_supports_flags(self):
        profile = self.library.linkage_profiles["boroxine_trimerization"]
        baseline = (
            profile.supports_binary_bridge_pair_generation,
            profile.supports_atomistic_realization,
            profile.supports_ring_forming_generation,
            profile.supports_topology_guided_generation,
        )
        # Flipping every informational field (including the
        # distinct-participant annotation) must change nothing: these fields
        # are documented as not enforced.
        rewritten = dataclasses.replace(
            profile,
            require_distinct_participant_copies=False,
            geometry_profile_id=None,
            validation_profile_id=None,
            library_layout="something_else",
        )
        self.assertEqual(
            (
                rewritten.supports_binary_bridge_pair_generation,
                rewritten.supports_atomistic_realization,
                rewritten.supports_ring_forming_generation,
                rewritten.supports_topology_guided_generation,
            ),
            baseline,
        )

    def test_ring_support_is_gated_by_coordination_not_the_annotation(self):
        profile = self.library.linkage_profiles["boroxine_trimerization"]
        # require_distinct_participant_copies stays True here: what actually
        # gates ring-forming support is the coordination/motif metadata.
        broken = dataclasses.replace(profile, ring_event_coordination=4)
        self.assertFalse(broken.supports_ring_forming_generation)
        missing_kind = dataclasses.replace(profile, ring_participant_motif_kind=None)
        self.assertFalse(missing_kind.supports_ring_forming_generation)

    def test_unprofiled_template_fallback_contract(self):
        ring_template = ReactionTemplate(
            id="custom_ring",
            arity=3,
            reactant_motif_kinds=("x", "x", "x"),
            product_name="custom",
            topology_role="ring",
        )
        self.assertEqual(bridge_target_distance(ring_template), 1.45)
        self.assertEqual(workflow_family(ring_template), "ring_forming")
        self.assertEqual(topology_assignment_mode(ring_template), "virtual_node_topology")
        self.assertFalse(supports_binary_bridge_pair_generation(ring_template))

        bridge_template = ReactionTemplate(
            id="custom_bridge",
            arity=2,
            reactant_motif_kinds=("x", "y"),
            product_name="custom",
        )
        self.assertEqual(bridge_target_distance(bridge_template), DEFAULT_BRIDGE_TARGET_DISTANCE)
        self.assertEqual(workflow_family(bridge_template), "template_driven")
        self.assertEqual(topology_assignment_mode(bridge_template), "topology_guided")


if __name__ == "__main__":
    unittest.main()
