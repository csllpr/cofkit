import inspect
import math
import unittest

import cofkit.linkage_geometry as linkage_geometry
import cofkit.reaction_realization as reaction_realization
from cofkit.linkage_geometry import (
    derived_origin_retraction_fraction,
    required_bridge_span,
)
from cofkit.reactions import bridge_geometry_priors


class RequiredBridgeSpanTests(unittest.TestCase):
    def test_equal_angle_closed_form_matches_hand_derivation(self):
        # Hand derivation (tau = 60 deg): tan a = 1.3 sin60 / (2.9 + 0.65),
        # d = 2.9 cos a + 1.3 cos(a - 60 deg) = 3.7242448899072143.
        span = required_bridge_span(1.45, 1.3, 1.45, 120.0, 120.0)
        self.assertIsNotNone(span)
        assert span is not None
        self.assertAlmostEqual(span, 3.7242448899072143, places=12)

    def test_span_closes_the_chain_at_both_angle_priors(self):
        # Independently rebuild the chain from the solved direction and
        # verify lateral closure, opposite-side placement, and both interior
        # angles against an asymmetric-priors case.
        outer_a, bridge, outer_b = 1.45, 1.3, 1.40
        angle_a, angle_b = 118.0, 122.0
        span = required_bridge_span(outer_a, bridge, outer_b, angle_a, angle_b)
        self.assertIsNotNone(span)
        assert span is not None
        tau_a = math.radians(180.0 - angle_a)
        tau_b = math.radians(180.0 - angle_b)
        s_coeff = outer_a + bridge * math.cos(tau_a) + outer_b * math.cos(tau_a - tau_b)
        d_coeff = bridge * math.sin(tau_a) + outer_b * math.sin(tau_a - tau_b)
        alpha = math.atan2(d_coeff, s_coeff)
        x_point = (outer_a * math.cos(alpha), outer_a * math.sin(alpha))
        y_point = (
            x_point[0] + bridge * math.cos(alpha - tau_a),
            x_point[1] + bridge * math.sin(alpha - tau_a),
        )
        b_point = (
            y_point[0] + outer_b * math.cos(alpha - tau_a + tau_b),
            y_point[1] + outer_b * math.sin(alpha - tau_a + tau_b),
        )
        self.assertGreater(x_point[1], 0.0)
        self.assertLess(y_point[1], 0.0)
        self.assertAlmostEqual(b_point[0], span, places=9)
        self.assertAlmostEqual(b_point[1], 0.0, places=9)

        def angle_deg(p1, p2, p3):
            v1 = (p1[0] - p2[0], p1[1] - p2[1])
            v2 = (p3[0] - p2[0], p3[1] - p2[1])
            n1 = math.hypot(*v1)
            n2 = math.hypot(*v2)
            return math.degrees(math.acos(max(-1.0, min(1.0, (v1[0] * v2[0] + v1[1] * v2[1]) / (n1 * n2)))))

        self.assertAlmostEqual(angle_deg((0.0, 0.0), x_point, y_point), angle_a, places=7)
        self.assertAlmostEqual(angle_deg(x_point, y_point, b_point), angle_b, places=7)

    def test_no_zigzag_solution_returns_none(self):
        # 180-degree priors are the collinear limit: no zig-zag closure.
        self.assertIsNone(required_bridge_span(1.45, 1.3, 1.45, 180.0, 180.0))
        # Degenerate lengths have no quadrilateral at all.
        self.assertIsNone(required_bridge_span(0.0, 1.3, 1.45, 120.0, 120.0))


class DerivedRetractionFractionTests(unittest.TestCase):
    def test_imine_fraction_matches_hand_derivation(self):
        # f = (1.3 + 2*1.45 - 3.7242448899072143) / (2*1.45) = 0.16405...
        fraction = derived_origin_retraction_fraction("imine_bridge", 1.45)
        self.assertAlmostEqual(fraction, 0.1640534862388917, places=12)

    def test_azine_fraction_matches_hand_derivation(self):
        # r* = sqrt(1.45^2 + 1.3^2 + 1.45*1.3) = 2.38276...;
        # f = (1.3 + 1.45 - r*) / (2*1.45) = 0.12663...
        expected = (1.3 + 1.45 - math.sqrt(1.45**2 + 1.3**2 + 1.45 * 1.3)) / (2.0 * 1.45)
        fraction = derived_origin_retraction_fraction("azine_bridge", 1.45)
        self.assertAlmostEqual(fraction, expected, places=12)
        self.assertAlmostEqual(fraction, 0.12663775465904448, places=12)

    def test_templates_without_priors_retract_nothing(self):
        self.assertEqual(derived_origin_retraction_fraction("hydrazone_bridge", 1.45), 0.0)
        self.assertEqual(derived_origin_retraction_fraction("keto_enamine_bridge", 1.45), 0.0)
        self.assertEqual(derived_origin_retraction_fraction("vinylene_bridge", 1.45), 0.0)
        self.assertEqual(derived_origin_retraction_fraction(None, 1.45), 0.0)

    def test_fraction_is_clamped_and_length_aware(self):
        self.assertEqual(derived_origin_retraction_fraction("imine_bridge", 0.0), 0.0)
        # Shorter side lengths need a larger fraction (the span deficit is a
        # larger share of the naive span).
        shorter = derived_origin_retraction_fraction("imine_bridge", 1.0)
        longer = derived_origin_retraction_fraction("imine_bridge", 1.45)
        self.assertGreater(shorter, longer)
        self.assertLessEqual(shorter, 0.9)


class PriorsOwnershipTests(unittest.TestCase):
    def test_profile_carries_the_geometry_priors(self):
        imine = bridge_geometry_priors("imine_bridge")
        self.assertIsNotNone(imine)
        assert imine is not None
        self.assertEqual(imine.carbon_angle_deg, 120.0)
        self.assertEqual(imine.nitrogen_angle_deg, 120.0)
        self.assertIsNone(imine.nn_target_distance)
        azine = bridge_geometry_priors("azine_bridge")
        self.assertIsNotNone(azine)
        assert azine is not None
        self.assertEqual(azine.carbon_angle_deg, 120.0)
        self.assertEqual(azine.nitrogen_angle_deg, 120.0)
        self.assertEqual(azine.nn_target_distance, 1.408)

    def test_no_hand_tuned_retraction_constants_remain(self):
        # Grep gate: the retired magic constants and the old imine scan
        # objective must not survive anywhere in the two owning modules.
        geometry_source = inspect.getsource(linkage_geometry)
        self.assertNotIn("IMINE_EFFECTIVE_ORIGIN_RETRACTION_FRACTION", geometry_source)
        self.assertNotIn("AZINE_EFFECTIVE_ORIGIN_RETRACTION_FRACTION", geometry_source)
        origin_source = inspect.getsource(linkage_geometry.effective_motif_origin)
        self.assertIn("derived_origin_retraction_fraction", origin_source)

        imine_fit_source = inspect.getsource(
            reaction_realization.ReactionRealizer._fit_imine_bridge_positions
        )
        self.assertNotIn("range(721)", imine_fit_source)
        self.assertNotIn("0.12 *", imine_fit_source)
        azine_endpoint_source = inspect.getsource(
            reaction_realization.ReactionRealizer._fit_azine_endpoint_carbon_position
        )
        self.assertNotIn("range(721)", azine_endpoint_source)


if __name__ == "__main__":
    unittest.main()
