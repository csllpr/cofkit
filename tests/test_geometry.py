import contextlib
import io
import math
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace


from cofkit.geometry import (
    Frame,
    angle_degrees,
    classify_2d_cell,
    classify_2d_cell_parameters,
    frame_axes,
    layer_normal_axis,
    matmul_vec,
    measure_layer_z_span,
    norm,
    orthogonal_component,
    planar_arrangement_mismatch,
    rotation_from_frame_to_axes,
    safe_normalize,
    smallest_covariance_axis,
    sub,
)


class GeometryTests(unittest.TestCase):
    def test_rotation_from_frame_to_axes_maps_orthonormalized_axes(self):
        source_frame = Frame(origin=(0.0, 0.0, 0.0), primary=(1.0, 2.0, 0.1), normal=(0.2, -0.1, 1.0))
        target_frame = Frame(origin=(0.0, 0.0, 0.0), primary=(0.0, 1.0, 0.2), normal=(0.0, 0.0, 1.0))

        rotation = rotation_from_frame_to_axes(
            source_frame,
            target_frame.primary,
            target_frame.normal,
        )
        source_axes = frame_axes(source_frame)
        target_axes = frame_axes(target_frame)

        self.assertLess(norm(sub(matmul_vec(rotation, source_axes[0]), target_axes[0])), 1e-8)
        self.assertLess(norm(sub(matmul_vec(rotation, source_axes[1]), target_axes[1])), 1e-8)
        self.assertLess(norm(sub(matmul_vec(rotation, source_axes[2]), target_axes[2])), 1e-8)


def _pose(translation=(0.0, 0.0, 0.0)):
    return SimpleNamespace(
        translation=translation,
        rotation_matrix=((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
    )


def _monomer(atom_positions):
    return SimpleNamespace(atom_positions=atom_positions)


class LayerSpanMeasurementTests(unittest.TestCase):
    """Regression tests for audit A7 (stacking-plan W7 tests 1-2 in spirit)."""

    def test_pose_translation_contributes_to_span(self):
        # Two planar instances of the same monomer at different heights:
        # omitting pose.translation (the C1 root cause) would report 0.0.
        monomer = _monomer([(0.0, 0.0, 0.0), (1.0, 1.0, 0.0)])
        report = measure_layer_z_span(
            poses={"i1": _pose((0.0, 0.0, 0.0)), "i2": _pose((0.0, 0.0, 2.0))},
            monomer_specs={"m": monomer},
            instance_to_monomer={"i1": "m", "i2": "m"},
            axis=(0.0, 0.0, 1.0),
            axis_label="c_hat",
        )
        self.assertAlmostEqual(report.span, 2.0)
        self.assertEqual(report.mode, "precursor_coordinates")
        self.assertEqual(report.axis, "c_hat")
        self.assertEqual(report.n_atoms, 4)
        self.assertTrue(report.z_translation_included)

    def test_tilted_c_cell_measures_along_layer_normal(self):
        # Tilted c (≈2.9° off the ab normal): an in-plane periodic image move
        # changes the c_hat coordinate but not the layer-normal coordinate.
        a_vec = (10.0, 0.0, 0.0)
        b_vec = (0.0, 10.0, 0.0)
        c_vec = (0.5, 0.0, 10.0)
        axis, label = layer_normal_axis(a_vec, b_vec, fallback_axis=c_vec)
        self.assertEqual(label, "layer_normal")
        self.assertAlmostEqual(axis[2], 1.0)

        c_hat = tuple(v / norm(c_vec) for v in c_vec)
        # The z=3 atom is represented in the (1, 0, 0) in-plane image — an
        # equivalent periodic representation of the same layer.
        base_monomer = _monomer([(0.0, 0.0, 0.0)])
        top_monomer = _monomer([(0.0, 0.0, 3.0)])
        poses = {"i1": _pose((0.0, 0.0, 0.0)), "i2": _pose(a_vec)}
        specs = {"m1": base_monomer, "m2": top_monomer}
        i2m = {"i1": "m1", "i2": "m2"}

        normal_span = measure_layer_z_span(
            poses=poses, monomer_specs=specs, instance_to_monomer=i2m,
            axis=axis, axis_label=label,
        )
        c_hat_span = measure_layer_z_span(
            poses=poses, monomer_specs=specs, instance_to_monomer=i2m,
            axis=c_hat, axis_label="c_hat",
        )
        # The layer-normal span is the true 3.0 thickness, invariant to the
        # in-plane image choice...
        self.assertAlmostEqual(normal_span.span, 3.0)
        # ...while the c_hat span is inflated by the image move.
        self.assertGreater(c_hat_span.span, 3.4)

    def test_genuinely_planar_layer_measures_zero_span(self):
        # A measured 0.0 span (planar layer) must be distinguishable from an
        # unmeasured fallback.
        monomer = _monomer([(0.0, 0.0, 0.0), (1.5, 0.5, 0.0), (2.0, 1.0, 0.0)])
        report = measure_layer_z_span(
            poses={"i1": _pose()}, monomer_specs={"m": monomer},
            instance_to_monomer={"i1": "m"}, axis=(0.0, 0.0, 1.0), axis_label="c_hat",
        )
        self.assertEqual(report.span, 0.0)
        self.assertEqual(report.mode, "precursor_coordinates")
        self.assertGreater(report.n_atoms, 0)

    def test_empty_coverage_reports_unavailable_not_zero_measurement(self):
        report = measure_layer_z_span(
            poses={}, monomer_specs={}, instance_to_monomer={},
            axis=(0.0, 0.0, 1.0), axis_label="c_hat",
        )
        self.assertEqual(report.span, 0.0)
        self.assertEqual(report.mode, "unavailable")
        self.assertEqual(report.n_atoms, 0)

    def test_degenerate_in_plane_cell_warns_and_falls_back(self):
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            axis, label = layer_normal_axis((1.0, 0.0, 0.0), (2.0, 0.0, 0.0), fallback_axis=(0.0, 0.0, 6.8))
        self.assertIn("warning:", stderr.getvalue())
        self.assertEqual(label, "c_hat")
        self.assertAlmostEqual(axis[2], 1.0)


class CellClassificationTests(unittest.TestCase):
    """Regression tests for audit A8 (stacking-plan W7 test 4)."""

    def test_fitted_60_09_degree_cell_is_hexagonal_60deg(self):
        gamma = 60.09
        gamma_rad = gamma * 3.141592653589793 / 180.0
        cell = (
            (15.1, 0.0, 0.0),
            (15.1 * math.cos(gamma_rad), 15.1 * math.sin(gamma_rad), 0.0),
            (0.0, 0.0, 7.0),
        )
        self.assertEqual(classify_2d_cell(cell), ("hexagonal", "60deg"))

    def test_classify_2d_cell_parameters_adapter(self):
        self.assertEqual(
            classify_2d_cell_parameters((15.1, 15.1, 3.6, 90.0, 90.0, 120.0)),
            ("hexagonal", "120deg"),
        )
        self.assertEqual(classify_2d_cell_parameters((10.0, 10.0, 90.0)), ("square", ""))
        self.assertEqual(classify_2d_cell_parameters((10.0, 12.0, 90.0)), ("orthogonal", ""))
        self.assertEqual(classify_2d_cell_parameters((10.0, 12.0, 80.0)), ("oblique", ""))
        self.assertEqual(classify_2d_cell_parameters((10.0,)), ("unknown", ""))

    def test_cross_module_classification_consistency(self):
        # The same definition cell must classify identically through the
        # migrated single-node and indexed-topology paths and the vector
        # classifier itself.
        from cofkit.indexed_topology_layouts import _metric_family as indexed_metric_family
        from cofkit.single_node_topologies import _metric_family as single_node_metric_family

        parameters = (20.0, 20.0, 3.4, 90.0, 90.0, 120.0)
        definition = SimpleNamespace(dimensionality="2D", metadata={"cell_parameters": parameters})
        gamma_rad = math.radians(120.0)
        vectors = (
            (20.0, 0.0, 0.0),
            (20.0 * math.cos(gamma_rad), 20.0 * math.sin(gamma_rad), 0.0),
            (0.0, 0.0, 3.4),
        )
        kind_from_vectors, setting = classify_2d_cell(vectors)
        self.assertEqual(single_node_metric_family(parameters), kind_from_vectors)
        self.assertEqual(indexed_metric_family(definition), kind_from_vectors)
        self.assertEqual((kind_from_vectors, setting), ("hexagonal", "120deg"))


class SharedPrimitiveTests(unittest.TestCase):
    """Audit A10: one implementation, explicit fallback/degenerate policies."""

    def test_safe_normalize_returns_unit_vector(self):
        result = safe_normalize((0.0, 3.0, 4.0))
        self.assertAlmostEqual(result[0], 0.0)
        self.assertAlmostEqual(result[1], 0.6)
        self.assertAlmostEqual(result[2], 0.8)

    def test_safe_normalize_degenerate_returns_fallback_verbatim(self):
        self.assertEqual(safe_normalize((0.0, 0.0, 0.0)), (0.0, 0.0, 1.0))
        self.assertEqual(safe_normalize((0.0, 0.0, 0.0), fallback=(1.0, 0.0, 0.0)), (1.0, 0.0, 0.0))

    def test_safe_normalize_warns_only_when_context_given(self):
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            safe_normalize((0.0, 0.0, 0.0))
        self.assertEqual(stderr.getvalue(), "")
        with contextlib.redirect_stderr(stderr):
            safe_normalize((0.0, 0.0, 0.0), warn_context="test context")
        self.assertIn("warning: test context", stderr.getvalue())

    def test_orthogonal_component_matches_for_unit_and_nonunit_axis(self):
        vector = (1.0, 2.0, 3.0)
        axis = (0.0, 0.0, 5.0)
        expected = orthogonal_component(vector, axis, axis_is_unit=False)
        self.assertEqual(expected, orthogonal_component(vector, (0.0, 0.0, 1.0), axis_is_unit=True))
        self.assertAlmostEqual(expected[2], 0.0)

    def test_orthogonal_component_degenerate_axis_returns_vector(self):
        self.assertEqual(orthogonal_component((1.0, 2.0, 3.0), (0.0, 0.0, 0.0)), (1.0, 2.0, 3.0))

    def test_angle_degrees_3d_and_2d(self):
        self.assertAlmostEqual(
            angle_degrees((1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 1.0, 0.0), on_degenerate=180.0),
            90.0,
        )
        self.assertAlmostEqual(
            angle_degrees((1.0, 0.0), (0.0, 0.0), (0.0, 1.0), on_degenerate=180.0),
            90.0,
        )

    def test_angle_degrees_degenerate_policy_is_caller_chosen(self):
        self.assertEqual(
            angle_degrees((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), on_degenerate=180.0),
            180.0,
        )
        self.assertIsNone(
            angle_degrees((0.0, 0.0), (0.0, 0.0), (1.0, 0.0), on_degenerate=None)
        )


class PlanarArrangementMismatchTests(unittest.TestCase):
    def test_regular_planar_polygon_has_zero_mismatch(self):
        points = tuple(
            (math.cos(2.0 * math.pi * index / 3), math.sin(2.0 * math.pi * index / 3), 0.0)
            for index in range(3)
        )
        mismatch = planar_arrangement_mismatch(points)
        self.assertIsNotNone(mismatch)
        self.assertLess(mismatch.angular_max_radians, 1e-6)
        self.assertLess(mismatch.radial_spread, 1e-6)
        self.assertLess(mismatch.planarity_rms, 1e-9)

    def test_skewed_arrangement_reports_angular_and_radial_terms(self):
        # One arm rotated off 120 degrees and shortened: the bug-report
        # monomer's lowest-energy conformer had exactly this failure mode.
        points = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (-0.45, -0.45, 0.0))
        mismatch = planar_arrangement_mismatch(points)
        self.assertIsNotNone(mismatch)
        self.assertGreater(mismatch.angular_max_radians, 0.05)
        self.assertGreater(mismatch.radial_spread, 0.05)
        self.assertLess(mismatch.planarity_rms, 1e-9)

    def test_out_of_plane_arrangement_reports_planarity_term(self):
        # Four points: three in-plane plus one lifted out of plane (any three
        # points are exactly coplanar, so planarity needs n >= 4).
        points = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.5), (-1.0, 0.0, 0.0), (0.0, -1.0, 0.5))
        mismatch = planar_arrangement_mismatch(points)
        self.assertIsNotNone(mismatch)
        self.assertGreater(mismatch.planarity_rms, 0.1)

    def test_degenerate_inputs_return_none(self):
        self.assertIsNone(planar_arrangement_mismatch(((1.0, 0.0, 0.0), (0.0, 1.0, 0.0))))
        self.assertIsNone(planar_arrangement_mismatch(((1.0, 0.0, 0.0),) * 3))

    def test_smallest_covariance_axis_finds_plane_normal(self):
        points = tuple(
            (math.cos(angle), math.sin(angle), 0.01 * index)
            for index, angle in enumerate((0.0, 2.1, 4.2))
        )
        axis = smallest_covariance_axis(points)
        self.assertIsNotNone(axis)
        self.assertAlmostEqual(abs(axis[2]), 1.0, places=3)


if __name__ == "__main__":
    unittest.main()
