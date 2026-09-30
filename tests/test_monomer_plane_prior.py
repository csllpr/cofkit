"""Regression tests for impact-review claim T1-8 (action A06).

The ditopic/monotopic plane prior used to be the world z axis of whatever
orientation the RDKit embedding happened to produce
(``chem/rdkit.py::_plane_normal`` returned ``(0, 0, 1)`` for fewer than three
motif anchors). The plane is now fitted from the conformer's own heavy-atom
covariance; genuinely non-planar, collinear, and degenerate conformers carry
no plane prior (zero motif frame normals, ``plane_normal=None`` metadata)
instead of borrowing physical meaning from a world axis, and the scoring /
optimizer consumers skip the plane-dependent terms for such motifs.

The key invariance property: rigidly rotating (and translating) the same
embedded conformer must leave the canonical monomer geometry — and therefore
chemical scores and assembled geometry — unchanged within numerical
tolerance. This is tested by rotating conformer coordinates, never by
comparing different RDKit seeds.
"""

import math
import unittest
from unittest.mock import patch

from cofkit import (
    AssignmentOutcome,
    AssignmentPlan,
    AssemblyState,
    CandidateScorer,
    ContinuousOptimizer,
    MonomerInstance,
    MotifRef,
    NetPlan,
    Pose,
    ReactionEvent,
    ReactionLibrary,
)
from cofkit.chem import rdkit as rdkit_module
from cofkit.chem.rdkit import build_rdkit_monomer
from cofkit.geometry import add, dot, matmul_vec, norm, sub
from cofkit.linkage_geometry import effective_motif_origin

try:
    from rdkit import Chem  # noqa: F401
except ImportError:  # pragma: no cover - environment-dependent
    Chem = None


def _rotation_about_axis(axis, angle):
    x, y, z = axis
    length = math.sqrt(x * x + y * y + z * z)
    x, y, z = x / length, y / length, z / length
    c = math.cos(angle)
    s = math.sin(angle)
    one_minus_c = 1.0 - c
    return (
        (c + x * x * one_minus_c, x * y * one_minus_c - z * s, x * z * one_minus_c + y * s),
        (y * x * one_minus_c + z * s, c + y * y * one_minus_c, y * z * one_minus_c - x * s),
        (z * x * one_minus_c - y * s, z * y * one_minus_c + x * s, c + z * z * one_minus_c),
    )


# A fixed proper rotation plus translation used for every invariance test.
_ROTATION = _rotation_about_axis((1.0, 2.0, 3.0), math.radians(37.0))
_TRANSLATION = (3.0, -2.0, 5.0)
_IDENTITY = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))


def _build_with_rigid_transform(monomer_id, smiles, motif_kind, rotation, translation):
    """Build a monomer from a conformer rigidly transformed after embedding.

    ``_embed_conformers`` is wrapped so the identical embedded conformer is
    rotated/translated before geometry construction; force-field minimization
    is stubbed out in both the reference and the transformed build so the two
    conformers differ by exactly the rigid transform, isolating the frame
    canonicalization under test.
    """
    original_embed = rdkit_module._embed_conformers

    def embed_then_transform(molecule, *, num_conformers, random_seed):
        result = original_embed(molecule, num_conformers=num_conformers, random_seed=random_seed)
        for conf_id in result.conformer_ids:
            conformer = molecule.GetConformer(conf_id)
            for atom_index in range(molecule.GetNumAtoms()):
                position = conformer.GetAtomPosition(atom_index)
                rotated = matmul_vec(rotation, (float(position.x), float(position.y), float(position.z)))
                conformer.SetAtomPosition(atom_index, add(rotated, translation))
        return result

    def skip_optimization(_molecule, conformer_ids, **_kwargs):
        return rdkit_module._ConformerSelectionResult(conformer_ids[0], None, "none", "skipped")

    with (
        patch.object(rdkit_module, "_embed_conformers", side_effect=embed_then_transform),
        patch.object(rdkit_module, "_optimize_conformers", side_effect=skip_optimization),
    ):
        return build_rdkit_monomer(monomer_id, monomer_id, smiles, motif_kind, num_conformers=1)


def _assert_vec_close(test_case, actual, expected, delta=1e-8):
    test_case.assertEqual(len(actual), len(expected))
    for a, e in zip(actual, expected):
        test_case.assertAlmostEqual(a, e, delta=delta)


def _assert_geometry_equivalent(test_case, rotated, reference, rotation):
    """The rotated-conformer build must reproduce the reference local geometry."""
    test_case.assertEqual(rotated.atom_symbols, reference.atom_symbols)
    test_case.assertEqual(rotated.bonds, reference.bonds)
    test_case.assertEqual(len(rotated.atom_positions), len(reference.atom_positions))
    for rotated_point, reference_point in zip(rotated.atom_positions, reference.atom_positions):
        _assert_vec_close(test_case, rotated_point, reference_point)
    test_case.assertEqual(len(rotated.motifs), len(reference.motifs))
    for rotated_motif, reference_motif in zip(rotated.motifs, reference.motifs):
        test_case.assertEqual(rotated_motif.id, reference_motif.id)
        test_case.assertEqual(rotated_motif.atom_ids, reference_motif.atom_ids)
        _assert_vec_close(test_case, rotated_motif.frame.origin, reference_motif.frame.origin)
        _assert_vec_close(test_case, rotated_motif.frame.primary, reference_motif.frame.primary)
        _assert_vec_close(test_case, rotated_motif.frame.normal, reference_motif.frame.normal)
    test_case.assertEqual(rotated.metadata["plane_status"], reference.metadata["plane_status"])
    reference_normal = reference.metadata["plane_normal"]
    rotated_normal = rotated.metadata["plane_normal"]
    if reference_normal is None:
        test_case.assertIsNone(rotated_normal)
    else:
        # The world-frame plane normal rotates covariantly with the conformer.
        test_case.assertAlmostEqual(dot(matmul_vec(rotation, reference_normal), rotated_normal), 1.0, delta=1e-8)


def _single_bridge_case(amine_spec, aldehyde_spec):
    """One imine bridge event between the first motifs of two real specs."""
    template = ReactionLibrary.builtin().get("imine_bridge")
    amine_motif = amine_spec.motifs[0]
    aldehyde_motif = aldehyde_spec.motifs[0]
    assignment_plan = AssignmentPlan(
        net_plan=NetPlan(topology=None, monomer_ids=(amine_spec.id, aldehyde_spec.id), reaction_ids=(template.id,)),
        slot_to_monomer={"slot1": amine_spec.id, "slot2": aldehyde_spec.id},
    )
    outcome = AssignmentOutcome(
        assignment_plan=assignment_plan,
        monomer_instances=(
            MonomerInstance(id="m1", monomer_id=amine_spec.id),
            MonomerInstance(id="m2", monomer_id=aldehyde_spec.id),
        ),
        events=(
            ReactionEvent(
                id="rxn1",
                template_id=template.id,
                participants=(
                    MotifRef(monomer_instance_id="m1", monomer_id=amine_spec.id, motif_id=amine_motif.id),
                    MotifRef(monomer_instance_id="m2", monomer_id=aldehyde_spec.id, motif_id=aldehyde_motif.id),
                ),
            ),
        ),
        unreacted_motifs=(),
        consumed_count=2,
    )
    state = AssemblyState(
        cell=((20.0, 0.0, 0.0), (0.0, 20.0, 0.0), (0.0, 0.0, 8.0)),
        monomer_poses={
            "m1": Pose(),
            "m2": Pose(
                translation=(4.0, 0.3, 0.2),
                rotation_matrix=_rotation_about_axis((0.0, 0.0, 1.0), math.pi),
            ),
        },
        stacking_state="disabled",
    )
    specs = {amine_spec.id: amine_spec, aldehyde_spec.id: aldehyde_spec}
    templates = {template.id: template}
    return specs, templates, outcome, state


def _bridge_report_args(case):
    """Reorder a ``_single_bridge_case`` tuple for ``bridge_geometry_report``."""
    specs, templates, outcome, state = case
    return outcome, state, specs, templates


@unittest.skipIf(Chem is None, "RDKit is not available")
class MolecularPlaneFitTests(unittest.TestCase):
    """Degeneracy handling of the heavy-atom plane fit itself."""

    def test_planar_point_set_yields_plane(self):
        points = ((-1.0, -1.0, 0.0), (1.0, -1.0, 0.0), (1.0, 1.0, 0.0), (-1.0, 1.0, 0.01))
        fit = rdkit_module._fit_molecular_plane(points)
        self.assertEqual(fit.status, "planar")
        self.assertIsNotNone(fit.normal)
        self.assertAlmostEqual(abs(fit.normal[2]), 1.0, delta=1e-3)
        self.assertLess(fit.rms_deviation, rdkit_module._MONOMER_PLANARITY_RMS_TOLERANCE_ANGSTROM)

    def test_nonplanar_point_set_is_uncertain(self):
        tetrahedron = ((1.0, 1.0, 1.0), (1.0, -1.0, -1.0), (-1.0, 1.0, -1.0), (-1.0, -1.0, 1.0))
        fit = rdkit_module._fit_molecular_plane(tetrahedron)
        self.assertEqual(fit.status, "nonplanar")
        self.assertIsNone(fit.normal)
        self.assertIsNotNone(fit.canonical_axis)

    def test_collinear_point_set_is_uncertain(self):
        line = ((0.0, 0.0, 0.0), (1.4, 0.0, 0.0), (2.8, 0.0, 0.0), (4.2, 0.0, 0.0))
        fit = rdkit_module._fit_molecular_plane(line)
        self.assertEqual(fit.status, "collinear")
        self.assertIsNone(fit.normal)

    def test_degenerate_point_sets_are_uncertain(self):
        for points in ((), ((1.0, 2.0, 3.0),), ((1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0))):
            fit = rdkit_module._fit_molecular_plane(points)
            self.assertEqual(fit.status, "degenerate")
            self.assertIsNone(fit.normal)
            self.assertIsNone(fit.canonical_axis)

    def test_plane_normal_rotates_covariantly(self):
        points = ((0.0, 0.0, 0.0), (2.1, 0.0, 0.0), (0.0, 1.3, 0.0), (2.1, 1.3, 0.01), (0.7, 0.4, -0.01))
        reference = rdkit_module._fit_molecular_plane(points)
        center = rdkit_module.centroid(points)
        rotated_points = tuple(
            add(matmul_vec(_ROTATION, rdkit_module.sub(point, center)), _TRANSLATION)
            for point in points
        )
        rotated = rdkit_module._fit_molecular_plane(rotated_points)
        self.assertEqual(rotated.status, reference.status)
        self.assertAlmostEqual(
            dot(matmul_vec(_ROTATION, reference.normal), rotated.normal),
            1.0,
            delta=1e-8,
        )
        self.assertAlmostEqual(rotated.rms_deviation, reference.rms_deviation, delta=1e-8)


@unittest.skipIf(Chem is None, "RDKit is not available")
class MonomerPlanePriorTests(unittest.TestCase):
    """The plane prior is molecular geometry, never the embedding's world z."""

    def test_planar_ditopic_monomer_gets_molecular_plane(self):
        monomer = build_rdkit_monomer("ppd", "ppd", "Nc1ccc(N)cc1", "amine", num_conformers=2)
        self.assertEqual(monomer.metadata["plane_status"], "planar")
        self.assertIsNotNone(monomer.metadata["plane_normal"])
        self.assertLessEqual(
            monomer.metadata["plane_rms_deviation"],
            rdkit_module._MONOMER_PLANARITY_RMS_TOLERANCE_ANGSTROM,
        )
        for motif in monomer.motifs:
            self.assertEqual(motif.frame.normal, (0.0, 0.0, 1.0))
        # The local frame really is the molecular plane: every heavy atom
        # lies near z=0 and the motif origins lie in the plane.
        heavy_z = [
            abs(position[2])
            for symbol, position in zip(monomer.atom_symbols, monomer.atom_positions)
            if symbol != "H"
        ]
        self.assertLess(max(heavy_z), 0.5)

    def test_planar_motif_axis_ditopic_monomer_gets_molecular_plane(self):
        # Terephthalonitrile's motif reactive/anchor points are collinear, so
        # an anchor-based plane cannot exist — but the molecule is planar and
        # the heavy-atom fit must find that plane.
        monomer = build_rdkit_monomer("tpn", "tpn", "N#Cc1ccc(C#N)cc1", "nitrile", num_conformers=2)
        self.assertEqual(monomer.metadata["plane_status"], "planar")
        for motif in monomer.motifs:
            self.assertEqual(motif.frame.normal, (0.0, 0.0, 1.0))

    def test_nonplanar_ditopic_monomer_is_uncertain(self):
        monomer = build_rdkit_monomer("dach", "dach", "NC1CCC(N)CC1", "amine", num_conformers=2)
        self.assertEqual(monomer.metadata["plane_status"], "nonplanar")
        self.assertIsNone(monomer.metadata["plane_normal"])
        for motif in monomer.motifs:
            self.assertEqual(motif.frame.normal, (0.0, 0.0, 0.0))

    def test_collinear_ditopic_monomer_is_uncertain(self):
        monomer = build_rdkit_monomer("dca", "dca", "N#CC#CC#N", "nitrile", num_conformers=2)
        self.assertEqual(monomer.metadata["plane_status"], "collinear")
        self.assertIsNone(monomer.metadata["plane_normal"])
        for motif in monomer.motifs:
            self.assertEqual(motif.frame.normal, (0.0, 0.0, 0.0))

    def test_planar_ditopic_geometry_invariant_under_rigid_rotation(self):
        reference = _build_with_rigid_transform("ppd_ref", "Nc1ccc(N)cc1", "amine", _IDENTITY, (0.0, 0.0, 0.0))
        rotated = _build_with_rigid_transform("ppd_rot", "Nc1ccc(N)cc1", "amine", _ROTATION, _TRANSLATION)
        _assert_geometry_equivalent(self, rotated, reference, _ROTATION)

    def test_nonplanar_ditopic_geometry_invariant_under_rigid_rotation(self):
        reference = _build_with_rigid_transform("dach_ref", "NC1CCC(N)CC1", "amine", _IDENTITY, (0.0, 0.0, 0.0))
        rotated = _build_with_rigid_transform("dach_rot", "NC1CCC(N)CC1", "amine", _ROTATION, _TRANSLATION)
        _assert_geometry_equivalent(self, rotated, reference, _ROTATION)

    def test_planar_motif_axis_ditopic_geometry_invariant_under_rigid_rotation(self):
        reference = _build_with_rigid_transform("tpn_ref", "N#Cc1ccc(C#N)cc1", "nitrile", _IDENTITY, (0.0, 0.0, 0.0))
        rotated = _build_with_rigid_transform("tpn_rot", "N#Cc1ccc(C#N)cc1", "nitrile", _ROTATION, _TRANSLATION)
        _assert_geometry_equivalent(self, rotated, reference, _ROTATION)

    def test_bridge_scores_invariant_under_rigid_rotation(self):
        # Chemical scores must not depend on the world orientation the
        # embedding produced: identical local geometry gives identical
        # per-event residuals through the full scoring path.
        reference_amine = _build_with_rigid_transform("ppd_ref", "Nc1ccc(N)cc1", "amine", _IDENTITY, (0.0, 0.0, 0.0))
        reference_aldehyde = _build_with_rigid_transform("tal_ref", "O=Cc1ccc(C=O)cc1", "aldehyde", _IDENTITY, (0.0, 0.0, 0.0))
        rotated_amine = _build_with_rigid_transform("ppd_rot", "Nc1ccc(N)cc1", "amine", _ROTATION, _TRANSLATION)
        rotated_aldehyde = _build_with_rigid_transform("tal_rot", "O=Cc1ccc(C=O)cc1", "aldehyde", _ROTATION, _TRANSLATION)
        scorer = CandidateScorer()
        reference_report = scorer.bridge_geometry_report(
            *_bridge_report_args(_single_bridge_case(reference_amine, reference_aldehyde))
        )
        rotated_report = scorer.bridge_geometry_report(
            *_bridge_report_args(_single_bridge_case(rotated_amine, rotated_aldehyde))
        )
        self.assertAlmostEqual(rotated_report.total_residual, reference_report.total_residual, delta=1e-8)
        for rotated_metrics, reference_metrics in zip(rotated_report.event_metrics, reference_report.event_metrics):
            self.assertAlmostEqual(rotated_metrics.distance_residual, reference_metrics.distance_residual, delta=1e-8)
            self.assertAlmostEqual(rotated_metrics.planarity_residual, reference_metrics.planarity_residual, delta=1e-8)
            self.assertAlmostEqual(rotated_metrics.alignment_residual, reference_metrics.alignment_residual, delta=1e-8)
            self.assertAlmostEqual(
                rotated_metrics.normal_misalignment_residual,
                reference_metrics.normal_misalignment_residual,
                delta=1e-8,
            )
            self.assertEqual(rotated_metrics.plane_prior_coverage, reference_metrics.plane_prior_coverage)


@unittest.skipIf(Chem is None, "RDKit is not available")
class UncertainPlaneConsumerTests(unittest.TestCase):
    """Consumers skip plane-dependent terms for plane-less (uncertain) motifs."""

    @classmethod
    def setUpClass(cls):
        cls.nonplanar_amine = build_rdkit_monomer("dach", "dach", "NC1CCC(N)CC1", "amine", num_conformers=2)
        cls.planar_aldehyde = build_rdkit_monomer("tal", "tal", "O=Cc1ccc(C=O)cc1", "aldehyde", num_conformers=2)
        cls.collinear_nitrile = build_rdkit_monomer("dca", "dca", "N#CC#CC#N", "nitrile", num_conformers=2)

    def test_scoring_skips_plane_terms_for_uncertain_motif(self):
        specs, templates, outcome, state = _single_bridge_case(self.nonplanar_amine, self.planar_aldehyde)
        report = CandidateScorer().bridge_geometry_report(outcome, state, specs, templates)
        metrics = report.event_metrics[0]
        self.assertEqual(metrics.plane_prior_coverage, "second")
        self.assertIsNone(metrics.normal_alignment)
        self.assertEqual(metrics.normal_misalignment_residual, 0.0)
        # Only the aldehyde side contributes a planarity term.
        motif1 = self.nonplanar_amine.motifs[0]
        motif2 = self.planar_aldehyde.motifs[0]
        origin1 = effective_motif_origin("imine_bridge", self.nonplanar_amine, motif1)
        origin2 = add(
            matmul_vec(
                state.monomer_poses["m2"].rotation_matrix,
                effective_motif_origin("imine_bridge", self.planar_aldehyde, motif2),
            ),
            state.monomer_poses["m2"].translation,
        )
        unit_vector = tuple(component / norm(sub(origin2, origin1)) for component in sub(origin2, origin1))
        expected = abs(dot(unit_vector, matmul_vec(state.monomer_poses["m2"].rotation_matrix, (0.0, 0.0, 1.0))))
        self.assertAlmostEqual(metrics.planarity_residual, expected, delta=1e-9)

    def test_scoring_reports_no_plane_coverage_when_both_motifs_uncertain(self):
        specs, templates, outcome, state = _single_bridge_case(self.nonplanar_amine, self.collinear_nitrile)
        report = CandidateScorer().bridge_geometry_report(outcome, state, specs, templates)
        metrics = report.event_metrics[0]
        self.assertEqual(metrics.plane_prior_coverage, "none")
        self.assertIsNone(metrics.normal_alignment)
        self.assertEqual(metrics.planarity_residual, 0.0)
        self.assertEqual(metrics.normal_misalignment_residual, 0.0)

    def test_optimizer_applies_distance_only_correction_without_plane_priors(self):
        specs, templates, outcome, state = _single_bridge_case(self.nonplanar_amine, self.collinear_nitrile)
        optimizer = ContinuousOptimizer(scorer=CandidateScorer())
        refined = optimizer._refine_translations(state, outcome, specs, templates)

        # With no plane prior on either side the update is purely the
        # along-bridge distance correction — no fabricated-plane term.
        motif1 = specs["dach"].motifs[0]
        motif2 = specs["dca"].motifs[0]
        origin1 = effective_motif_origin("imine_bridge", specs["dach"], motif1)
        origin2 = add(
            matmul_vec(
                state.monomer_poses["m2"].rotation_matrix,
                effective_motif_origin("imine_bridge", specs["dca"], motif2),
            ),
            state.monomer_poses["m2"].translation,
        )
        delta = sub(origin2, origin1)
        separation = norm(delta)
        direction = tuple(component / separation for component in delta)
        target = optimizer.scorer._target_distance(templates["imine_bridge"])
        expected_step = 0.5 * (separation - target) * optimizer.config.translation_step
        update1 = refined.monomer_poses["m1"].translation
        for component, unit in zip(update1, direction):
            self.assertAlmostEqual(component, unit * expected_step, delta=1e-9)

        # The full greedy loop must also tolerate plane-less motifs.
        result = optimizer.optimize(outcome, state, specs, templates)
        self.assertIn("final_residual", result.metrics)


if __name__ == "__main__":
    unittest.main()
