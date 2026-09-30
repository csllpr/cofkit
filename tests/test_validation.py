import contextlib
import io
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


from cofkit.validation import CoarseStructureValidator, CoarseValidationThresholds, classify_batch_output

try:
    import gemmi  # noqa: F401
except ImportError:  # pragma: no cover - environment-dependent
    gemmi = None


def _write_test_cif(path: Path, *, atoms: list[tuple[str, str, float, float, float]], bonds: list[tuple], bond_distance: float = 1.5) -> None:
    lines = [
        f"data_{path.stem}",
        "_audit_creation_method 'cofkit test'",
        "_space_group_name_H-M_alt 'P 1'",
        "_space_group_IT_number 1",
        "_cell_length_a 10.0",
        "_cell_length_b 10.0",
        "_cell_length_c 10.0",
        "_cell_angle_alpha 90.0",
        "_cell_angle_beta 90.0",
        "_cell_angle_gamma 90.0",
        "",
        "loop_",
        "_space_group_symop_operation_xyz",
        "'x,y,z'",
        "",
        "loop_",
        "_atom_site_label",
        "_atom_site_type_symbol",
        "_atom_site_fract_x",
        "_atom_site_fract_y",
        "_atom_site_fract_z",
        "_atom_site_occupancy",
    ]
    for label, symbol, fx, fy, fz in atoms:
        lines.append(f"{label} {symbol} {fx:.6f} {fy:.6f} {fz:.6f} 1.00")
    if bonds:
        lines.extend(
            [
                "",
                "loop_",
                "_geom_bond_atom_site_label_1",
                "_geom_bond_atom_site_label_2",
                "_geom_bond_site_symmetry_1",
                "_geom_bond_site_symmetry_2",
                "_geom_bond_distance",
            ]
        )
        for bond in bonds:
            left, right = bond[0], bond[1]
            sym1 = bond[2] if len(bond) > 2 else "."
            sym2 = bond[3] if len(bond) > 3 else "."
            distance = bond[4] if len(bond) > 4 else bond_distance
            lines.append(f"{left} {right} {sym1} {sym2} {distance:.6f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_boundary_bond_cif(path: Path) -> None:
    """Single-instance CIF with a 0.8 A bond crossing the x cell boundary
    (symmetry 1_455) plus a genuinely nonbonded pair at 2.9 A — above the
    0.75 vdW-ratio cutoff for C...C (2.55 A) and the 2.2 A severe-overlap
    floor, so it must not be flagged."""
    _write_test_cif(
        path,
        atoms=[
            ("mol_a_C1", "C", 0.05, 0.10, 0.10),
            ("mol_a_C2", "C", 0.97, 0.10, 0.10),
            ("mol_a_C3", "C", 0.40, 0.50, 0.50),
            ("mol_a_C4", "C", 0.40, 0.50, 0.79),
        ],
        bonds=[("mol_a_C1", "mol_a_C2", ".", "1_455", 0.8)],
    )


def _summary_record(
    cif_path: Path,
    *,
    structure_id: str,
    topology_id: str = "hcb",
    distance_residual: float = 0.1,
    actual_distance: float = 1.4,
    target_distance: float = 1.3,
    template_id: str | None = None,
    include_bridge_event: bool = True,
) -> dict[str, object]:
    bridge_event_metrics = (
        [
            {
                "distance_residual": distance_residual,
                "actual_distance": actual_distance,
                "target_distance": target_distance,
            }
        ]
        if include_bridge_event
        else []
    )
    return {
        "structure_id": structure_id,
        "pair_id": structure_id,
        "pair_mode": "3+2-node-linker",
        "status": "ok",
        "amine_record_id": "amine",
        "aldehyde_record_id": "aldehyde",
        "amine_connectivity": 3,
        "aldehyde_connectivity": 2,
        "topology_id": topology_id,
        "score": 1.0,
        "flags": [],
        "cif_path": str(cif_path),
        "metadata": {
            "graph_summary": {
                "n_monomer_instances": 2,
                "n_reaction_events": 1 if include_bridge_event else 0,
                "reaction_templates": {"imine_bridge": 1} if include_bridge_event else {},
            },
            **({"template_id": template_id} if template_id is not None else {}),
            "score_metadata": {
                "n_unreacted_motifs": 0,
                "bridge_event_metrics": bridge_event_metrics,
            },
        },
    }


@unittest.skipIf(gemmi is None, "gemmi is not available")
class CoarseValidationTests(unittest.TestCase):
    def test_validator_marks_bridge_distance_metadata_as_needs_optimization(self):
        # The verdict comes from the measured final coordinates: the realized
        # inter-monomer bond is 2.4 A long against the 1.3 A imine target.
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "valid.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("a1_C1", "C", 0.1, 0.1, 0.1),
                    ("b1_C1", "C", 0.34, 0.1, 0.1),
                ],
                bonds=[("a1_C1", "b1_C1")],
                bond_distance=2.4,
            )
            record = _summary_record(
                cif_path,
                structure_id="too_far",
                distance_residual=1.2,
                actual_distance=2.4,
            )

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "needs_optimization")
        self.assertFalse(report.passes_hard_validation)
        self.assertIn("bridge_distance_residual_max_hard", report.needs_optimization_reasons)
        self.assertIn("bridge_distance_too_long", report.needs_optimization_reasons)
        self.assertEqual(report.hard_invalid_reasons, ())
        self.assertEqual(report.metrics["bridge_metrics_source"], "final_cif")
        self.assertEqual(report.coverage["bridge_geometry"], "measured")

    def test_distorted_final_geometry_fails_despite_clean_seed_metrics(self):
        """A04 regression: a deliberately distorted repaired linkage must not
        pass because the seed/assembly metrics look clean. The verdict follows
        the measured final coordinates."""
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "distorted.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("a1_C1", "C", 0.1, 0.1, 0.1),
                    ("b1_C1", "C", 0.34, 0.1, 0.1),
                ],
                bonds=[("a1_C1", "b1_C1")],
                bond_distance=2.4,
            )
            record = _summary_record(
                cif_path,
                structure_id="distorted",
                distance_residual=0.0,
                actual_distance=1.3,
                target_distance=1.3,
            )

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "needs_optimization")
        self.assertFalse(report.is_valid)
        self.assertIn("bridge_distance_residual_max_hard", report.needs_optimization_reasons)
        # Seed metrics survive as clearly-labeled informational assembly metrics.
        self.assertEqual(report.metrics["seed_max_bridge_distance_residual"], 0.0)
        self.assertAlmostEqual(report.metrics["max_bridge_distance_residual"], 1.1, places=3)
        self.assertEqual(report.metrics["bridge_metrics_source"], "final_cif")

    def test_clean_seed_metrics_survive_when_final_geometry_matches(self):
        """Complementary control: undistorted export with clean seed metrics."""
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "clean.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("a1_C1", "C", 0.1, 0.1, 0.1),
                    ("b1_C1", "C", 0.23, 0.1, 0.1),
                ],
                bonds=[("a1_C1", "b1_C1")],
                bond_distance=1.3,
            )
            record = _summary_record(
                cif_path,
                structure_id="clean",
                distance_residual=0.0,
                actual_distance=1.3,
                target_distance=1.3,
            )

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "valid")
        self.assertEqual(report.coverage["bridge_geometry"], "measured")
        self.assertEqual(report.unmeasured_required_checks, ())

    def test_validator_marks_overlong_bridge_as_hard_hard_invalid(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "hard_hard.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("a1_C1", "C", 0.1, 0.1, 0.1),
                    ("b1_C1", "C", 0.36, 0.1, 0.1),
                ],
                bonds=[("a1_C1", "b1_C1")],
                bond_distance=2.6,
            )
            record = _summary_record(
                cif_path,
                structure_id="hard_hard",
                distance_residual=1.3,
                actual_distance=2.5,
            )

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "hard_hard_invalid")
        self.assertFalse(report.passes_hard_validation)
        self.assertTrue(report.blocks_cif_export)
        self.assertIn("bridge_distance_exceeds_cif_export_limit", report.hard_hard_invalid_reasons)

    def test_validator_marks_moderate_bridge_drift_as_warning(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "warning.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("a1_C1", "C", 0.1, 0.1, 0.1),
                    ("b1_C1", "C", 0.29, 0.1, 0.1),
                ],
                bonds=[("a1_C1", "b1_C1")],
                bond_distance=1.9,
            )
            record = _summary_record(
                cif_path,
                structure_id="warning",
                distance_residual=0.6,
                actual_distance=1.9,
            )

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "warning")
        self.assertTrue(report.passes_hard_validation)
        self.assertIn("bridge_distance_residual_mean", report.warning_reasons)
        self.assertIn("bridge_distance_residual_fraction", report.warning_reasons)

    def test_validator_surfaces_monomer_geometry_degradation_as_warning(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "degraded.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("a1_C1", "C", 0.1, 0.1, 0.1),
                    ("b1_C1", "C", 0.25, 0.1, 0.1),
                ],
                bonds=[("a1_C1", "b1_C1")],
            )
            record = _summary_record(cif_path, structure_id="degraded")
            record["metadata"]["monomer_geometry_warnings"] = [
                "monomer_geometry_degraded:tapb:conformer embedding fell back to a 2D planar depiction (rdkit-2d)",
                "monomer_geometry_degraded:tapb:conformer is not force-field minimized",
            ]

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "warning")
        self.assertTrue(report.passes_hard_validation)
        self.assertIn("monomer_geometry_degraded", report.warning_reasons)
        self.assertEqual(
            tuple(report.metrics["monomer_geometry_degraded_details"]),
            (
                "monomer_geometry_degraded:tapb:conformer embedding fell back to a 2D planar depiction (rdkit-2d)",
                "monomer_geometry_degraded:tapb:conformer is not force-field minimized",
            ),
        )

    def test_validator_without_monomer_geometry_warnings_stays_valid(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "clean.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("a1_C1", "C", 0.1, 0.1, 0.1),
                    ("b1_C1", "C", 0.25, 0.1, 0.1),
                ],
                bonds=[("a1_C1", "b1_C1")],
            )
            record = _summary_record(cif_path, structure_id="clean")

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "valid")
        self.assertNotIn("monomer_geometry_degraded", report.warning_reasons)
        self.assertNotIn("monomer_geometry_degraded_details", report.metrics)

    def test_validator_rejects_disconnected_instance_graph(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "disconnected.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("a1_C1", "C", 0.1, 0.1, 0.1),
                    ("b1_C1", "C", 0.7, 0.7, 0.7),
                ],
                bonds=[],
            )
            record = _summary_record(cif_path, structure_id="disconnected")

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "hard_invalid")
        self.assertIn("disconnected_instance_graph", report.hard_invalid_reasons)
        self.assertEqual(report.metrics["n_instance_components"], 2)

    def test_validator_allows_disconnected_bilayer_when_stacking_metadata_matches_components(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "stacked_bilayer.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("a1L0_C1", "C", 0.10, 0.10, 0.10),
                    ("b1L0_C1", "C", 0.25, 0.10, 0.10),
                    ("a1L1_C1", "C", 0.10, 0.10, 0.60),
                    ("b1L1_C1", "C", 0.25, 0.10, 0.60),
                ],
                bonds=[("a1L0_C1", "b1L0_C1"), ("a1L1_C1", "b1L1_C1")],
            )
            record = _summary_record(cif_path, structure_id="stacked_bilayer")
            record["metadata"]["graph_summary"] = {
                "n_monomer_instances": 4,
                "n_reaction_events": 2,
                "reaction_templates": {"imine_bridge": 2},
            }
            record["metadata"]["stacking"] = {"id": "AA", "layer_count": 2}

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "valid")
        self.assertNotIn("disconnected_instance_graph", report.hard_invalid_reasons)
        self.assertEqual(report.metrics["n_instance_components"], 2)
        self.assertEqual(tuple(report.metrics["instance_component_sizes"]), (2, 2))
        self.assertTrue(report.metrics["disconnected_instance_graph_allowed"])

    def test_validator_rejects_heavy_atom_clash(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "clash.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("x1_C1", "C", 0.1, 0.1, 0.1),
                    ("x1_C2", "C", 0.15, 0.1, 0.1),
                ],
                bonds=[],
            )
            record = _summary_record(cif_path, structure_id="clash")
            record["metadata"]["graph_summary"] = {"n_monomer_instances": 1, "n_reaction_events": 1}

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "hard_invalid")
        self.assertIn("heavy_atom_clash", report.hard_invalid_reasons)
        self.assertLess(report.metrics["min_nonbonded_heavy_distance"], 1.05)

    def test_validator_flags_two_angstrom_nonbonded_heavy_pair(self):
        """Plan W7 test 6: a 2.0 A heavy-heavy nonbonded pair must be flagged.

        2.0 A / (1.7 + 1.7) A = 0.588 < 0.75 (vdW ratio criterion) and below
        the 2.2 A severe-overlap floor; the old 1.05 A backstop missed it.
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "interpenetrated.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("x1_C1", "C", 0.40, 0.50, 0.50),
                    ("x1_C2", "C", 0.60, 0.50, 0.50),
                ],
                bonds=[],
            )
            record = _summary_record(cif_path, structure_id="interpenetrated")
            record["metadata"]["graph_summary"] = {"n_monomer_instances": 1, "n_reaction_events": 1}

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "hard_invalid")
        self.assertIn("heavy_atom_clash", report.hard_invalid_reasons)
        self.assertAlmostEqual(report.metrics["min_nonbonded_heavy_distance"], 2.0, places=3)
        self.assertAlmostEqual(report.metrics["min_nonbonded_heavy_vdw_ratio"], 2.0 / 3.4, places=3)
        self.assertEqual(report.metrics["n_heavy_atom_clash_pairs"], 1)
        detail = report.metrics["heavy_atom_clash_min"]
        self.assertEqual(tuple(detail["labels"]), ("x1_C1", "x1_C2"))
        self.assertEqual(tuple(detail["image"]), (0, 0, 0))

    def test_ideal_benzene_ring_is_not_flagged(self):
        """Negative control: ideal benzene (1.4 A bonds) has 1-3 C...C at 2.42 A,
        below the naive 0.75 vdW-ratio cutoff (2.55 A); with bond-graph 1-3/1-4
        exclusions nothing may be flagged."""
        cell = 10.0
        bond_cc = 1.397
        bond_ch = 1.09
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "benzene.cif"
            atoms = []
            bonds = []
            for k in range(6):
                theta = math.radians(60 * k)
                cx = 0.5 + bond_cc * math.cos(theta) / cell
                cy = 0.5 + bond_cc * math.sin(theta) / cell
                atoms.append((f"bz_C{k + 1}", "C", cx, cy, 0.5))
                hx = 0.5 + (bond_cc + bond_ch) * math.cos(theta) / cell
                hy = 0.5 + (bond_cc + bond_ch) * math.sin(theta) / cell
                atoms.append((f"bz_H{k + 1}", "H", hx, hy, 0.5))
            for k in range(6):
                bonds.append((f"bz_C{k + 1}", f"bz_C{(k + 1) % 6 + 1}", ".", ".", bond_cc))
                bonds.append((f"bz_C{k + 1}", f"bz_H{k + 1}", ".", ".", bond_ch))
            _write_test_cif(cif_path, atoms=atoms, bonds=bonds)
            record = _summary_record(cif_path, structure_id="benzene", include_bridge_event=False)

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "valid")
        self.assertNotIn("heavy_atom_clash", report.reasons)
        self.assertNotIn("hydrogen_atom_clash", report.reasons)
        self.assertNotIn("excluded_pair_severe_overlap", report.reasons)
        # Every heavy pair is bonded / 1-3 / 1-4: no nonbonded contact remains.
        self.assertIsNone(report.metrics["min_nonbonded_heavy_distance"])

    def test_ordinary_angle_neighbors_are_not_flagged(self):
        """Negative control: a C-C-C angle at 109.5 deg puts the 1-3 pair at
        2.51 A — below the naive 2.55 A ratio cutoff — yet must not flag."""
        bond = 1.54
        half_angle = math.radians(109.5) / 2.0
        ax, ay = 3.46, 5.0
        bx, by = 5.0, 5.0
        # BC makes a 109.5 deg angle with BA = (-1, 0)
        bc_dir = (math.cos(math.radians(70.5)), math.sin(math.radians(70.5)))
        cx = bx + bond * bc_dir[0]
        cy = by + bond * bc_dir[1]
        self.assertAlmostEqual(math.dist((ax, ay), (cx, cy)), 2.0 * bond * math.sin(half_angle), places=3)
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "angle_chain.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("ch_C1", "C", ax / 10.0, ay / 10.0, 0.5),
                    ("ch_C2", "C", bx / 10.0, by / 10.0, 0.5),
                    ("ch_C3", "C", cx / 10.0, cy / 10.0, 0.5),
                ],
                bonds=[("ch_C1", "ch_C2", ".", ".", bond), ("ch_C2", "ch_C3", ".", ".", bond)],
            )
            record = _summary_record(cif_path, structure_id="angle_chain", include_bridge_event=False)

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "valid")
        self.assertNotIn("heavy_atom_clash", report.reasons)
        self.assertIsNone(report.metrics["min_nonbonded_heavy_distance"])

    def _write_cis_chain_cif(self, path: Path, angle_deg: float) -> float:
        """Planar cis C-C-C-C chain (1.45 A bonds); returns the 1-4 distance."""
        bond = 1.45
        external = math.radians(180.0 - angle_deg)
        a = (4.0, 4.0)
        b = (a[0] + bond, a[1])
        c = (b[0] + bond * math.cos(external), b[1] + bond * math.sin(external))
        d = (c[0] + bond * math.cos(2.0 * external), c[1] + bond * math.sin(2.0 * external))
        _write_test_cif(
            path,
            atoms=[
                ("ch_C1", "C", a[0] / 10.0, a[1] / 10.0, 0.5),
                ("ch_C2", "C", b[0] / 10.0, b[1] / 10.0, 0.5),
                ("ch_C3", "C", c[0] / 10.0, c[1] / 10.0, 0.5),
                ("ch_C4", "C", d[0] / 10.0, d[1] / 10.0, 0.5),
            ],
            bonds=[
                ("ch_C1", "ch_C2", ".", ".", bond),
                ("ch_C2", "ch_C3", ".", ".", bond),
                ("ch_C3", "ch_C4", ".", ".", bond),
            ],
        )
        return math.dist(a, d)

    def test_cis_1_4_contact_below_ratio_cutoff_is_not_flagged(self):
        """1-4 policy: cis torsion contacts legitimately sit below
        0.75 * sum(r_vdw); excluded from the ratio check, above the 2.2 A
        severe-overlap floor -> not flagged."""
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "cis_chain.cif"
            d14 = self._write_cis_chain_cif(cif_path, 110.0)
            self.assertGreater(d14, 2.2)  # above the severe-overlap floor
            self.assertLess(d14, 0.75 * 3.4)  # below the naive ratio cutoff
            record = _summary_record(cif_path, structure_id="cis_chain", include_bridge_event=False)

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "valid")
        self.assertNotIn("heavy_atom_clash", report.reasons)
        self.assertNotIn("excluded_pair_severe_overlap", report.reasons)
        self.assertIsNone(report.metrics["min_nonbonded_heavy_distance"])

    def test_1_4_pair_below_severe_overlap_floor_is_flagged(self):
        """1-4 policy: the exclusion never hides a severe overlap — a heavy 1-4
        pair below the 2.2 A plain-distance floor means a broken structure."""
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "collapsed_chain.cif"
            d14 = self._write_cis_chain_cif(cif_path, 100.0)
            self.assertLess(d14, 2.2)
            record = _summary_record(cif_path, structure_id="collapsed_chain")

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "hard_invalid")
        self.assertIn("excluded_pair_severe_overlap", report.hard_invalid_reasons)
        self.assertEqual(report.metrics["n_excluded_pair_severe_overlaps"], 1)
        detail = report.metrics["excluded_pair_severe_overlap_min"]
        self.assertEqual(detail["exclusion"], "torsion_1-4")
        self.assertAlmostEqual(detail["distance"], d14, places=3)

    def test_bonded_pair_below_severe_overlap_floor_is_flagged(self):
        """Even a bonded pair at fused-nuclei distance is a broken structure."""
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "fused_bond.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("f_C1", "C", 0.40, 0.50, 0.50),
                    ("f_C2", "C", 0.45, 0.50, 0.50),
                ],
                bonds=[("f_C1", "f_C2", ".", ".", 0.5)],
            )
            record = _summary_record(cif_path, structure_id="fused_bond")

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "hard_invalid")
        self.assertIn("excluded_pair_severe_overlap", report.hard_invalid_reasons)
        self.assertNotIn("heavy_atom_clash", report.hard_invalid_reasons)
        self.assertEqual(report.metrics["excluded_pair_severe_overlap_min"]["exclusion"], "bonded_1-2")

    def test_hydrogen_contacts_reported_in_hydrogen_channel_not_heavy_clash(self):
        """An H...H contact below the heavy-atom thresholds lands in the
        hydrogen metric channel as a warning, never as heavy_atom_clash."""
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "hydrogen_contact.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("h2_H1", "H", 0.40, 0.50, 0.50),
                    ("h2_H2", "H", 0.57, 0.50, 0.50),
                ],
                bonds=[],
            )
            record = _summary_record(cif_path, structure_id="hydrogen_contact", include_bridge_event=False)
            record["metadata"]["graph_summary"] = {"n_monomer_instances": 1, "n_reaction_events": 1}

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "warning")
        self.assertIn("hydrogen_atom_clash", report.warning_reasons)
        self.assertNotIn("heavy_atom_clash", report.reasons)
        self.assertIsNone(report.metrics["min_nonbonded_heavy_distance"])
        self.assertAlmostEqual(report.metrics["min_nonbonded_hydrogen_distance"], 1.7, places=3)
        self.assertAlmostEqual(report.metrics["min_nonbonded_hydrogen_vdw_ratio"], 1.7 / 2.4, places=3)

    def test_hydrogen_contacts_skip_the_heavy_plain_distance_floor(self):
        """A C...H contact at 2.0 A is below the 2.2 A heavy floor but that
        floor is heavy-heavy only; and an H...H contact at 2.1 A (ratio 0.875)
        is no clash at all."""
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "ch_contact.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("h2_C1", "C", 0.30, 0.50, 0.50),
                    ("h2_H1", "H", 0.50, 0.50, 0.50),
                    ("h2_H2", "H", 0.50, 0.71, 0.50),
                ],
                bonds=[],
            )
            record = _summary_record(cif_path, structure_id="ch_contact", include_bridge_event=False)
            record["metadata"]["graph_summary"] = {"n_monomer_instances": 1, "n_reaction_events": 1}

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "warning")
        self.assertIn("hydrogen_atom_clash", report.warning_reasons)
        self.assertNotIn("heavy_atom_clash", report.reasons)
        # Min H contact is the 2.0 A C...H pair (ratio 0.69), not the 2.1 A H...H.
        self.assertAlmostEqual(report.metrics["min_nonbonded_hydrogen_distance"], 2.0, places=3)

    def test_hydrogen_contacts_above_ratio_cutoff_do_not_warn(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "hh_ok.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("h2_H1", "H", 0.40, 0.50, 0.50),
                    ("h2_H2", "H", 0.61, 0.50, 0.50),
                ],
                bonds=[],
            )
            record = _summary_record(cif_path, structure_id="hh_ok", include_bridge_event=False)
            record["metadata"]["graph_summary"] = {"n_monomer_instances": 1, "n_reaction_events": 1}

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "valid")
        self.assertAlmostEqual(report.metrics["min_nonbonded_hydrogen_distance"], 2.1, places=3)

    def test_unsupported_element_uses_fallback_radius_with_stderr_warning(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "unsupported_element.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("z1_Zn1", "Zn", 0.40, 0.50, 0.50),
                    ("z1_Zn2", "Zn", 0.64, 0.50, 0.50),
                ],
                bonds=[],
            )
            record = _summary_record(cif_path, structure_id="unsupported_element")
            record["metadata"]["graph_summary"] = {"n_monomer_instances": 1, "n_reaction_events": 1}

            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                report = CoarseStructureValidator().validate_manifest_record(record)

        # Fallback radius 1.70 A: 2.4 / 3.4 = 0.706 < 0.75 -> heavy clash.
        self.assertIn("heavy_atom_clash", report.hard_invalid_reasons)
        self.assertIn("warning: vdw: no Bondi vdW radius for element 'Zn'", stderr.getvalue())


        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "stretched_boronate.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("m1_B1", "B", 0.1, 0.1, 0.1),
                    ("m2_O1", "O", 0.29, 0.1, 0.1),
                ],
                bonds=[("m1_B1", "m2_O1")],
                bond_distance=1.9,
            )
            record = _summary_record(
                cif_path,
                structure_id="stretched_boronate",
                distance_residual=0.1,
                actual_distance=0.55,
                target_distance=0.55,
                template_id="boronate_ester_bridge",
            )

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "needs_optimization")
        self.assertIn("realized_bridge_bond_distance", report.needs_optimization_reasons)
        self.assertEqual(report.hard_invalid_reasons, ())
        self.assertEqual(report.metrics["realized_bridge_bond_count"], 1)
        self.assertAlmostEqual(report.metrics["realized_bridge_bond_distance_max"], 1.9)

    def test_validator_accepts_in_window_realized_boronate_bond(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "closed_boronate.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("m1_B1", "B", 0.1, 0.1, 0.1),
                    ("m2_O1", "O", 0.247, 0.1, 0.1),
                ],
                bonds=[("m1_B1", "m2_O1")],
                bond_distance=1.47,
            )
            record = _summary_record(
                cif_path,
                structure_id="closed_boronate",
                distance_residual=0.1,
                actual_distance=0.55,
                target_distance=0.55,
                template_id="boronate_ester_bridge",
            )

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "valid")
        self.assertEqual(report.metrics["realized_bridge_bond_count"], 1)

    def test_validator_parses_instance_ids_containing_underscores(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "underscore_instances.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("bex_node_B1", "B", 0.10, 0.10, 0.10),
                    ("bex_node_B2", "B", 0.10, 0.40, 0.10),
                    ("bex_link_O1", "O", 0.247, 0.10, 0.10),
                    ("bex_link_O2", "O", 0.247, 0.40, 0.10),
                ],
                bonds=[("bex_node_B1", "bex_link_O1"), ("bex_node_B2", "bex_link_O2")],
                bond_distance=1.47,
            )
            record = _summary_record(
                cif_path,
                structure_id="underscore_instances",
                distance_residual=0.1,
                actual_distance=0.55,
                target_distance=0.55,
                template_id="boronate_ester_bridge",
            )

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.metrics["n_instance_nodes"], 2)
        self.assertGreater(report.metrics["n_inter_instance_edges"], 0)
        self.assertEqual(report.metrics["n_instance_components"], 1)
        self.assertEqual(report.metrics["realized_bridge_bond_count"], 2)
        self.assertEqual(report.classification, "valid")

    def test_validator_excludes_bonded_pair_across_cell_boundary_from_clash_check(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "boundary_bond.cif"
            _write_boundary_bond_cif(cif_path)
            record = _summary_record(cif_path, structure_id="boundary_bond", include_bridge_event=False)

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "valid")
        self.assertNotIn("heavy_atom_clash", report.hard_invalid_reasons)
        self.assertAlmostEqual(report.metrics["min_nonbonded_heavy_distance"], 2.9, places=3)

    def test_soft_relax_excludes_bonded_pair_across_cell_boundary_from_clash_check(self):
        from cofkit import soft_relax

        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "boundary_bond.cif"
            _write_boundary_bond_cif(cif_path)
            block = gemmi.cif.read_file(str(cif_path)).sole_block()
            system, _warnings = soft_relax._parse_system(block, soft_relax.SoftRelaxConfig())
            minimum, clashes = soft_relax._min_heavy_distance_and_clashes(system, 1.05)

        self.assertIsNotNone(minimum)
        self.assertAlmostEqual(minimum, 2.9, places=3)
        self.assertEqual(clashes, 0)

    def test_validator_runs_clash_checks_even_when_metadata_invalid(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "metadata_bad_clash.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("x1_C1", "C", 0.1, 0.1, 0.1),
                    ("x1_C2", "C", 0.15, 0.1, 0.1),
                ],
                bonds=[],
            )
            record = _summary_record(
                cif_path,
                structure_id="metadata_bad_clash",
                include_bridge_event=False,
            )
            record["metadata"]["score_metadata"]["n_unreacted_motifs"] = 1

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "hard_invalid")
        self.assertIn("unreacted_motifs", report.hard_invalid_reasons)
        self.assertIn("heavy_atom_clash", report.hard_invalid_reasons)
        self.assertNotIn("cif_checks_skipped", report.metrics)
        self.assertLess(report.metrics["min_nonbonded_heavy_distance"], 1.05)

    def test_validator_skips_only_metadata_dependent_cif_checks_when_metadata_invalid(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "metadata_bad_boronate.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("m1_B1", "B", 0.1, 0.1, 0.1),
                    ("m2_O1", "O", 0.29, 0.1, 0.1),
                ],
                bonds=[("m1_B1", "m2_O1")],
                bond_distance=1.9,
            )
            record = _summary_record(
                cif_path,
                structure_id="metadata_bad_boronate",
                distance_residual=1.2,
                actual_distance=2.4,
                template_id="boronate_ester_bridge",
            )
            record["metadata"]["score_metadata"]["n_unreacted_motifs"] = 1

            skipped = CoarseStructureValidator().validate_manifest_record(record)
            enforced = CoarseStructureValidator(
                thresholds=CoarseValidationThresholds(skip_cif_checks_when_metadata_invalid=False)
            ).validate_manifest_record(record)

        self.assertEqual(skipped.classification, "hard_invalid")
        self.assertIn("unreacted_motifs", skipped.hard_invalid_reasons)
        self.assertTrue(skipped.metrics["cif_metadata_checks_skipped"])
        self.assertNotIn("realized_bridge_bond_distance", skipped.reasons)
        self.assertEqual(skipped.coverage["bridge_geometry"], "skipped")
        self.assertIn("min_nonbonded_heavy_distance", skipped.metrics)
        self.assertEqual(enforced.classification, "hard_invalid")
        self.assertIn("realized_bridge_bond_distance", enforced.hard_invalid_reasons)
        self.assertEqual(enforced.coverage["bridge_geometry"], "measured")

    def test_validator_does_not_treat_missing_bridge_distance_data_as_zero_residual(self):
        # One of two seed bridge events carries no distance data; the seed
        # aggregates are informational (seed_*), and the verdict comes from
        # the measured final coordinates (bond at 2.1 A vs the 1.3 A target).
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "missing_bridge_data.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("a1_C1", "C", 0.1, 0.1, 0.1),
                    ("b1_C1", "C", 0.31, 0.1, 0.1),
                ],
                bonds=[("a1_C1", "b1_C1")],
                bond_distance=2.1,
            )
            record = _summary_record(cif_path, structure_id="missing_bridge_data")
            record["metadata"]["score_metadata"]["bridge_event_metrics"] = [
                {},
                {"distance_residual": 0.8, "actual_distance": 2.1, "target_distance": 1.3},
            ]

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.metrics["n_bridge_events"], 2)
        self.assertEqual(report.metrics["n_bridge_events_missing_distance_data"], 1)
        self.assertAlmostEqual(report.metrics["seed_max_bridge_distance_residual"], 0.8)
        self.assertAlmostEqual(report.metrics["seed_mean_bridge_distance_residual"], 0.8)
        self.assertAlmostEqual(report.metrics["max_bridge_distance_residual"], 0.8)
        self.assertAlmostEqual(report.metrics["mean_bridge_distance_residual"], 0.8)
        self.assertEqual(report.metrics["bridge_metrics_source"], "final_cif")
        self.assertIn("bridge_distance_residual_max", report.warning_reasons)

    def test_missing_bond_loop_marks_bridge_geometry_unmeasured_not_valid(self):
        """A04: bridge events are claimed but the CIF has no bond loop, so the
        linkage geometry is required-but-unmeasured: the record must land in
        an explicit unvalidated state, never the valid bucket."""
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "no_bond_loop.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("a1_C1", "C", 0.1, 0.1, 0.1),
                    ("a1_C2", "C", 0.55, 0.1, 0.1),
                ],
                bonds=[],
            )
            record = _summary_record(cif_path, structure_id="no_bond_loop")

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "unvalidated")
        self.assertIsNone(report.is_valid)
        self.assertIsNone(report.passes_hard_validation)
        self.assertEqual(report.coverage["bridge_geometry"], "missing_data")
        self.assertIn("bridge_geometry", report.unmeasured_required_checks)
        self.assertIn("unmeasured_required_checks", report.reasons)
        self.assertEqual(report.metrics["bridge_metrics_source"], None)

    def test_completed_contact_scan_with_zero_neighbors_is_not_missing_data(self):
        """A04/T2-14: a completed bounded contact search that finds no
        assessable neighbor is `no_contacts`, distinct from `missing_data`."""
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "isolated_pair.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("m1_C1", "C", 0.40, 0.50, 0.50),
                    ("m1_C2", "C", 0.55, 0.50, 0.50),
                ],
                bonds=[("m1_C1", "m1_C2", ".", ".", 1.5)],
            )
            record = _summary_record(cif_path, structure_id="isolated_pair", include_bridge_event=False)

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "valid")
        self.assertEqual(report.coverage["contact_scan"], "no_contacts")
        self.assertEqual(report.metrics["n_nonbonded_pairs_scanned"], 0)
        self.assertIsNone(report.metrics["min_nonbonded_heavy_distance"])
        self.assertEqual(report.coverage["bridge_geometry"], "not_applicable")

    def test_clash_contact_scan_reports_measured_coverage(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "clash.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("x1_C1", "C", 0.1, 0.1, 0.1),
                    ("x1_C2", "C", 0.15, 0.1, 0.1),
                ],
                bonds=[],
            )
            record = _summary_record(cif_path, structure_id="clash", include_bridge_event=False)

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "hard_invalid")
        self.assertEqual(report.coverage["contact_scan"], "measured")
        self.assertGreater(report.metrics["n_nonbonded_pairs_scanned"], 0)

    def test_topology_dimensionality_warns_on_unexpected_lookup_failure(self):
        validator = CoarseStructureValidator()
        self.assertIsNone(validator._topology_dimensionality("definitely-unknown-topology"))
        stderr = io.StringIO()
        with mock.patch("cofkit.validation.get_topology_hint", side_effect=RuntimeError("boom")):
            with contextlib.redirect_stderr(stderr):
                result = validator._topology_dimensionality("hcb")
        self.assertIsNone(result)
        self.assertIn("warning: topology dimensionality lookup failed", stderr.getvalue())

    def test_classifier_sorts_valid_and_invalid_outputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source = root / "batch_out"
            cifs = source / "cifs"
            cifs.mkdir(parents=True)
            valid_cif = cifs / "valid.cif"
            warning_cif = cifs / "warning.cif"
            invalid_cif = cifs / "invalid.cif"
            needs_optimization_cif = cifs / "needs_optimization.cif"
            hard_hard_cif = cifs / "hard_hard.cif"
            unvalidated_cif = cifs / "unvalidated.cif"
            _write_test_cif(
                valid_cif,
                atoms=[
                    ("a1_C1", "C", 0.1, 0.1, 0.1),
                    ("b1_C1", "C", 0.25, 0.1, 0.1),
                ],
                bonds=[("a1_C1", "b1_C1")],
            )
            _write_test_cif(
                warning_cif,
                atoms=[
                    ("w1_C1", "C", 0.1, 0.1, 0.1),
                    ("w2_C1", "C", 0.29, 0.1, 0.1),
                ],
                bonds=[("w1_C1", "w2_C1")],
                bond_distance=1.9,
            )
            _write_test_cif(
                invalid_cif,
                atoms=[
                    ("x1_C1", "C", 0.1, 0.1, 0.1),
                    ("x1_C2", "C", 0.15, 0.1, 0.1),
                ],
                bonds=[],
            )
            _write_test_cif(
                needs_optimization_cif,
                atoms=[
                    ("n1_C1", "C", 0.1, 0.1, 0.1),
                    ("n2_C1", "C", 0.34, 0.1, 0.1),
                ],
                bonds=[("n1_C1", "n2_C1")],
                bond_distance=2.4,
            )
            _write_test_cif(
                hard_hard_cif,
                atoms=[
                    ("h1_C1", "C", 0.1, 0.1, 0.1),
                    ("h2_C1", "C", 0.36, 0.1, 0.1),
                ],
                bonds=[("h1_C1", "h2_C1")],
                bond_distance=2.6,
            )
            # Bridge events claimed but no bond loop: required bridge-geometry
            # coverage is unmeasured, so this must not land in a valid bucket.
            _write_test_cif(
                unvalidated_cif,
                atoms=[
                    ("u1_C1", "C", 0.1, 0.1, 0.1),
                    ("u1_C2", "C", 0.55, 0.1, 0.1),
                ],
                bonds=[],
            )
            manifest_rows = [
                _summary_record(valid_cif, structure_id="valid"),
                _summary_record(
                    warning_cif,
                    structure_id="warning",
                    distance_residual=0.6,
                    actual_distance=1.9,
                ),
                _summary_record(
                    hard_hard_cif,
                    structure_id="hard_hard",
                    distance_residual=1.3,
                    actual_distance=2.5,
                ),
                _summary_record(
                    needs_optimization_cif,
                    structure_id="needs_optimization",
                    distance_residual=1.2,
                    actual_distance=2.4,
                ),
                _summary_record(invalid_cif, structure_id="invalid"),
                _summary_record(unvalidated_cif, structure_id="unvalidated"),
            ]
            (source / "manifest.jsonl").write_text(
                "".join(json.dumps(row, sort_keys=True) + "\n" for row in manifest_rows),
                encoding="utf-8",
            )

            output = root / "classified"
            summary = classify_batch_output(source, output, link_mode="copy", max_workers=2)

            self.assertEqual(summary.total_structures, 6)
            self.assertEqual(summary.valid_structures, 1)
            self.assertEqual(summary.warning_structures, 1)
            self.assertEqual(summary.needs_optimization_structures, 1)
            self.assertEqual(summary.hard_hard_invalid_structures, 1)
            self.assertEqual(summary.hard_invalid_structures, 1)
            self.assertEqual(summary.unvalidated_structures, 1)
            self.assertEqual(
                summary.total_structures,
                summary.valid_structures
                + summary.warning_structures
                + summary.needs_optimization_structures
                + summary.hard_hard_invalid_structures
                + summary.hard_invalid_structures
                + summary.unvalidated_structures,
            )
            self.assertTrue((output / "valid" / "cifs" / "valid.cif").is_file())
            self.assertTrue((output / "warning" / "cifs" / "warning.cif").is_file())
            self.assertTrue((output / "warning" / "reasons" / "bridge_distance_residual_mean" / "warning.cif").is_file())
            self.assertTrue((output / "needs_optimization" / "cifs" / "needs_optimization.cif").is_file())
            self.assertTrue(
                (
                    output
                    / "needs_optimization"
                    / "reasons"
                    / "bridge_distance_residual_max_hard"
                    / "needs_optimization.cif"
                ).is_file()
            )
            self.assertTrue((output / "hard_hard_invalid" / "cifs" / "hard_hard.cif").is_file())
            self.assertTrue(
                (
                    output
                    / "hard_hard_invalid"
                    / "reasons"
                    / "bridge_distance_exceeds_cif_export_limit"
                    / "hard_hard.cif"
                ).is_file()
            )
            self.assertTrue((output / "hard_invalid" / "cifs" / "invalid.cif").is_file())
            self.assertTrue((output / "hard_invalid" / "reasons" / "heavy_atom_clash" / "invalid.cif").is_file())
            self.assertTrue((output / "unvalidated" / "cifs" / "unvalidated.cif").is_file())
            self.assertTrue((output / "unvalidated" / "manifest.jsonl").is_file())

            # Serialized output distinguishes missing data from a completed
            # contact search that found no neighbors.
            classified_rows = [
                json.loads(line)
                for line in (output / "classification_manifest.jsonl").read_text(encoding="utf-8").splitlines()
                if line.strip()
            ]
            by_id = {row["structure_id"]: row for row in classified_rows}
            unvalidated_row = by_id["unvalidated"]["validation"]
            self.assertEqual(unvalidated_row["classification"], "unvalidated")
            self.assertIsNone(unvalidated_row["is_valid"])
            self.assertEqual(unvalidated_row["coverage"]["bridge_geometry"], "missing_data")
            self.assertIn("bridge_geometry", unvalidated_row["unmeasured_required_checks"])
            valid_row = by_id["valid"]["validation"]
            self.assertEqual(valid_row["coverage"]["bridge_geometry"], "measured")
            self.assertEqual(valid_row["coverage"]["contact_scan"], "no_contacts")
            self.assertIsNone(valid_row["metrics"]["min_nonbonded_heavy_distance"])


if __name__ == "__main__":
    unittest.main()
