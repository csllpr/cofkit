import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit.validation import CoarseStructureValidator, classify_batch_output

try:
    import gemmi  # noqa: F401
except ImportError:  # pragma: no cover - environment-dependent
    gemmi = None


def _write_test_cif(path: Path, *, atoms: list[tuple[str, str, float, float, float]], bonds: list[tuple[str, str]]) -> None:
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
        for left, right in bonds:
            lines.append(f"{left} {right} . . 1.500000")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _summary_record(
    cif_path: Path,
    *,
    structure_id: str,
    topology_id: str = "hcb",
    distance_residual: float = 0.1,
    actual_distance: float = 1.4,
    target_distance: float = 1.3,
) -> dict[str, object]:
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
            "graph_summary": {"n_monomer_instances": 2, "n_reaction_events": 1, "reaction_templates": {"imine_bridge": 1}},
            "score_metadata": {
                "n_unreacted_motifs": 0,
                "bridge_event_metrics": [
                    {
                        "distance_residual": distance_residual,
                        "actual_distance": actual_distance,
                        "target_distance": target_distance,
                    }
                ],
            },
        },
    }


@unittest.skipIf(gemmi is None, "gemmi is not available")
class CoarseValidationTests(unittest.TestCase):
    def test_validator_rejects_bridge_distance_metadata(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "valid.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("a1_C1", "C", 0.1, 0.1, 0.1),
                    ("b1_C1", "C", 0.25, 0.1, 0.1),
                ],
                bonds=[("a1_C1", "b1_C1")],
            )
            record = _summary_record(
                cif_path,
                structure_id="too_far",
                distance_residual=1.2,
                actual_distance=2.4,
            )

            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "hard_invalid")
        self.assertFalse(report.passes_hard_validation)
        self.assertIn("bridge_distance_residual_max_hard", report.hard_invalid_reasons)
        self.assertIn("bridge_distance_too_long", report.hard_invalid_reasons)

    def test_validator_marks_overlong_bridge_as_hard_hard_invalid(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cif_path = root / "hard_hard.cif"
            _write_test_cif(
                cif_path,
                atoms=[
                    ("a1_C1", "C", 0.1, 0.1, 0.1),
                    ("b1_C1", "C", 0.25, 0.1, 0.1),
                ],
                bonds=[("a1_C1", "b1_C1")],
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
                    ("b1_C1", "C", 0.25, 0.1, 0.1),
                ],
                bonds=[("a1_C1", "b1_C1")],
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

    def test_classifier_sorts_valid_and_invalid_outputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source = root / "batch_out"
            cifs = source / "cifs"
            cifs.mkdir(parents=True)
            valid_cif = cifs / "valid.cif"
            warning_cif = cifs / "warning.cif"
            invalid_cif = cifs / "invalid.cif"
            hard_hard_cif = cifs / "hard_hard.cif"
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
                    ("w2_C1", "C", 0.25, 0.1, 0.1),
                ],
                bonds=[("w1_C1", "w2_C1")],
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
                hard_hard_cif,
                atoms=[
                    ("h1_C1", "C", 0.1, 0.1, 0.1),
                    ("h2_C1", "C", 0.25, 0.1, 0.1),
                ],
                bonds=[("h1_C1", "h2_C1")],
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
                _summary_record(invalid_cif, structure_id="invalid"),
            ]
            (source / "manifest.jsonl").write_text(
                "".join(json.dumps(row, sort_keys=True) + "\n" for row in manifest_rows),
                encoding="utf-8",
            )

            output = root / "classified"
            summary = classify_batch_output(source, output, link_mode="copy", max_workers=2)

            self.assertEqual(summary.total_structures, 4)
            self.assertEqual(summary.valid_structures, 1)
            self.assertEqual(summary.warning_structures, 1)
            self.assertEqual(summary.hard_hard_invalid_structures, 1)
            self.assertEqual(summary.hard_invalid_structures, 1)
            self.assertTrue((output / "valid" / "cifs" / "valid.cif").is_file())
            self.assertTrue((output / "warning" / "cifs" / "warning.cif").is_file())
            self.assertTrue((output / "warning" / "reasons" / "bridge_distance_residual_mean" / "warning.cif").is_file())
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


if __name__ == "__main__":
    unittest.main()
