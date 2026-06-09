import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit import BatchGenerationConfig, BatchMonomerRecord, BatchStructureGenerator
from cofkit.cli import main as cli_main
from cofkit.cofid import COFidMonomer, canonicalize_smiles, serialize_cofid
from cofkit.validate import validate_cif_against_cofid


TAPB = "C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N"
TFB = "C1=C(C=C(C=C1C=O)C=O)C=O"
PDA = "O=Cc1ccc(C=O)cc1"


def _generator() -> BatchStructureGenerator:
    return BatchStructureGenerator(
        BatchGenerationConfig(
            allowed_reactions=("imine_bridge",),
            single_node_topology_ids=("hcb",),
            topology_ids=("hcb",),
            enumerate_all_topologies=False,
            write_cif=True,
            rdkit_num_conformers=2,
        )
    )


def _record(
    record_id: str,
    smiles: str,
    motif_kind: str,
    connectivity: int,
) -> BatchMonomerRecord:
    return BatchMonomerRecord(
        id=record_id,
        name=record_id,
        smiles=smiles,
        motif_kind=motif_kind,
        expected_connectivity=connectivity,
    )


def _generated_hcb_imine(temp_path: Path):
    return _generator().generate_pair_candidate(
        _record("tapb", TAPB, "amine", 3),
        _record("tfb", TFB, "aldehyde", 3),
        out_dir=temp_path,
        write_cif=True,
    )


def _wrong_imine_cofid() -> str:
    return serialize_cofid(
        monomers=(
            COFidMonomer(3, "amine", canonicalize_smiles(TAPB)),
            COFidMonomer(2, "aldehyde", canonicalize_smiles(PDA)),
        ),
        topology="hcb",
        linkage="imine",
    )


class ValidateCliTests(unittest.TestCase):
    def test_validate_simple_forces_distance_inferred_decomposition(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generated_hcb_imine(temp_path)

            result = validate_cif_against_cofid(summary.metadata["cofid"], summary.cif_path)

        self.assertTrue(result.ok, result.reason)
        self.assertFalse(result.topology_compared)
        self.assertTrue(result.monomers_match)
        self.assertTrue(result.linkage_match)
        self.assertEqual(result.decomposition.metadata["bond_source"], "distance_inferred")
        self.assertGreater(result.decomposition.metadata["n_distance_inferred_bonds"], 0)

    def test_validate_simple_cli_prints_match_summary(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generated_hcb_imine(temp_path)
            buffer = io.StringIO()

            with contextlib.redirect_stdout(buffer):
                cli_main(["validate", "simple", summary.metadata["cofid"], str(summary.cif_path)])

        output = buffer.getvalue()
        self.assertIn("status: match", output)
        self.assertIn("topology_compared: false", output)
        self.assertIn("bond_source: distance_inferred", output)

    def test_validate_simple_cli_exits_nonzero_on_mismatch_after_printing_result(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generated_hcb_imine(temp_path)
            buffer = io.StringIO()

            with self.assertRaises(SystemExit) as raised, contextlib.redirect_stdout(buffer):
                cli_main(["validate", "simple", _wrong_imine_cofid(), str(summary.cif_path)])

        self.assertEqual(raised.exception.code, 1)
        output = buffer.getvalue()
        self.assertIn("status: mismatch", output)
        self.assertIn("monomers_match: false", output)
        self.assertIn("linkage_match: true", output)

    def test_validate_optimize_cli_uses_default_lammps_pipeline_then_distance_decomposes(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generated_hcb_imine(temp_path)
            fake_lammps_result = SimpleNamespace(
                optimized_cif=str(summary.cif_path),
                output_dir=str(temp_path / "lammps_out"),
                to_dict=lambda: {
                    "optimized_cif": str(summary.cif_path),
                    "output_dir": str(temp_path / "lammps_out"),
                },
            )
            buffer = io.StringIO()

            with patch("cofkit.validate.optimize_cif_with_lammps", return_value=fake_lammps_result) as optimized:
                with contextlib.redirect_stdout(buffer):
                    cli_main(
                        [
                            "validate",
                            "optimize",
                            summary.metadata["cofid"],
                            str(summary.cif_path),
                            "--output-dir",
                            str(temp_path / "lammps_out"),
                            "--json",
                        ]
                    )

        optimized.assert_called_once()
        self.assertNotIn("settings", optimized.call_args.kwargs)
        report = json.loads(buffer.getvalue())
        self.assertEqual(report["status"], "match")
        self.assertEqual(report["optimized_cif"], str(summary.cif_path))
        self.assertEqual(report["decomposition"]["metadata"]["bond_source"], "distance_inferred")


if __name__ == "__main__":
    unittest.main()
