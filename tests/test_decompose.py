import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit import BatchGenerationConfig, BatchMonomerRecord, BatchStructureGenerator
from cofkit.cli import main as cli_main
from cofkit.decompose import decompose_cif_to_cofid


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


def _without_cofid_comment(source_path: str | Path, target_path: Path) -> Path:
    lines = Path(source_path).read_text(encoding="utf-8").splitlines()
    if lines and lines[0].startswith("# COFid: "):
        lines = lines[1:]
    target_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target_path


class DecomposeRoundTripTests(unittest.TestCase):
    def test_decompose_generated_hcb_three_three_imine_returns_original_inputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generator().generate_pair_candidate(
                _record("tapb", TAPB, "amine", 3),
                _record("tfb", TFB, "aldehyde", 3),
                out_dir=temp_path,
                write_cif=True,
            )
            input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped.cif")

            result = decompose_cif_to_cofid(input_cif, topology="hcb")

        self.assertTrue(result.ok, result.reason)
        self.assertEqual(result.cofid, summary.metadata["cofid"])
        self.assertEqual(result.metadata["n_imine_linkage_bonds"], 3)
        self.assertEqual(result.metadata["n_unique_monomers"], 2)

    def test_decompose_generated_hcb_three_two_imine_returns_original_inputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generator().generate_pair_candidate(
                _record("tapb", TAPB, "amine", 3),
                _record("pda", PDA, "aldehyde", 2),
                out_dir=temp_path,
                write_cif=True,
            )
            input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped.cif")

            result = decompose_cif_to_cofid(input_cif, topology="hcb")

        self.assertTrue(result.ok, result.reason)
        self.assertEqual(result.cofid, summary.metadata["cofid"])
        self.assertEqual(result.metadata["n_imine_linkage_bonds"], 6)
        self.assertEqual(result.metadata["n_unique_monomers"], 2)

    def test_analyze_decompose_cli_prints_recovered_cofid(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generator().generate_pair_candidate(
                _record("tapb", TAPB, "amine", 3),
                _record("tfb", TFB, "aldehyde", 3),
                out_dir=temp_path,
                write_cif=True,
            )
            input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped.cif")

            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                cli_main(["analyze", "decompose", str(input_cif), "--topology", "hcb"])

        self.assertEqual(buffer.getvalue().strip(), summary.metadata["cofid"])


if __name__ == "__main__":
    unittest.main()
