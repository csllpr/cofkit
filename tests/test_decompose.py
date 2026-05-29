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
from cofkit.reactions import ReactionLibrary


TAPB = "C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N"
TFB = "C1=C(C=C(C=C1C=O)C=O)C=O"
PDA = "O=Cc1ccc(C=O)cc1"
COF42_HYDRAZIDE = "CCOc1cc(C(=O)NN)cc(C(=O)NN)c1OCC"
HYDRAZINE = "NN"
TP = "O=Cc1c(O)c(C=O)c(O)c(C=O)c1O"
PPD = "Nc1ccc(N)cc1"
BDBA = "OB(O)c1ccc(B(O)O)cc1"
HHTP = "OC1=C(O)C=C2C(=C1)C1=CC(O)=C(O)C=C1C1=CC(O)=C(O)C=C21"
TMT = "Cc1nc(C)nc(C)n1"


def _generator(template_id: str = "imine_bridge") -> BatchStructureGenerator:
    return BatchStructureGenerator(
        BatchGenerationConfig(
            allowed_reactions=(template_id,),
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


ALL_BINARY_LINKAGE_ROUND_TRIP_CASES = (
    (
        "hydrazone_bridge",
        "hydrazone",
        _record("hdz", COF42_HYDRAZIDE, "hydrazide", 2),
        _record("tfb", TFB, "aldehyde", 3),
        6,
    ),
    (
        "azine_bridge",
        "azine",
        _record("hyd", HYDRAZINE, "hydrazine", 2),
        _record("tfb", TFB, "aldehyde", 3),
        6,
    ),
    (
        "keto_enamine_bridge",
        "bken",
        _record("tp", TP, "keto_aldehyde", 3),
        _record("ppd", PPD, "amine", 2),
        6,
    ),
    (
        "boronate_ester_bridge",
        "boest",
        _record("bdba", BDBA, "boronic_acid", 2),
        _record("hhtp", HHTP, "catechol", 3),
        12,
    ),
    (
        "vinylene_bridge",
        "vinylene",
        _record("tmt", TMT, "activated_methylene", 3),
        _record("pda", PDA, "aldehyde", 2),
        6,
    ),
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

    def test_decompose_generated_hcb_cofs_for_all_build_binary_linkages(self):
        for template_id, linkage, first, second, expected_linkage_bonds in ALL_BINARY_LINKAGE_ROUND_TRIP_CASES:
            with self.subTest(template_id=template_id), tempfile.TemporaryDirectory() as temp_dir:
                temp_path = Path(temp_dir)
                summary, _candidate = _generator(template_id).generate_pair_candidate(
                    first,
                    second,
                    out_dir=temp_path,
                    write_cif=True,
                )
                input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped.cif")

                result = decompose_cif_to_cofid(input_cif, topology="hcb", linkage=linkage)

                self.assertTrue(result.ok, result.reason)
                self.assertEqual(result.cofid, summary.metadata["cofid"])
                self.assertEqual(result.metadata[f"n_{linkage}_linkage_bonds"], expected_linkage_bonds)
                self.assertEqual(result.metadata["n_unique_monomers"], 2)

    def test_generated_decompose_cases_cover_buildable_binary_linkages(self):
        library = ReactionLibrary.builtin()
        buildable = {
            template_id
            for template_id, profile in library.linkage_profiles.items()
            if profile.supports_binary_bridge_pair_generation and profile.supports_atomistic_realization
        }
        covered = {"imine_bridge"} | {
            template_id for template_id, _linkage, _first, _second, _expected in ALL_BINARY_LINKAGE_ROUND_TRIP_CASES
        }

        self.assertEqual(covered, buildable)

    def test_decompose_accepts_template_id_linkage_alias(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generator("hydrazone_bridge").generate_pair_candidate(
                _record("hdz", COF42_HYDRAZIDE, "hydrazide", 2),
                _record("tfb", TFB, "aldehyde", 3),
                out_dir=temp_path,
                write_cif=True,
            )
            input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped.cif")

            result = decompose_cif_to_cofid(input_cif, topology="hcb", linkage="hydrazone_bridge")

        self.assertTrue(result.ok, result.reason)
        self.assertEqual(result.linkage, "hydrazone")
        self.assertEqual(result.cofid, summary.metadata["cofid"])


if __name__ == "__main__":
    unittest.main()
