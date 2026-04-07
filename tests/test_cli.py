import contextlib
import io
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit import __version__ as package_version
from cofkit.batch_models import BatchPairSummary
from cofkit.cli import main as cli_main

try:
    from rdkit import Chem  # noqa: F401
except ImportError:  # pragma: no cover - environment-dependent
    Chem = None


def _canonical(smiles: str) -> str:
    assert Chem is not None
    molecule = Chem.MolFromSmiles(smiles)
    assert molecule is not None
    return str(Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False))


class BatchSummaryCompatibilityTests(unittest.TestCase):
    def test_legacy_aliases_resolve_from_role_metadata(self):
        summary = BatchPairSummary(
            structure_id="pair__hcb",
            pair_id="pair",
            pair_mode="3+2-node-linker",
            status="ok",
            reactant_a_record_id="tri_amine",
            reactant_b_record_id="di_aldehyde",
            reactant_a_connectivity=3,
            reactant_b_connectivity=2,
            metadata={
                "reactant_roles": ("amine", "aldehyde"),
                "reactant_record_ids": {"amine": "tri_amine", "aldehyde": "di_aldehyde"},
                "reactant_connectivities": {"amine": 3, "aldehyde": 2},
            },
        )

        self.assertEqual(summary.amine_record_id, "tri_amine")
        self.assertEqual(summary.aldehyde_record_id, "di_aldehyde")
        self.assertEqual(summary.amine_connectivity, 3)
        self.assertEqual(summary.aldehyde_connectivity, 2)


class CliTests(unittest.TestCase):
    def test_root_help_lists_grouped_namespaces(self):
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            cli_main([])

        help_text = buffer.getvalue()
        self.assertIn("build", help_text)
        self.assertIn("analyze", help_text)
        self.assertIn("calculate", help_text)

    def test_package_version_uses_calver(self):
        self.assertRegex(package_version, r"^\d{4}\.\d{1,2}\.\d{1,2}(?:\.post\d+)?$")

    def test_root_version_reports_package_version(self):
        buffer = io.StringIO()
        with self.assertRaises(SystemExit) as raised, contextlib.redirect_stdout(buffer):
            cli_main(["--version"])

        self.assertEqual(raised.exception.code, 0)
        self.assertEqual(buffer.getvalue().strip(), f"cofkit {package_version}")

    def test_cli_loads_dotenv_from_parent_directory(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            (temp_path / ".env").write_text("COFKIT_TEST_SENTINEL=from_dotenv\n", encoding="utf-8")
            nested_dir = temp_path / "nested" / "child"
            nested_dir.mkdir(parents=True)
            previous_cwd = Path.cwd()
            try:
                with patch.dict(os.environ, {}, clear=True):
                    os.chdir(nested_dir)
                    buffer = io.StringIO()
                    with contextlib.redirect_stdout(buffer):
                        cli_main(["build", "list-templates", "--json"])
                    self.assertEqual(os.environ.get("COFKIT_TEST_SENTINEL"), "from_dotenv")
            finally:
                os.chdir(previous_cwd)

    def test_cli_dotenv_does_not_override_existing_environment(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            (temp_path / ".env").write_text("COFKIT_TEST_SENTINEL=from_dotenv\n", encoding="utf-8")
            previous_cwd = Path.cwd()
            try:
                with patch.dict(os.environ, {"COFKIT_TEST_SENTINEL": "explicit"}, clear=True):
                    os.chdir(temp_path)
                    buffer = io.StringIO()
                    with contextlib.redirect_stdout(buffer):
                        cli_main(["build", "list-templates", "--json"])
                    self.assertEqual(os.environ.get("COFKIT_TEST_SENTINEL"), "explicit")
            finally:
                os.chdir(previous_cwd)

    def test_list_templates_json_reports_execution_capabilities(self):
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            cli_main(["build", "list-templates", "--json"])
        rows = json.loads(buffer.getvalue())

        azine = next(row for row in rows if row["template_id"] == "azine_bridge")
        imine = next(row for row in rows if row["template_id"] == "imine_bridge")
        boroxine = next(row for row in rows if row["template_id"] == "boroxine_trimerization")

        self.assertTrue(azine["supports_pair_generation"])
        self.assertTrue(azine["supports_atomistic_realization"])
        self.assertEqual(azine["workflow_family"], "binary_bridge")
        self.assertTrue(imine["supports_pair_generation"])
        self.assertTrue(imine["supports_atomistic_realization"])
        self.assertEqual(imine["workflow_family"], "binary_bridge")
        self.assertFalse(boroxine["supports_pair_generation"])
        self.assertEqual(boroxine["workflow_family"], "ring_forming")

    def test_legacy_list_templates_alias_still_works(self):
        stdout_buffer = io.StringIO()
        stderr_buffer = io.StringIO()
        with contextlib.redirect_stdout(stdout_buffer), contextlib.redirect_stderr(stderr_buffer):
            cli_main(["list-templates", "--json"])

        rows = json.loads(stdout_buffer.getvalue())
        self.assertTrue(any(row["template_id"] == "imine_bridge" for row in rows))
        self.assertIn("deprecated", stderr_buffer.getvalue())
        self.assertIn("cofkit build list-templates", stderr_buffer.getvalue())

    @unittest.skipIf(Chem is None, "RDKit is not available")
    def test_single_pair_cli_writes_cif(self):
        tapb = "C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N"
        tfb = "C1=C(C=C(C=C1C=O)C=O)C=O"
        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "single_pair"
            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                cli_main(
                    [
                        "build",
                        "single-pair",
                        "--template-id",
                        "imine_bridge",
                        "--first-smiles",
                        tapb,
                        "--second-smiles",
                        tfb,
                        "--first-id",
                        "tapb",
                        "--second-id",
                        "tfb",
                        "--first-motif-kind",
                        "amine",
                        "--second-motif-kind",
                        "aldehyde",
                        "--topology",
                        "hcb",
                        "--output-dir",
                        str(output_dir),
                    ]
                )

            summary_path = output_dir / "summary.json"
            self.assertTrue(summary_path.exists())
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            self.assertEqual(summary["template_id"], "imine_bridge")
            self.assertEqual(summary["successful_structures"], 1)
            cif_path = Path(summary["results"][0]["cif_path"])
            self.assertTrue(cif_path.exists())
            self.assertEqual(cif_path.suffix, ".cif")
            cif_lines = cif_path.read_text(encoding="utf-8").splitlines()
            self.assertEqual(cif_lines[0], f"# COFid: {summary['results'][0]['metadata']['cofid']}")

    @unittest.skipIf(Chem is None, "RDKit is not available")
    def test_single_pair_cli_accepts_cofid(self):
        tapb = "C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N"
        tfb = "C1=C(C=C(C=C1C=O)C=O)C=O"
        monomers = (
            (3, "amine", _canonical(tapb)),
            (3, "aldehyde", _canonical(tfb)),
        )
        cofid = ".".join(
            f"{connectivity}:{reactive_group}:{smiles}"
            for connectivity, reactive_group, smiles in sorted(
                monomers,
                key=lambda item: (-item[0], item[2]),
            )
        ) + "&&hcb&&imine"

        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "single_pair_cofid"
            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                cli_main(
                    [
                        "build",
                        "single-pair",
                        "--cofid",
                        cofid,
                        "--output-dir",
                        str(output_dir),
                    ]
                )

            summary = json.loads((output_dir / "summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["requested_cofid"], cofid)
            self.assertEqual(summary["template_id"], "imine_bridge")
            self.assertEqual(summary["successful_structures"], 1)
            self.assertEqual(summary["results"][0]["metadata"]["cofid"], cofid)

    def test_single_pair_cli_rejects_internal_post_build_conversion_flag(self):
        stderr_buffer = io.StringIO()
        with self.assertRaises(SystemExit) as raised, contextlib.redirect_stderr(stderr_buffer):
            cli_main(
                [
                    "build",
                    "single-pair",
                    "--template-id",
                    "imine_bridge",
                    "--first-smiles",
                    "NC1=CC=C(C=C1)N",
                    "--second-smiles",
                    "O=CC1=CC=C(C=C1)C=O",
                    "--annotate-post-build-conversion",
                    "sulfur_enabled_imine_conversion",
                ]
            )

        self.assertEqual(raised.exception.code, 2)
        self.assertIn("unrecognized arguments", stderr_buffer.getvalue())
        self.assertIn("--annotate-post-build-conversion", stderr_buffer.getvalue())


if __name__ == "__main__":
    unittest.main()
