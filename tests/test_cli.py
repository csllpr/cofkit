import contextlib
import io
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit import __version__ as package_version
from cofkit.batch_models import BatchPairSummary
from cofkit.cli import main as cli_main
from cofkit import cli_build

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


class DefaultRepoPathTests(unittest.TestCase):
    def test_source_checkout_resolves_against_repo_root(self):
        resolved = cli_build._default_repo_path("out", "single_pair_generation")

        expected_root = Path(cli_build.__file__).resolve().parents[2]
        self.assertEqual(resolved, str(expected_root / "out" / "single_pair_generation"))

    def test_installed_package_falls_back_to_cwd_relative(self):
        fake_installed_file = "/venv/lib/python3.11/site-packages/cofkit/cli_build.py"
        with patch.object(cli_build, "__file__", fake_installed_file):
            resolved = cli_build._default_repo_path("out", "single_pair_generation")

        self.assertEqual(resolved, str(Path("out") / "single_pair_generation"))


@unittest.skipIf(Chem is None, "RDKit is not available")
class DefaultLibraryOutputGuardTests(unittest.TestCase):
    def _make_input_dir(self, root: Path) -> Path:
        input_dir = root / "input"
        input_dir.mkdir()
        (input_dir / "amines_count_2.txt").write_text(
            "smiles\nNc1ccc(N)cc1\n",
            encoding="utf-8",
        )
        return input_dir

    def _run(self, input_dir: Path, output_dir: Path, *, force: bool = False) -> None:
        args = SimpleNamespace(
            input_dir=str(input_dir),
            output_dir=str(output_dir),
            num_conformers=1,
            force=force,
        )
        with contextlib.redirect_stdout(io.StringIO()):
            cli_build._run_default_library(args)

    def test_refuses_to_delete_foreign_non_empty_directory(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            input_dir = self._make_input_dir(root)
            output_dir = root / "output"
            output_dir.mkdir()
            sentinel = output_dir / "important.txt"
            sentinel.write_text("do not delete", encoding="utf-8")

            with self.assertRaises(SystemExit):
                self._run(input_dir, output_dir)

            self.assertTrue(sentinel.is_file())

    def test_force_replaces_foreign_directory_and_writes_marker(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            input_dir = self._make_input_dir(root)
            output_dir = root / "output"
            output_dir.mkdir()
            sentinel = output_dir / "important.txt"
            sentinel.write_text("delete me", encoding="utf-8")

            self._run(input_dir, output_dir, force=True)

            self.assertFalse(sentinel.exists())
            self.assertTrue((output_dir / cli_build._DEFAULT_LIBRARY_MARKER).is_file())

    def test_rerun_without_force_succeeds_when_marker_present(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            input_dir = self._make_input_dir(root)
            output_dir = root / "output"

            self._run(input_dir, output_dir)
            self.assertTrue((output_dir / cli_build._DEFAULT_LIBRARY_MARKER).is_file())
            # Second run must not raise even without --force.
            self._run(input_dir, output_dir)

            self.assertTrue((output_dir / "registry.jsonl").is_file())

    def test_missing_input_dir_fails_with_clear_message(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            args = SimpleNamespace(
                input_dir=str(root / "does-not-exist"),
                output_dir=str(root / "output"),
                num_conformers=1,
                force=False,
            )
            with self.assertRaises(SystemExit) as ctx:
                cli_build._run_default_library(args)

            self.assertIn("input directory does not exist", str(ctx.exception))


class CliTests(unittest.TestCase):
    def test_root_help_lists_grouped_namespaces(self):
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            cli_main([])

        help_text = buffer.getvalue()
        self.assertIn("build", help_text)
        self.assertIn("analyze", help_text)
        self.assertIn("calculate", help_text)
        self.assertIn("validate", help_text)

    def test_batch_build_help_lists_repair_geometry(self):
        buffer = io.StringIO()
        with self.assertRaises(SystemExit) as raised, contextlib.redirect_stdout(buffer):
            cli_main(["build", "batch-binary-bridge", "--help"])

        self.assertEqual(raised.exception.code, 0)
        self.assertIn("--repair-geometry", buffer.getvalue())

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
        self.assertTrue(boroxine["supports_ring_generation"])
        self.assertTrue(boroxine["supports_atomistic_realization"])
        self.assertTrue(boroxine["supports_topology_guided_generation"])
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

    def test_single_pair_cli_writes_stacking_suffix_in_cif_comment(self):
        tapb = "C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N"
        tfb = "C1=C(C=C(C=C1C=O)C=O)C=O"
        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "single_pair_stacked"
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
                        "--stacking",
                        "AA",
                        "--output-dir",
                        str(output_dir),
                    ]
                )

            summary = json.loads((output_dir / "summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["successful_structures"], 1)
            self.assertEqual(summary["results"][0]["structure_id"], "tapb__tfb__hcb__AA")
            self.assertEqual(summary["results"][0]["metadata"]["stacking"]["id"], "AA")
            cif_path = Path(summary["results"][0]["cif_path"])
            cif_lines = cif_path.read_text(encoding="utf-8").splitlines()
            self.assertEqual(cif_lines[0], f"# COFid: {summary['results'][0]['metadata']['cofid']} stacking=AA")

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
