"""CLI exceptions, stderr and result exit status form one public contract."""

import json
import subprocess
import sys
import sysconfig
from pathlib import Path
from types import SimpleNamespace

import pytest

from cofkit import cli, cli_analyze, cli_build, cli_validate
from cofkit.decompose import CifDecompositionResult


@pytest.mark.parametrize("args, message", [
    (["build", "single-pair", "--cofid", "nonsense"], "COFid must have exactly 3 blocks"),
    (["build", "ring-forming", "--cofid", "nonsense"], "COFid must have exactly 3 blocks"),
    (["analyze", "classify-output", "missing"], "expected manifest"),
    (["build", "batch-binary-bridge", "--input-dir", "missing"], "input directory does not exist"),
])
def test_expected_errors_have_numeric_exit_and_stderr(tmp_path, monkeypatch, capsys, args, message):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit) as error:
        cli.main(args)
    assert error.value.code == 1
    output = capsys.readouterr()
    assert output.out == ""
    assert message in output.err
    assert output.err.startswith("error: ")
    assert "Traceback" not in output.err


@pytest.mark.parametrize("position", ["root", "group", "command", "legacy"])
def test_debug_preserves_original_wrapped_exception(monkeypatch, position):
    def fail(args):
        try:
            raise ValueError("original cause")
        except ValueError as exc:
            raise SystemExit("error: wrapper") from exc

    monkeypatch.setattr(cli_build, "_run_list_templates", fail)
    arguments = {
        "root": ["--debug", "build", "list-templates"],
        "group": ["build", "--debug", "list-templates"],
        "command": ["build", "list-templates", "--debug"],
        "legacy": ["--debug", "list-templates"],
    }[position]
    with pytest.raises(ValueError, match="original cause"):
        cli.main(arguments)


def test_unexpected_error_is_explicit_and_debug_reraises(monkeypatch, capsys):
    def fail(args):
        raise KeyError("unexpected")

    monkeypatch.setattr(cli_build, "_run_list_templates", fail)
    with pytest.raises(SystemExit) as error:
        cli.main(["build", "list-templates"])
    assert error.value.code == 1
    output = capsys.readouterr()
    assert "KeyError" in output.err and "--debug" in output.err
    with pytest.raises(KeyError):
        cli.main(["--debug", "build", "list-templates"])


def test_dotenv_failure_is_also_handled(monkeypatch, capsys):
    def fail():
        raise PermissionError("unreadable .env")

    monkeypatch.setattr(cli, "_load_dotenv", fail)
    with pytest.raises(SystemExit) as error:
        cli.main(["--version"])
    assert error.value.code == 1
    assert "PermissionError: unreadable .env" in capsys.readouterr().err


@pytest.mark.parametrize("status", ["success", "skipped", "ambiguous", "error"])
@pytest.mark.parametrize("json_output", [False, True])
def test_decomposition_status_is_independent_of_output_format(status, json_output, monkeypatch, capsys):
    result = CifDecompositionResult(status=status, input_cif="input.cif", topology="hcb", linkage="imine",
                                    cofid="result" if status == "success" else None, reason="diagnostic")
    monkeypatch.setattr(cli_analyze, "decompose_cif_to_cofid", lambda *a, **kw: result)
    args = ["analyze", "decompose", "input.cif"] + (["--json"] if json_output else [])
    if result.ok:
        cli.main(args)
    else:
        with pytest.raises(SystemExit) as error:
            cli.main(args)
        assert error.value.code == 1
    output = capsys.readouterr()
    if json_output:
        assert json.loads(output.out)["status"] == status
    elif result.ok:
        assert output.out.strip() == "result"
    else:
        assert output.out == "" and "diagnostic" in output.err


@pytest.mark.parametrize("json_output", [False, True])
def test_failed_validation_uses_command_failure_status(json_output, monkeypatch, capsys):
    monkeypatch.setattr(cli_validate, "_print_human_result", lambda result: print("failed"))
    result = SimpleNamespace(ok=False, status="failed", to_dict=lambda: {"status": "failed"})
    with pytest.raises(SystemExit) as error:
        cli_validate._emit_result(result, json_output=json_output)
    assert error.value.code == 1
    assert capsys.readouterr().out


def test_repair_warning_goes_to_stderr(capsys):
    summary = SimpleNamespace(structure_id="broken", metadata={"validation": {"geometry_repair": {"status": "failed", "error": "cause"}}})
    cli_build._print_single_pair_geometry_repair_warning(summary)
    output = capsys.readouterr()
    assert output.out == ""
    assert "geometry repair failed" in output.err


@pytest.mark.parametrize("module", ["cofkit", "cofkit.cli"])
def test_module_entrypoints_match(module, tmp_path):
    result = subprocess.run([sys.executable, "-m", module, "--version"], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0
    assert result.stdout.startswith("cofkit ") and not result.stderr
    bad = subprocess.run([sys.executable, "-m", module, "build", "single-pair", "--cofid", "bad"], cwd=tmp_path, capture_output=True, text=True)
    assert bad.returncode == 1
    assert "COFid must" in bad.stderr and "Traceback" not in bad.stderr
    debug = subprocess.run([sys.executable, "-m", module, "--debug", "build", "single-pair", "--cofid", "bad"], cwd=tmp_path, capture_output=True, text=True)
    assert debug.returncode == 1
    assert "Traceback" in debug.stderr and "ValueError" in debug.stderr


def test_console_script_uses_same_boundary(tmp_path):
    executable = Path(sysconfig.get_path("scripts")) / ("cofkit.exe" if sys.platform == "win32" else "cofkit")
    result = subprocess.run([str(executable), "build", "single-pair", "--cofid", "bad"], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 1
    assert "COFid must" in result.stderr and "Traceback" not in result.stderr


def test_usage_error_remains_argparse_status(capsys):
    with pytest.raises(SystemExit) as error:
        cli.main(["build", "single-pair", "--unknown"])
    assert error.value.code == 2
    assert "unrecognized arguments" in capsys.readouterr().err


def test_interrupt_uses_shell_status(monkeypatch, capsys):
    def interrupt(args):
        raise KeyboardInterrupt

    monkeypatch.setattr(cli_build, "_run_list_templates", interrupt)
    with pytest.raises(SystemExit) as error:
        cli.main(["build", "list-templates"])
    assert error.value.code == 130
    assert capsys.readouterr().err == "error: interrupted\n"
