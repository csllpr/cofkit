"""Scientific input identity and singular-geometry checks for calculation adapters."""

import json
from pathlib import Path

import pytest

from cofkit.cif_checks import charge_input_geometry_report, read_ordered_structure
from cofkit.graspa import EqeqExecutionError, assign_eqeq_charges_to_cif
from cofkit import lammps as lmp
import test_lammps as fixtures


@pytest.mark.parametrize("duplicate", ["a1 a2 . . 1.0 D", "a2 a1 . . 1.0 D", "a1 a2 1_655 1_655 1.0 D"])
def test_redundant_periodic_edges_are_counted_once(tmp_path, duplicate):
    path = tmp_path / "input.cif"
    path.write_text(fixtures.LammpsTests()._example_cif_text() + duplicate + "\n")
    parsed = lmp._parse_explicit_bond_cif(path)
    assert len(parsed.bonds) == 2
    assert len(parsed.angles) == 1
    assert [b.bond_id for b in parsed.bonds] == [1, 2]
    assert parsed.warnings == ("normalized 1 redundant CIF bond rows representing the same periodic edges",)


def test_translated_reverse_edge_is_one_bond(tmp_path):
    path = tmp_path / "input.cif"
    path.write_text(fixtures.LammpsTests()._periodic_example_cif_text() + "a2 a1 . 1_455 1.0 S\n")
    parsed = lmp._parse_explicit_bond_cif(path)
    assert len(parsed.bonds) == 1
    assert lmp._sub_shift(parsed.bonds[0].shift_2, parsed.bonds[0].shift_1) == (1, 0, 0)
    assert parsed.bonds[0].equilibrium_distance == pytest.approx(1)
    lmp._build_optimization_model(parsed, settings=lmp.LammpsOptimizationSettings())


def test_duplicate_edge_with_conflicting_order_is_rejected(tmp_path):
    path = tmp_path / "input.cif"
    path.write_text(fixtures.LammpsTests()._example_cif_text() + "a2 a1 . . 1.0 S\n")
    with pytest.raises(lmp.LammpsInputError, match="Conflicting bond orders"):
        lmp._parse_explicit_bond_cif(path)


@pytest.mark.parametrize("known_first", [True, False])
def test_unambiguous_duplicate_supplies_a_missing_order(tmp_path, known_first):
    path = tmp_path / "input.cif"
    text = fixtures.LammpsTests()._example_cif_text()
    if known_first:
        text += "a2 a1 . . 1.0 .\n"
    else:
        text = text.replace("a1 a2 . . 1.000000 D", "a1 a2 . . 1.000000 .") + "a2 a1 . . 1.0 D\n"
    path.write_text(text)
    parsed = lmp._parse_explicit_bond_cif(path)
    assert len(parsed.bonds) == 2
    assert parsed.bonds[0].bond_order == 2


@pytest.mark.parametrize("extra", ["a1 a2 . 1_655 11.0 D", "a1 a1 . 1_655 10.0 S"])
def test_distinct_periodic_edges_still_require_a_supercell(tmp_path, extra):
    path = tmp_path / "input.cif"
    path.write_text(fixtures.LammpsTests()._example_cif_text() + extra + "\n")
    with pytest.raises(lmp.LammpsInputError, match="explicit larger supercell"):
        lmp._parse_explicit_bond_cif(path)


@pytest.mark.parametrize("identity", [".", "1_555"])
def test_identity_bond_is_not_replaced_by_a_shorter_periodic_copy(tmp_path, identity):
    path = tmp_path / "input.cif"
    path.write_text(fixtures.LammpsTests()._periodic_example_cif_text().replace(". 1_655 1.000000 S", f". {identity} 9.000000 S"))
    parsed = lmp._parse_explicit_bond_cif(path)
    assert parsed.bonds[0].equilibrium_distance == pytest.approx(9)
    assert parsed.bonds[0].shift_2 == (0, 0, 0)
    with pytest.raises(lmp.LammpsInputError, match="intended edge was not reassigned"):
        lmp._build_optimization_model(parsed, settings=lmp.LammpsOptimizationSettings())


def test_equally_close_bonded_images_require_an_unambiguous_supercell(tmp_path):
    path = tmp_path / "input.cif"
    path.write_text(fixtures.LammpsTests()._periodic_example_cif_text().replace("0.950000", "0.550000").replace("1_655", "."))
    parsed = lmp._parse_explicit_bond_cif(path)
    assert parsed.bonds[0].equilibrium_distance == pytest.approx(5)
    with pytest.raises(lmp.LammpsInputError, match="unambiguous LAMMPS bonded images"):
        lmp._build_optimization_model(parsed, settings=lmp.LammpsOptimizationSettings())


def test_unwrapped_atom_representatives_do_not_limit_image_search_to_one_cell(tmp_path):
    path = tmp_path / "input.cif"
    path.write_text(fixtures.LammpsTests()._periodic_example_cif_text().replace("0.950000", "1.950000").replace("1_655", "1_755"))
    parsed = lmp._parse_explicit_bond_cif(path)
    assert parsed.bonds[0].equilibrium_distance == pytest.approx(1)
    lmp._build_optimization_model(parsed, settings=lmp.LammpsOptimizationSettings())


@pytest.mark.parametrize("position", ["0.100000", "1.100000"])
def test_eqeq_preflight_reports_coincident_atoms_without_running_backend(tmp_path, position):
    path = tmp_path / "input.cif"
    text = fixtures.LammpsTests()._example_cif_text().replace("a2 C 0.200000", f"a2 C {position}")
    path.write_text(text)
    binary = tmp_path / "eqeq"
    binary.write_text("#!/bin/sh\ntouch backend_ran\nexit 1\n")
    binary.chmod(0o755)
    output = tmp_path / "run"
    with pytest.raises(EqeqExecutionError, match="coincident atom sites 'a1' and 'a2'"):
        assign_eqeq_charges_to_cif(path, output_dir=output, eqeq_path=binary)
    report = json.loads((output / "input_geometry_validation.json").read_text())
    assert report["status"] == "invalid"
    assert report["first_collision"]["labels"] == ["a1", "a2"]
    assert report["first_collision"]["distance_angstrom"] < 1e-12
    assert not (output / "backend_ran").exists()
    assert not (output / "execution.json").exists()
    assert (output / path.name).read_bytes() == path.read_bytes()


def test_nonperiodic_charge_preflight_distinguishes_separate_cell_copies(tmp_path):
    path = tmp_path / "input.cif"
    path.write_text(fixtures.LammpsTests()._example_cif_text().replace("a2 C 0.200000", "a2 C 1.100000"))
    _, small = read_ordered_structure(path)
    assert charge_input_geometry_report(small)["status"] == "invalid"
    assert charge_input_geometry_report(small, periodic=False)["status"] == "valid"


def test_coincidence_guard_is_not_a_general_clash_filter(tmp_path):
    path = tmp_path / "input.cif"
    path.write_text(fixtures.LammpsTests()._example_cif_text().replace("a2 C 0.200000", "a2 C 0.101000"))
    _, small = read_ordered_structure(path)
    assert charge_input_geometry_report(small)["status"] == "valid"


def test_successful_eqeq_exposes_its_input_geometry_check(tmp_path):
    path = tmp_path / "input.cif"
    path.write_text(fixtures.LammpsTests()._example_cif_text())
    binary = fixtures.LammpsTests()._write_fake_eqeq_binary(tmp_path / "eqeq")
    result = assign_eqeq_charges_to_cif(path, output_dir=tmp_path / "run", eqeq_path=binary)
    report_path = Path(result.input_geometry_validation_path)
    assert result.to_dict()["input_geometry_validation_path"] == str(report_path)
    assert json.loads(report_path.read_text())["status"] == "valid"


@pytest.mark.parametrize("error_stream", ["stdout", "log"])
def test_lammps_failure_exposes_the_engine_error_and_preserves_logs(tmp_path, error_stream):
    binary = tmp_path / "lmp"
    binary.write_text(
        "#!/usr/bin/env python3\nimport sys\nfrom pathlib import Path\n"
        "error='ERROR on proc 0: Bond atoms 1 2 missing at step 10'\n"
        + ("print(error,flush=True)\n" if error_stream == "stdout" else "Path(sys.argv[sys.argv.index('-log')+1]).write_text(error+'\\n')\n")
        + "sys.stderr.write('MPI_ABORT was invoked\\n')\nraise SystemExit(1)\n"
    )
    binary.chmod(0o755)
    script = tmp_path / "input.in"
    script.write_text("units real\n")
    with pytest.raises(lmp.LammpsExecutionError, match="Bond atoms 1 2 missing at step 10"):
        lmp._run_lammps(binary=binary, input_script_path=script, log_path=tmp_path / "lammps.log",
                        stdout_log_path=tmp_path / "stdout.log", stderr_log_path=tmp_path / "stderr.log",
                        timeout_seconds=30, enable_omp=False)
    assert "MPI_ABORT" in (tmp_path / "stderr.log").read_text()
    command = json.loads((tmp_path / "execution.json").read_text())["command"]
    assert "-nonbuf" in command
    assert "-screen" not in command
