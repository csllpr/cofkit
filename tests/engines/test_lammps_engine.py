"""Real-engine checks; the dedicated CI job requires the configured executable."""

import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest


@pytest.fixture
def lmp():
    binary = os.environ.get("COFKIT_TEST_LMP")
    if not binary:
        if os.environ.get("COFKIT_REQUIRE_LAMMPS") == "1":
            pytest.fail("LAMMPS qualification requested but COFKIT_TEST_LMP is missing")
        pytest.skip("Real LAMMPS evidence requires COFKIT_TEST_LMP")
    resolved = shutil.which(binary)
    assert resolved, f"Required LAMMPS executable is unavailable: {binary}"
    return resolved


def test_lennard_jones_minimum(lmp, tmp_path):
    script = tmp_path / "pair.in"
    script.write_text(f"""units real
atom_style atomic
boundary p p p
region cell block 0 20 0 20 0 20
create_box 1 cell
create_atoms 1 single 5 5 5
create_atoms 1 single {5 + 2 ** (1 / 6):.16g} 5 5
mass 1 12
pair_style lj/cut 2.5
pair_coeff 1 1 1 1
run 0
variable energy equal pe
print "COFKIT_PAIR_ENERGY ${{energy}}"
""")
    result = subprocess.run(
        [lmp, "-in", str(script)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    values = re.findall(
        r"^COFKIT_PAIR_ENERGY ([-+\deE.]+)$", result.stdout, re.MULTILINE
    )
    assert len(values) == 1
    assert float(values[0]) == pytest.approx(-1, abs=1e-10)


def test_optimizer_adapter_convergence(lmp, tmp_path):
    from cofkit.lammps import LammpsOptimizationSettings, optimize_cif_with_lammps

    cif = tmp_path / "input.cif"
    cif.write_text(
        (
            Path(__file__).resolve().parents[1] / "fixtures" / "engine_molecule.cif"
        ).read_text()
    )
    result = optimize_cif_with_lammps(
        cif,
        lmp_path=lmp,
        settings=LammpsOptimizationSettings(
            enable_omp=False,
            charge_model="none",
            pre_minimization_steps=0,
            relax_cell=False,
            max_iterations=10000,
            max_evaluations=100000,
        ),
    )
    assert result.convergence["converged"]
    assert result.convergence["stages"][-1]["force_norm_kcal_mol_angstrom"] <= 1e-6


def test_dreiding_hbond_optimizer_run(lmp, tmp_path):
    import math

    from cofkit.lammps import LammpsOptimizationSettings, optimize_cif_with_lammps

    cif = tmp_path / "input.cif"
    cif.write_text(
        (
            Path(__file__).resolve().parents[1] / "fixtures" / "engine_methylamine.cif"
        ).read_text()
    )
    result = optimize_cif_with_lammps(
        cif,
        lmp_path=lmp,
        settings=LammpsOptimizationSettings(
            enable_omp=False,
            charge_model="none",
            pre_minimization_steps=0,
            relax_cell=False,
            force_tolerance=1e-4,
            max_iterations=10000,
            max_evaluations=100000,
        ),
    )
    assert result.convergence["converged"]
    assert "H__HB" in set(result.atom_type_symbols.values())
    script_text = Path(result.lammps_input_script_path).read_text()
    assert "hbond/dreiding/lj" in script_text
    log_text = Path(result.lammps_log_path).read_text()
    energy_blocks = re.findall(
        r"Energy initial, next-to-last, final =\s*\n\s*([-+\deE. ]+)\n", log_text
    )
    assert energy_blocks
    energies = [float(value) for block in energy_blocks for value in block.split()]
    assert energies
    assert all(math.isfinite(value) for value in energies)


def test_optimizer_rejects_real_iteration_exhaustion(lmp, tmp_path):
    from cofkit.lammps import (
        LammpsExecutionError,
        LammpsOptimizationSettings,
        optimize_cif_with_lammps,
    )

    cif = tmp_path / "input.cif"
    cif.write_text(
        (
            Path(__file__).resolve().parents[1] / "fixtures" / "engine_molecule.cif"
        ).read_text()
    )
    with pytest.raises(LammpsExecutionError, match="unconverged"):
        optimize_cif_with_lammps(
            cif,
            lmp_path=lmp,
            settings=LammpsOptimizationSettings(
                enable_omp=False,
                charge_model="none",
                pre_minimization_steps=0,
                relax_cell=False,
                two_stage_protocol=False,
                position_restraint_force_constant=0,
                max_iterations=1,
                max_evaluations=10,
                force_tolerance=1e-12,
            ),
        )
    assert not list(tmp_path.rglob("*_optimized.cif"))
