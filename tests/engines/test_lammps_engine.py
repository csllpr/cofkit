"""Real-engine checks; the dedicated CI job requires the configured executable."""

import json
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


def test_uff_periodic_angle_energy_matches_reference(lmp, tmp_path):
    import math

    from cofkit.lammps import (
        _CifAngleRecord,
        _compute_uff_angle_coefficients,
        _load_uff_parameters,
    )

    parameters = _load_uff_parameters()
    r_ab = 1.40
    r_bc = 1.45
    theta_degrees = 100.0
    for center_type, flank_type in (("C_2", "C_R"), ("Bi3+3", "F_")):
        coefficients = _compute_uff_angle_coefficients(
            angle=_CifAngleRecord(
                angle_id=1,
                atom_id_1=1,
                atom_id_2=2,
                atom_id_3=3,
                equilibrium_degrees=parameters[center_type].theta0,
            ),
            atom_type_by_atom_id={1: flank_type, 2: center_type, 3: flank_type},
            parameters=parameters,
            bond_reference={(1, 2): (0.0, r_ab), (2, 3): (0.0, r_bc)},
        )
        style, c_coeff, b_coeff, n_coeff = coefficients
        assert style == "cosine/periodic"
        theta = math.radians(theta_degrees)
        x3 = r_bc * math.cos(theta)
        y3 = r_bc * math.sin(theta)
        data_file = tmp_path / f"angle_{center_type.replace('+', 'p')}.data"
        data_file.write_text(f"""cofkit UFF periodic angle check

3 atoms
1 angles

1 atom types
1 angle types

-15.0 15.0 xlo xhi
-15.0 15.0 ylo yhi
-15.0 15.0 zlo zhi

Masses

1 12.011

Angle Coeffs

1 {c_coeff:.16g} {b_coeff} {n_coeff}

Atoms # molecular

1 1 1 {r_ab:.16g} 0.0 0.0
2 1 1 0.0 0.0 0.0
3 1 1 {x3:.16g} {y3:.16g} 0.0

Angles

1 1 1 2 3
""")
        script = tmp_path / f"angle_{center_type.replace('+', 'p')}.in"
        script.write_text(f"""units real
atom_style molecular
boundary f f f
pair_style none
angle_style cosine/periodic
read_data {data_file}
run 0
variable energy equal eangle
print "COFKIT_ANGLE_ENERGY ${{energy}}"
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
            r"^COFKIT_ANGLE_ENERGY ([-+\deE.]+)$", result.stdout, re.MULTILINE
        )
        assert len(values) == 1
        # Independent reference expression: E = ka / n^2 *
        # (1 - cos(n*theta0) * cos(n*theta)) with ka from the UFF paper
        # (Rappe et al., JACS 1992, 114, 10024-10035, eq. 10).
        theta0 = math.radians(parameters[center_type].theta0)
        cos_theta0 = math.cos(theta0)
        r_ac = math.sqrt(r_ab**2 + r_bc**2 - 2.0 * r_ab * r_bc * cos_theta0)
        beta = 664.12 / (r_ab * r_bc)
        ka = beta * (parameters[flank_type].z1**2 / (r_ac**5.0)) * r_ab * r_bc
        ka *= 3.0 * r_ab * r_bc * (1.0 - cos_theta0 * cos_theta0) - r_ac * r_ac * cos_theta0
        expected = (ka / n_coeff**2) * (
            1.0 - math.cos(n_coeff * theta0) * math.cos(math.radians(n_coeff * theta_degrees))
        )
        assert float(values[0]) == pytest.approx(expected, rel=1e-8, abs=1e-10)


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
    assert result.convergence["stages"][-1]["force_norm_kcal_mol_angstrom"] <= 1e-4


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


def test_optimizer_warns_on_real_iteration_exhaustion(lmp, tmp_path):
    from cofkit.lammps import (
        LammpsOptimizationSettings,
        optimize_cif_with_lammps,
    )

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
            two_stage_protocol=False,
            position_restraint_force_constant=0,
            max_iterations=1,
            max_evaluations=10,
            force_tolerance=1e-12,
        ),
    )
    assert not result.convergence["converged"]
    assert any("unconverged" in warning for warning in result.warnings)
    optimized = Path(result.optimized_cif)
    assert optimized.is_file()
    assert "treat this structure as unconverged" in optimized.read_text()
    report = json.loads(Path(result.report_path).read_text())
    assert report["convergence"]["converged"] is False
    assert any("unconverged" in warning for warning in report["warnings"])
    convergence_json = json.loads((Path(result.output_dir) / "convergence.json").read_text())
    assert convergence_json["converged"] is False


def test_lammps_geometry_repair_reports_final_coordinate_measurements(lmp, tmp_path):
    """A04: the LAMMPS geometry-repair route must report linkage measurements
    recomputed from the repaired CIF, not the pre-repair seed metrics."""
    pytest.importorskip("rdkit")
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from test_batch import (  # noqa: E402
        TAPB,
        TEREPHTHALALDEHYDE,
        _DistortingCifWriter,
    )
    from cofkit import (  # noqa: E402
        BatchGenerationConfig,
        BatchMonomerRecord,
        BatchStructureGenerator,
        CoarseValidationThresholds,
    )

    generator = BatchStructureGenerator(
        BatchGenerationConfig(
            rdkit_num_conformers=1,
            retain_top_results=1,
            single_node_topology_ids=("hcb",),
            hard_hard_max_bridge_distance=10.0,
            repair_geometry=True,
            repair_geometry_lmp_path=lmp,
            validation_thresholds=CoarseValidationThresholds(hard_hard_max_bridge_distance=10.0),
        )
    )
    # Distort the exported linkage so the seed metrics are clean but the
    # final coordinates are not; repair then runs for real.
    generator.cif_writer = _DistortingCifWriter(generator.cif_writer, {})
    amine = BatchMonomerRecord(
        id="tapb", name="tapb", smiles=TAPB, motif_kind="amine", expected_connectivity=3
    )
    aldehyde = BatchMonomerRecord(
        id="tpal",
        name="tpal",
        smiles=TEREPHTHALALDEHYDE,
        motif_kind="aldehyde",
        expected_connectivity=2,
    )

    summary, _candidate = generator.generate_pair_candidate(
        amine, aldehyde, out_dir=tmp_path / "run", write_cif=True
    )

    assert summary.status == "ok"
    validation = summary.metadata["validation"]
    assert validation["classification"] == "needs_optimization"
    assert validation["metrics"]["bridge_metrics_source"] == "final_cif"
    repair = validation["geometry_repair"]
    assert repair["status"] == "ok"
    post = repair["post_repair_validation"]
    assert post["scope"] == "optimized_cif_final_geometry_checks"
    assert post["metrics"]["bridge_metrics_source"] == "final_cif"
    assert post["metrics"]["n_measured_bridge_bonds"] > 0
    assert "max_bridge_distance_residual" in post["metrics"]
    assert post["coverage"]["bridge_geometry"] == "measured"


def _fixture_text(name):
    return (Path(__file__).resolve().parents[1] / "fixtures" / name).read_text()


def test_soft_pre_minimization_skips_clash_free_input(lmp, tmp_path):
    from cofkit.lammps import LammpsOptimizationSettings, optimize_cif_with_lammps

    cif = tmp_path / "input.cif"
    cif.write_text(_fixture_text("engine_molecule.cif"))
    result = optimize_cif_with_lammps(
        cif,
        lmp_path=lmp,
        settings=LammpsOptimizationSettings(
            enable_omp=False,
            charge_model="none",
            pre_minimization_mode="soft",
            relax_cell=False,
        ),
    )
    soft = result.convergence["soft_pre_minimization"]
    assert soft["initial_contact_energy_kcal_mol"] == 0
    assert soft["ramp_ran"] is False
    assert result.convergence["converged"]


def test_soft_pre_minimization_separates_clashing_molecules(lmp, tmp_path):
    from cofkit.lammps import LammpsOptimizationSettings, optimize_cif_with_lammps

    # Two methylamine copies whose carbons sit 0.84 angstrom apart: far
    # inside the 2.55 angstrom C...C clash cutoff.
    text = _fixture_text("engine_methylamine.cif")
    atoms, bonds = text.split("\nloop_\n_geom_bond_atom_site_label_1", 1)
    copies = []
    for line in atoms.splitlines():
        parts = line.split()
        if len(parts) == 6 and parts[1] in {"C", "N", "H"}:
            parts[0] = parts[0] + "b"
            parts[4] = f"{float(parts[4]) + 0.07:.6f}"
            copies.append(" ".join(parts))
    bond_rows = [line for line in bonds.splitlines() if line.endswith(" S")]
    copy_bonds = [
        " ".join([row.split()[0] + "b", row.split()[1] + "b", *row.split()[2:]]) for row in bond_rows
    ]
    cif = tmp_path / "clash.cif"
    cif.write_text(
        atoms + "\n" + "\n".join(copies) + "\n\nloop_\n_geom_bond_atom_site_label_1" + bonds + "\n".join(copy_bonds) + "\n"
    )
    result = optimize_cif_with_lammps(
        cif,
        lmp_path=lmp,
        settings=LammpsOptimizationSettings(
            enable_omp=False,
            charge_model="none",
            pre_minimization_mode="soft",
            relax_cell=False,
            min_style="cg",
        ),
    )
    soft = result.convergence["soft_pre_minimization"]
    assert soft["ramp_ran"] is True
    assert soft["initial_contact_energy_kcal_mol"] > 0
    assert soft["residual_contact_energy_kcal_mol"] < 1e-2 * soft["initial_contact_energy_kcal_mol"]
    assert result.convergence["converged"]


def test_cell_relaxation_cycles_run_in_real_lammps(lmp, tmp_path):
    from cofkit.lammps import LammpsOptimizationSettings, optimize_cif_with_lammps

    cif = tmp_path / "input.cif"
    cif.write_text(_fixture_text("engine_methylamine.cif"))
    result = optimize_cif_with_lammps(
        cif,
        lmp_path=lmp,
        settings=LammpsOptimizationSettings(
            enable_omp=False,
            charge_model="none",
            pre_minimization_steps=0,
            cell_relaxation_max_cycles=2,
        ),
    )
    cell = result.convergence["cell_relaxation"]
    assert 1 <= cell["cycles_run"] <= 2
    labels = [stage["stage"] for stage in result.convergence["stages"]]
    assert labels[-1] == f"fixed_cell_{cell['cycles_run']}"
    assert labels.count("box_relax_1") == 1
    # A cell that was stable after one cycle must have left the loop.
    if cell["cycles"][0]["max_box_change_angstrom"] <= cell["stability_tolerance_angstrom"]:
        assert cell["cycles_run"] == 1


def test_dsf_coulomb_runs_charged_optimization(lmp, tmp_path):
    import math

    from cofkit.lammps import (
        LammpsOptimizationSettings,
        _convergence_diagnostics_available,
        optimize_cif_with_lammps,
    )

    text = _fixture_text("engine_methylamine.cif").replace(
        "_atom_site_occupancy\n", "_atom_site_occupancy\n_atom_site_charge\n"
    )
    charges = {"c1": -0.2, "n1": -0.5, "h1": 0.2, "h2": 0.2, "h3": 0.1, "h4": 0.1, "h5": 0.1}
    lines = []
    for line in text.splitlines():
        parts = line.split()
        if len(parts) == 6 and parts[0] in charges and parts[5] == "1.00":
            line = f"{line} {charges[parts[0]]:.6f}"
        lines.append(line)
    cif = tmp_path / "charged.cif"
    cif.write_text("\n".join(lines) + "\n")
    result = optimize_cif_with_lammps(
        cif,
        lmp_path=lmp,
        settings=LammpsOptimizationSettings(
            enable_omp=False,
            charge_model="none",
            pre_minimization_steps=0,
            relax_cell=False,
        ),
    )
    script_text = Path(result.lammps_input_script_path).read_text()
    assert "lj/cut/coul/dsf" in script_text
    assert "kspace_style" not in script_text
    # An isolated charged molecule has flat modes that stall both Coulomb
    # methods; this checks that the DSF export runs and reports finite,
    # readable diagnostics.
    assert _convergence_diagnostics_available(result.convergence)
    assert all(
        math.isfinite(stage["force_norm_kcal_mol_angstrom"]) for stage in result.convergence["stages"]
    )
