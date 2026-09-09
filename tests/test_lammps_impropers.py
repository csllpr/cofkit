"""Umbrella regressions; set COFKIT_LMP_PATH to run the LAMMPS checks."""

import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pytest


import cofkit.lammps as lammps_module


@pytest.fixture
def formaldehyde(tmp_path):
    # An asymmetric trigonal site catches a duplicated umbrella axis. The
    # central carbon deliberately has atom ID 2, between unlike outer atoms.
    path = tmp_path / "formaldehyde.cif"
    path.write_text(
        """data_formaldehyde
_space_group_name_H-M_alt 'P 1'
_space_group_IT_number 1
_cell_length_a 10
_cell_length_b 10
_cell_length_c 10
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
O1 O 0.620 0.500 0.500 1
C1 C 0.500 0.500 0.500 1
H1 H 0.445 0.595 0.500 1
H2 H 0.435 0.410 0.500 1
loop_
_geom_bond_atom_site_label_1
_geom_bond_atom_site_label_2
_geom_bond_site_symmetry_1
_geom_bond_site_symmetry_2
_geom_bond_distance
_ccdc_geom_bond_type
C1 O1 . . 1.200 D
C1 H1 . . 1.098 S
C1 H2 . . 1.110 S
""",
        encoding="utf-8",
    )
    return lammps_module._parse_explicit_bond_cif(path)


def _section_rows(prepared, name):
    _, sections = lammps_module._split_lammps_data_text(prepared.data_text)
    return [row.split() for row in lammps_module._lammps_section_rows(sections, name)]


def test_dreiding_export_preserves_three_center_first_umbrella_axes(formaldehyde):
    prepared = lammps_module._prepare_dreiding_lammps_system(formaldehyde)
    rows = [[int(value) for value in row] for row in _section_rows(prepared, "Impropers")]
    assert rows == [
        [1, 1, 2, 3, 4, 1],
        [2, 1, 2, 1, 4, 3],
        [3, 1, 2, 1, 3, 4],
    ]
    assert prepared.improper_styles == ("umbrella",)
    assert prepared.n_impropers == 3
    assert prepared.n_improper_types == 1
    coefficients = _section_rows(prepared, "Improper Coeffs")
    assert len(coefficients) == 1
    assert [float(value) for value in coefficients[0]] == pytest.approx([1, 40 / 3, 0])


def test_shared_and_uff_impropers_keep_center_second(formaldehyde):
    expected = [[1, 2, 3, 4], [3, 2, 1, 4], [4, 2, 1, 3]]
    prepared = lammps_module._prepare_uff_lammps_system(formaldehyde)
    assert prepared.improper_styles == ("fourier",)
    assert [
        [int(value) for value in row[2:]] for row in _section_rows(prepared, "Impropers")
    ] == expected
    assert [
        [term.atom_id_1, term.atom_id_2, term.atom_id_3, term.atom_id_4]
        for term in formaldehyde.impropers
    ] == expected


@pytest.fixture(scope="module")
def lammps_binary():
    configured = os.environ.get(lammps_module.COFKIT_LMP_ENV_VAR)
    candidate = configured or shutil.which("lmp") or shutil.which("lmp_mpi")
    if candidate is None:
        pytest.skip("Set COFKIT_LMP_PATH to a LAMMPS binary with the MOLECULE package")
    # An explicitly configured but invalid binary must fail, not silently skip.
    return lammps_module.resolve_lammps_binary(candidate)


def _umbrella_energy(positions):
    # Independent geometry reference: average the three inversion energies
    # around carbon (index 1), using each neighbor as the out-of-plane axis.
    # E = K(1 - cos(omega)), omega0 = 0:
    # https://docs.lammps.org/improper_umbrella.html
    bonds = positions[[0, 2, 3]] - positions[1]
    energy = 0.0
    for axis in range(3):
        plane = [index for index in range(3) if index != axis]
        normal = np.cross(bonds[plane[0]], bonds[plane[1]])
        sin_omega = np.dot(normal, bonds[axis]) / (
            np.linalg.norm(normal) * np.linalg.norm(bonds[axis])
        )
        energy += (40 / 3) * (1 - np.sqrt(1 - sin_omega**2))
    return energy


def _reference_forces(positions):
    forces = np.zeros_like(positions)
    step = 1e-5
    for atom in range(len(positions)):
        for component in range(3):
            plus, minus = positions.copy(), positions.copy()
            plus[atom, component] += step
            minus[atom, component] -= step
            forces[atom, component] = -(_umbrella_energy(plus) - _umbrella_energy(minus)) / (2 * step)
    return forces


def _run_improper_only(prepared, tmp_path, binary, displaced_atom=2, displacement=0.0):
    (tmp_path / "system.data").write_text(prepared.data_text, encoding="utf-8")
    coefficients = _section_rows(prepared, "Improper Coeffs")
    script = [
        "units real",
        "atom_style molecular",
        "boundary f f f",
        "pair_style zero 4.0",
        "bond_style harmonic",
        "angle_style harmonic",
        "improper_style umbrella",
        # Keep the exported topology and improper coefficients; zero all other
        # interactions so every measured energy and force is an improper term.
        "read_data system.data nocoeff",
        "pair_coeff * *",
        "bond_coeff * 0.0 1.0",
        "angle_coeff * 0.0 120.0",
        *("improper_coeff " + " ".join(row) for row in coefficients),
        f"group displaced id {displaced_atom}",
        f"displace_atoms displaced move 0 0 {displacement} units box",
        "thermo_style custom step pe eimp",
        "run 0",
        'print "IMPROPER_ENERGY $(eimp:%.16g)"',
        "write_dump all custom forces.dump id x y z fx fy fz modify sort id format float %.16g",
    ]
    (tmp_path / "run.in").write_text("\n".join(script) + "\n", encoding="utf-8")
    result = subprocess.run(
        [str(binary), "-in", "run.in", "-log", "none"],
        cwd=tmp_path,
        env={**os.environ, "OMP_NUM_THREADS": "1"},
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    energy = next(
        float(line.split()[1]) for line in result.stdout.splitlines() if line.startswith("IMPROPER_ENERGY ")
    )
    atoms = np.loadtxt(tmp_path / "forces.dump", skiprows=9)
    np.testing.assert_array_equal(atoms[:, 0], [1, 2, 3, 4])
    return energy, atoms[:, 1:4], atoms[:, 4:7]


def test_dreiding_planar_improper_energy_and_forces(formaldehyde, tmp_path, lammps_binary):
    prepared = lammps_module._prepare_dreiding_lammps_system(formaldehyde)
    energy, positions, forces = _run_improper_only(prepared, tmp_path, lammps_binary)
    assert energy == pytest.approx(0.0, abs=1e-10)
    np.testing.assert_allclose(forces, 0.0, atol=1e-10)
    assert energy == pytest.approx(_umbrella_energy(positions), abs=1e-10)


@pytest.mark.parametrize("displaced_atom", [1, 2, 3, 4])
@pytest.mark.parametrize("displacement", [-0.15, 0.15])
def test_dreiding_out_of_plane_improper_energy_and_forces(
    formaldehyde, tmp_path, lammps_binary, displaced_atom, displacement
):
    prepared = lammps_module._prepare_dreiding_lammps_system(formaldehyde)
    energy, positions, forces = _run_improper_only(
        prepared, tmp_path, lammps_binary, displaced_atom, displacement
    )
    assert energy > 0.0
    assert energy == pytest.approx(_umbrella_energy(positions), rel=1e-7, abs=1e-9)
    np.testing.assert_allclose(forces, _reference_forces(positions), rtol=1e-6, atol=1e-7)
    assert forces[displaced_atom - 1, 2] * displacement < 0.0
    np.testing.assert_allclose(forces.sum(axis=0), 0.0, atol=1e-10)
