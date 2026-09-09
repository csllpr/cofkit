import math
from unittest.mock import Mock, patch

import pytest

from cofkit.chem.rdkit import _optimize_conformers
from cofkit.graspa import (
    GraspaIsothermSettings,
    GraspaMixtureSettings,
    GraspaWidomSettings,
    _render_isotherm_simulation_input,
)
from cofkit.lammps import (
    LammpsOptimizationSettings,
    LammpsParseError,
    _build_minimization_stages,
    _check_minimization_convergence,
    _parse_lammps_dump_last_frame,
    _validate_settings,
)


def test_sampling_defaults_distinguish_adsorption_and_widom():
    assert not GraspaIsothermSettings().use_max_step
    assert not GraspaMixtureSettings(components=()).use_max_step
    assert GraspaWidomSettings().use_max_step
    for backend in ("graspa", "raspa2"):
        settings = GraspaIsothermSettings(backend=backend)
        text = _render_isotherm_simulation_input(
            settings, unit_cells=(2, 2, 2), pressure=1e6
        )
        if backend == "graspa":
            assert "UseMaxStep  no" in text
            assert "FugacityCoefficient      PR-EOS" in text
        else:
            assert (
                "FugacityCoefficient" not in text
            )  # native PR-EOS, explicit coefficients override
        zero = _render_isotherm_simulation_input(
            settings, unit_cells=(2, 2, 2), pressure=0
        )
        assert "FugacityCoefficient      1" in zero


def test_zero_tolerances_survive_stage_overrides():
    settings = LammpsOptimizationSettings(
        energy_tolerance=1e-5, stage2_energy_tolerance=0, box_relax_energy_tolerance=0
    )
    _validate_settings(settings)
    assert [s.energy_tolerance for s in _build_minimization_stages(settings)] == [
        1e-5,
        0,
        0,
    ]


@pytest.mark.parametrize(
    "reason,force",
    [
        ("max iterations", 0),
        ("linesearch alpha is zero", 0),
        ("energy tolerance", 1),
        ("force tolerance", math.nan),
    ],
)
def test_inadequate_minimization_evidence_is_not_convergence(tmp_path, reason, force):
    path = tmp_path / "log"
    path.write_text(
        f"COFKIT_STAGE stage1\nStopping criterion = {reason}\nCOFKIT_MINIMUM stage1 {force} {force} 0 0 0 0 0 0\n"
    )
    result = _check_minimization_convergence(
        path,
        LammpsOptimizationSettings(two_stage_protocol=False, relax_cell=False),
        None,
    )
    assert not result["converged"]


def test_missing_convergence_and_excess_stress_are_rejected(tmp_path):
    path = tmp_path / "log"
    path.write_text("")
    assert not _check_minimization_convergence(
        path, LammpsOptimizationSettings(), None
    )["converged"]
    path.write_text(
        "COFKIT_STAGE box_relax\nStopping criterion = force tolerance\nCOFKIT_MINIMUM box_relax 0 0 10 10 10 0 0 0\n"
    )
    assert not _check_minimization_convergence(
        path, LammpsOptimizationSettings(box_relax_mode="iso"), None
    )["converged"]


def test_rdkit_nonconverged_conformer_is_never_optimized():
    field = Mock()
    field.Minimize.return_value = 1
    field.CalcEnergy.return_value = -10
    with (
        patch("cofkit.chem.rdkit.AllChem.MMFFHasAllMoleculeParams", return_value=False),
        patch("cofkit.chem.rdkit.AllChem.UFFGetMoleculeForceField", return_value=field),
    ):
        with pytest.raises(ValueError, match="did not converge"):
            _optimize_conformers(Mock(), (0,))


def test_streaming_dump_rejects_truncated_final_frame(tmp_path):
    path = tmp_path / "dump"
    frame = "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\nITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\nITEM: ATOMS id x y z ix iy iz\n1 1 1 1 0 0 0\n"
    path.write_text(frame * 1000)
    result = _parse_lammps_dump_last_frame(path, expected_atoms=1)
    assert result.cartesian_positions == {1: (1, 1, 1)}
    with path.open("a") as handle:
        handle.write("ITEM: TIMESTEP\n")
    with pytest.raises(LammpsParseError, match="Truncated"):
        _parse_lammps_dump_last_frame(path, expected_atoms=1)
