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
    _convergence_diagnostics_available,
    _parse_lammps_dump_last_frame,
    _unconverged_optimization_warning,
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
        energy_tolerance=1e-5,
        stage2_energy_tolerance=0,
        box_relax_energy_tolerance=0,
        cell_relaxation_max_cycles=1,
    )
    _validate_settings(settings)
    assert [s.energy_tolerance for s in _build_minimization_stages(settings)] == [
        1e-5,
        0,
        0,
        0,
    ]


def test_force_tolerance_defaults_loosen_for_the_box_relax_stage():
    settings = LammpsOptimizationSettings()
    assert settings.force_tolerance == 1e-4
    assert settings.box_relax_force_tolerance == 1e-3
    _validate_settings(settings)
    stages = _build_minimization_stages(settings)
    cycles = settings.cell_relaxation_max_cycles
    assert [s.label for s in stages] == ["stage1", "stage2"] + [
        label for cycle in range(1, cycles + 1) for label in (f"box_relax_{cycle}", f"fixed_cell_{cycle}")
    ]
    # Each cycle's fixed-cell minimization keeps the final-stage criteria.
    assert [s.force_tolerance for s in stages] == [1e-4, 1e-4] + [1e-3, 1e-4] * cycles
    assert [s.relax_cell for s in stages] == [False, False] + [True, False] * cycles
    inherited = LammpsOptimizationSettings(box_relax_force_tolerance=None, cell_relaxation_max_cycles=1)
    assert [s.force_tolerance for s in _build_minimization_stages(inherited)] == [
        1e-4,
        1e-4,
        1e-4,
        1e-4,
    ]


def _stage_log(label, reason="force tolerance", force=0, stress=(0, 0, 0, 0, 0, 0)):
    return (
        f"COFKIT_STAGE {label}\nStopping criterion = {reason}\n"
        f"COFKIT_MINIMUM {label} {force} {force} {' '.join(map(str, stress))}\n"
    )


def _cell_row(cycle, start, end):
    return f"COFKIT_CELL_CYCLE {cycle} {' '.join(map(str, start))} {' '.join(map(str, end))}\n"


_BOX = (20.0, 20.0, 7.0, 0.0, 0.0, 0.0)


def _cell_relax_log(*cycle_ends, box_reason="linesearch alpha is zero", final_stress=(0, 0, 0, 0, 0, 0)):
    text = _stage_log("stage1") + _stage_log("stage2")
    start = _BOX
    for cycle, end in enumerate(cycle_ends, start=1):
        text += _stage_log(f"box_relax_{cycle}", reason=box_reason, force=5)
        stress = final_stress if cycle == len(cycle_ends) else (0, 0, 0, 0, 0, 0)
        text += _stage_log(f"fixed_cell_{cycle}", stress=stress)
        text += _cell_row(cycle, start, end)
        start = end
    return text


def test_cell_relaxation_accepts_a_stable_cell_despite_stalled_box_relax(tmp_path):
    path = tmp_path / "log"
    moved = (19.5, 20.2, 7.2, 0.0, 0.0, 0.0)
    settled = (19.505, 20.2, 7.2, 0.0, 0.0, 0.0)
    path.write_text(_cell_relax_log(moved, settled))
    result = _check_minimization_convergence(path, LammpsOptimizationSettings(box_relax_mode="aniso"), None)
    assert result["converged"], result
    assert result["final_stage"] == "fixed_cell_2"
    assert [s["stage"] for s in result["stages"]] == [
        "stage1", "stage2", "box_relax_1", "fixed_cell_1", "box_relax_2", "fixed_cell_2"
    ]
    assert result["stages"][2]["stopping_criterion"] == "linesearch alpha is zero"
    cell = result["cell_relaxation"]
    assert cell["cycles_run"] == 2
    assert [round(c["max_box_change_angstrom"], 6) for c in cell["cycles"]] == [0.5, 0.005]
    assert cell["cycles"][-1]["box_angstrom"]["lx"] == 19.505


def test_cell_relaxation_rejects_a_cell_still_changing_after_the_last_cycle(tmp_path):
    path = tmp_path / "log"
    ends = [(20.0 - 0.1 * k, 20.0, 7.0, 0.0, 0.0, 0.0) for k in range(1, 4)]
    path.write_text(_cell_relax_log(*ends))
    result = _check_minimization_convergence(
        path, LammpsOptimizationSettings(box_relax_mode="aniso", cell_relaxation_max_cycles=3), None
    )
    assert not result["converged"]
    assert result["criteria"] == {
        "final_stage_termination_and_force": True,
        "cell_stable": False,
        "residual_stress": True,
    }
    assert _convergence_diagnostics_available(result)
    warning = _unconverged_optimization_warning(result)
    assert "cell still changing after 3 of 3 cell-relaxation cycles" in warning


def test_cell_relaxation_rejects_excess_residual_stress(tmp_path):
    path = tmp_path / "log"
    path.write_text(_cell_relax_log(_BOX, final_stress=(0, 0, 0, 0, 0, 0)))
    settings = LammpsOptimizationSettings(box_relax_mode="tri")
    assert _check_minimization_convergence(path, settings, None)["converged"]
    tolerance = settings.pressure_tolerance
    path.write_text(_cell_relax_log(_BOX, final_stress=(0, 0, 0, 0, 0, 2 * tolerance)))
    result = _check_minimization_convergence(path, settings, None)
    assert not result["converged"]
    assert result["criteria"]["residual_stress"] is False
    assert result["cell_relaxation"]["max_residual_stress_atm"] == 2 * tolerance
    assert "residual stress" in _unconverged_optimization_warning(result)
    # aniso does not relax shear, so shear stress is not part of acceptance.
    aniso = _check_minimization_convergence(path, LammpsOptimizationSettings(box_relax_mode="aniso"), None)
    assert aniso["converged"]


def test_cell_relaxation_without_cell_rows_is_unreadable(tmp_path):
    path = tmp_path / "log"
    text = _stage_log("stage1") + _stage_log("stage2") + _stage_log("box_relax_1") + _stage_log("fixed_cell_1")
    path.write_text(text)
    result = _check_minimization_convergence(path, LammpsOptimizationSettings(box_relax_mode="aniso"), None)
    assert not result["converged"]
    assert result["cell_relaxation"]["reason"] == "missing_cell_metrics"
    assert not _convergence_diagnostics_available(result)


def test_cell_relaxation_requires_the_final_fixed_cell_stage_to_converge(tmp_path):
    path = tmp_path / "log"
    text = (
        _stage_log("stage1") + _stage_log("stage2") + _stage_log("box_relax_1", force=5)
        + _stage_log("fixed_cell_1", reason="max iterations", force=1) + _cell_row(1, _BOX, _BOX)
    )
    path.write_text(text)
    result = _check_minimization_convergence(path, LammpsOptimizationSettings(box_relax_mode="aniso"), None)
    assert not result["converged"]
    assert result["criteria"]["final_stage_termination_and_force"] is False
    assert "final stage fixed_cell_1 stopped on 'max iterations'" in _unconverged_optimization_warning(result)


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
    path.write_text(_cell_relax_log(_BOX, final_stress=(10, 10, 10, 0, 0, 0)))
    assert not _check_minimization_convergence(
        path, LammpsOptimizationSettings(box_relax_mode="iso", pressure_tolerance=1.0), None
    )["converged"]


@pytest.mark.parametrize(
    "reason",
    ["energy tolerance", "max iterations", "linesearch alpha is zero"],
)
def test_unconverged_with_readable_diagnostics_downgrades_to_warning(tmp_path, reason):
    path = tmp_path / "log"
    path.write_text(
        f"COFKIT_STAGE stage1\nStopping criterion = {reason}\nCOFKIT_MINIMUM stage1 5 5 0 0 0 0 0 0\n"
    )
    convergence = _check_minimization_convergence(
        path,
        LammpsOptimizationSettings(two_stage_protocol=False, relax_cell=False),
        None,
    )
    assert not convergence["converged"]
    assert _convergence_diagnostics_available(convergence)
    warning = _unconverged_optimization_warning(convergence)
    assert "unconverged" in warning
    assert "force norm" in warning
    assert repr(reason) in warning


@pytest.mark.parametrize(
    "log_text",
    [
        "",
        "COFKIT_STAGE stage1\nStopping criterion = force tolerance\nCOFKIT_MINIMUM stage1 nan nan 0 0 0 0 0 0\n",
        "COFKIT_STAGE stage1\nno metrics here\n",
    ],
)
def test_unreadable_convergence_diagnostics_keep_hard_failure(tmp_path, log_text):
    path = tmp_path / "log"
    path.write_text(log_text)
    convergence = _check_minimization_convergence(
        path,
        LammpsOptimizationSettings(two_stage_protocol=False, relax_cell=False),
        None,
    )
    assert not convergence["converged"]
    assert not _convergence_diagnostics_available(convergence)


def test_rdkit_nonconverged_conformer_downgrades_to_unconverged():
    field = Mock()
    field.Minimize.return_value = 1
    field.CalcEnergy.return_value = -10
    with (
        patch("cofkit.chem.rdkit.AllChem.MMFFHasAllMoleculeParams", return_value=False),
        patch("cofkit.chem.rdkit.AllChem.UFFGetMoleculeForceField", return_value=field),
    ):
        result = _optimize_conformers(Mock(), (0,))
    assert result.optimization_status == "unconverged"
    assert result.conformer_id == 0
    assert result.energy == -10
    assert result.forcefield == "UFF"
    assert any("did not converge" in diagnostic for diagnostic in result.diagnostics)


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
