"""Real gRASPA engine checks; gated the same way as tests/engines/test_lammps_engine.py.

Set COFKIT_TEST_GRASPA (plus COFKIT_TEST_EQEQ / COFKIT_EQEQ_PATH for the
charge-staging leg, and COFKIT_TEST_LMP for the full hybrid leg) to run these
against real executables. COFKIT_REQUIRE_GRASPA=1 turns a missing gRASPA
binary into a failure instead of a skip.
"""

import os
import shutil
from pathlib import Path

import pytest


def _resolve_binary(env_var: str) -> str | None:
    binary = os.environ.get(env_var)
    if not binary:
        return None
    resolved = shutil.which(binary)
    if resolved is None and Path(binary).is_file():
        resolved = binary
    return resolved


@pytest.fixture
def graspa():
    resolved = _resolve_binary("COFKIT_TEST_GRASPA")
    if not resolved:
        if os.environ.get("COFKIT_REQUIRE_GRASPA") == "1":
            pytest.fail("gRASPA qualification requested but COFKIT_TEST_GRASPA is missing")
        pytest.skip("Real gRASPA evidence requires COFKIT_TEST_GRASPA")
    return resolved


@pytest.fixture
def eqeq():
    resolved = _resolve_binary("COFKIT_TEST_EQEQ") or _resolve_binary("COFKIT_EQEQ_PATH")
    if not resolved:
        if os.environ.get("COFKIT_REQUIRE_GRASPA") == "1":
            pytest.fail("gRASPA qualification requested but no EQeq executable (COFKIT_TEST_EQEQ)")
        pytest.skip("gRASPA workflows stage framework charges with EQeq; set COFKIT_TEST_EQEQ")
    return resolved


@pytest.fixture
def lmp():
    resolved = _resolve_binary("COFKIT_TEST_LMP")
    if not resolved:
        pytest.skip("The full hybrid MD/MC handoff requires COFKIT_TEST_LMP")
    return resolved


def _framework_cif(tmp_path) -> Path:
    cif = tmp_path / "framework.cif"
    cif.write_text(
        (Path(__file__).resolve().parents[1] / "fixtures" / "engine_molecule.cif").read_text()
    )
    return cif


def test_graspa_snapshot_discovery_and_restart_handoff(graspa, eqeq, tmp_path):
    from cofkit.graspa import GraspaIsothermSettings, run_graspa_isotherm_workflow
    from cofkit.guest_restart import (
        build_lammps_guest_restart_state_from_gcmc_result,
        find_latest_gcmc_movie_snapshot,
        write_graspa_restart_file,
    )

    result = run_graspa_isotherm_workflow(
        _framework_cif(tmp_path),
        output_dir=tmp_path / "isotherm",
        eqeq_path=eqeq,
        graspa_path=graspa,
        graspa_timeout_seconds=600.0,
        isotherm_settings=GraspaIsothermSettings(
            component="CO2_DREIDING",
            pressures=(100000.0,),
            fugacity_coefficient=1.0,
            initialization_cycles=20,
            equilibration_cycles=20,
            production_cycles=40,
            number_of_blocks=1,
            movies_every=39,
            create_number_of_molecules=2,
        ),
    )
    simulation_input = Path(result.point_results[0].simulation_input_path).read_text()
    assert "MoviesEvery 39" in simulation_input

    snapshot_path = find_latest_gcmc_movie_snapshot(result)
    assert snapshot_path is not None, "gRASPA did not produce a discoverable Movies/System_0 snapshot"

    # The supported multisite CO2 model must parse from the real snapshot;
    # framework atoms in the movie are excluded but explicitly accounted.
    state = build_lammps_guest_restart_state_from_gcmc_result(result, components=("CO2_DREIDING",))
    assert state.n_atoms > 0
    assert state.n_atoms % 3 == 0
    assert state.components == ("CO2_DREIDING",)
    assert state.skipped_molecules == ()
    assert state.snapshot_cell is not None
    framework_rows = sum(count for _label, count in state.skipped_unknown_site_counts)
    assert framework_rows > 0  # real gRASPA movies include framework atoms
    assert any("counts by label" in warning for warning in state.warnings)

    # The MD->MC boundary artifact: the restart file the hybrid workflow
    # writes after the MD leg must be accepted by the real backend.
    restart_result = write_graspa_restart_file(
        state,
        tmp_path / "md_to_gcmc_restartfile",
        cell=state.snapshot_cell,
        component_order=("CO2_DREIDING",),
    )
    assert restart_result.n_adsorbate_molecules == state.n_atoms // 3
    assert restart_result.n_adsorbate_atoms == state.n_atoms

    restarted = run_graspa_isotherm_workflow(
        _framework_cif(tmp_path),
        output_dir=tmp_path / "isotherm_restarted",
        eqeq_path=eqeq,
        graspa_path=graspa,
        graspa_timeout_seconds=600.0,
        initial_restart_file=restart_result.restart_file_path,
        isotherm_settings=GraspaIsothermSettings(
            component="CO2_DREIDING",
            pressures=(100000.0,),
            fugacity_coefficient=1.0,
            initialization_cycles=0,
            equilibration_cycles=0,
            production_cycles=20,
            number_of_blocks=1,
            movies_every=19,
            restart_file=True,
        ),
    )
    restarted_state = build_lammps_guest_restart_state_from_gcmc_result(
        restarted, components=("CO2_DREIDING",)
    )
    # The restarted run completed, so the real backend accepted the handoff
    # restart file. The parsed population must remain consistent: complete
    # molecules only, and an empty population is reported explicitly rather
    # than silently.
    assert restarted_state.skipped_molecules == ()
    assert restarted_state.n_atoms % 3 == 0
    if restarted_state.n_atoms == 0:
        assert any(
            "adsorbate population is being treated as empty" in warning
            for warning in restarted_state.warnings
        )


def test_hybrid_mdmc_guest_restart_with_real_engines(graspa, eqeq, lmp, tmp_path):
    from cofkit.graspa import EqeqChargeSettings, GraspaMixtureComponentSettings
    from cofkit.hybrid_mdmc import HybridMdMcSettings, run_hybrid_mdmc_workflow
    from cofkit.lammps import LammpsMdSettings

    result = run_hybrid_mdmc_workflow(
        _framework_cif(tmp_path),
        output_dir=tmp_path / "hybrid",
        lmp_path=lmp,
        eqeq_path=eqeq,
        graspa_path=graspa,
        lammps_timeout_seconds=600.0,
        graspa_timeout_seconds=600.0,
        settings=HybridMdMcSettings(
            cycles=2,
            exchange_mode="guest_restart",
            pressure=100000.0,
            components=(
                GraspaMixtureComponentSettings(
                    component="CO2_DREIDING",
                    mol_fraction=1.0,
                    fugacity_coefficient=1.0,
                    create_number_of_molecules=2,
                ),
            ),
            initialization_cycles=10,
            equilibration_cycles=10,
            production_cycles=20,
        ),
        lammps_md_settings=LammpsMdSettings(
            forcefield="uff",
            charge_model="none",
            steps=20,
            enable_omp=False,
        ),
        lammps_eqeq_settings=EqeqChargeSettings(),
        raspa_eqeq_settings=EqeqChargeSettings(),
    )

    assert len(result.cycle_results) == 2
    first, second = result.cycle_results
    assert first.n_output_guest_atoms > 0
    # Guest population is conserved across the MC->MD and MD->MC handoffs.
    assert second.n_input_guest_atoms == first.n_output_guest_atoms
    assert second.n_md_output_guest_atoms == second.n_input_guest_atoms
    assert second.gcmc_initial_restart_file_path is not None
    assert Path(second.gcmc_initial_restart_file_path).is_file()
    assert first.gcmc_initial_restart_file_path is None
    cycle2_input = Path(second.gcmc_result.point_results[0].simulation_input_path).read_text()
    assert "MoviesEvery 19" in cycle2_input
