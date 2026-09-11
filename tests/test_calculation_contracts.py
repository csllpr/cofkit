"""Regression tests for scientific result acceptance, independent of engines."""

import json
import math
from dataclasses import replace
from pathlib import Path

import pytest

from cofkit.batch import BatchStructureGenerator
from cofkit.batch_models import BatchMonomerRecord
from cofkit.cif_checks import validate_charge_assignment
from cofkit.graspa import (
    EqeqExecutionError,
    GraspaConfigurationError,
    GraspaIsothermSettings,
    GraspaWidomSettings,
    _validate_isotherm_settings,
    _validate_widom_settings,
    _render_widom_simulation_input,
    assign_eqeq_charges_to_cif,
)
from cofkit.lammps import (
    LammpsInputError,
    LammpsParseError,
    LammpsOptimizationSettings,
    _validate_settings,
    optimize_cif_with_lammps,
)
import test_lammps as lammps_fixtures


@pytest.mark.parametrize("value", [math.nan, math.inf, -math.inf])
def test_nonfinite_settings_fail_before_rendering(value):
    with pytest.raises(ValueError, match="finite"):
        _validate_isotherm_settings(GraspaIsothermSettings(temperature=value))
    with pytest.raises(ValueError, match="finite"):
        _validate_isotherm_settings(GraspaIsothermSettings(pressures=(value,)))
    with pytest.raises(ValueError, match="finite"):
        _validate_settings(LammpsOptimizationSettings(pair_cutoff=value))


def test_invalid_sampling_contracts():
    for changes in [
        dict(production_cycles=1.5),
        dict(max_step_per_cycle=0),
        dict(swap_probability=0),
        dict(framework_name="../escape"),
    ]:
        with pytest.raises(ValueError):
            _validate_isotherm_settings(GraspaIsothermSettings(**changes))


def test_unsupported_hydrogen_fails_and_defaults_remain_usable():
    _validate_widom_settings(GraspaWidomSettings())
    with pytest.raises(GraspaConfigurationError, match="Feynman-Hibbs"):
        _validate_widom_settings(GraspaWidomSettings(components=("H2_DREIDING",)))


def test_raspa2_seed_is_effective():
    settings = GraspaWidomSettings(backend="raspa2", random_seed=456)
    first = _render_widom_simulation_input(settings, unit_cells=(1, 1, 1))
    second = _render_widom_simulation_input(
        replace(settings, random_seed=789), unit_cells=(1, 1, 1)
    )
    assert first != second
    assert "RandomSeed                    456" in first


def test_cache_uses_chemistry_and_preserves_source_identity():
    generator = BatchStructureGenerator()
    first = BatchMonomerRecord("same", "same", "CN", "amine", 1, source_path="a")
    second = replace(first, smiles="CCN", source_path="b")
    a = generator.build_monomer(first)
    b = generator.build_monomer(second)
    assert a.ok and b.ok
    assert len(a.monomer.atom_symbols) != len(b.monomer.atom_symbols)
    assert b.record == second
    assert (
        generator.build_monomer(replace(first, source_path="c")).record.source_path
        == "c"
    )


def test_spatial_cache_distinguishes_geometry_with_the_same_id():
    from cofkit.geometry import Frame
    from cofkit.model import MonomerSpec, ReactiveMotif

    targets = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0))
    motifs = tuple(
        ReactiveMotif(str(i), "amine", (), Frame(v, v, (0.0, 0.0, 1.0)))
        for i, v in enumerate(targets)
    )
    first = MonomerSpec("same", "same", motifs)
    second = replace(
        first,
        motifs=tuple(
            replace(
                m,
                frame=replace(
                    m.frame,
                    origin=tuple(-v for v in m.frame.origin),
                    primary=tuple(-v for v in m.frame.primary),
                ),
            )
            for m in motifs
        ),
    )
    generator = BatchStructureGenerator()
    a, _, _ = generator._rotation_for_spatial_motifs(first, targets)
    b, _, _ = generator._rotation_for_spatial_motifs(second, targets)
    assert a[0][0] == pytest.approx(1)
    assert b[0][0] == pytest.approx(-1)


def test_child_seeds_are_reproducible_distinct_and_engine_compatible():
    from cofkit.calculation_io import derive_seed, MAX_ENGINE_SEED

    seeds = [derive_seed(246813, "md", cycle) for cycle in range(100)]
    assert seeds == [derive_seed(246813, "md", cycle) for cycle in range(100)]
    assert len(set(seeds)) == len(seeds)
    assert all(0 < seed <= MAX_ENGINE_SEED for seed in seeds)
    assert derive_seed(246813, "md", 1) != derive_seed(246813, "mc", 1)


def test_rerun_cannot_reuse_lammps_dump(tmp_path):
    fixture = lammps_fixtures.LammpsTests()
    cif = tmp_path / "input.cif"
    cif.write_text(fixture._example_cif_text())
    binary = fixture._write_fake_lammps_binary(tmp_path / "lmp")
    root = tmp_path / "runs"
    result = optimize_cif_with_lammps(
        cif,
        output_dir=root,
        lmp_path=binary,
        settings=LammpsOptimizationSettings(forcefield="uff", charge_model="none"),
    )
    old = Path(result.optimized_cif).read_bytes()
    binary.write_text("#!/bin/sh\nexit 0\n")
    with pytest.raises(LammpsParseError):
        optimize_cif_with_lammps(
            cif,
            output_dir=root,
            lmp_path=binary,
            settings=LammpsOptimizationSettings(forcefield="uff", charge_model="none"),
        )
    assert Path(result.optimized_cif).read_bytes() == old
    manifests = [json.loads(p.read_text()) for p in root.glob("**/attempt.json")]
    assert {m["status"] for m in manifests} == {"completed", "incomplete"}


def test_eqeq_rejects_stale_file(tmp_path):
    cif = tmp_path / "input.cif"
    cif.write_text(lammps_fixtures.LammpsTests()._example_cif_text())
    root = tmp_path / "runs"
    root.mkdir()
    stale = root / "input.cif_EQeq_ewald_1.20_-2.00.cif"
    stale.write_bytes(cif.read_bytes())
    binary = tmp_path / "eqeq"
    binary.write_text("#!/bin/sh\nexit 0\n")
    binary.chmod(0o755)
    with pytest.raises(EqeqExecutionError, match="without writing"):
        assign_eqeq_charges_to_cif(cif, output_dir=root, eqeq_path=binary)
    assert stale.is_file()


def test_fractional_occupancy_rejected_before_charge_assignment(tmp_path):
    cif = tmp_path / "input.cif"
    cif.write_text(
        lammps_fixtures.LammpsTests()._example_cif_text().replace("1.00", "0.50")
    )
    binary = tmp_path / "unused"
    binary.write_text("#!/bin/sh\nexit 99\n")
    binary.chmod(0o755)
    with pytest.raises(LammpsInputError, match="occupancy"):
        optimize_cif_with_lammps(cif, lmp_path=binary)


@pytest.mark.parametrize("charges", [(0.1, -0.1), (math.nan, 0), (0.1, 0.1)])
def test_charge_bijection_and_total(tmp_path, charges):
    header = "data_x\n_cell_length_a 20\n_cell_length_b 20\n_cell_length_c 20\n_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\nloop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
    source = tmp_path / "source.cif"
    source.write_text(header + "A C 0 0 0\nB N 0.1 0 0\n")
    output = tmp_path / "charged.cif"
    output.write_text(
        header
        + f"_atom_site_charge\nrenamed N 0.1 0 0 {charges[0]}\nother C 0 0 0 {charges[1]}\n"
    )
    if charges == (0.1, -0.1):
        checks = validate_charge_assignment(
            source, output, target_charge=0, tolerance=1e-3
        )
        assert checks["net_charge"] == 0
        assert "B N" in output.read_text()
    else:
        with pytest.raises(ValueError):
            validate_charge_assignment(source, output, target_charge=0, tolerance=1e-3)


def test_cell_serialization_rounding_accepted(tmp_path):
    header = (
        "data_x\n"
        "_cell_length_a 14.502860\n_cell_length_b 14.502969\n"
        "_cell_length_c 6.800283\n_cell_angle_alpha 90.000000\n"
        "_cell_angle_beta 90.000000\n_cell_angle_gamma 119.999752\n"
        "loop_\n_atom_site_label\n_atom_site_type_symbol\n"
        "_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
    )
    source = tmp_path / "source.cif"
    source.write_text(header + "A C 0 0 0\nB N 0.1 0 0\n")
    output = tmp_path / "charged.cif"
    output.write_text(
        "data_x\n"
        "_cell_length_a 14.50286\n_cell_length_b 14.50297\n"
        "_cell_length_c 6.80028\n_cell_angle_alpha 90.00000\n"
        "_cell_angle_beta 90.00000\n_cell_angle_gamma 119.99975\n"
        "loop_\n_atom_site_label\n_atom_site_type_symbol\n"
        "_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
        "_atom_site_charge\nA C 0 0 0 0.1\nB N 0.1 0 0 -0.1\n"
    )
    checks = validate_charge_assignment(source, output, target_charge=0, tolerance=1e-3)
    assert checks["net_charge"] == 0
    output.write_text(output.read_text().replace("14.50286", "14.50386"))
    with pytest.raises(ValueError, match="unit cell"):
        validate_charge_assignment(source, output, target_charge=0, tolerance=1e-3)


def test_fractional_rounding_accepted_in_skewed_cell(tmp_path):
    header = (
        "data_x\n"
        "_cell_length_a 14.502860\n_cell_length_b 14.502969\n"
        "_cell_length_c 6.800283\n_cell_angle_alpha 90.000000\n"
        "_cell_angle_beta 90.000000\n_cell_angle_gamma 119.999752\n"
        "loop_\n_atom_site_label\n_atom_site_type_symbol\n"
        "_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
    )
    source = tmp_path / "source.cif"
    source.write_text(header + "A C 0.279585 0.559184 0.750004\n")
    output = tmp_path / "charged.cif"
    # Five-decimal EQeq serialization of the same position; in this skewed
    # cell the Cartesian shift exceeds 1e-4 angstrom.
    output.write_text(
        header + "_atom_site_charge\nA C 0.27959 0.55918 0.75000 0.0\n"
    )
    checks = validate_charge_assignment(source, output, target_charge=0, tolerance=1e-3)
    assert checks["net_charge"] == 0
    output.write_text(
        header + "_atom_site_charge\nA C 0.28059 0.55918 0.75000 0.0\n"
    )
    with pytest.raises(ValueError, match="atom mapping"):
        validate_charge_assignment(source, output, target_charge=0, tolerance=1e-3)
