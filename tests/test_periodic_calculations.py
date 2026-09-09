import math

import gemmi
import pytest

from cofkit.graspa import _compute_unit_cells_from_cif
from cofkit.lammps import LammpsInputError, _parse_explicit_bond_cif
from cofkit.periodic_geometry import images_within
from cofkit.validation import CoarseStructureValidator
from test_validation import _write_test_cif, _summary_record


def test_oblique_cell_supercell_uses_face_widths(tmp_path):
    cif = tmp_path / "cell.cif"
    _write_test_cif(cif, atoms=[("a1_C", "C", 0.1, 0.1, 0.1)], bonds=[])
    cif.write_text(
        cif.read_text()
        .replace("10.0", "28.0")
        .replace("_cell_angle_gamma 90.0", "_cell_angle_gamma 120.0")
    )
    assert _compute_unit_cells_from_cif(cif, cutoff=12.8) == (2, 2, 1)


@pytest.mark.parametrize("origin", [0, 0.23, 1.0])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("coincident", [False, True])
def test_periodic_and_coincident_clashes_are_invariant(
    tmp_path, origin, reverse, coincident
):
    atoms = [
        ("a1_C", "C", 0.01 + origin, 0.1, 0.1),
        ("a1_N", "N", (0.01 if coincident else 0.99) + origin, 0.1, 0.1),
    ]
    cif = tmp_path / "clash.cif"
    _write_test_cif(cif, atoms=atoms[::-1] if reverse else atoms, bonds=[])
    result = CoarseStructureValidator().validate_manifest_record(
        _summary_record(cif, structure_id="x")
    )
    assert "heavy_atom_clash" in result.reasons
    assert result.metrics["min_nonbonded_heavy_distance"] == pytest.approx(
        0 if coincident else 0.2
    )


def test_self_images_and_specific_bond_exclusion():
    small = gemmi.SmallStructure()
    small.cell = gemmi.UnitCell(10, 10, 0.5, 90, 90, 90)
    site = gemmi.SmallStructure.Site()
    site.label, site.element, site.fract = (
        "C",
        gemmi.Element("C"),
        gemmi.Fractional(0.1, 0.1, 0.1),
    )
    small.add_site(site)
    validator = CoarseStructureValidator()
    assert validator._min_nonbonded_heavy_distance_below_cutoff(
        small, set()
    ) == pytest.approx(0.5)
    # Excluding the nearest bonded image must not hide the next periodic contact.
    exclusions = {("C", "C", (0, 0, -1)), ("C", "C", (0, 0, 1))}
    assert validator._min_nonbonded_heavy_distance_below_cutoff(
        small, exclusions
    ) == pytest.approx(1)


def test_image_enumeration_matches_replicated_cell():
    cell = gemmi.UnitCell(10, 10, 0.5, 90, 90, 90)
    reference = gemmi.Fractional(0.1, 0.1, 0.1)
    distances = [
        d
        for shift, d in images_within(cell, reference, reference, 1.05)
        if shift != (0, 0, 0)
    ]
    assert sorted(distances) == pytest.approx([0.5, 0.5, 1, 1])


def test_periodic_multiedges_rejected_before_export(tmp_path):
    cif = tmp_path / "chain.cif"
    cif.write_text("""data_chain
_cell_length_a 2.8
_cell_length_b 20
_cell_length_c 20
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C1 C 0 0 0
C2 C 0.5 0 0
loop_
_geom_bond_atom_site_label_1
_geom_bond_atom_site_label_2
_geom_bond_site_symmetry_1
_geom_bond_site_symmetry_2
_geom_bond_distance
_ccdc_geom_bond_type
C1 C2 . . 1.4 S
C1 C2 . 1_455 1.4 S
""")
    with pytest.raises(LammpsInputError, match="supercell"):
        _parse_explicit_bond_cif(cif)


def test_circular_rotation_near_branch_cut():
    from cofkit.batch import BatchStructureGenerator
    from cofkit.geometry import Frame
    from cofkit.model import MonomerSpec, ReactiveMotif

    spec = MonomerSpec(
        "m",
        "m",
        tuple(
            ReactiveMotif(str(i), "amine", (), Frame((x, y, 0), (x, y, 0), (0, 0, 1)))
            for i, (x, y) in enumerate([(1.0, 0.0), (0.0, 1.0)])
        ),
    )
    directions = tuple(
        (math.cos(a), math.sin(a), 0.0) for a in (math.radians(179), math.radians(-89))
    )
    generator = BatchStructureGenerator()
    rotation, _, _ = generator._rotation_for_planar_motifs(spec, directions)
    assert rotation[0][0] == pytest.approx(-1)
