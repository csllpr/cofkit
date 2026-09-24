import math
from pathlib import Path

import gemmi

from cofkit.batch import BatchStructureGenerator
from cofkit.batch import BatchGenerationConfig
from cofkit.batch_models import BatchMonomerRecord
from cofkit.soft_relax import SoftRelaxConfig, relax_cif_clashes
from cofkit.validation import CoarseStructureValidator

TAPB = "C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N"
TEREPHTHALALDEHYDE = "O=Cc1ccc(C=O)cc1"


def _write_clashing_cif(path: Path) -> None:
    # Fragment F1 has a periodic bond across the x boundary; fragment F2
    # clashes with F1_C1 at 0.8 A (non-bonded) while its own bond is clean.
    path.write_text(
        """data_soft_relax_test
_cell_length_a 15.0
_cell_length_b 15.0
_cell_length_c 15.0
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 90.0

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
F1_C1 C 0.050000 0.500000 0.500000 1.00
F1_C2 C 0.950000 0.500000 0.500000 1.00
F2_C1 C 0.103333 0.500000 0.500000 1.00
F2_C2 C 0.103333 0.500000 0.600000 1.00

loop_
_geom_bond_atom_site_label_1
_geom_bond_atom_site_label_2
_geom_bond_distance
_geom_bond_site_symmetry_1
_geom_bond_site_symmetry_2
_ccdc_geom_bond_type
F1_C1 F1_C2 1.500 . 1_455 S
F2_C1 F2_C2 1.500 . . S
"""
    )


def _bond_distance(block, label_a, label_b):
    small = gemmi.make_small_structure_from_block(block)
    labels1 = block.find_loop("_geom_bond_atom_site_label_1")
    labels2 = block.find_loop("_geom_bond_atom_site_label_2")
    sites = {str(s.label): s for s in small.sites}
    for k in range(len(labels1)):
        if {labels1[k], labels2[k]} == {label_a, label_b}:
            pos_a = small.cell.orthogonalize(sites[label_a].fract)
            best = None
            # minimum over the ±x images of the second site
            for nx in (-1, 0, 1):
                shifted = gemmi.Fractional(
                    sites[label_b].fract.x + nx,
                    sites[label_b].fract.y,
                    sites[label_b].fract.z,
                )
                d = pos_a.dist(small.cell.orthogonalize(shifted))
                if best is None or d < best:
                    best = d
            return best
    raise AssertionError(f"bond {label_a}-{label_b} not found")


def test_soft_relax_relieves_clash_and_preserves_bonds(tmp_path):
    cif_path = tmp_path / "clashing.cif"
    _write_clashing_cif(cif_path)

    validator = CoarseStructureValidator()
    _metrics, reasons_before = validator._validate_cif(
        cif_path, topology_id=None, template_id=None
    )
    assert "heavy_atom_clash" in reasons_before

    report = relax_cif_clashes(cif_path, tmp_path / "relaxed.cif")

    assert report.clashes_before == 1
    assert report.clashes_after == 0
    assert report.max_bond_drift < 0.1

    _metrics, reasons_after = validator._validate_cif(
        Path(report.output_path), topology_id=None, template_id=None
    )
    assert "heavy_atom_clash" not in reasons_after

    block = gemmi.cif.read_file(report.output_path).sole_block()
    assert len(block.find_loop("_geom_bond_atom_site_label_1")) == 2
    for label_a, label_b in (("F1_C1", "F1_C2"), ("F2_C1", "F2_C2")):
        assert abs(_bond_distance(block, label_a, label_b) - 1.5) < 0.1


def test_soft_relax_leaves_clean_structure_essentially_untouched(tmp_path):
    cif_path = tmp_path / "clean.cif"
    _write_clashing_cif(cif_path)
    # Move F2 far away so nothing clashes.
    text = cif_path.read_text().replace("0.103333", "0.400000")
    cif_path.write_text(text)

    report = relax_cif_clashes(cif_path, tmp_path / "relaxed_clean.cif")

    assert report.clashes_before == 0
    assert report.clashes_after == 0
    assert report.max_bond_drift < 0.05
    block = gemmi.cif.read_file(report.output_path).sole_block()
    for label_a, label_b in (("F1_C1", "F1_C2"), ("F2_C1", "F2_C2")):
        assert abs(_bond_distance(block, label_a, label_b) - 1.5) < 0.05


def test_build_pipeline_soft_relax_flag_records_diagnostics(tmp_path):
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

    off_generator = BatchStructureGenerator(
        BatchGenerationConfig(
            rdkit_num_conformers=1, retain_top_results=1, single_node_topology_ids=("hcb",)
        )
    )
    on_generator = BatchStructureGenerator(
        BatchGenerationConfig(
            rdkit_num_conformers=1,
            retain_top_results=1,
            single_node_topology_ids=("hcb",),
            soft_relax=True,
        )
    )

    summary_off, _candidate = off_generator.generate_pair_candidate(
        amine, aldehyde, out_dir=tmp_path / "off", write_cif=True
    )
    summary_on, _candidate = on_generator.generate_pair_candidate(
        amine, aldehyde, out_dir=tmp_path / "on", write_cif=True
    )

    assert summary_off.status == "ok" and summary_on.status == "ok"
    validation_off = summary_off.metadata.get("validation", {})
    validation_on = summary_on.metadata.get("validation", {})
    assert "soft_relax" not in validation_off
    soft = validation_on["soft_relax"]
    assert soft["applied"] is True
    assert soft["clashes_after"] == 0
    assert soft["max_bond_drift"] < 0.1

    # The relaxed export still parses cleanly and keeps its bond loop.
    block = gemmi.cif.read_file(summary_on.cif_path).sole_block()
    assert len(block.find_loop("_geom_bond_atom_site_label_1")) > 0


def _write_cif(path: Path, atoms, bonds, cell=(15.0, 15.0, 15.0)) -> None:
    # atoms: (label, symbol, (x, y, z)) cartesian; bonds: (label1, label2, sym2)
    a, b, c = cell
    lines = [
        "data_soft_relax_fixture",
        f"_cell_length_a {a}",
        f"_cell_length_b {b}",
        f"_cell_length_c {c}",
        "_cell_angle_alpha 90.0",
        "_cell_angle_beta 90.0",
        "_cell_angle_gamma 90.0",
        "",
        "loop_",
        "_atom_site_label",
        "_atom_site_type_symbol",
        "_atom_site_fract_x",
        "_atom_site_fract_y",
        "_atom_site_fract_z",
        "_atom_site_occupancy",
    ]
    for label, symbol, (x, y, z) in atoms:
        lines.append(f"{label} {symbol} {x / a:.6f} {y / b:.6f} {z / c:.6f} 1.00")
    lines += [
        "",
        "loop_",
        "_geom_bond_atom_site_label_1",
        "_geom_bond_atom_site_label_2",
        "_geom_bond_distance",
        "_geom_bond_site_symmetry_1",
        "_geom_bond_site_symmetry_2",
        "_ccdc_geom_bond_type",
    ]
    for label1, label2, sym2 in bonds:
        lines.append(f"{label1} {label2} 1.500 . {sym2} S")
    path.write_text("\n".join(lines) + "\n")


def _distance(path: Path, label_a: str, label_b: str) -> float:
    """Minimum-image distance between two labeled sites (bonded or not)."""
    block = gemmi.cif.read_file(str(path)).sole_block()
    small = gemmi.make_small_structure_from_block(block)
    sites = {str(s.label): s for s in small.sites}
    pa = small.cell.orthogonalize(sites[label_a].fract)
    best = None
    for nx in (-1, 0, 1):
        for ny in (-1, 0, 1):
            for nz in (-1, 0, 1):
                shifted = gemmi.Fractional(
                    sites[label_b].fract.x + nx,
                    sites[label_b].fract.y + ny,
                    sites[label_b].fract.z + nz,
                )
                d = pa.dist(small.cell.orthogonalize(shifted))
                if best is None or d < best:
                    best = d
    return best


def test_soft_relax_preserves_hydrogen_bond(tmp_path):
    # O-H...N contact at 1.8 A in an otherwise clean structure.
    cif_path = tmp_path / "hbond.cif"
    _write_cif(
        cif_path,
        atoms=[
            ("Ca", "C", (5.0, 7.5, 7.5)),
            ("O1", "O", (6.43, 7.5, 7.5)),
            ("H1", "H", (6.43, 7.5, 8.5)),
            ("N1", "N", (6.43, 7.5, 10.3)),
            ("Cb", "C", (6.43, 7.5, 11.77)),
        ],
        bonds=[("Ca", "O1", "."), ("O1", "H1", "."), ("N1", "Cb", ".")],
    )
    assert abs(_distance(cif_path, "H1", "N1") - 1.8) < 1e-3

    report = relax_cif_clashes(cif_path, tmp_path / "hbond_relaxed.cif")

    assert 1.5 < _distance(Path(report.output_path), "H1", "N1") < 2.2


def test_soft_relax_repairs_pathological_one_three_contact(tmp_path):
    # X-A-Y with the angle at A collapsed so the 1-3 pair X...Y sits at 0.6 A:
    # a Urey-Bradley restraint at the initial distance would hold the clash.
    cif_path = tmp_path / "one_three.cif"
    _write_cif(
        cif_path,
        atoms=[
            ("A1", "C", (7.5, 7.5, 7.5)),
            ("X1", "C", (9.0, 7.5, 7.5)),
            ("Y1", "C", (8.9, 8.083, 7.5)),
        ],
        bonds=[("A1", "X1", "."), ("A1", "Y1", ".")],
    )
    assert _distance(cif_path, "X1", "Y1") < 0.7

    report = relax_cif_clashes(cif_path, tmp_path / "one_three_relaxed.cif")

    assert _distance(Path(report.output_path), "X1", "Y1") > 1.05
    for pair in (("A1", "X1"), ("A1", "Y1")):
        assert abs(_distance(Path(report.output_path), *pair) - 1.5) < 0.15


def test_soft_relax_leaves_stacked_layers_in_place(tmp_path):
    # Two eclipsed benzene-like layers at 3.4 A; layer L1 carries an
    # out-of-plane N-H pointing into the gap (H...C above at 2.5 A), like a
    # pyramidal amine in a physically stacked cell. Nothing here is a clash,
    # so the layers must not slide.
    import math as _math

    atoms = []
    bonds = []
    radius = 1.4
    for layer, z in (("L1", 6.0), ("L2", 9.4)):
        for k in range(6):
            angle = _math.radians(60 * k)
            x, y = 7.5 + radius * _math.cos(angle), 7.5 + radius * _math.sin(angle)
            symbol = "N" if (layer == "L1" and k == 0) else "C"
            atoms.append((f"{layer}_C{k}", symbol, (x, y, z)))
            bonds.append((f"{layer}_C{k}", f"{layer}_C{(k + 1) % 6}", "."))
            hx, hy = 7.5 + (radius + 1.08) * _math.cos(angle), 7.5 + (radius + 1.08) * _math.sin(angle)
            atoms.append((f"{layer}_H{k}", "H", (hx, hy, z)))
            bonds.append((f"{layer}_C{k}", f"{layer}_H{k}", "."))
    atoms.append(("L1_HN", "H", (7.5 + radius, 7.5, 6.9)))
    bonds.append(("L1_C0", "L1_HN", "."))
    cif_path = tmp_path / "stacked.cif"
    _write_cif(cif_path, atoms, bonds)

    report = relax_cif_clashes(cif_path, tmp_path / "stacked_relaxed.cif")

    block = gemmi.cif.read_file(report.output_path).sole_block()
    small = gemmi.make_small_structure_from_block(block)
    pos = {str(s.label): small.cell.orthogonalize(s.fract) for s in small.sites}

    def centroid(prefix, z_target):
        pts = [p for label, p in pos.items() if label.startswith(f"{prefix}_C")]
        return (
            sum(p.x for p in pts) / len(pts),
            sum(p.y for p in pts) / len(pts),
        )

    before = (centroid("L1", 6.0), centroid("L2", 9.4))
    dx = (before[1][0] - before[0][0]) - 0.0
    dy = (before[1][1] - before[0][1]) - 0.0
    assert math.hypot(dx, dy) < 0.05  # eclipsed alignment preserved
    assert report.max_bond_drift < 0.1
