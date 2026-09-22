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
