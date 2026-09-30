import contextlib
import io
import json
import math
import re
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import gemmi

from cofkit import (
    AromaticityRestoreError,
    BatchGenerationConfig,
    BatchMonomerRecord,
    BatchStructureGenerator,
    CoarseValidationThresholds,
    Frame,
    MonomerRoleResolver,
    MonomerSpec,
    ReactiveMotif,
    ReactionRealizer,
    build_rdkit_monomer,
)
from cofkit.geometry import add, dot, matmul_vec, normalize, scale, sub
from cofkit.chem import rdkit as rdkit_module
from cofkit.validation import CoarseStructureValidator

try:
    from rdkit import Chem  # noqa: F401
except ImportError:  # pragma: no cover - environment-dependent
    Chem = None


TAPB = "C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N"
TFB = "C1=C(C=C(C=C1C=O)C=O)C=O"
TEREPHTHALALDEHYDE = "O=Cc1ccc(C=O)cc1"
PPD = "Nc1ccc(N)cc1"
ASYMMETRIC_TRIAMINE = (
    "COC(=O)[C@@H]1C[C@H](OCc2c(-c3ccc(N)cc3)cc(-c3ccc(N)cc3)cc2-c2ccc(N)cc2)CN1C(=O)OC(C)(C)C"
)
LONG_DIALDEHYDE = "O=Cc1ccc(C#Cc2c3ccccc3c(C#Cc3ccc(C=O)cc3)c3ccccc23)cc1"
ASYMMETRIC_TRIALDEHYDE = "O=Cc1ccccc1C#Cc1cc(C#Cc2ccccc2C=O)cc(C#Cc2ccccc2C=O)c1"
TETRA_AMINE = "Nc1ccc(C(c2ccc(N)cc2)(c2ccc(N)cc2)c2ccc(N)cc2)cc1"
TETRA_ALDEHYDE = "O=Cc1ccc(C(c2ccc(C=O)cc2)(c2ccc(C=O)cc2)c2ccc(C=O)cc2)cc1"
D2H_TETRA_AMINE = "Nc1cc(N)c(N)cc1N"
BEX_D2H_ALDEHYDE = (
    "C1C=C(N(C2C=CC(C3C4C(=NON=4)C(C4C=CC(N(C5C=CC(C=O)=CC=5)C5C=CC(C=O)=CC=5)=CC=4)=CC=3)=CC=2)"
    "C2C=CC(C=O)=CC=2)C=CC=1C=O"
)
BEX_D2H_AMINE = (
    "Nc1ccc(-c2ccc3c(c2)C(=C2c4cc(-c5ccc(N)cc5)ccc4-c4ccc(-c5ccc(N)cc5)cc42)c2cc(-c4ccc(N)cc4)ccc2-3)cc1"
)
HEXA_AMINE = "Nc1cc2c3cc(N)c(N)cc3c3cc(N)c(N)cc3c2cc1N"
HEXA_ALDEHYDE = "O=Cc1ccc(-c2cc3c(cc2-c2ccc(C=O)cc2)C2c4cc(-c5ccc(C=O)cc5)c(-c5ccc(C=O)cc5)cc4C3c3cc(-c4ccc(C=O)cc4)c(-c4ccc(C=O)cc4)cc32)cc1"
TP = "O=Cc1c(O)c(C=O)c(O)c(C=O)c1O"
PHENYLENEDIACETONITRILE = "N#CCc1ccc(CC#N)cc1"
COF42_HYDRAZIDE = "CCOc1cc(C(=O)NN)cc(C(=O)NN)c1OCC"


def _distort_first_bridge_bond(cif_path: Path, shift_angstrom: float = 2.1) -> str:
    """Deliberately distort one realized linkage bond in a written CIF.

    Moves the second atom of the first inter-instance bond perpendicular to
    the bond direction (sign chosen to maximize its clearance from other
    atoms) and rewrites its coordinate row. Returns the pristine CIF text so
    a fake repair step can restore it. Used to prove that validation verdicts
    follow the final coordinates, not the seed/assembly metrics.
    """
    from cofkit.cif_checks import cif_value_str
    from cofkit.periodic_geometry import p1_shift

    pristine = cif_path.read_text(encoding="utf-8")
    block = gemmi.cif.read_file(str(cif_path)).sole_block()
    small = gemmi.make_small_structure_from_block(block)
    positions = {str(site.label): (site.fract.x, site.fract.y, site.fract.z) for site in small.sites}
    labels1 = block.find_loop("_geom_bond_atom_site_label_1")
    labels2 = block.find_loop("_geom_bond_atom_site_label_2")
    syms1 = block.find_loop("_geom_bond_site_symmetry_1")
    syms2 = block.find_loop("_geom_bond_site_symmetry_2")
    moved_label = None
    moved_fract = None
    for index in range(len(labels1)):
        label_a = cif_value_str(labels1[index])
        label_b = cif_value_str(labels2[index])
        if label_a.rsplit("_", 1)[0] == label_b.rsplit("_", 1)[0]:
            continue
        shift = tuple(
            b - a
            for a, b in zip(
                p1_shift(cif_value_str(syms1[index]) if len(syms1) else "."),
                p1_shift(cif_value_str(syms2[index]) if len(syms2) else "."),
            )
        )
        cart_a = small.cell.orthogonalize(gemmi.Fractional(*positions[label_a]))
        cart_b = small.cell.orthogonalize(
            gemmi.Fractional(
                positions[label_b][0] + shift[0],
                positions[label_b][1] + shift[1],
                positions[label_b][2] + shift[2],
            )
        )
        bond_vec = (cart_b.x - cart_a.x, cart_b.y - cart_a.y, cart_b.z - cart_a.z)
        norm = math.sqrt(sum(component * component for component in bond_vec))
        bond_vec = tuple(component / norm for component in bond_vec)
        # A perpendicular shift leaves the moved atom's other bonds stretched
        # rather than compressed, avoiding fabricated severe overlaps.
        perp = (bond_vec[1], -bond_vec[0], 0.0)
        perp_norm = math.sqrt(sum(component * component for component in perp))
        if perp_norm < 1e-6:
            perp = (0.0, 0.0, 1.0)
            perp_norm = 1.0
        perp = tuple(component / perp_norm for component in perp)
        others = [
            small.cell.orthogonalize(gemmi.Fractional(*fract))
            for label, fract in positions.items()
            if label != label_b
        ]

        def clearance(sign: float):
            moved = gemmi.Position(
                cart_b.x + sign * shift_angstrom * perp[0],
                cart_b.y + sign * shift_angstrom * perp[1],
                cart_b.z + sign * shift_angstrom * perp[2],
            )
            return min(moved.dist(other) for other in others), moved

        forward, backward = clearance(1.0), clearance(-1.0)
        _best_clearance, moved_position = max(forward, backward, key=lambda item: item[0])
        moved_fract = small.cell.fractionalize(moved_position)
        moved_label = label_b
        break
    assert moved_label is not None and moved_fract is not None
    rewritten_lines = []
    for line in pristine.splitlines():
        tokens = line.split()
        if tokens and cif_value_str(tokens[0]) == moved_label and len(tokens) >= 6:
            tokens[2] = f"{moved_fract.x:.6f}"
            tokens[3] = f"{moved_fract.y:.6f}"
            tokens[4] = f"{moved_fract.z:.6f}"
            line = " ".join(tokens)
        rewritten_lines.append(line)
    cif_path.write_text("\n".join(rewritten_lines) + "\n", encoding="utf-8")
    return pristine


class _DistortingCifWriter:
    """Test wrapper: delegate to a real CIFWriter, then distort one bridge
    bond in the written file so the exported geometry no longer matches the
    clean seed/assembly metrics. Records the pristine text by file name."""

    def __init__(self, delegate, pristine_by_name: dict[str, str]):
        self._delegate = delegate
        self._pristine_by_name = pristine_by_name

    def write_candidate(self, path, *args, **kwargs):
        result = self._delegate.write_candidate(path, *args, **kwargs)
        path = Path(path)
        self._pristine_by_name[path.name] = _distort_first_bridge_bond(path)
        return result

    def __getattr__(self, name):
        return getattr(self._delegate, name)


@unittest.skipIf(Chem is None, "RDKit is not available")
class BatchStructureGeneratorTests(unittest.TestCase):
    def _world_position(self, candidate, instance_id, local_position, *, image=(0, 0, 0)):
        pose = candidate.state.monomer_poses[instance_id]
        cell = candidate.state.cell
        return add(
            add(pose.translation, matmul_vec(pose.rotation_matrix, local_position)),
            add(
                add(scale(cell[0], image[0]), scale(cell[1], image[1])),
                scale(cell[2], image[2]),
            ),
        )

    def _angle(self, point_a, point_b, point_c) -> float:
        left = normalize(sub(point_a, point_b))
        right = normalize(sub(point_c, point_b))
        cosine = max(-1.0, min(1.0, dot(left, right)))
        import math

        return math.degrees(math.acos(cosine))

    def setUp(self):
        self.generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
                single_node_topology_ids=("hcb",),
            )
        )

    def test_load_smiles_library_parses_fixture_format(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "amines_count_2.txt"
            path.write_text("smiles\nNc1ccc(N)cc1\nNc1ccnc(N)c1\n", encoding="utf-8")

            records = self.generator.load_smiles_library(
                path,
                motif_kind="amine",
                expected_connectivity=2,
            )

        self.assertEqual([record.id for record in records], ["amines_count_2_0001", "amines_count_2_0002"])
        self.assertEqual(records[0].smiles, "Nc1ccc(N)cc1")
        self.assertEqual(records[0].expected_connectivity, 2)

    def test_infer_monomer_record_detects_amine(self):
        record = self.generator.infer_monomer_record(
            PPD,
            record_id="generic_0001",
        )

        self.assertEqual(record.motif_kind, "amine")
        self.assertEqual(record.expected_connectivity, 2)
        self.assertTrue(record.metadata["auto_detected"])
        self.assertEqual(record.metadata["detected_motif_kinds"], ("amine",))

    def test_infer_monomer_record_prefers_keto_aldehyde_over_generic_aldehyde(self):
        record = self.generator.infer_monomer_record(
            TP,
            record_id="generic_0001",
        )

        self.assertEqual(record.motif_kind, "keto_aldehyde")
        self.assertEqual(record.expected_connectivity, 3)
        self.assertEqual(record.metadata["detected_motif_kinds"], ("aldehyde", "keto_aldehyde"))
        self.assertEqual(record.metadata["detected_connectivities"]["keto_aldehyde"], 3)

    def test_forced_kind_warnings_flag_keto_aldehyde_behind_imine_assignment(self):
        resolver = MonomerRoleResolver.builtin()

        warnings = resolver.forced_kind_warnings(
            TP,
            assigned_kind="aldehyde",
            monomer_id="tp",
            template_id="imine_bridge",
        )

        self.assertEqual(len(warnings), 1)
        self.assertIn("keto_aldehyde", warnings[0])
        self.assertIn("keto_enamine_bridge", warnings[0])

    def test_forced_kind_warnings_silent_for_keto_enamine_assignment(self):
        resolver = MonomerRoleResolver.builtin()

        warnings = resolver.forced_kind_warnings(
            TP,
            assigned_kind="keto_aldehyde",
            monomer_id="tp",
            template_id="keto_enamine_bridge",
        )

        self.assertEqual(warnings, ())

    def test_forced_kind_warnings_silent_for_plain_aldehyde(self):
        resolver = MonomerRoleResolver.builtin()

        warnings = resolver.forced_kind_warnings(
            TFB,
            assigned_kind="aldehyde",
            monomer_id="tfb",
            template_id="imine_bridge",
        )

        self.assertEqual(warnings, ())

    def test_forced_kind_warnings_silent_for_nitrile_activated_methylene_overlap(self):
        resolver = MonomerRoleResolver.builtin()

        warnings = resolver.forced_kind_warnings(
            PHENYLENEDIACETONITRILE,
            assigned_kind="activated_methylene",
            monomer_id="pda",
        )

        self.assertEqual(warnings, ())

    def test_infer_monomer_record_with_restricted_kinds_warns_about_keto_aldehyde(self):
        record = self.generator.infer_monomer_record(
            TP,
            record_id="generic_0001",
            allowed_motif_kinds=("amine", "aldehyde"),
        )

        self.assertEqual(record.motif_kind, "aldehyde")
        warnings = record.metadata["overlap_warnings"]
        self.assertEqual(len(warnings), 1)
        self.assertIn("keto_enamine_bridge", warnings[0])

    def test_infer_monomer_record_plain_aldehyde_has_no_overlap_warnings(self):
        record = self.generator.infer_monomer_record(
            TFB,
            record_id="generic_0001",
            allowed_motif_kinds=("amine", "aldehyde"),
        )

        self.assertEqual(record.motif_kind, "aldehyde")
        self.assertEqual(record.metadata["overlap_warnings"], ())

    def test_load_smiles_library_attaches_overlap_warnings(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "aldehydes_count_3.txt"
            path.write_text(f"smiles\n{TP}\n{TFB}\n", encoding="utf-8")

            records = self.generator.load_smiles_library(
                path,
                motif_kind="aldehyde",
                expected_connectivity=3,
            )

        keto_warnings = records[0].metadata["overlap_warnings"]
        self.assertEqual(len(keto_warnings), 1)
        self.assertIn("keto_enamine_bridge", keto_warnings[0])
        self.assertEqual(records[1].metadata["overlap_warnings"], ())

    def test_load_binary_bridge_test_set_auto_groups_generic_files_by_detected_role_and_connectivity(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                allowed_reactions=("keto_enamine_bridge",),
                rdkit_num_conformers=2,
                retain_top_results=5,
                single_node_topology_ids=("hcb",),
            )
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            (temp_path / "set_a.txt").write_text(f"smiles\n{PPD}\n", encoding="utf-8")
            (temp_path / "set_b.txt").write_text(f"smiles\n{TP}\n", encoding="utf-8")

            libraries = generator.load_binary_bridge_test_set_auto(
                temp_path,
                template_id="keto_enamine_bridge",
            )

        self.assertEqual(tuple(sorted(libraries)), ("amines_count_2", "keto_aldehydes_count_3"))
        self.assertEqual(libraries["amines_count_2"][0].motif_kind, "amine")
        self.assertEqual(libraries["keto_aldehydes_count_3"][0].motif_kind, "keto_aldehyde")

    def test_available_binary_bridge_template_ids_discovers_supported_templates(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            (temp_path / "amines_count_3.txt").write_text(f"smiles\n{TAPB}\n", encoding="utf-8")
            (temp_path / "aldehydes_count_3.txt").write_text(f"smiles\n{TFB}\n", encoding="utf-8")
            (temp_path / "hydrazines_count_2.txt").write_text("smiles\nNN\n", encoding="utf-8")
            (temp_path / "hydrazides_count_2.txt").write_text(f"smiles\n{COF42_HYDRAZIDE}\n", encoding="utf-8")
            (temp_path / "keto_aldehydes_count_3.txt").write_text(f"smiles\n{TP}\n", encoding="utf-8")

            template_ids = self.generator.available_binary_bridge_template_ids(temp_path)

        self.assertEqual(template_ids, ("azine_bridge", "hydrazone_bridge", "imine_bridge", "keto_enamine_bridge"))

    def test_custom_smiles_monomer_builder_can_be_injected(self):
        calls = []

        def fake_builder(monomer_id, name, smiles, motif_kind, *, num_conformers, random_seed):
            del name, smiles
            calls.append((monomer_id, motif_kind, num_conformers, random_seed))
            return MonomerSpec(
                id=monomer_id,
                name=monomer_id,
                motifs=(
                    ReactiveMotif(id=f"{motif_kind}1", kind=motif_kind, atom_ids=(1,), frame=Frame.xy()),
                    ReactiveMotif(id=f"{motif_kind}2", kind=motif_kind, atom_ids=(2,), frame=Frame.yz()),
                ),
            )

        generator = BatchStructureGenerator(
            BatchGenerationConfig(rdkit_num_conformers=3, rdkit_random_seed=123),
            smiles_monomer_builder=fake_builder,
        )
        record = BatchMonomerRecord(
            id="fake_amine",
            name="fake_amine",
            smiles="N",
            motif_kind="amine",
            expected_connectivity=2,
        )

        built = generator.build_monomer(record)

        self.assertTrue(built.ok)
        self.assertEqual(calls, [("fake_amine", "amine", 3, 123)])

    def test_three_plus_three_pair_builds_hcb_candidate(self):
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tfb",
            name="tfb",
            smiles=TFB,
            motif_kind="aldehyde",
            expected_connectivity=3,
        )

        summary, candidate = self.generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "3+3-node-node")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "hcb")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "single-node-bipartite")
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 3)

    def test_top_rung_monomers_produce_no_geometry_degradation_warning(self):
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tfb",
            name="tfb",
            smiles=TFB,
            motif_kind="aldehyde",
            expected_connectivity=3,
        )
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            summary, candidate = self.generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertIsNotNone(candidate)
        self.assertEqual(summary.metadata["monomer_geometry_warnings"], ())
        self.assertNotIn("monomer_geometry_degraded", stderr.getvalue())

    def test_degraded_monomer_geometry_surfaces_in_summary_metadata_and_validation(self):
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tfb",
            name="tfb",
            smiles=TFB,
            motif_kind="aldehyde",
            expected_connectivity=3,
        )
        stderr = io.StringIO()
        with tempfile.TemporaryDirectory() as temp_dir:
            with patch.object(rdkit_module.AllChem, "MMFFHasAllMoleculeParams", return_value=False):
                with contextlib.redirect_stderr(stderr):
                    summary, candidate = self.generator.generate_pair_candidate(
                        amine,
                        aldehyde,
                        out_dir=temp_dir,
                        write_cif=True,
                    )

        self.assertEqual(summary.status, "ok")
        self.assertIsNotNone(candidate)
        warnings = summary.metadata["monomer_geometry_warnings"]
        self.assertTrue(
            any("monomer_geometry_degraded:tapb:" in warning and "UFF" in warning for warning in warnings)
        )
        self.assertTrue(
            any("monomer_geometry_degraded:tfb:" in warning and "UFF" in warning for warning in warnings)
        )
        self.assertIn("warning: monomer_geometry_degraded:tapb:", stderr.getvalue())
        self.assertIn("warning: monomer_geometry_degraded:tfb:", stderr.getvalue())
        validation = summary.metadata["validation"]
        self.assertIn("monomer_geometry_degraded", validation["warning_reasons"])
        details = validation["metrics"]["monomer_geometry_degraded_details"]
        self.assertTrue(any("monomer_geometry_degraded:tapb:" in detail for detail in details))

    def test_three_plus_three_pair_can_enumerate_and_export_stacked_variant(self):
        base_generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
                single_node_topology_ids=("hcb",),
            )
        )
        stacked_generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
                single_node_topology_ids=("hcb",),
                stacking_ids=("AA",),
            )
        )
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tfb",
            name="tfb",
            smiles=TFB,
            motif_kind="aldehyde",
            expected_connectivity=3,
        )

        base_summary, base_candidate = base_generator.generate_pair_candidate(amine, aldehyde)
        self.assertEqual(base_summary.status, "ok")
        self.assertIsNotNone(base_candidate)
        assert base_candidate is not None

        with tempfile.TemporaryDirectory() as temp_dir:
            summaries, candidates, attempted_structures = stacked_generator.generate_pair_candidates(
                amine,
                aldehyde,
                out_dir=temp_dir,
                write_cif=True,
            )

            self.assertEqual(attempted_structures, 1)
            self.assertEqual(len(summaries), 1)
            self.assertEqual(len(candidates), 1)
            summary = summaries[0]
            candidate = candidates[0]
            self.assertEqual(summary.status, "ok")
            self.assertEqual(summary.structure_id, "tapb__tfb__hcb__AA")
            self.assertEqual(candidate.id, "tapb__tfb__hcb__single_node_node_node__AA")
            self.assertEqual(summary.metadata["stacking"]["id"], "AA")
            self.assertEqual(summary.metadata["stacking"]["comment_suffix"], "stacking=AA")
            self.assertEqual(candidate.state.stacking_state, "AA")
            self.assertTrue(candidate.metadata["embedding"]["stacking_enabled"])
            self.assertIn("stacked_2d", candidate.flags)
            self.assertEqual(
                candidate.metadata["graph_summary"]["n_monomer_instances"],
                2 * base_candidate.metadata["graph_summary"]["n_monomer_instances"],
            )
            self.assertEqual(
                candidate.metadata["graph_summary"]["n_reaction_events"],
                2 * base_candidate.metadata["graph_summary"]["n_reaction_events"],
            )
            self.assertIsNotNone(summary.cif_path)
            cif_text = Path(summary.cif_path).read_text(encoding="utf-8")

        self.assertEqual(cif_text.splitlines()[0], f"# COFid: {summary.metadata['cofid']} stacking=AA")

    def _min_interlayer_contact_from_cif(self, cif_path: Path) -> float:
        """Atomistic minimum-image contact between the two layers of a stacked
        export, measured independently of the stacking metadata (layers are
        identified by the L0/L1 atom-label suffix)."""
        doc = gemmi.cif.read_file(str(cif_path))
        block = doc.sole_block()
        cell = gemmi.UnitCell(
            float(block.find_value("_cell_length_a")),
            float(block.find_value("_cell_length_b")),
            float(block.find_value("_cell_length_c")),
            float(block.find_value("_cell_angle_alpha")),
            float(block.find_value("_cell_angle_beta")),
            float(block.find_value("_cell_angle_gamma")),
        )
        labels = [str(value) for value in block.find_values("_atom_site_label")]
        fract_x = [float(value) for value in block.find_values("_atom_site_fract_x")]
        fract_y = [float(value) for value in block.find_values("_atom_site_fract_y")]
        fract_z = [float(value) for value in block.find_values("_atom_site_fract_z")]
        layer_zero: list[gemmi.Position] = []
        layer_one: list[gemmi.Position] = []
        for label, fx, fy, fz in zip(labels, fract_x, fract_y, fract_z):
            position = cell.orthogonalize(gemmi.Fractional(fx, fy, fz))
            if "L0" in label:
                layer_zero.append(position)
            elif "L1" in label:
                layer_one.append(position)
        self.assertTrue(layer_zero)
        self.assertTrue(layer_one)
        return min(
            cell.find_nearest_image(first, second).dist()
            for first in layer_zero
            for second in layer_one
        )

    def test_tapb_tfb_stacked_exports_carry_honest_span_and_cell_derivation(self):
        """Plan W7 test 1 / definition-of-done criterion 1: canonical TAPB/TFB
        imine hcb AA+AB+slipped end-to-end regression. The originally shipped
        exports recorded layer_z_span 0.0 with c = 6.8 A and a 2.02 A
        interlayer C...C contact that validated as `valid`."""
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                single_node_topology_ids=("hcb",),
                stacking_ids=("AA", "AB", "slipped"),
            )
        )
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tfb",
            name="tfb",
            smiles=TFB,
            motif_kind="aldehyde",
            expected_connectivity=3,
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            summaries, candidates, attempted_structures = generator.generate_pair_candidates(
                amine,
                aldehyde,
                out_dir=temp_dir,
                write_cif=True,
            )

            self.assertEqual(attempted_structures, 3)
            self.assertEqual(len(summaries), 3)
            self.assertEqual(len(candidates), 3)
            summaries_by_registry = {}
            for summary in summaries:
                self.assertEqual(summary.status, "ok")
                stacking = summary.metadata["stacking"]
                registry = stacking["id"]
                summaries_by_registry[registry] = summary
                # The span must be a real atomistic measurement along the
                # layer normal, well above the shipped 0.0 fallback.
                self.assertGreater(stacking["layer_z_span"], 3.0)
                self.assertEqual(stacking["layer_z_span_mode"], "atomistic_product")
                self.assertEqual(stacking["layer_z_span_axis"], "layer_normal")
                # c2c = interlayer_clearance + layer_z_span, and the stacked
                # cell is built along the layer normal with c = 2*c2c.
                self.assertAlmostEqual(
                    stacking["center_to_center_distance"],
                    stacking["interlayer_clearance"] + stacking["layer_z_span"],
                    places=9,
                )
                self.assertEqual(stacking["layer_count"], 2)
                self.assertEqual(stacking["c_axis_basis"], "layer_normal")
                # Atomistic interlayer contact self-check ran and found no
                # contact below the cutoff (None is a lower bound, not a
                # missing value); the shipped AA export measured 2.02 A.
                self.assertIsNone(stacking["min_interlayer_contact"])
                self.assertEqual(stacking["min_interlayer_contact_mode"], "atomistic_product")
                self.assertNotIn("stacking_clash", summary.flags)
                self.assertIsNotNone(summary.cif_path)
                cif_text = Path(summary.cif_path).read_text(encoding="utf-8")
                cell_c = float(re.search(r"_cell_length_c (\S+)", cif_text).group(1))
                self.assertAlmostEqual(
                    cell_c,
                    stacking["layer_count"] * stacking["center_to_center_distance"],
                    places=4,
                )
                self.assertIn("# stacking-geometry:", cif_text)
                self.assertIn(f"registry={registry}", cif_text)
                self.assertIn("# c-axis-semantics: periodic_bilayer", cif_text)
                self.assertEqual(summary.metadata["embedding"]["c_axis_semantics"], "periodic_bilayer")
                self.assertEqual(
                    summary.metadata["validation"]["classification"],
                    "valid",
                )

            self.assertEqual(set(summaries_by_registry), {"AA", "AB", "slipped"})
            for registry, summary in summaries_by_registry.items():
                contact = self._min_interlayer_contact_from_cif(Path(summary.cif_path))
                self.assertGreaterEqual(
                    contact,
                    2.5,
                    f"{registry} min interlayer contact {contact:.3f} A regressed below 2.5 A",
                )

            # W7 test 3 (60-degree half): the hexagonal AB shift maps a vertex
            # onto the pore center, so its cartesian magnitude is a/sqrt(3).
            ab_candidate = next(
                candidate
                for candidate in candidates
                if candidate.metadata["stacking"]["id"] == "AB"
            )
            ab_stacking = ab_candidate.metadata["stacking"]
            self.assertEqual(ab_stacking["cell_classification"]["setting"], "60deg")
            cell = ab_candidate.state.cell
            shift = ab_stacking["registry_shift_fractional"]
            shift_cartesian = [
                shift[0] * cell[0][axis] + shift[1] * cell[1][axis]
                for axis in range(3)
            ]
            a_length = math.sqrt(sum(component**2 for component in cell[0]))
            shift_length = math.sqrt(sum(component**2 for component in shift_cartesian))
            # The fitted cell is only approximately hexagonal, so compare
            # with a 1% relative tolerance rather than exactly.
            self.assertAlmostEqual(
                shift_length,
                a_length / math.sqrt(3.0),
                delta=0.01 * a_length / math.sqrt(3.0),
            )

    def test_validator_flags_c_squashed_interpenetrated_variant(self):
        """Squash the fresh AA export's c axis to the shipped 6.8 A while
        keeping cartesian coordinates (recreating the original interpenetrated
        geometry) and confirm the validator now flags it."""
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                single_node_topology_ids=("hcb",),
                stacking_ids=("AA",),
            )
        )
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tfb",
            name="tfb",
            smiles=TFB,
            motif_kind="aldehyde",
            expected_connectivity=3,
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            summaries, _candidates, _attempted = generator.generate_pair_candidates(
                amine,
                aldehyde,
                out_dir=temp_dir,
                write_cif=True,
            )
            summary = summaries[0]
            self.assertEqual(summary.status, "ok")
            cif_text = Path(summary.cif_path).read_text(encoding="utf-8")
            original_c = float(re.search(r"_cell_length_c (\S+)", cif_text).group(1))
            squashed_c = 6.8
            self.assertGreater(original_c, squashed_c)

            # Keep cartesian coordinates fixed: fractional z scales by
            # original_c / squashed_c, pushing each layer through the
            # periodic image of the other (the shipped failure mode).
            atom_site_pattern = re.compile(
                r"^(\S+) (\S+) (\S+) (\S+) (\S+) (\S+)$"
            )
            squashed_lines = []
            for line in cif_text.splitlines():
                if line.startswith("_cell_length_c"):
                    squashed_lines.append(f"_cell_length_c {squashed_c:.6f}")
                    continue
                match = atom_site_pattern.match(line)
                if match and "L" in match.group(1):
                    label, symbol, fx, fy, fz, occupancy = match.groups()
                    scaled_z = float(fz) * original_c / squashed_c
                    squashed_lines.append(f"{label} {symbol} {fx} {fy} {scaled_z:.6f} {occupancy}")
                else:
                    squashed_lines.append(line)
            squashed_path = Path(temp_dir) / "squashed_AA.cif"
            squashed_path.write_text("\n".join(squashed_lines) + "\n", encoding="utf-8")

            record = {
                "structure_id": "squashed_AA",
                "pair_id": "squashed_AA",
                "status": "ok",
                "topology_id": "hcb",
                "score": None,
                "flags": [],
                "cif_path": str(squashed_path),
                "metadata": summary.metadata,
            }
            report = CoarseStructureValidator().validate_manifest_record(record)

        self.assertEqual(report.classification, "hard_invalid")
        self.assertIn("heavy_atom_clash", report.hard_invalid_reasons)
        self.assertGreater(report.metrics["n_heavy_atom_clash_pairs"], 0)
        self.assertLess(report.metrics["min_nonbonded_heavy_distance"], 2.5)

    def test_generate_monomer_pair_candidate_accepts_direct_monomer_specs(self):
        amine_record = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde_record = BatchMonomerRecord(
            id="tfb",
            name="tfb",
            smiles=TFB,
            motif_kind="aldehyde",
            expected_connectivity=3,
        )

        built_amine = self.generator.build_monomer(amine_record)
        built_aldehyde = self.generator.build_monomer(aldehyde_record)

        assert built_amine.monomer is not None
        assert built_aldehyde.monomer is not None
        summary, candidate = self.generator.generate_monomer_pair_candidate(
            built_aldehyde.monomer,
            built_amine.monomer,
        )

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_id, "tapb__tfb")
        self.assertEqual(summary.amine_record_id, "tapb")
        self.assertEqual(summary.aldehyde_record_id, "tfb")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "hcb")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "single-node-bipartite")

    def test_three_plus_two_node_linker_pair_builds_hcb_candidate(self):
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        summary, candidate = self.generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "3+2-node-linker")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "hcb")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "node-linker-single-node")

    def test_experimental_azine_pair_builds_hcb_candidate(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                allowed_reactions=("azine_bridge",),
                rdkit_num_conformers=2,
                retain_top_results=5,
                single_node_topology_ids=("hcb",),
            )
        )
        hydrazine = BatchMonomerRecord(
            id="hydrazine",
            name="hydrazine",
            smiles="NN",
            motif_kind="hydrazine",
            expected_connectivity=2,
        )
        tfb = BatchMonomerRecord(
            id="tfb",
            name="tfb",
            smiles=TFB,
            motif_kind="aldehyde",
            expected_connectivity=3,
        )

        summary, candidate = generator.generate_pair_candidate(hydrazine, tfb)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "3+2-node-linker")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "hcb")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "node-linker-single-node")
        self.assertGreater(candidate.metadata["graph_summary"]["n_reaction_events"], 0)
        self.assertEqual(candidate.metadata["graph_summary"]["n_monomer_instances"], 5)
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 6)
        self.assertNotIn("unreacted_motifs:1", candidate.flags)
        self.assertLess(self._max_bridge_distance_error(candidate), 0.02)
        self.assertLess(candidate.metadata["score_metadata"]["bridge_geometry_residual"], 0.2)

    def test_reverse_polarity_two_plus_three_pair_builds_hcb_candidate(self):
        amine = BatchMonomerRecord(
            id="ppd",
            name="ppd",
            smiles=PPD,
            motif_kind="amine",
            expected_connectivity=2,
        )
        aldehyde = BatchMonomerRecord(
            id="tfb",
            name="tfb",
            smiles=TFB,
            motif_kind="aldehyde",
            expected_connectivity=3,
        )

        summary, candidate = self.generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "3+2-node-linker")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "hcb")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "node-linker-single-node")
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 6)
        self.assertLess(self._max_bridge_distance_error(candidate), 0.02)
        self.assertLess(candidate.metadata["score_metadata"]["bridge_geometry_residual"], 0.2)

    def test_three_plus_two_pair_handles_asymmetric_tritopic_node(self):
        amine = BatchMonomerRecord(
            id="asymmetric_triamine",
            name="asymmetric_triamine",
            smiles=ASYMMETRIC_TRIAMINE,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="long_dialdehyde",
            name="long_dialdehyde",
            smiles=LONG_DIALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        summary, candidate = self.generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "3+2-node-linker")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["embedding"]["cell_kind"], "oblique")
        self.assertLess(self._max_bridge_distance_error(candidate), 0.02)
        self.assertLess(self._max_linker_translation_mismatch(candidate), 0.02)

    def test_three_plus_three_pair_handles_asymmetric_tritopic_nodes(self):
        amine = BatchMonomerRecord(
            id="asymmetric_triamine",
            name="asymmetric_triamine",
            smiles=ASYMMETRIC_TRIAMINE,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="asymmetric_trialdehyde",
            name="asymmetric_trialdehyde",
            smiles=ASYMMETRIC_TRIALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=3,
        )

        summary, candidate = self.generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "3+3-node-node")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["embedding"]["cell_kind"], "oblique")
        self.assertLess(self._max_bridge_distance_error(candidate), 0.02)

    def test_three_plus_two_pair_can_export_cif(self):
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            summary, candidate = self.generator.generate_pair_candidate(
                amine,
                aldehyde,
                out_dir=temp_dir,
                write_cif=True,
            )

            self.assertEqual(summary.status, "ok")
            self.assertIsNotNone(candidate)
            self.assertIsNotNone(summary.cif_path)
            self.assertIn(f"{Path(temp_dir) / 'valid'}", summary.cif_path)
            cif_text = Path(summary.cif_path).read_text(encoding="utf-8")

        self.assertIn("cofid", summary.metadata)
        self.assertEqual(cif_text.splitlines()[0], f"# COFid: {summary.metadata['cofid']}")
        self.assertIn("_geom_bond_atom_site_label_1", cif_text)
        self.assertIn("_atom_site_label", cif_text)

    def test_three_plus_two_pair_exports_cif_by_default_when_output_dir_is_provided(self):
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            summary, candidate = self.generator.generate_pair_candidate(
                amine,
                aldehyde,
                out_dir=temp_dir,
            )

            self.assertEqual(summary.status, "ok")
            self.assertIsNotNone(candidate)
            self.assertIsNotNone(summary.cif_path)
            self.assertTrue(Path(summary.cif_path).is_file())
            self.assertIn(f"{Path(temp_dir) / 'valid'}", summary.cif_path)

    def test_hard_hard_invalid_bridge_blocks_cif_export(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
                single_node_topology_ids=("hcb",),
                hard_hard_max_bridge_distance=1.0,
            )
        )
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            summary, candidate = generator.generate_pair_candidate(
                amine,
                aldehyde,
                out_dir=temp_dir,
                write_cif=True,
            )

        self.assertEqual(summary.status, "ok")
        self.assertIsNotNone(candidate)
        self.assertIsNone(summary.cif_path)
        self.assertEqual(
            tuple(summary.metadata["hard_hard_invalid_reasons"]),
            ("bridge_distance_exceeds_cif_export_limit",),
        )
        self.assertTrue(summary.metadata["cif_export_blocked"])

    def test_repairable_bridge_geometry_exports_to_needs_optimization_bucket(self):
        # A04: the verdict must come from the final exported coordinates.
        # The exported CIF is deliberately distorted after writing; the clean
        # seed metrics must not rescue it.
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
                single_node_topology_ids=("hcb",),
                hard_hard_max_bridge_distance=10.0,
                validation_thresholds=CoarseValidationThresholds(
                    hard_hard_max_bridge_distance=10.0,
                ),
            )
        )
        pristine_by_name: dict[str, str] = {}
        generator.cif_writer = _DistortingCifWriter(generator.cif_writer, pristine_by_name)
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            summary, candidate = generator.generate_pair_candidate(
                amine,
                aldehyde,
                out_dir=temp_dir,
                write_cif=True,
            )

        self.assertEqual(summary.status, "ok")
        self.assertIsNotNone(candidate)
        self.assertIsNotNone(summary.cif_path)
        assert summary.cif_path is not None
        self.assertIn(f"{Path(temp_dir) / 'needs_optimization'}", summary.cif_path)
        validation = summary.metadata["validation"]
        self.assertEqual(validation["classification"], "needs_optimization")
        self.assertIn("optimization_hint", validation)
        self.assertIn("bridge_distance_residual_max_hard", validation["needs_optimization_reasons"])
        # The verdict came from the measured final coordinates: the seed
        # assembly metrics were clean and are informational only.
        self.assertEqual(validation["metrics"]["bridge_metrics_source"], "final_cif")
        self.assertLessEqual(validation["metrics"]["seed_max_bridge_distance_residual"], 0.05)
        self.assertGreater(validation["metrics"]["max_bridge_distance_residual"], 1.0)
        self.assertEqual(validation["coverage"]["bridge_geometry"], "measured")
        self.assertEqual(summary.validation_classification, "needs_optimization")
        self.assertEqual(summary.validation_coverage.get("bridge_geometry"), "measured")
        self.assertEqual(summary.unmeasured_required_checks, ())

    def test_repair_geometry_runs_lammps_for_needs_optimization_candidate(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
                single_node_topology_ids=("hcb",),
                hard_hard_max_bridge_distance=10.0,
                repair_geometry=True,
                validation_thresholds=CoarseValidationThresholds(
                    hard_hard_max_bridge_distance=10.0,
                ),
            )
        )
        pristine_by_name: dict[str, str] = {}
        generator.cif_writer = _DistortingCifWriter(generator.cif_writer, pristine_by_name)
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        def fake_optimize(cif_path, *, output_dir, lmp_path, settings, timeout_seconds, eqeq_timeout_seconds):
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            optimized_cif = output_path / "optimized.cif"
            # The fake repair restores the pristine (pre-distortion) export.
            optimized_cif.write_text(pristine_by_name[Path(cif_path).name], encoding="utf-8")
            return SimpleNamespace(
                optimized_cif=str(optimized_cif),
                output_dir=str(output_dir),
                report_path=str(Path(output_dir) / "lammps_report.json"),
                settings=settings,
                warnings=(),
            )

        with tempfile.TemporaryDirectory() as temp_dir, patch("cofkit.batch.optimize_cif_with_lammps", side_effect=fake_optimize) as optimized:
            summary, candidate = generator.generate_pair_candidate(
                amine,
                aldehyde,
                out_dir=temp_dir,
                write_cif=True,
            )

        self.assertEqual(summary.status, "ok")
        self.assertIsNotNone(candidate)
        optimized.assert_called_once()
        validation = summary.metadata["validation"]
        self.assertEqual(validation["classification"], "needs_optimization")
        self.assertEqual(validation["metrics"]["bridge_metrics_source"], "final_cif")
        repair = validation["geometry_repair"]
        self.assertEqual(repair["status"], "ok")
        self.assertEqual(repair["engine"], "lammps")
        self.assertEqual(repair["settings"]["forcefield"], "dreiding")
        self.assertEqual(repair["settings"]["charge_model"], "none")
        self.assertEqual(repair["settings"]["pre_minimization_mode"], "soft")
        post_repair = repair["post_repair_validation"]
        self.assertEqual(post_repair["classification"], "valid")
        self.assertTrue(post_repair["became_valid"])
        # The repaired CIF is re-measured from its own coordinates.
        self.assertEqual(post_repair["scope"], "optimized_cif_final_geometry_checks")
        self.assertEqual(post_repair["metrics"]["bridge_metrics_source"], "final_cif")
        self.assertLess(post_repair["metrics"]["max_bridge_distance_residual"], 0.1)
        self.assertEqual(post_repair["coverage"]["bridge_geometry"], "measured")

    def test_batch_geometry_repair_failures_are_counted_and_recorded(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
                single_node_topology_ids=("hcb",),
                enumerate_all_topologies=False,
                hard_hard_max_bridge_distance=10.0,
                repair_geometry=True,
                max_workers=1,
                validation_thresholds=CoarseValidationThresholds(
                    hard_hard_max_bridge_distance=10.0,
                    hard_max_bridge_distance_residual=0.0,
                    hard_mean_bridge_distance_residual=0.0,
                ),
            )
        )

        with tempfile.TemporaryDirectory() as temp_dir, patch(
            "cofkit.batch.optimize_cif_with_lammps",
            side_effect=RuntimeError("LAMMPS unavailable"),
        ):
            temp_root = Path(temp_dir)
            input_dir = temp_root / "input"
            output_dir = temp_root / "out"
            input_dir.mkdir()
            (input_dir / "amines_count_3.txt").write_text(f"smiles\n{TAPB}\n", encoding="utf-8")
            (input_dir / "aldehydes_count_2.txt").write_text(f"smiles\n{TEREPHTHALALDEHYDE}\n", encoding="utf-8")

            summary = generator.run_imine_batch(input_dir, output_dir, max_pairs=1, write_cif=True)

            assert summary.geometry_repair_failed_records_path is not None
            failure_path = Path(summary.geometry_repair_failed_records_path)
            failure_rows = [json.loads(line) for line in failure_path.read_text(encoding="utf-8").splitlines() if line]

        self.assertEqual(summary.geometry_repair_counts["failed"], 1)
        self.assertEqual(len(failure_rows), 1)
        repair = failure_rows[0]["metadata"]["validation"]["geometry_repair"]
        self.assertEqual(repair["status"], "failed")
        self.assertIn("RuntimeError: LAMMPS unavailable", repair["error"])

    def test_three_plus_two_pair_can_target_fes_single_node_topology(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
                single_node_topology_ids=("fes",),
            )
        )
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "fes")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "node-linker-single-node-expanded")
        self.assertEqual(candidate.metadata["graph_summary"]["n_monomer_instances"], 10)
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 12)

    def test_three_plus_two_pair_can_target_hca_single_node_topology(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
                single_node_topology_ids=("hca",),
            )
        )
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "hca")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "node-linker-single-node-expanded")
        self.assertEqual(candidate.metadata["graph_summary"]["n_monomer_instances"], 15)
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 18)

    def test_three_plus_two_pair_can_target_fxt_single_node_topology(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
                single_node_topology_ids=("fxt",),
            )
        )
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "fxt")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "node-linker-single-node-expanded")
        self.assertEqual(candidate.metadata["graph_summary"]["n_monomer_instances"], 30)
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 36)

    def test_three_plus_three_pair_can_target_fes_single_node_topology(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
                single_node_topology_ids=("fes",),
            )
        )
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tfb",
            name="tfb",
            smiles=TFB,
            motif_kind="aldehyde",
            expected_connectivity=3,
        )

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "fes")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "single-node-expanded-node-node")
        self.assertEqual(candidate.metadata["graph_summary"]["n_monomer_instances"], 16)
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 24)

    def test_three_plus_three_pair_can_target_fxt_single_node_topology(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
                single_node_topology_ids=("fxt",),
            )
        )
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tfb",
            name="tfb",
            smiles=TFB,
            motif_kind="aldehyde",
            expected_connectivity=3,
        )

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "fxt")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "single-node-expanded-node-node")
        self.assertEqual(candidate.metadata["graph_summary"]["n_monomer_instances"], 12)
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 18)

    def test_three_plus_three_pair_rejects_non_bipartite_single_node_topologies(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
                single_node_topology_ids=("hca",),
            )
        )
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tfb",
            name="tfb",
            smiles=TFB,
            motif_kind="aldehyde",
            expected_connectivity=3,
        )

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "generation-failed")
        self.assertIsNone(candidate)
        self.assertIn("not bipartite", summary.metadata["failed_topologies"]["hca"])

    def test_default_topology_selection_respects_single_node_category_rules(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
            )
        )
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )

        three_plus_two_summary, _ = generator.generate_pair_candidate(
            amine,
            BatchMonomerRecord(
                id="tpal",
                name="tpal",
                smiles=TEREPHTHALALDEHYDE,
                motif_kind="aldehyde",
                expected_connectivity=2,
            ),
        )
        three_plus_three_summary, _ = generator.generate_pair_candidate(
            amine,
            BatchMonomerRecord(
                id="tfb",
                name="tfb",
                smiles=TFB,
                motif_kind="aldehyde",
                expected_connectivity=3,
            ),
        )

        self.assertEqual(
            three_plus_two_summary.metadata["topology_selection"]["available_topologies"],
            ("hcb", "hca", "fes", "fxt", "srs"),
        )
        self.assertEqual(
            three_plus_three_summary.metadata["topology_selection"]["available_topologies"],
            ("hcb", "fes", "fxt", "srs"),
        )

    def test_generate_pair_candidates_returns_all_supported_topologies(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
            )
        )
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        summaries, candidates, attempted_structures = generator.generate_pair_candidates(amine, aldehyde)

        self.assertEqual(attempted_structures, 5)
        self.assertEqual(len(summaries), 5)
        self.assertEqual(len(candidates), 5)
        self.assertEqual(
            {summary.topology_id for summary in summaries},
            {"hcb", "hca", "fes", "fxt", "srs"},
        )
        instance_counts = {
            candidate.metadata["net_plan"]["topology"]: candidate.metadata["graph_summary"]["n_monomer_instances"]
            for candidate in candidates
        }
        self.assertEqual(instance_counts, {"hcb": 5, "hca": 15, "fes": 10, "fxt": 30, "srs": 20})
        self.assertTrue(all(summary.cif_path is None for summary in summaries))

    def test_stacking_expansion_failure_isolates_failing_candidate(self):
        from cofkit import batch as batch_module

        real_expand = batch_module.enumerate_candidate_stackings

        def flaky_expand(candidate, *, registry_ids=(), monomer_specs=None):
            net_plan = candidate.metadata.get("net_plan", {})
            if isinstance(net_plan, dict) and net_plan.get("topology") == "hcb":
                raise RuntimeError("boom")
            return real_expand(candidate, registry_ids=registry_ids, monomer_specs=monomer_specs)

        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
                stacking_ids=("AA",),
            )
        )
        amine = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        stderr = io.StringIO()
        with (
            patch.object(batch_module, "enumerate_candidate_stackings", flaky_expand),
            contextlib.redirect_stderr(stderr),
        ):
            summaries, candidates, attempted_structures = generator.generate_pair_candidates(amine, aldehyde)

        self.assertEqual(attempted_structures, 5)
        self.assertEqual(len(summaries), 4)
        self.assertEqual(len(candidates), 4)
        self.assertEqual(
            {summary.topology_id for summary in summaries},
            {"hca", "fes", "fxt", "srs"},
        )
        self.assertTrue(all(summary.status == "ok" for summary in summaries))
        for summary in summaries:
            errors = summary.metadata["stacking_expansion_errors"]
            self.assertEqual(len(errors), 1)
            self.assertEqual(list(errors.values()), ["RuntimeError: boom"])
        self.assertIn("warning: stacking expansion failed for candidate", stderr.getvalue())
        self.assertIn("RuntimeError: boom", stderr.getvalue())

    def test_four_plus_two_pair_uses_expanded_default_topology_pool(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=5,
            )
        )
        amine = BatchMonomerRecord(
            id="tetra_amine",
            name="tetra_amine",
            smiles=TETRA_AMINE,
            motif_kind="amine",
            expected_connectivity=4,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "4+2-node-linker")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(
            set(summary.metadata["topology_selection"]["available_topologies"]),
            {"dia", "pts", "lon", "qtz"},
        )
        self.assertIn(candidate.metadata["net_plan"]["topology"], {"dia", "pts", "lon", "qtz"})
        self.assertGreater(candidate.metadata["graph_summary"]["n_reaction_events"], 0)

    def test_four_plus_two_pair_returns_explicit_2d_topologies_when_requested(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=5,
                single_node_topology_ids=("sql", "kgm", "htb"),
                shape_aware_topology_filter=False,
            )
        )
        amine = BatchMonomerRecord(
            id="tetra_amine",
            name="tetra_amine",
            smiles=TETRA_AMINE,
            motif_kind="amine",
            expected_connectivity=4,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        summaries, candidates, attempted_structures = generator.generate_pair_candidates(amine, aldehyde)

        self.assertEqual(attempted_structures, 3)
        self.assertEqual({summary.topology_id for summary in summaries}, {"sql", "kgm", "htb"})
        self.assertTrue(all(summary.pair_mode == "4+2-node-linker" for summary in summaries))
        instance_counts = {
            candidate.metadata["net_plan"]["topology"]: candidate.metadata["graph_summary"]["n_monomer_instances"]
            for candidate in candidates
        }
        self.assertEqual(instance_counts, {"sql": 3, "kgm": 9, "htb": 18})

    def _d2h_pair_records(self):
        return (
            BatchMonomerRecord(
                id="d2h_tetra_amine",
                name="d2h_tetra_amine",
                smiles=D2H_TETRA_AMINE,
                motif_kind="amine",
                expected_connectivity=4,
            ),
            BatchMonomerRecord(
                id="tpal",
                name="tpal",
                smiles=TEREPHTHALALDEHYDE,
                motif_kind="aldehyde",
                expected_connectivity=2,
            ),
        )

    def test_explicit_topology_shape_conflict_is_a_warning_not_a_skip(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=1,
                topology_ids=("sql",),
                write_cif=False,
            )
        )
        amine, aldehyde = self._d2h_pair_records()

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.metadata["reactant_node_shapes"]["amine"], "rectangular")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "sql")
        warnings = summary.metadata["shape_warnings"]
        self.assertEqual(len(warnings), 1)
        self.assertIn("requested topology 'sql'", warnings[0])
        self.assertIn("rectangular", warnings[0])

    def test_enumerated_topologies_still_filter_shape_conflicts(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(rdkit_num_conformers=1, retain_top_results=1, write_cif=False)
        )
        amine, aldehyde = self._d2h_pair_records()

        summary, _candidate = generator.generate_pair_candidate(amine, aldehyde)

        available = set(summary.metadata["topology_selection"]["available_topologies"])
        self.assertNotIn("sql", available)
        self.assertIn("kgm", available)
        self.assertEqual(summary.metadata["shape_warnings"], ())

    def test_shape_aware_topology_filter_off_enumerates_conflicting_topologies(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=1,
                write_cif=False,
                shape_aware_topology_filter=False,
            )
        )
        amine, aldehyde = self._d2h_pair_records()

        summary, _candidate = generator.generate_pair_candidate(amine, aldehyde)

        available = set(summary.metadata["topology_selection"]["available_topologies"])
        self.assertIn("sql", available)
        self.assertEqual(summary.metadata["shape_warnings"], ())

    def test_shape_filtered_internal_enumeration_drops_conflicting_topology(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=1,
                topology_ids=("sql",),
                write_cif=False,
                shape_filter_explicit_topologies=True,
            )
        )
        amine, aldehyde = self._d2h_pair_records()

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "generation-failed")
        self.assertIsNone(candidate)
        self.assertEqual(summary.metadata["shape_warnings"], ())
        self.assertIn("rectangular", summary.metadata["failed_topologies"]["sql"])

    def test_explicit_decorated_bex_request_warns_about_non_rectangular_node(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=1,
                topology_ids=("bex",),
                write_cif=False,
            )
        )
        amine = BatchMonomerRecord(
            id="tetra_amine",
            name="tetra_amine",
            smiles=TETRA_AMINE,
            motif_kind="amine",
            expected_connectivity=4,
        )
        aldehyde = BatchMonomerRecord(
            id="tetra_aldehyde",
            name="tetra_aldehyde",
            smiles=TETRA_ALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=4,
        )

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "decorated-bex-node-node")
        warnings = summary.metadata["shape_warnings"]
        self.assertEqual(len(warnings), 1)
        self.assertIn("decorated bex", warnings[0])
        self.assertIn("tetrahedral", warnings[0])

    def test_four_plus_two_dia_imine_batch_build_keeps_bent_bridge_geometry(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=1,
                target_dimensionality="3D",
                topology_ids=("dia",),
                write_cif=False,
            )
        )
        amine_record = BatchMonomerRecord(
            id="amines_count_2_0188",
            name="amines_count_2_0188",
            smiles="C1=CC2=C(C=C1N)C(=O)C3=C(C=C(C=C3)N)C2=O",
            motif_kind="amine",
            expected_connectivity=2,
        )
        aldehyde_record = BatchMonomerRecord(
            id="aldehydes_count_4_0008",
            name="aldehydes_count_4_0008",
            smiles="C1=C(C=CC(=C1)C(C2=CC=C(C=C2)C=O)(C3=CC=C(C=C3)C=O)C4=CC=C(C=C4)C=O)C=O",
            motif_kind="aldehyde",
            expected_connectivity=4,
        )

        built_amine = generator.build_monomer(amine_record)
        built_aldehyde = generator.build_monomer(aldehyde_record)
        assert built_amine.monomer is not None
        assert built_aldehyde.monomer is not None

        summary, candidate = generator.generate_monomer_pair_candidate(
            built_amine.monomer,
            built_aldehyde.monomer,
        )

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.topology_id, "dia")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "node-linker-single-node-3d")

        realizer = ReactionRealizer()
        realization = realizer.realize(
            candidate,
            {
                built_amine.monomer.id: built_amine.monomer,
                built_aldehyde.monomer.id: built_aldehyde.monomer,
            },
            candidate.metadata["instance_to_monomer"],
        )

        self.assertIsNotNone(realization)
        assert realization is not None
        realized_by_instance = {
            instance_id: {atom.atom_id: atom.local_position for atom in atoms}
            for instance_id, atoms in realization.atoms_by_instance.items()
        }

        min_anchor_c_n = float("inf")
        min_h_c_n = float("inf")
        min_c_n_anchor = float("inf")
        for event in candidate.events:
            first_ref, second_ref = event.participants
            if first_ref.monomer_id == built_amine.monomer.id:
                amine_ref, aldehyde_ref = first_ref, second_ref
            else:
                amine_ref, aldehyde_ref = second_ref, first_ref
            amine_motif = built_amine.monomer.motif_by_id(amine_ref.motif_id)
            aldehyde_motif = built_aldehyde.monomer.motif_by_id(aldehyde_ref.motif_id)
            carbon_atom_id = aldehyde_motif.metadata["reactive_atom_id"]
            carbon_anchor_atom_id = aldehyde_motif.metadata["anchor_atom_id"]
            nitrogen_atom_id = amine_motif.metadata["reactive_atom_id"]
            nitrogen_anchor_atom_id = amine_motif.metadata["anchor_atom_id"]
            retained_hydrogen_atom_id = realizer._attached_hydrogen_ids_for_parent(
                built_aldehyde.monomer,
                aldehyde_motif,
                carbon_atom_id,
                set(realization.removed_atom_ids_by_instance[aldehyde_ref.monomer_instance_id]),
            )[0]

            carbon_world = self._world_position(
                candidate,
                aldehyde_ref.monomer_instance_id,
                realized_by_instance[aldehyde_ref.monomer_instance_id][carbon_atom_id],
                image=aldehyde_ref.periodic_image,
            )
            carbon_anchor_world = self._world_position(
                candidate,
                aldehyde_ref.monomer_instance_id,
                realized_by_instance[aldehyde_ref.monomer_instance_id][carbon_anchor_atom_id],
                image=aldehyde_ref.periodic_image,
            )
            hydrogen_world = self._world_position(
                candidate,
                aldehyde_ref.monomer_instance_id,
                realized_by_instance[aldehyde_ref.monomer_instance_id][retained_hydrogen_atom_id],
                image=aldehyde_ref.periodic_image,
            )
            nitrogen_world = self._world_position(
                candidate,
                amine_ref.monomer_instance_id,
                realized_by_instance[amine_ref.monomer_instance_id][nitrogen_atom_id],
                image=amine_ref.periodic_image,
            )
            nitrogen_anchor_world = self._world_position(
                candidate,
                amine_ref.monomer_instance_id,
                realized_by_instance[amine_ref.monomer_instance_id][nitrogen_anchor_atom_id],
                image=amine_ref.periodic_image,
            )

            min_anchor_c_n = min(min_anchor_c_n, self._angle(carbon_anchor_world, carbon_world, nitrogen_world))
            min_h_c_n = min(min_h_c_n, self._angle(hydrogen_world, carbon_world, nitrogen_world))
            min_c_n_anchor = min(min_c_n_anchor, self._angle(carbon_world, nitrogen_world, nitrogen_anchor_world))

        # Post-rework the bridge is constructed at the 120-degree template
        # priors; the dia 3D placement lands ~0.03 angstrom off the
        # prior-consistent span, so the best-effort closure balances both
        # angle residuals a few degrees below/above the priors (measured
        # 117.0-118.5 across the eight events) instead of the old
        # displacement-dominated fit's ~131 degrees.
        self.assertGreater(min_anchor_c_n, 112.0)
        self.assertLess(min_anchor_c_n, 128.0)
        self.assertGreater(min_h_c_n, 108.0)
        self.assertLess(min_h_c_n, 125.0)
        self.assertGreater(min_c_n_anchor, 112.0)
        self.assertLess(min_c_n_anchor, 128.0)

    def test_reverse_polarity_two_plus_four_pair_returns_explicit_2d_topologies_when_requested(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=5,
                single_node_topology_ids=("sql", "kgm", "htb"),
                shape_aware_topology_filter=False,
            )
        )
        amine = BatchMonomerRecord(
            id="ppd",
            name="ppd",
            smiles=PPD,
            motif_kind="amine",
            expected_connectivity=2,
        )
        aldehyde = BatchMonomerRecord(
            id="tetra_aldehyde",
            name="tetra_aldehyde",
            smiles=TETRA_ALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=4,
        )

        summaries, candidates, attempted_structures = generator.generate_pair_candidates(amine, aldehyde)

        self.assertEqual(attempted_structures, 3)
        self.assertEqual({summary.topology_id for summary in summaries}, {"sql", "kgm", "htb"})
        self.assertTrue(all(summary.pair_mode == "4+2-node-linker" for summary in summaries))
        self.assertEqual({candidate.metadata["net_plan"]["topology"] for candidate in candidates}, {"sql", "kgm", "htb"})

    def test_four_plus_four_pair_uses_expanded_default_topology_pool(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=5,
            )
        )
        amine = BatchMonomerRecord(
            id="tetra_amine",
            name="tetra_amine",
            smiles=TETRA_AMINE,
            motif_kind="amine",
            expected_connectivity=4,
        )
        aldehyde = BatchMonomerRecord(
            id="tetra_aldehyde",
            name="tetra_aldehyde",
            smiles=TETRA_ALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=4,
        )

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "4+4-node-node")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(
            set(summary.metadata["topology_selection"]["available_topologies"]),
            {"dia", "pts", "lon", "qtz"},
        )
        self.assertIn(candidate.metadata["net_plan"]["topology"], {"dia", "pts", "lon", "qtz"})
        self.assertGreater(candidate.metadata["graph_summary"]["n_reaction_events"], 0)

    def test_four_plus_four_pair_can_target_sql_bipartite_topology(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=5,
                single_node_topology_ids=("sql",),
                shape_aware_topology_filter=False,
            )
        )
        amine = BatchMonomerRecord(
            id="tetra_amine",
            name="tetra_amine",
            smiles=TETRA_AMINE,
            motif_kind="amine",
            expected_connectivity=4,
        )
        aldehyde = BatchMonomerRecord(
            id="tetra_aldehyde",
            name="tetra_aldehyde",
            smiles=TETRA_ALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=4,
        )

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "4+4-node-node")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "sql")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "single-node-expanded-node-node")
        self.assertEqual(candidate.metadata["graph_summary"]["n_monomer_instances"], 4)
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 8)

    def test_four_plus_four_pair_can_target_decorated_bex_topology(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=5,
                topology_ids=("bex",),
                write_cif=False,
                shape_aware_topology_filter=False,
            )
        )
        amine = BatchMonomerRecord(
            id="tetra_amine",
            name="tetra_amine",
            smiles=TETRA_AMINE,
            motif_kind="amine",
            expected_connectivity=4,
        )
        aldehyde = BatchMonomerRecord(
            id="tetra_aldehyde",
            name="tetra_aldehyde",
            smiles=TETRA_ALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=4,
        )

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "4+4-node-node")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "bex")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "decorated-bex-node-node")
        self.assertEqual(candidate.metadata["graph_summary"]["n_monomer_instances"], 2)
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 4)
        self.assertEqual(
            set(candidate.metadata["assignment"]),
            {"bex_c4_node", "bex_c3_pair"},
        )
        self.assertEqual(
            sorted(event.participants[1].periodic_image for event in candidate.events),
            [(-1, -1, 0), (-1, 0, 0), (0, -1, 0), (0, 0, 0)],
        )
        self.assertIn("batch_decorated_bex_node_node", candidate.flags)

    def test_reported_d2h_bex_smiles_build_decorated_topology(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=5,
                topology_ids=("bex",),
                write_cif=False,
                shape_aware_topology_filter=False,
            )
        )
        amine = BatchMonomerRecord(
            id="bex_d2h_amine",
            name="bex_d2h_amine",
            smiles=BEX_D2H_AMINE,
            motif_kind="amine",
            expected_connectivity=4,
        )
        aldehyde = BatchMonomerRecord(
            id="bex_d2h_aldehyde",
            name="bex_d2h_aldehyde",
            smiles=BEX_D2H_ALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=4,
        )

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "bex")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "decorated-bex-node-node")
        self.assertEqual(
            candidate.metadata["net_plan"]["metadata"]["decorated_topology_unit"],
            "c3-edge-c3",
        )
        layer_z_span = candidate.metadata["embedding"]["layer_z_span"]
        self.assertGreater(layer_z_span, 0.0)
        # A7: span provenance is reported honestly (shared measurer, layer
        # normal axis).
        self.assertEqual(candidate.metadata["embedding"]["layer_z_span_axis"], "layer_normal")
        self.assertIn(
            candidate.metadata["embedding"]["layer_z_span_mode"],
            ("atomistic_product", "precursor_coordinates"),
        )
        self.assertAlmostEqual(
            candidate.state.cell[2][2],
            layer_z_span + generator.config.embedding_config.default_layer_spacing,
            places=5,
        )
        self.assertEqual(candidate.state.cell[2][0], 0.0)
        self.assertEqual(candidate.state.cell[2][1], 0.0)

    def test_decorated_bex_pair_can_enumerate_and_export_stacked_variant(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=5,
                topology_ids=("bex",),
                stacking_ids=("AA",),
                shape_aware_topology_filter=False,
                write_cif=True,
                hard_hard_max_bridge_distance=10.0,
                validation_thresholds=CoarseValidationThresholds(hard_hard_max_bridge_distance=10.0),
            )
        )
        amine = BatchMonomerRecord(
            id="bex_d2h_amine",
            name="bex_d2h_amine",
            smiles=BEX_D2H_AMINE,
            motif_kind="amine",
            expected_connectivity=4,
        )
        aldehyde = BatchMonomerRecord(
            id="bex_d2h_aldehyde",
            name="bex_d2h_aldehyde",
            smiles=BEX_D2H_ALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=4,
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            summaries, candidates, attempted_structures = generator.generate_pair_candidates(
                amine,
                aldehyde,
                out_dir=temp_dir,
                write_cif=True,
            )

            self.assertEqual(attempted_structures, 1)
            self.assertEqual(len(summaries), 1)
            self.assertEqual(len(candidates), 1)
            summary = summaries[0]
            candidate = candidates[0]
            self.assertEqual(summary.status, "ok")
            self.assertEqual(summary.structure_id, "bex_d2h_amine__bex_d2h_aldehyde__bex__AA")
            self.assertEqual(candidate.state.stacking_state, "AA")
            self.assertEqual(summary.metadata["stacking"]["id"], "AA")
            self.assertEqual(summary.metadata["stacking"]["comment_suffix"], "stacking=AA")
            self.assertTrue(candidate.metadata["embedding"]["stacking_enabled"])
            self.assertIn("stacked_2d", candidate.flags)
            self.assertEqual(candidate.metadata["graph_summary"]["n_monomer_instances"], 4)
            self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 8)
            layer_z_span = float(summary.metadata["stacking"]["layer_z_span"])
            self.assertGreater(layer_z_span, 0.0)
            self.assertEqual(summary.metadata["stacking"]["layer_z_span_axis"], "layer_normal")
            self.assertAlmostEqual(
                candidate.state.cell[2][2],
                2.0 * (layer_z_span + summary.metadata["stacking"]["interlayer_clearance"]),
                places=5,
            )
            self.assertIsNotNone(summary.cif_path)
            cif_text = Path(summary.cif_path).read_text(encoding="utf-8")

        self.assertIn("stacking=AA", cif_text.splitlines()[0])

    def test_six_plus_two_pair_uses_expanded_default_topology_pool(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=5,
            )
        )
        amine = BatchMonomerRecord(
            id="hexa_amine",
            name="hexa_amine",
            smiles=HEXA_AMINE,
            motif_kind="amine",
            expected_connectivity=6,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "6+2-node-linker")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(
            set(summary.metadata["topology_selection"]["available_topologies"]),
            {"pcu", "hxl", "acs"},
        )
        self.assertIn(candidate.metadata["net_plan"]["topology"], {"pcu", "hxl", "acs"})
        self.assertGreater(candidate.metadata["graph_summary"]["n_reaction_events"], 0)

    def test_six_plus_two_pair_can_target_hxl_when_requested(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=5,
                single_node_topology_ids=("hxl",),
            )
        )
        amine = BatchMonomerRecord(
            id="hexa_amine",
            name="hexa_amine",
            smiles=HEXA_AMINE,
            motif_kind="amine",
            expected_connectivity=6,
        )
        aldehyde = BatchMonomerRecord(
            id="tpal",
            name="tpal",
            smiles=TEREPHTHALALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=2,
        )

        summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertEqual(summary.pair_mode, "6+2-node-linker")
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate.metadata["net_plan"]["topology"], "hxl")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "node-linker-single-node-expanded")
        self.assertEqual(candidate.metadata["graph_summary"]["n_monomer_instances"], 4)
        self.assertEqual(candidate.metadata["graph_summary"]["n_reaction_events"], 6)

    def test_extended_topology_selection_respects_single_node_category_rules(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=5,
            )
        )
        amine_4 = BatchMonomerRecord(
            id="tetra_amine",
            name="tetra_amine",
            smiles=TETRA_AMINE,
            motif_kind="amine",
            expected_connectivity=4,
        )
        aldehyde_4 = BatchMonomerRecord(
            id="tetra_aldehyde",
            name="tetra_aldehyde",
            smiles=TETRA_ALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=4,
        )
        amine_6 = BatchMonomerRecord(
            id="hexa_amine",
            name="hexa_amine",
            smiles=HEXA_AMINE,
            motif_kind="amine",
            expected_connectivity=6,
        )

        four_plus_two_summary, _ = generator.generate_pair_candidate(
            amine_4,
            BatchMonomerRecord(
                id="tpal",
                name="tpal",
                smiles=TEREPHTHALALDEHYDE,
                motif_kind="aldehyde",
                expected_connectivity=2,
            ),
        )
        four_plus_four_summary, _ = generator.generate_pair_candidate(amine_4, aldehyde_4)
        six_plus_two_summary, _ = generator.generate_pair_candidate(
            amine_6,
            BatchMonomerRecord(
                id="tpal",
                name="tpal",
                smiles=TEREPHTHALALDEHYDE,
                motif_kind="aldehyde",
                expected_connectivity=2,
            ),
        )

        self.assertEqual(
            set(four_plus_two_summary.metadata["topology_selection"]["available_topologies"]),
            {"dia", "pts", "lon", "qtz"},
        )
        self.assertEqual(
            set(four_plus_four_summary.metadata["topology_selection"]["available_topologies"]),
            {"dia", "pts", "lon", "qtz"},
        )
        self.assertEqual(
            set(six_plus_two_summary.metadata["topology_selection"]["available_topologies"]),
            {"pcu", "hxl", "acs"},
        )

    def test_default_selector_includes_curated_compatible_topologies(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=5,
            )
        )

        three_plus_two = generator._topology_ids_for_pair(connectivities=(3, 2), pair_mode="node-linker")
        three_plus_three = generator._topology_ids_for_pair(connectivities=(3, 3), pair_mode="node-node")
        three_plus_four = generator._topology_ids_for_pair(connectivities=(3, 4), pair_mode="node-node")
        three_plus_six = generator._topology_ids_for_pair(connectivities=(3, 6), pair_mode="node-node")
        four_plus_two = generator._topology_ids_for_pair(connectivities=(4, 2), pair_mode="node-linker")
        four_plus_four = generator._topology_ids_for_pair(connectivities=(4, 4), pair_mode="node-node")
        six_plus_two = generator._topology_ids_for_pair(connectivities=(6, 2), pair_mode="node-linker")
        six_plus_six = generator._topology_ids_for_pair(connectivities=(6, 6), pair_mode="node-node")

        self.assertIn("srs", three_plus_two)
        self.assertIn("srs", three_plus_three)
        self.assertIn("ctn", three_plus_four)
        self.assertIn("bor", three_plus_four)
        self.assertIn("tbo", three_plus_four)
        self.assertIn("kgd", three_plus_six)
        self.assertIn("sql", four_plus_two)
        self.assertIn("kgm", four_plus_two)
        self.assertIn("pts", four_plus_two)
        self.assertIn("lon", four_plus_two)
        self.assertIn("qtz", four_plus_two)
        self.assertIn("sql", four_plus_four)
        self.assertIn("pts", four_plus_four)
        self.assertIn("lon", four_plus_four)
        self.assertIn("qtz", four_plus_four)
        self.assertIn("hxl", six_plus_two)
        self.assertIn("acs", six_plus_two)
        self.assertIn("acs", six_plus_six)

    def test_file_driven_batch_runner_writes_manifest(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_root = Path(temp_dir)
            input_dir = temp_root / "input"
            output_dir = temp_root / "out"
            input_dir.mkdir()
            (input_dir / "amines_count_2.txt").write_text("smiles\nNc1ccc(N)cc1\n", encoding="utf-8")
            (input_dir / "amines_count_3.txt").write_text(f"smiles\n{TAPB}\n", encoding="utf-8")
            (input_dir / "aldehydes_count_2.txt").write_text(f"smiles\n{TEREPHTHALALDEHYDE}\n", encoding="utf-8")
            (input_dir / "aldehydes_count_3.txt").write_text(f"smiles\n{TFB}\n", encoding="utf-8")

            summary = self.generator.run_imine_batch(
                input_dir,
                output_dir,
                max_pairs=3,
                write_cif=False,
            )

            self.assertEqual(summary.attempted_pairs, 3)
            self.assertEqual(summary.successful_pairs, 3)
            self.assertTrue((output_dir / "manifest.jsonl").exists())
            self.assertTrue((output_dir / "summary.md").exists())
            self.assertEqual(sum(1 for _ in (output_dir / "manifest.jsonl").open()), 3)

    def test_file_driven_batch_runner_parallelizes_pair_generation(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
                max_workers=2,
            )
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_root = Path(temp_dir)
            input_dir = temp_root / "input"
            output_dir = temp_root / "out"
            input_dir.mkdir()
            (input_dir / "amines_count_2.txt").write_text("smiles\nNc1ccc(N)cc1\n", encoding="utf-8")
            (input_dir / "amines_count_3.txt").write_text(f"smiles\n{TAPB}\n", encoding="utf-8")
            (input_dir / "aldehydes_count_2.txt").write_text(f"smiles\n{TEREPHTHALALDEHYDE}\n", encoding="utf-8")
            (input_dir / "aldehydes_count_3.txt").write_text(f"smiles\n{TFB}\n", encoding="utf-8")

            summary = generator.run_imine_batch(
                input_dir,
                output_dir,
                max_pairs=4,
                write_cif=False,
            )

            self.assertEqual(summary.attempted_pairs, 3)
            self.assertEqual(summary.successful_pairs, 3)
            self.assertTrue((output_dir / "manifest.jsonl").exists())
            self.assertEqual(
                len(
                    [
                        line
                        for line in (output_dir / "manifest.jsonl").read_text(encoding="utf-8").splitlines()
                        if line
                    ]
                ),
                summary.successful_structures,
            )

    def test_file_driven_batch_runner_enumerates_all_topologies_by_default(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=2,
                retain_top_results=5,
            )
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_root = Path(temp_dir)
            input_dir = temp_root / "input"
            output_dir = temp_root / "out"
            input_dir.mkdir()
            (input_dir / "amines_count_2.txt").write_text("smiles\n", encoding="utf-8")
            (input_dir / "amines_count_3.txt").write_text(f"smiles\n{TAPB}\n", encoding="utf-8")
            (input_dir / "aldehydes_count_2.txt").write_text(f"smiles\n{TEREPHTHALALDEHYDE}\n", encoding="utf-8")
            (input_dir / "aldehydes_count_3.txt").write_text("smiles\n", encoding="utf-8")

            summary = generator.run_imine_batch(
                input_dir,
                output_dir,
            )

            manifest_rows = [line for line in (output_dir / "manifest.jsonl").read_text(encoding="utf-8").splitlines() if line]
            cif_count = sum(1 for _ in (output_dir / "cifs").rglob("*.cif"))

            self.assertEqual(summary.attempted_pairs, 1)
            self.assertEqual(summary.successful_pairs, 1)
            self.assertEqual(summary.attempted_structures, 5)
            self.assertEqual(summary.successful_structures, 5)
            self.assertEqual(len(manifest_rows), 5)
            self.assertEqual(cif_count, 5)
            self.assertTrue((output_dir / "cifs" / "valid").is_dir())

    def test_file_driven_batch_runner_discovers_four_and_six_connectivity_libraries(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                retain_top_results=5,
            )
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_root = Path(temp_dir)
            input_dir = temp_root / "input"
            output_dir = temp_root / "out"
            input_dir.mkdir()
            (input_dir / "amines_count_2.txt").write_text("smiles\nNc1ccc(N)cc1\n", encoding="utf-8")
            (input_dir / "amines_count_4.txt").write_text(f"smiles\n{TETRA_AMINE}\n", encoding="utf-8")
            (input_dir / "amines_count_6.txt").write_text(f"smiles\n{HEXA_AMINE}\n", encoding="utf-8")
            (input_dir / "aldehydes_count_2.txt").write_text(f"smiles\n{TEREPHTHALALDEHYDE}\n", encoding="utf-8")
            (input_dir / "aldehydes_count_4.txt").write_text(f"smiles\n{TETRA_ALDEHYDE}\n", encoding="utf-8")
            (input_dir / "aldehydes_count_6.txt").write_text(f"smiles\n{HEXA_ALDEHYDE}\n", encoding="utf-8")

            summary = generator.run_imine_batch(
                input_dir,
                output_dir,
                write_cif=False,
            )

            self.assertEqual(summary.attempted_pairs, 8)
            self.assertEqual(summary.successful_pairs, 8)
            self.assertGreaterEqual(summary.attempted_structures, 8)
            self.assertGreater(summary.attempted_structures, summary.successful_structures)
            self.assertGreater(summary.successful_structures, 8)
            self.assertIn("4+4-node-node", summary.mode_counts)
            self.assertIn("4+6-node-node", summary.mode_counts)
            self.assertIn("4+2-node-linker", summary.mode_counts)
            self.assertIn("6+2-node-linker", summary.mode_counts)
            self.assertIn("dia", summary.topology_counts)
            self.assertIn("pcu", summary.topology_counts)
            self.assertIn("cor", summary.topology_counts)

    def _max_bridge_distance_error(self, candidate) -> float:
        return max(
            abs(float(metric["actual_distance"]) - float(metric["target_distance"]))
            for metric in candidate.metadata["score_metadata"]["bridge_event_metrics"]
        )

    def _max_linker_translation_mismatch(self, candidate) -> float:
        return max(
            float(details["translation_mismatch"])
            for details in candidate.metadata["embedding"]["poses"].values()
            if details.get("role") == "linker"
        )


@unittest.skipIf(Chem is None, "RDKit is not available")
class ScoringModeTests(unittest.TestCase):
    """Candidates carry score=None and are ranked by geometry residual."""

    def _generate(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                write_cif=False,
            )
        )
        first = generator.infer_monomer_record(TAPB, record_id="tapb")
        second = generator.infer_monomer_record(TFB, record_id="tfb")
        return generator.generate_pair_candidates(first, second, write_cif=False)

    def test_candidates_are_unscored_and_rank_best_geometry_first(self):
        summaries, candidates, _attempted = self._generate()

        self.assertTrue(candidates)
        for candidate in candidates:
            self.assertIsNone(candidate.score)
        # TAPB+TFB: hcb has by far the best bridge geometry (mean per-event
        # residual ~0.006) but the fewest reaction events; the removed legacy
        # score ranked it last because it scaled with unit-cell event count.
        self.assertEqual(summaries[0].topology_id, "hcb")
        self.assertIsNone(summaries[0].score)
        mean_residuals = [
            candidate.metadata["score_metadata"]["bridge_geometry_residual"]
            / len(candidate.metadata["score_metadata"]["bridge_event_metrics"])
            for candidate in candidates
        ]
        self.assertEqual(mean_residuals, sorted(mean_residuals))


@unittest.skipIf(Chem is None, "RDKit is not available")
class BatchProcessPoolFallbackTests(unittest.TestCase):
    def _write_tiny_library(self, root: Path) -> None:
        root.mkdir(parents=True, exist_ok=True)
        (root / "amines_count_3.txt").write_text(f"smiles\n{TAPB}\n{TAPB}\n", encoding="utf-8")
        (root / "aldehydes_count_3.txt").write_text(f"smiles\n{TFB}\n{TFB}\n", encoding="utf-8")

    def _make_generator(self) -> BatchStructureGenerator:
        return BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                write_cif=False,
                max_workers=2,
                single_node_topology_ids=("hcb",),
            )
        )

    def _manifest_structure_ids(self, output_dir: Path) -> list[str]:
        manifest_path = output_dir / "manifest.jsonl"
        return [
            json.loads(line)["structure_id"]
            for line in manifest_path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]

    def test_process_pool_construction_failure_falls_back_to_threads(self):
        import cofkit.batch as batch_module

        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            self._write_tiny_library(root / "input")
            generator = self._make_generator()

            with patch(
                "cofkit.batch.ProcessPoolExecutor",
                side_effect=batch_module.BrokenProcessPool("boom"),
            ):
                summary = generator.run_binary_bridge_batch(
                    root / "input",
                    root / "output",
                    write_cif=False,
                )

            self.assertEqual(summary.attempted_pairs, 4)
            self.assertEqual(summary.successful_pairs, 4)
            structure_ids = self._manifest_structure_ids(root / "output")
            self.assertEqual(len(structure_ids), len(set(structure_ids)))

    def test_mid_run_process_pool_failure_does_not_duplicate_results(self):
        import cofkit.batch as batch_module

        class _ExplodingExecutor:
            def __init__(self, *args, **kwargs):
                pass

            def __enter__(self):
                return self

            def __exit__(self, *exc_info):
                return False

            def map(self, fn, tasks):
                tasks = tuple(tasks)
                yield fn(tasks[0])
                raise OSError("simulated mid-run pool failure")

        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            self._write_tiny_library(root / "input")
            generator = self._make_generator()

            batch_module._init_batch_process_worker(generator.config, generator.reaction_library)
            try:
                with patch("cofkit.batch.ProcessPoolExecutor", _ExplodingExecutor):
                    summary = generator.run_binary_bridge_batch(
                        root / "input",
                        root / "output",
                        write_cif=False,
                    )
            finally:
                batch_module._PROCESS_BATCH_GENERATOR = None

            self.assertEqual(summary.attempted_pairs, 4)
            self.assertEqual(summary.successful_pairs, 4)
            structure_ids = self._manifest_structure_ids(root / "output")
            self.assertEqual(len(structure_ids), 4)
            self.assertEqual(len(structure_ids), len(set(structure_ids)))


@unittest.skipIf(Chem is None, "RDKit is not available")
class BatchMonomerBuildFailureTests(unittest.TestCase):
    """A monomer build error is recorded against its record, never raised."""

    def _write_tiny_library(self, root: Path) -> None:
        root.mkdir(parents=True, exist_ok=True)
        (root / "amines_count_3.txt").write_text(f"smiles\n{TAPB}\n{TAPB}\n", encoding="utf-8")
        (root / "aldehydes_count_3.txt").write_text(f"smiles\n{TFB}\n{TFB}\n", encoding="utf-8")

    def test_aromaticity_restore_error_is_recorded_per_record(self):
        failing_id = "amines_count_3_0001"

        def builder(monomer_id, name, smiles, motif_kind, *, num_conformers, random_seed):
            if monomer_id == failing_id:
                raise AromaticityRestoreError(
                    f"monomer {monomer_id!r}: MMFF aromaticity unresolvable"
                )
            return build_rdkit_monomer(
                monomer_id,
                name,
                smiles,
                motif_kind,
                num_conformers=num_conformers,
                random_seed=random_seed,
            )

        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            self._write_tiny_library(root / "input")
            generator = BatchStructureGenerator(
                BatchGenerationConfig(
                    rdkit_num_conformers=1,
                    write_cif=False,
                    max_workers=2,
                    single_node_topology_ids=("hcb",),
                ),
                smiles_monomer_builder=builder,
            )

            summary = generator.run_binary_bridge_batch(
                root / "input",
                root / "output",
                write_cif=False,
            )
            summary_text = (root / "output" / "summary.md").read_text(encoding="utf-8")

        self.assertEqual(summary.failed_monomers, 1)
        self.assertEqual(set(summary.build_failures), {failing_id})
        self.assertIn("AromaticityRestoreError", summary.build_failures[failing_id])
        self.assertGreater(summary.successful_pairs, 0)
        self.assertIn("## Build failures", summary_text)
        self.assertIn("AromaticityRestoreError", summary_text)


class DegenerateVectorPolicyTests(unittest.TestCase):
    """Audit A10 follow-up: degenerate vectors raise instead of fabricating +x.

    Probe evidence (2026-09-26): the old silent (1,0,0) fallback never fired
    on any live path in optimizer/embedding/batch across real builds and the
    full test suite, so degenerate geometry is now a visible ValueError that
    the per-topology build try/except records against the failed record.
    """

    def test_batch_safe_normalize_raises_on_degenerate_vector(self):
        generator = BatchStructureGenerator(BatchGenerationConfig(write_cif=False))
        with self.assertRaises(ValueError) as ctx:
            generator._safe_normalize((0.0, 0.0, 0.0))
        self.assertIn("BatchStructureGenerator", str(ctx.exception))
        self.assertIn("degenerate", str(ctx.exception))


@unittest.skipIf(Chem is None, "RDKit is not available")
class ShapeAwareConformerGateTests(unittest.TestCase):
    """Universal activation of shape-aware conformer selection, gated to
    trigonal (3-connecting) monomers by A/B evidence on the default monomer
    library (2026-09-29): flexible tritopic monomers build clash-free hcb
    structures where energy-only selection produced hydrogen-clash warnings,
    while tetrahedral 4-connecting monomers (dia) got mildly worse because no
    planar conformer exists for a tetrahedral node."""

    TAM = "Nc1ccc(C(c2ccc(N)cc2)(c2ccc(N)cc2)c2ccc(N)cc2)cc1"

    def _build(self, smiles, kind, connectivity, **config_kwargs):
        generator = BatchStructureGenerator(BatchGenerationConfig(write_cif=False, **config_kwargs))
        record = BatchMonomerRecord(
            id="probe",
            name="probe",
            smiles=smiles,
            motif_kind=kind,
            expected_connectivity=connectivity,
        )
        built = generator.build_monomer(record)
        self.assertTrue(built.ok, getattr(built, "error", None))
        return built.monomer

    def test_tritopic_monomer_gets_shape_selected_conformer(self):
        monomer = self._build("Nc1cc(N)cc(N)c1", "amine", 3)
        self.assertEqual(monomer.metadata["conformer_selection"], "motif_shape")
        self.assertIn("selected_conformer_shape_score", monomer.metadata)

    def test_tetratopic_monomer_keeps_energy_selection(self):
        # Tetrahedral TAM: the planar regular-polygon target is wrong for a
        # tetrahedral node, so 4-connecting monomers stay on energy selection.
        monomer = self._build(self.TAM, "amine", 4)
        self.assertEqual(monomer.metadata["conformer_selection"], "energy")
        self.assertNotIn("selected_conformer_shape_score", monomer.metadata)

    def test_ditopic_monomer_keeps_energy_selection(self):
        monomer = self._build("Nc1ccc(N)cc1", "amine", 2)
        self.assertEqual(monomer.metadata["conformer_selection"], "energy")

    def test_opt_out_keeps_energy_selection_for_tritopic(self):
        monomer = self._build("Nc1cc(N)cc(N)c1", "amine", 3, shape_aware_conformer=False)
        self.assertEqual(monomer.metadata["conformer_selection"], "energy")
        self.assertNotIn("selected_conformer_shape_score", monomer.metadata)


if __name__ == "__main__":
    unittest.main()
