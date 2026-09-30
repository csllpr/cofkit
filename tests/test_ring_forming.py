import contextlib
import io
import json
import sys
import tempfile
import unittest
from collections import Counter
from dataclasses import replace
from math import acos, cos, degrees, pi, sin, sqrt
from pathlib import Path


import gemmi

from cofkit import COFEngine, COFProject, CoarseStructureValidator
from cofkit.build_workflows.ring_forming import RingFormationConfig, RingFormingStructureGenerator
from cofkit.chem.rdkit import SHAPE_SELECTION_MIN_CONFORMERS
from cofkit.chem.rdkit import build_rdkit_monomer
from cofkit.cif import CIFWriter
from cofkit.cli import main as cli_main
from cofkit.cofid import cofid_to_build_request, generate_cofid
from cofkit.geometry import Frame
from cofkit.model import AssemblyState, MonomerSpec, MotifRef, Pose, ReactionEvent, ReactiveMotif
from cofkit.reaction_realization import ReactionRealizer
from cofkit.ring_geometry import (
    RING_ATTACHMENT_IDEAL_ANGLE_DEGREES,
    ring_attachment_report,
    ring_geometry_profile,
    validate_ring_geometry,
)
from cofkit.stacking import DEFAULT_INTERLAYER_CLEARANCE_ANGSTROM


class RingFormingWorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cof1_precursor = build_rdkit_monomer(
            "cof1_precursor",
            "benzene-1,4-diboronic acid",
            "OB(O)c1ccc(B(O)O)cc1",
            "boronic_acid",
            num_conformers=2,
        )
        cls.ctf1_precursor = build_rdkit_monomer(
            "ctf1_precursor",
            "terephthalonitrile",
            "N#Cc1ccc(C#N)cc1",
            "nitrile",
            num_conformers=2,
        )

    def test_hcb_events_use_correct_periodic_incidence(self):
        result = RingFormingStructureGenerator().build(
            self.cof1_precursor,
            "boroxine_trimerization",
        )

        self.assertEqual(len(result.outcome.monomer_instances), 3)
        self.assertEqual(len(result.outcome.events), 2)
        self.assertEqual(result.outcome.consumed_count, 6)
        self.assertEqual(result.graph.unreacted_motifs({self.cof1_precursor.id: self.cof1_precursor}), ())
        for event in result.outcome.events:
            self.assertEqual(len(event.participants), 3)
            self.assertEqual(len({ref.monomer_instance_id for ref in event.participants}), 3)
        self.assertEqual(result.candidate.metadata["net_plan"]["topology"], "hcb")
        self.assertEqual(result.candidate.metadata["ring_validation"]["classification"], "accepted")
        self.assertLess(
            result.candidate.metadata["score_metadata"]["ring_geometry"]["total_residual"],
            1e-3,
        )

    def test_ring_geometry_validator_rejects_displaced_precursor(self):
        candidate = RingFormingStructureGenerator().generate(
            self.ctf1_precursor,
            "triazine_trimerization",
        )
        poses = dict(candidate.state.monomer_poses)
        original = poses["p1"]
        poses["p1"] = replace(
            original,
            translation=(original.translation[0] + 1.0, original.translation[1], original.translation[2]),
        )
        displaced_state = replace(candidate.state, monomer_poses=poses)

        validation = validate_ring_geometry(
            candidate.events,
            displaced_state,
            {self.ctf1_precursor.id: self.ctf1_precursor},
        )

        self.assertEqual(validation.classification, "rejected")
        self.assertTrue(any("radial residual" in reason for reason in validation.reasons))

    def test_boroxine_realization_has_product_formula_and_water_loss(self):
        candidate = RingFormingStructureGenerator().generate(
            self.cof1_precursor,
            "boroxine_trimerization",
        )
        realization = ReactionRealizer().realize(
            candidate,
            {self.cof1_precursor.id: self.cof1_precursor},
            candidate.metadata["instance_to_monomer"],
        )
        self.assertIsNotNone(realization)
        counts = Counter(
            atom.symbol
            for atoms in realization.atoms_by_instance.values()
            for atom in atoms
        )

        self.assertEqual(counts, Counter({"C": 18, "H": 12, "B": 6, "O": 6}))
        self.assertEqual(realization.metadata["removed_atom_symbols"], {"O": 6, "H": 12})
        self.assertEqual(realization.metadata["applied_event_count"], 2)
        self.assertEqual(len(realization.bonds), 6)
        self.assertTrue(all(abs(bond.distance - 1.38) < 1e-3 for bond in realization.bonds))

        export = CIFWriter().export_candidate(
            candidate,
            {self.cof1_precursor.id: self.cof1_precursor},
        )
        block = gemmi.cif.read_string(export.text).sole_block()
        structure = gemmi.make_small_structure_from_block(block)
        self.assertEqual(Counter(site.element.name for site in structure.sites), counts)
        self.assertEqual(export.n_sites, 42)

    def test_triazine_realization_rewrites_nitrile_bond_orders(self):
        candidate = RingFormingStructureGenerator().generate(
            self.ctf1_precursor,
            "triazine_trimerization",
        )
        export = CIFWriter().export_candidate(
            candidate,
            {self.ctf1_precursor.id: self.ctf1_precursor},
        )
        block = gemmi.cif.read_string(export.text).sole_block()
        structure = gemmi.make_small_structure_from_block(block)
        bond_types = Counter(str(value) for value in block.find_values("_ccdc_geom_bond_type"))

        self.assertEqual(
            Counter(site.element.name for site in structure.sites),
            Counter({"C": 24, "H": 12, "N": 6}),
        )
        self.assertEqual(bond_types["D"], 6)
        self.assertEqual(bond_types["T"], 0)
        realization = export.metadata["reaction_realization"]
        self.assertEqual(realization["removed_atom_count"], 0)
        self.assertEqual(realization["applied_event_count"], 2)

    def test_ring_cofid_round_trip(self):
        candidate = RingFormingStructureGenerator().generate(
            self.ctf1_precursor,
            "triazine_trimerization",
        )
        cofid = generate_cofid(candidate, {self.ctf1_precursor.id: self.ctf1_precursor})
        request = cofid_to_build_request(cofid)

        self.assertEqual(cofid, "2:nitrile:N#Cc1ccc(C#N)cc1&&hcb&&triazine")
        self.assertEqual(request.template_id, "triazine_trimerization")
        self.assertEqual(len(request.monomers), 1)
        self.assertEqual(request.monomers[0].motif_kind, "nitrile")

    def test_boroxine_aa_stacking_duplicates_layer_graph_and_atomistic_product(self):
        base = RingFormingStructureGenerator().generate(
            self.cof1_precursor,
            "boroxine_trimerization",
        )
        stacked = RingFormingStructureGenerator(
            RingFormationConfig(stacking_ids=("AA",))
        ).generate(
            self.cof1_precursor,
            "boroxine_trimerization",
        )
        validation = validate_ring_geometry(
            stacked.events,
            stacked.state,
            {self.cof1_precursor.id: self.cof1_precursor},
        )
        cofid = generate_cofid(stacked, {self.cof1_precursor.id: self.cof1_precursor})
        export = CIFWriter().export_candidate(
            stacked,
            {self.cof1_precursor.id: self.cof1_precursor},
            cofid=cofid,
            cofid_comment_suffix=stacked.metadata["stacking"]["comment_suffix"],
        )
        block = gemmi.cif.read_string(export.text).sole_block()
        structure = gemmi.make_small_structure_from_block(block)

        self.assertEqual(stacked.id, "ring-candidate-1__AA")
        self.assertEqual(stacked.state.stacking_state, "AA")
        self.assertEqual(stacked.metadata["stacking"]["registry_shift_fractional"], (0.0, 0.0))
        self.assertEqual(len(stacked.state.monomer_poses), 2 * len(base.state.monomer_poses))
        self.assertEqual(len(stacked.events), 2 * len(base.events))
        self.assertEqual(stacked.metadata["graph_summary"]["n_reaction_events"], 4)
        self.assertEqual(len(stacked.metadata["ring_validation"]["metrics"]["events"]), 4)
        self.assertEqual(validation.classification, "accepted")
        self.assertEqual(
            Counter(site.element.name for site in structure.sites),
            Counter({"C": 36, "H": 24, "B": 12, "O": 12}),
        )
        self.assertEqual(export.metadata["reaction_realization"]["applied_event_count"], 4)
        self.assertAlmostEqual(
            structure.cell.c,
            2.0
            * (
                stacked.metadata["stacking"]["interlayer_clearance"]
                + stacked.metadata["stacking"]["layer_z_span"]
            ),
            places=5,
        )
        self.assertEqual(export.text.splitlines()[0], f"# COFid: {cofid} stacking=AA")

    def test_boroxine_aa_stacked_export_renders_stacking_geometry_comment(self):
        """W3.2 wiring: real stacked ring-forming exports carry the
        stacking-geometry derivation comment (registry/shift/clearance/span/
        c2c/atomistic contact) next to the periodic_bilayer c-axis-semantics
        label; the monolayer export carries neither the comment nor the
        bilayer label."""
        monolayer = RingFormingStructureGenerator().generate(
            self.cof1_precursor,
            "boroxine_trimerization",
        )
        monolayer_export = CIFWriter().export_candidate(
            monolayer,
            {self.cof1_precursor.id: self.cof1_precursor},
        )
        stacked = RingFormingStructureGenerator(
            RingFormationConfig(stacking_ids=("AA",))
        ).generate(
            self.cof1_precursor,
            "boroxine_trimerization",
        )
        stacked_export = CIFWriter().export_candidate(
            stacked,
            {self.cof1_precursor.id: self.cof1_precursor},
        )

        comment_lines = [
            line for line in stacked_export.text.splitlines() if line.startswith("# stacking-geometry:")
        ]
        self.assertEqual(len(comment_lines), 1)
        comment = comment_lines[0]
        stacking_meta = stacked.metadata["stacking"]
        self.assertIn("registry=AA", comment)
        self.assertIn("shift_frac=0,0", comment)
        self.assertIn(f"interlayer_clearance={DEFAULT_INTERLAYER_CLEARANCE_ANGSTROM:.3g}", comment)
        self.assertIn(f"layer_z_span={stacking_meta['layer_z_span']:.4g}", comment)
        self.assertIn(f"c2c={stacking_meta['center_to_center_distance']:.4g}", comment)
        self.assertIn("c=2*c2c", comment)
        # The boroxine monolayer is near-flat, so the interlayer contact equals
        # the clearance (3.5 Å) and sits exactly on the 3.5 Å search-cutoff
        # boundary; float noise decides whether it renders as the measured
        # contact or the honest ">cutoff" lower bound. Accept either.
        self.assertIn("atomistic contact, mode=atomistic_product", comment)
        if stacking_meta.get("min_interlayer_contact") is None:
            self.assertIn(
                f"min_interlayer_contact>{stacking_meta['min_interlayer_contact_cutoff']:.3g}",
                comment,
            )
            self.assertIn("none below cutoff", comment)
        else:
            self.assertIn(f"min_interlayer_contact={stacking_meta['min_interlayer_contact']:.3g}", comment)
            self.assertIn("atoms=", comment)
            self.assertIn("involves_hydrogen=", comment)
        self.assertIn("# c-axis-semantics: periodic_bilayer", stacked_export.text)

        self.assertNotIn("# stacking-geometry:", monolayer_export.text)
        self.assertNotIn("stacking", monolayer_export.metadata)
        self.assertIn("# c-axis-semantics: vacuum_slab", monolayer_export.text)

    def test_stacking_span_provenance_and_fitted_cell_classification(self):
        """Audit A7/A8: the span is measured along the layer normal with
        honest mode/axis provenance, and the embedding cell_kind classifies
        the fitted cell (declared RCSR family kept as provenance only)."""
        stacked = RingFormingStructureGenerator(
            RingFormationConfig(stacking_ids=("AA",))
        ).generate(
            self.cof1_precursor,
            "boroxine_trimerization",
        )
        stacking_meta = stacked.metadata["stacking"]
        embedding = stacked.metadata["embedding"]

        self.assertEqual(embedding["layer_z_span_mode"], "atomistic_product")
        self.assertEqual(embedding["layer_z_span_axis"], "layer_normal")
        self.assertGreaterEqual(embedding["layer_z_span"], 0.0)
        self.assertEqual(embedding["cell_kind"], "hexagonal")
        self.assertEqual(embedding["declared_metric_family"], "hexagonal")

        # Stacking metadata honestly propagates the measurement provenance.
        self.assertEqual(stacking_meta["layer_z_span_mode"], "atomistic_product")
        self.assertEqual(stacking_meta["layer_z_span_axis"], "layer_normal")
        self.assertAlmostEqual(stacking_meta["layer_z_span"], embedding["layer_z_span"], places=6)

    def test_triazine_ab_stacking_uses_hexagonal_registry_and_passes_coarse_validation(self):
        stacked = RingFormingStructureGenerator(
            RingFormationConfig(stacking_ids=("AB",))
        ).generate(
            self.ctf1_precursor,
            "triazine_trimerization",
        )
        stacking_meta = stacked.metadata["stacking"]
        cell_setting = stacking_meta["cell_classification"].get("setting")
        # The hexagonal AB shift is cell-setting-aware (W2.2): 60° cells use
        # (1/3, 1/3); 120° cells (built by the ring-forming path) use (1/3, 2/3).
        self.assertEqual(cell_setting, "120deg")
        self.assertEqual(
            stacking_meta["registry_shift_fractional"],
            (1.0 / 3.0, 2.0 / 3.0),
        )
        # W7 test 3 (120-degree half): the shift maps a vertex onto the pore
        # center, so its cartesian magnitude is a/sqrt(3).
        cell = stacked.state.cell
        shift = stacking_meta["registry_shift_fractional"]
        shift_cartesian = [
            shift[0] * cell[0][axis] + shift[1] * cell[1][axis]
            for axis in range(3)
        ]
        a_length = sqrt(sum(component**2 for component in cell[0]))
        shift_length = sqrt(sum(component**2 for component in shift_cartesian))
        # The fitted cell is only approximately hexagonal, so compare with a
        # 1% relative tolerance rather than exactly.
        self.assertAlmostEqual(
            shift_length,
            a_length / sqrt(3.0),
            delta=0.01 * a_length / sqrt(3.0),
        )
        with tempfile.TemporaryDirectory() as temporary_dir:
            cif_path = Path(temporary_dir) / "ctf1_ab.cif"
            cofid = generate_cofid(stacked, {self.ctf1_precursor.id: self.ctf1_precursor})
            export = CIFWriter().write_candidate(
                cif_path,
                stacked,
                {self.ctf1_precursor.id: self.ctf1_precursor},
                cofid=cofid,
                cofid_comment_suffix=stacked.metadata["stacking"]["comment_suffix"],
            )
            report = CoarseStructureValidator().validate_manifest_record(
                {
                    "topology_id": "hcb",
                    "cif_path": str(cif_path),
                    "flags": list(stacked.flags),
                    "metadata": dict(stacked.metadata),
                }
            )

        self.assertEqual(export.n_sites, 84)
        self.assertEqual(report.classification, "valid")
        self.assertEqual(report.metrics["n_instance_components"], 2)
        self.assertTrue(report.metrics["disconnected_instance_graph_allowed"])

    def test_monolayer_export_uses_unified_vacuum_slab_c_axis(self):
        """W1.4: ring-forming monolayers export c = 8.0 (vacuum-slab padding,
        unified with the embedding path) and label the semantics explicitly."""
        candidate = RingFormingStructureGenerator().generate(
            self.cof1_precursor,
            "boroxine_trimerization",
        )
        export = CIFWriter().export_candidate(
            candidate,
            {self.cof1_precursor.id: self.cof1_precursor},
        )
        block = gemmi.cif.read_string(export.text).sole_block()
        structure = gemmi.make_small_structure_from_block(block)

        self.assertAlmostEqual(candidate.state.cell[2][2], 8.0)
        self.assertAlmostEqual(structure.cell.c, 8.0)
        self.assertEqual(candidate.metadata["embedding"]["c_axis_semantics"], "vacuum_slab")
        self.assertEqual(export.metadata["c_axis_semantics"], "vacuum_slab")
        self.assertIn("# c-axis-semantics: vacuum_slab", export.text)

    def test_user_set_layer_spacing_is_still_vacuum_slab_padding(self):
        """A user-overridden layer_spacing changes the padding thickness, not
        the semantics: the c axis is still vacuum-slab padding, not a physical
        stacking repeat."""
        candidate = RingFormingStructureGenerator(
            RingFormationConfig(layer_spacing=12.0)
        ).generate(
            self.cof1_precursor,
            "boroxine_trimerization",
        )
        export = CIFWriter().export_candidate(
            candidate,
            {self.cof1_precursor.id: self.cof1_precursor},
        )

        self.assertAlmostEqual(candidate.state.cell[2][2], 12.0)
        self.assertEqual(export.metadata["c_axis_semantics"], "vacuum_slab")
        self.assertIn("# c-axis-semantics: vacuum_slab", export.text)

    def test_stacked_export_labels_periodic_bilayer_not_vacuum_slab(self):
        stacked = RingFormingStructureGenerator(
            RingFormationConfig(stacking_ids=("AA",))
        ).generate(
            self.cof1_precursor,
            "boroxine_trimerization",
        )
        export = CIFWriter().export_candidate(
            stacked,
            {self.cof1_precursor.id: self.cof1_precursor},
        )

        self.assertEqual(stacked.metadata["embedding"]["c_axis_semantics"], "periodic_bilayer")
        self.assertEqual(export.metadata["c_axis_semantics"], "periodic_bilayer")
        self.assertIn("# c-axis-semantics: periodic_bilayer", export.text)
        self.assertNotIn("vacuum_slab", export.text)

    def test_engine_enumerates_requested_ring_stackings(self):
        ensemble = COFEngine().run(
            COFProject(
                monomers=(self.ctf1_precursor,),
                allowed_reactions=("triazine_trimerization",),
                target_dimensionality="2D",
                target_topologies=("hcb",),
                stacking_ids=("AA", "AB"),
            )
        )

        self.assertEqual(len(ensemble.candidates), 2)
        self.assertEqual(
            {candidate.state.stacking_state for candidate in ensemble.candidates},
            {"AA", "AB"},
        )

    def test_indexed_kgd_supports_hexatopic_precursor_and_periodic_copy_identity(self):
        count = 6
        positions = tuple(
            (4.0 * cos(2.0 * pi * index / count), 4.0 * sin(2.0 * pi * index / count), 0.0)
            for index in range(count)
        )
        precursor = MonomerSpec(
            id="hexatopic",
            name="idealized hexatopic precursor",
            motifs=tuple(
                ReactiveMotif(
                    id=f"b{index}",
                    kind="boronic_acid",
                    atom_ids=(index,),
                    frame=Frame(
                        origin=positions[index],
                        primary=positions[index],
                        normal=(0.0, 0.0, 1.0),
                    ),
                )
                for index in range(count)
            ),
            atom_symbols=("B",) * count,
            atom_positions=positions,
        )
        result = RingFormingStructureGenerator().build(
            precursor,
            "boroxine_trimerization",
            topology_id="kgd",
        )

        self.assertEqual(len(result.outcome.monomer_instances), 1)
        self.assertEqual(len(result.outcome.events), 2)
        self.assertEqual(result.candidate.metadata["ring_validation"]["classification"], "accepted")
        for event in result.outcome.events:
            physical_copies = {(ref.monomer_instance_id, ref.periodic_image) for ref in event.participants}
            self.assertEqual(len(physical_copies), 3)

    def test_tritopic_triazine_build_with_shape_selected_conformer_is_accepted(self):
        # Regression for the twisted-triazine-ring bug report: this precursor's
        # lowest-energy conformer has nitrile arms at 108/130/122-degree gaps,
        # which the rigid node placement cannot repair, so the exported ring
        # was an irregular hexagon (bonds 0.94-1.67 A). Shape-aware conformer
        # selection picks a topology-compatible conformer instead.
        precursor = build_rdkit_monomer(
            "tris_nitrile",
            "tris_nitrile",
            "N#Cc1ccc(NC(=O)c2cc(C(=O)Nc3ccc(C#N)cc3)cc(C(=O)Nc3ccc(C#N)cc3)c2)cc1",
            "nitrile",
            num_conformers=SHAPE_SELECTION_MIN_CONFORMERS,
            select_conformer_by_motif_shape=True,
        )
        self.assertEqual(precursor.metadata["conformer_selection"], "motif_shape")

        result = RingFormingStructureGenerator().build(precursor, "triazine_trimerization")

        self.assertEqual(result.candidate.metadata["ring_validation"]["classification"], "accepted")
        self.assertNotIn("ring_geometry_rejected", result.candidate.flags)
        validation = validate_ring_geometry(
            result.outcome.events,
            result.candidate.state,
            {precursor.id: precursor},
        )
        self.assertEqual(validation.classification, "accepted")
        self.assertEqual(validation.reasons, ())

    def test_ring_forming_cli_writes_summary_and_cif(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                cli_main(
                    [
                        "build",
                        "ring-forming",
                        "--template-id",
                        "triazine_trimerization",
                        "--smiles",
                        "N#Cc1ccc(C#N)cc1",
                        "--num-conformers",
                        "1",
                        "--output-dir",
                        temporary_dir,
                    ]
                )
            report = json.loads((Path(temporary_dir) / "summary.json").read_text())

            self.assertEqual(report["attempted_structures"], 1)
            self.assertEqual(report["successful_structures"], 1)
            self.assertEqual(report["cifs_written"], 1)
            self.assertEqual(report["precursor"]["geometry"]["embedding_method"], "etkdg-v3")
            self.assertFalse(report["precursor"]["geometry"]["fallback"])
            self.assertEqual(report["results"][0]["status"], "ok")
            self.assertEqual(report["results"][0]["reaction_realization_status"], "completed")
            self.assertIsNotNone(report["results"][0]["reaction_realization"])
            self.assertEqual(report["ring_validation"]["classification"], "accepted")
            self.assertEqual(report["graph_summary"]["n_reaction_events"], 2)
            self.assertEqual(report["cif_sites"], 42)
            self.assertTrue(Path(report["cif_path"]).is_file())

    def test_ring_forming_cli_no_cif_still_reports_success(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            with contextlib.redirect_stdout(io.StringIO()):
                cli_main(
                    [
                        "build",
                        "ring-forming",
                        "--cofid",
                        "2:nitrile:N#Cc1ccc(C#N)cc1&&hcb&&triazine",
                        "--num-conformers",
                        "1",
                        "--no-write-cif",
                        "--output-dir",
                        temporary_dir,
                    ]
                )
            report = json.loads((Path(temporary_dir) / "summary.json").read_text())

        self.assertEqual(report["attempted_structures"], 1)
        self.assertEqual(report["successful_structures"], 1)
        self.assertEqual(report["cifs_written"], 0)
        self.assertEqual(report["results"][0]["status"], "ok")
        self.assertEqual(report["results"][0]["reaction_realization_status"], "not_requested")
        self.assertIsNone(report["results"][0]["reaction_realization"])

    def test_long_ditopic_precursor_tolerates_topology_coordinate_rounding(self):
        precursor = build_rdkit_monomer(
            "ta_por",
            "TA-Por-sp2-COF precursor",
            "N#Cc1ccc(-c2c3nc(cc4ccc([nH]4)c(-c4ccc(C#N)cc4)c4nc(cc5ccc2[nH]5)C=C4)C=C3)cc1",
            "nitrile",
            num_conformers=1,
        )

        candidate = RingFormingStructureGenerator().generate(
            precursor,
            "triazine_trimerization",
        )

        placement = candidate.metadata["edge_placement"]
        self.assertLessEqual(placement["max_residual"], placement["tolerance"])
        self.assertEqual(candidate.metadata["ring_validation"]["classification"], "accepted")

    def test_ring_forming_cli_enumerates_multiple_stacking_registries(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            with contextlib.redirect_stdout(io.StringIO()):
                cli_main(
                    [
                        "build",
                        "ring-forming",
                        "--template-id",
                        "triazine_trimerization",
                        "--smiles",
                        "N#Cc1ccc(C#N)cc1",
                        "--stacking",
                        "AA",
                        "--stacking",
                        "AB",
                        "--num-conformers",
                        "1",
                        "--output-dir",
                        temporary_dir,
                    ]
                )
            report = json.loads((Path(temporary_dir) / "summary.json").read_text())

            self.assertEqual(report["stacking_requested"], ["AA", "AB"])
            self.assertEqual(len(report["results"]), 2)
            self.assertEqual({row["stacking"]["id"] for row in report["results"]}, {"AA", "AB"})
            for row in report["results"]:
                self.assertEqual(row["ring_validation"]["classification"], "accepted")
                self.assertEqual(row["cif_sites"], 84)
                cif_lines = Path(row["cif_path"]).read_text().splitlines()
                self.assertEqual(
                    cif_lines[0],
                    f"# COFid: {row['generated_cofid']} stacking={row['stacking']['id']}",
                )


class RingAttachmentGeometryTests(unittest.TestCase):
    """A05: exocyclic ring-monomer attachment geometry is measured and
    classified separately from the ring-participant arrangement verdict."""

    FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures" / "ring_attachment"

    @classmethod
    def setUpClass(cls):
        cls.para_precursor = build_rdkit_monomer(
            "para",
            "benzene-1,4-diboronic acid",
            "OB(O)c1ccc(B(O)O)cc1",
            "boronic_acid",
            num_conformers=1,
        )
        cls.ortho_precursor = build_rdkit_monomer(
            "ortho",
            "benzene-1,2-diboronic acid",
            "OB(O)c1ccccc1B(O)O",
            "boronic_acid",
            num_conformers=1,
        )

    def _attachment_angles_from_cif(self, cif_path: Path) -> list[float]:
        block = gemmi.cif.read_file(str(cif_path)).sole_block()
        structure = gemmi.make_small_structure_from_block(block)
        by_label = {str(site.label): site for site in structure.sites}
        neighbors: dict[str, set[str]] = {}
        labels1 = block.find_loop("_geom_bond_atom_site_label_1")
        labels2 = block.find_loop("_geom_bond_atom_site_label_2")
        for index in range(min(len(labels1), len(labels2))):
            first, second = str(labels1[index]), str(labels2[index])
            neighbors.setdefault(first, set()).add(second)
            neighbors.setdefault(second, set()).add(first)
        angles: list[float] = []
        for site in structure.sites:
            if site.element.name != "B":
                continue
            boron = structure.cell.orthogonalize(site.fract)
            for neighbor_label in neighbors.get(str(site.label), ()):
                neighbor = by_label[neighbor_label]
                if neighbor.element.name != "C":
                    continue
                for other_label in neighbors.get(str(site.label), ()):
                    other = by_label[other_label]
                    if other.element.name != "O":
                        continue
                    vectors = [
                        structure.cell.find_nearest_pbc_position(
                            boron, structure.cell.orthogonalize(partner.fract), 0
                        )
                        - boron
                        for partner in (neighbor, other)
                    ]
                    first, second = vectors
                    cosine = max(-1.0, min(1.0, first.dot(second) / first.length() / second.length()))
                    angles.append(degrees(acos(cosine)))
        return sorted(angles)

    def test_saved_para_ortho_fixtures_distinguish_attachment_angles(self):
        """The saved fixtures reproduce the review evidence: ordinary aromatic
        sp2 attachment at ~120 deg (para) versus the pathological ~64/175 deg
        case (ortho)."""
        profile = ring_geometry_profile("boroxine_trimerization")
        para_angles = self._attachment_angles_from_cif(self.FIXTURE_DIR / "para.cif")
        ortho_angles = self._attachment_angles_from_cif(self.FIXTURE_DIR / "ortho.cif")

        self.assertEqual(len(para_angles), 12)
        self.assertEqual(len(ortho_angles), 12)
        for angle in para_angles:
            self.assertLessEqual(
                abs(angle - RING_ATTACHMENT_IDEAL_ANGLE_DEGREES),
                profile.attachment_warning_deviation_degrees,
            )
        # Ortho pathology: ~64.3/~175.3 deg, i.e. ~55 deg off ideal — well
        # beyond the rejection tolerance.
        self.assertAlmostEqual(ortho_angles[0], 64.31, delta=0.5)
        self.assertAlmostEqual(ortho_angles[-1], 175.31, delta=0.5)
        self.assertGreater(
            min(abs(angle - RING_ATTACHMENT_IDEAL_ANGLE_DEGREES) for angle in ortho_angles[:6]),
            profile.attachment_rejection_deviation_degrees,
        )
        self.assertGreater(
            min(abs(angle - RING_ATTACHMENT_IDEAL_ANGLE_DEGREES) for angle in ortho_angles[6:]),
            profile.attachment_rejection_deviation_degrees,
        )

    def test_para_attachment_accepted_and_ortho_attachment_rejected(self):
        para = RingFormingStructureGenerator().build(self.para_precursor, "boroxine_trimerization")
        ortho = RingFormingStructureGenerator().build(self.ortho_precursor, "boroxine_trimerization")

        para_validation = para.candidate.metadata["ring_validation"]
        self.assertEqual(para_validation["classification"], "accepted")
        self.assertEqual(para_validation["arrangement_classification"], "accepted")
        self.assertEqual(para_validation["attachment_classification"], "accepted")
        self.assertEqual(para_validation["attachment_status"], "measured")
        self.assertNotIn("ring_attachment_rejected", para.candidate.flags)
        self.assertNotIn("ring_attachment_warning", para.candidate.flags)
        for measurement in para_validation["attachment"]["measurements"]:
            for angle in measurement["ring_angles_degrees"]:
                self.assertAlmostEqual(angle, RING_ATTACHMENT_IDEAL_ANGLE_DEGREES, delta=1.0)

        ortho_validation = ortho.candidate.metadata["ring_validation"]
        # The arrangement channel still passes (the review evidence: both
        # precursors pass the ring test); the attachment channel rejects.
        self.assertEqual(ortho_validation["arrangement_classification"], "accepted")
        self.assertEqual(ortho_validation["attachment_classification"], "rejected")
        self.assertEqual(ortho_validation["classification"], "rejected")
        self.assertEqual(ortho_validation["attachment_status"], "measured")
        self.assertIn("ring_attachment_rejected", ortho.candidate.flags)
        self.assertNotIn("ring_geometry_rejected", ortho.candidate.flags)
        ortho_angles = sorted(
            angle
            for measurement in ortho_validation["attachment"]["measurements"]
            for angle in measurement["ring_angles_degrees"]
        )
        self.assertAlmostEqual(ortho_angles[0], 64.31, delta=0.5)
        self.assertAlmostEqual(ortho_angles[-1], 175.31, delta=0.5)
        self.assertTrue(ortho_validation["reasons"])

    def test_coarse_validation_reports_attachment_channel_separately(self):
        """validate_manifest_record: the ortho pathology lands in the
        attachment channel (ring_attachment_invalid), not the arrangement
        channel (ring_geometry_invalid), and the coverage map records both
        channels as measured."""
        for name, precursor, expected_classification in (
            ("para", self.para_precursor, "valid"),
            ("ortho", self.ortho_precursor, "hard_invalid"),
        ):
            candidate = RingFormingStructureGenerator().generate(precursor, "boroxine_trimerization")
            with tempfile.TemporaryDirectory() as temporary_dir:
                cif_path = Path(temporary_dir) / f"{name}.cif"
                CIFWriter().write_candidate(cif_path, candidate, {precursor.id: precursor})
                report = CoarseStructureValidator().validate_manifest_record(
                    {
                        "topology_id": "hcb",
                        "cif_path": str(cif_path),
                        "flags": list(candidate.flags),
                        "metadata": dict(candidate.metadata),
                    }
                )
            self.assertEqual(report.classification, expected_classification)
            self.assertEqual(report.coverage["ring_geometry"], "measured")
            self.assertEqual(report.coverage["ring_attachment"], "measured")
            if name == "ortho":
                self.assertIn("ring_attachment_invalid", report.hard_invalid_reasons)
                self.assertNotIn("ring_geometry_invalid", report.hard_invalid_reasons)
                self.assertEqual(report.metrics["ring_attachment_classification"], "rejected")

    def test_legacy_ring_record_without_attachment_channel_reports_missing_data(self):
        """Records whose ring_validation metadata predates the attachment
        channel keep their arrangement verdict but honestly report the
        attachment channel as unmeasured."""
        candidate = RingFormingStructureGenerator().generate(self.para_precursor, "boroxine_trimerization")
        legacy_ring_validation = {
            "classification": candidate.metadata["ring_validation"]["classification"],
            "reasons": [],
            "metrics": dict(candidate.metadata["ring_validation"]["metrics"]),
        }
        metadata = dict(candidate.metadata)
        metadata["ring_validation"] = legacy_ring_validation
        with tempfile.TemporaryDirectory() as temporary_dir:
            cif_path = Path(temporary_dir) / "para.cif"
            CIFWriter().write_candidate(cif_path, candidate, {self.para_precursor.id: self.para_precursor})
            report = CoarseStructureValidator().validate_manifest_record(
                {
                    "topology_id": "hcb",
                    "cif_path": str(cif_path),
                    "flags": list(candidate.flags),
                    "metadata": metadata,
                }
            )
        self.assertEqual(report.coverage["ring_attachment"], "missing_data")
        self.assertEqual(report.coverage["ring_geometry"], "measured")
        # ring_attachment is not a required check: the record can still be
        # valid, but the missing channel is explicit.
        self.assertEqual(report.classification, "valid")
        self.assertNotIn("ring_attachment", report.unmeasured_required_checks)

    def test_arrangement_and_attachment_verdicts_are_independent(self):
        """An off-radius participant fails the arrangement channel while the
        attachment channel (radially outward exocyclic bond) stays accepted.

        The 0.20 A offset exceeds the arrangement radial tolerance (0.18 A)
        but stretches the exocyclic bond by only 0.20/1.55 = 12.9%, below the
        attachment warning fraction (15%)."""
        events, state, spec = _synthetic_ring_assembly(radial_offset=0.20)
        validation = validate_ring_geometry(events, state, {spec.id: spec})

        self.assertEqual(validation.arrangement_classification, "rejected")
        self.assertEqual(validation.attachment_classification, "accepted")
        self.assertEqual(validation.classification, "rejected")
        self.assertTrue(any("radial residual" in reason for reason in validation.reasons))

    def test_attachment_warning_tier_for_moderately_strained_geometry(self):
        events, state, spec = _synthetic_ring_assembly(deviation_degrees=25.0)
        profile = ring_geometry_profile("boroxine_trimerization")
        self.assertGreater(25.0, profile.attachment_warning_deviation_degrees)
        self.assertLess(25.0, profile.attachment_rejection_deviation_degrees)

        validation = validate_ring_geometry(events, state, {spec.id: spec})

        self.assertEqual(validation.arrangement_classification, "accepted")
        self.assertEqual(validation.attachment_classification, "warning")
        self.assertEqual(validation.classification, "warning")
        self.assertTrue(any("warning tolerance" in reason for reason in validation.reasons))

    def test_attachment_rejection_tier_for_severely_strained_geometry(self):
        events, state, spec = _synthetic_ring_assembly(deviation_degrees=55.0)

        validation = validate_ring_geometry(events, state, {spec.id: spec})

        self.assertEqual(validation.arrangement_classification, "accepted")
        self.assertEqual(validation.attachment_classification, "rejected")
        self.assertEqual(validation.classification, "rejected")

    def test_attachment_missing_atom_metadata_is_not_a_pass(self):
        events, state, spec = _synthetic_ring_assembly(with_atom_metadata=False)
        validation = validate_ring_geometry(events, state, {spec.id: spec})
        report = ring_attachment_report(events, state, {spec.id: spec})

        self.assertEqual(report.status, "missing_data")
        self.assertEqual(report.n_unmeasured_participants, 3)
        self.assertEqual(report.classification, "not_applicable")
        self.assertEqual(validation.attachment.status, "missing_data")
        # Arrangement acceptance is unaffected; the missing attachment channel
        # is reported through the status, not silently counted as passing.
        self.assertEqual(validation.arrangement_classification, "accepted")
        self.assertEqual(validation.classification, "accepted")


def _synthetic_ring_assembly(
    deviation_degrees: float = 0.0,
    radial_offset: float = 0.0,
    with_atom_metadata: bool = True,
) -> tuple[tuple[ReactionEvent, ...], AssemblyState, MonomerSpec]:
    """Three two-atom (B/C) monomer instances on a regular boroxine ring.

    The anchor (C) sits 1.55 A from the reactive atom (B), pointing radially
    outward plus `deviation_degrees` of in-plane rotation, so both exocyclic
    attachment angles deviate from the 120-degree ideal by that amount.
    `radial_offset` moves every participant off the ring radius to exercise
    the arrangement channel independently.
    """
    profile = ring_geometry_profile("boroxine_trimerization")
    radius = profile.ring_atom_bond_length
    metadata = {"reactive_atom_id": 0, "anchor_atom_id": 1} if with_atom_metadata else {}
    spec = MonomerSpec(
        id="synthetic",
        name="synthetic boronic acid",
        motifs=(
            ReactiveMotif(
                id="bor1",
                kind="boronic_acid",
                atom_ids=(0, 1),
                frame=Frame(origin=(0.0, 0.0, 0.0), primary=(1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)),
                metadata=metadata,
            ),
        ),
        atom_symbols=("B", "C"),
        atom_positions=((0.0, 0.0, 0.0), (1.55, 0.0, 0.0)),
    )
    deviation = deviation_degrees * pi / 180.0
    poses: dict[str, Pose] = {}
    refs: list[MotifRef] = []
    for index in range(3):
        ring_angle = 2.0 * pi * index / 3.0
        vertex = (
            (radius + radial_offset) * cos(ring_angle),
            (radius + radial_offset) * sin(ring_angle),
            0.0,
        )
        anchor_angle = ring_angle + deviation
        rotation = (
            (cos(anchor_angle), -sin(anchor_angle), 0.0),
            (sin(anchor_angle), cos(anchor_angle), 0.0),
            (0.0, 0.0, 1.0),
        )
        instance_id = f"p{index + 1}"
        poses[instance_id] = Pose(translation=vertex, rotation_matrix=rotation)
        refs.append(MotifRef(instance_id, spec.id, "bor1"))
    event = ReactionEvent(
        id="r1",
        template_id="boroxine_trimerization",
        participants=tuple(refs),
        product_state="boroxine",
        metadata={
            "ring_center_fractional": (0.0, 0.0, 0.0),
            "ring_normal": (0.0, 0.0, 1.0),
            "ring_atom_bond_length": radius,
        },
    )
    state = AssemblyState(
        cell=((30.0, 0.0, 0.0), (0.0, 30.0, 0.0), (0.0, 0.0, 8.0)),
        monomer_poses=poses,
        stacking_state="disabled",
    )
    return (event,), state, spec


if __name__ == "__main__":
    unittest.main()
