import contextlib
import io
import json
import sys
import tempfile
import unittest
from collections import Counter
from dataclasses import replace
from math import cos, pi, sin, sqrt
from pathlib import Path


import gemmi

from cofkit import COFEngine, COFProject, CoarseStructureValidator
from cofkit.build_workflows.ring_forming import RingFormationConfig, RingFormingStructureGenerator
from cofkit.chem.rdkit import build_rdkit_monomer
from cofkit.cif import CIFWriter
from cofkit.cli import main as cli_main
from cofkit.cofid import cofid_to_build_request, generate_cofid
from cofkit.geometry import Frame
from cofkit.model import MonomerSpec, ReactiveMotif
from cofkit.reaction_realization import ReactionRealizer
from cofkit.ring_geometry import validate_ring_geometry
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


if __name__ == "__main__":
    unittest.main()
