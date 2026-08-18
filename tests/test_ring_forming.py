import contextlib
import io
import json
import sys
import tempfile
import unittest
from collections import Counter
from dataclasses import replace
from math import cos, pi, sin
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

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
        self.assertEqual(stacked.metadata["stacking"]["lateral_shift_fractional"], (0.0, 0.0))
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
                stacked.metadata["stacking"]["interlayer_distance"]
                + stacked.metadata["stacking"]["layer_z_span"]
            ),
            places=5,
        )
        self.assertEqual(export.text.splitlines()[0], f"# COFid: {cofid} stacking=AA")

    def test_triazine_ab_stacking_uses_hexagonal_registry_and_passes_coarse_validation(self):
        stacked = RingFormingStructureGenerator(
            RingFormationConfig(stacking_ids=("AB",))
        ).generate(
            self.ctf1_precursor,
            "triazine_trimerization",
        )
        self.assertEqual(
            stacked.metadata["stacking"]["lateral_shift_fractional"],
            (1.0 / 3.0, 1.0 / 3.0),
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

            self.assertEqual(report["ring_validation"]["classification"], "accepted")
            self.assertEqual(report["graph_summary"]["n_reaction_events"], 2)
            self.assertEqual(report["cif_sites"], 42)
            self.assertTrue(Path(report["cif_path"]).is_file())

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
