"""Regression tests for A18 (impact-review claims T1-7 / T3.5 / T3.11):
effective conformer-construction settings — requested vs actual budgets,
seeds, selection mode, fallback status, and implementation identity — are
recorded in durable outputs, and CLI defaults reference their semantic
owners.
"""

import argparse
import contextlib
import inspect
import io
import json
import tempfile
import unittest
from pathlib import Path

from cofkit import __version__
from cofkit import cli_build
from cofkit.batch import BatchGenerationConfig, BatchStructureGenerator
from cofkit.batch_models import BatchMonomerRecord
from cofkit.chem.rdkit import (
    DEFAULT_RDKIT_RANDOM_SEED,
    SHAPE_SELECTION_MIN_CONFORMERS,
    build_rdkit_monomer,
)
from cofkit.cli import main as cli_main
from cofkit.monomer_library import AUTODETECT_MAX_CONFORMERS

try:
    from rdkit import Chem  # noqa: F401
except ImportError:  # pragma: no cover - environment-dependent
    Chem = None


TAPB = "C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N"
TFB = "C1=C(C=C(C=C1C=O)C=O)C=O"
PPD = "Nc1ccc(N)cc1"


@unittest.skipIf(Chem is None, "RDKit is not available")
class BuilderMetadataTests(unittest.TestCase):
    """The monomer builder records the effective settings it ran with."""

    def test_builder_records_budget_seed_selection_and_identity(self):
        monomer = build_rdkit_monomer("ppd", "ppd", PPD, "amine", num_conformers=3)
        metadata = monomer.metadata
        self.assertEqual(metadata["requested_num_conformers"], 3)
        self.assertGreaterEqual(metadata["n_conformers"], 1)
        self.assertEqual(metadata["random_seed"], DEFAULT_RDKIT_RANDOM_SEED)
        self.assertEqual(metadata["shape_aware_selection_requested"], False)
        self.assertEqual(metadata["conformer_selection"], "energy")
        self.assertNotIn("conformer_selection_note", metadata)
        self.assertEqual(metadata["cofkit_version"], __version__)
        self.assertTrue(metadata["rdkit_version"])
        self.assertEqual(metadata["embedding_fallback"], False)
        self.assertEqual(metadata["embedding_method"], "etkdg-v3")

    def test_shape_request_below_motif_gate_records_note(self):
        # Ditopic monomer: shape-aware selection was requested but the motif
        # count is below SHAPE_SELECTION_MIN_MOTIFS, so the record must say
        # why energy selection was kept instead of silently reading as a
        # plain energy-selected build.
        monomer = build_rdkit_monomer(
            "ppd",
            "ppd",
            PPD,
            "amine",
            num_conformers=2,
            select_conformer_by_motif_shape=True,
        )
        metadata = monomer.metadata
        self.assertEqual(metadata["shape_aware_selection_requested"], True)
        self.assertEqual(metadata["conformer_selection"], "energy")
        self.assertIn("below", metadata["conformer_selection_note"])


class CLIDefaultOwnershipTests(unittest.TestCase):
    """CLI argparse defaults reference the semantic owners instead of
    retyping literals, where they represent the same setting."""

    def _subparser(self, add_parser, name):
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers()
        add_parser(subparsers)
        return subparsers.choices[name]

    def test_num_conformers_defaults_reference_owner_constant(self):
        for add_parser, name in (
            (cli_build._add_single_pair_parser, "single-pair"),
            (cli_build._add_ring_forming_parser, "ring-forming"),
            (cli_build._add_batch_binary_bridge_parser, "batch-binary-bridge"),
        ):
            parser = self._subparser(add_parser, name)
            self.assertEqual(
                parser.get_default("num_conformers"),
                cli_build.DEFAULT_CLI_NUM_CONFORMERS,
                f"{name} --num-conformers default drifted from DEFAULT_CLI_NUM_CONFORMERS",
            )

    def test_embedding_seed_defaults_reference_owner_constant(self):
        # BatchGenerationConfig and the builder signatures share the chem.rdkit
        # owner constant for the default RDKit embedding seed.
        self.assertEqual(BatchGenerationConfig().rdkit_random_seed, DEFAULT_RDKIT_RANDOM_SEED)
        signature = inspect.signature(build_rdkit_monomer)
        self.assertEqual(
            signature.parameters["random_seed"].default,
            DEFAULT_RDKIT_RANDOM_SEED,
        )

    def test_autodetect_budget_references_owner_constant(self):
        # The two-conformer autodetection clamp is owned by monomer_library.
        generator = BatchStructureGenerator(BatchGenerationConfig(rdkit_num_conformers=8))
        self.assertEqual(generator._autodetect_num_conformers(), AUTODETECT_MAX_CONFORMERS)
        small = BatchStructureGenerator(BatchGenerationConfig(rdkit_num_conformers=1))
        self.assertEqual(small._autodetect_num_conformers(), 1)


@unittest.skipIf(Chem is None, "RDKit is not available")
class SinglePairSummaryProvenanceTests(unittest.TestCase):
    """A single-pair build's summary.json carries requested vs actual
    conformer budget, seed, selection mode, and fallback status."""

    def test_single_pair_summary_records_conformer_settings(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "single_pair"
            with contextlib.redirect_stdout(io.StringIO()):
                cli_main(
                    [
                        "build",
                        "single-pair",
                        "--template-id",
                        "imine_bridge",
                        "--first-smiles",
                        TAPB,
                        "--second-smiles",
                        TFB,
                        "--first-id",
                        "tapb",
                        "--second-id",
                        "tfb",
                        "--first-motif-kind",
                        "amine",
                        "--second-motif-kind",
                        "aldehyde",
                        "--topology",
                        "hcb",
                        "--num-conformers",
                        "2",
                        "--no-write-cif",
                        "--output-dir",
                        str(output_dir),
                    ]
                )
            summary = json.loads((output_dir / "summary.json").read_text(encoding="utf-8"))

        for key in ("first", "second"):
            conformers = summary[key]["conformers"]
            # Requested CLI budget vs the effective budget after the
            # shape-aware ensemble floor (both monomers are 3-connecting).
            self.assertEqual(conformers["requested_num_conformers"], 2)
            self.assertEqual(conformers["effective_num_conformers"], SHAPE_SELECTION_MIN_CONFORMERS)
            self.assertGreaterEqual(conformers["actual_num_conformers"], 1)
            self.assertLessEqual(conformers["actual_num_conformers"], conformers["effective_num_conformers"])
            self.assertEqual(conformers["random_seed"], DEFAULT_RDKIT_RANDOM_SEED)
            self.assertEqual(conformers["shape_aware_requested"], True)
            self.assertEqual(conformers["shape_aware_applied"], True)
            self.assertEqual(conformers["conformer_selection"], "motif_shape")
            self.assertEqual(conformers["embedding_fallback"], False)
            self.assertEqual(conformers["cofkit_version"], __version__)
            self.assertTrue(conformers["rdkit_version"])

        # Per-structure provenance on the result rows matches the monomers.
        result = summary["results"][0]
        provenance = result["metadata"]["reactant_conformer_provenance"]
        self.assertEqual(set(provenance), {"tapb", "tfb"})
        for entry in provenance.values():
            self.assertEqual(entry["effective_num_conformers"], SHAPE_SELECTION_MIN_CONFORMERS)
            self.assertEqual(entry["random_seed"], DEFAULT_RDKIT_RANDOM_SEED)
            self.assertEqual(entry["conformer_selection"], "motif_shape")
            self.assertEqual(entry["embedding_fallback"], False)


@unittest.skipIf(Chem is None, "RDKit is not available")
class BatchRunProvenanceTests(unittest.TestCase):
    """A batch run's durable records carry the effective settings per record
    and reconcile with the run summary."""

    def _write_tiny_library(self, root: Path) -> None:
        root.mkdir(parents=True, exist_ok=True)
        (root / "amines_count_3.txt").write_text(f"smiles\n{TAPB}\n", encoding="utf-8")
        (root / "aldehydes_count_3.txt").write_text(f"smiles\n{TFB}\n", encoding="utf-8")

    def test_batch_records_reconcile_with_summary(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            self._write_tiny_library(root / "input")
            generator = BatchStructureGenerator(
                BatchGenerationConfig(write_cif=False, rdkit_num_conformers=2, max_workers=1)
            )

            summary = generator.run_binary_bridge_batch(root / "input", root / "output")

            settings = summary.conformer_settings
            self.assertIsNotNone(settings)
            self.assertEqual(settings.requested_num_conformers, 2)
            self.assertEqual(settings.autodetect_num_conformers, AUTODETECT_MAX_CONFORMERS)
            self.assertEqual(settings.shape_aware_ensemble_floor, SHAPE_SELECTION_MIN_CONFORMERS)
            self.assertTrue(settings.shape_aware_conformer)
            self.assertEqual(settings.random_seed, DEFAULT_RDKIT_RANDOM_SEED)
            self.assertEqual(settings.cofkit_version, __version__)

            # Per-monomer ledger rows carry the same requested budget and the
            # shape-aware ensemble increase as the effective budget.
            monomer_rows = [
                json.loads(line)
                for line in Path(summary.monomer_records_path).read_text(encoding="utf-8").splitlines()
                if line.strip()
            ]
            self.assertEqual(len(monomer_rows), 2)
            for row in monomer_rows:
                provenance = row["conformer_provenance"]
                self.assertEqual(provenance["requested_num_conformers"], settings.requested_num_conformers)
                self.assertEqual(provenance["effective_num_conformers"], SHAPE_SELECTION_MIN_CONFORMERS)
                self.assertGreaterEqual(provenance["actual_num_conformers"], 1)
                self.assertEqual(provenance["random_seed"], settings.random_seed)
                self.assertEqual(provenance["shape_aware_requested"], True)
                self.assertEqual(provenance["shape_aware_applied"], True)
                self.assertEqual(provenance["conformer_selection"], "motif_shape")
                self.assertEqual(provenance["embedding_fallback"], False)
                self.assertEqual(provenance["builder"], "build_rdkit_monomer")
                self.assertEqual(provenance["cofkit_version"], __version__)

            # Per-structure manifest rows carry per-reactant selected-conformer
            # provenance consistent with the ledger.
            manifest_rows = [
                json.loads(line)
                for line in (root / "output" / "manifest.jsonl").read_text(encoding="utf-8").splitlines()
                if line.strip()
            ]
            self.assertTrue(manifest_rows)
            ledger_by_id = {row["id"]: row for row in monomer_rows}
            for row in manifest_rows:
                provenance = row["metadata"]["reactant_conformer_provenance"]
                self.assertEqual(set(provenance), set(ledger_by_id))
                for monomer_id, entry in provenance.items():
                    self.assertEqual(entry["effective_num_conformers"], SHAPE_SELECTION_MIN_CONFORMERS)
                    self.assertEqual(entry["random_seed"], settings.random_seed)
                    self.assertEqual(
                        entry["conformer_selection"],
                        ledger_by_id[monomer_id]["conformer_provenance"]["conformer_selection"],
                    )

            # The markdown summary explains the autodetection clamp, the
            # shape-aware ensemble increase, and where the per-monomer
            # provenance lives.
            summary_text = (root / "output" / "summary.md").read_text(encoding="utf-8")
            self.assertIn("## Construction settings (conformer provenance)", summary_text)
            self.assertIn("autodetection", summary_text)
            self.assertIn(str(SHAPE_SELECTION_MIN_CONFORMERS), summary_text)
            self.assertIn("selected", summary_text)
            self.assertIn("monomers.jsonl", summary_text)

    def test_shape_aware_opt_out_records_energy_selection(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(write_cif=False, rdkit_num_conformers=2, shape_aware_conformer=False)
        )
        record = BatchMonomerRecord(
            id="tapb",
            name="tapb",
            smiles=TAPB,
            motif_kind="amine",
            expected_connectivity=3,
        )
        built = generator.build_monomer(record)
        self.assertTrue(built.ok, built.error)
        provenance = built.conformer_provenance
        self.assertIsNotNone(provenance)
        # No ensemble increase when the shape-aware gate is off.
        self.assertEqual(provenance.effective_num_conformers, 2)
        self.assertEqual(provenance.requested_num_conformers, 2)
        self.assertEqual(provenance.shape_aware_requested, False)
        self.assertEqual(provenance.shape_aware_applied, False)
        self.assertEqual(provenance.conformer_selection, "energy")


if __name__ == "__main__":
    unittest.main()
