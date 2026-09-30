import sys
import tempfile
import unittest
import unittest.mock
from pathlib import Path


from cofkit import BatchGenerationConfig, BatchMonomerRecord, BatchStructureGenerator, default_motif_kind_registry
from cofkit.cofid import (
    cofid_comment_line,
    cofid_to_build_request,
    generate_cofid_with_cause,
    read_cofid_from_cif,
)

try:
    from rdkit import Chem
except ImportError:  # pragma: no cover - environment-dependent
    Chem = None


PPD = "Nc1ccc(N)cc1"
TP = "O=Cc1c(O)c(C=O)c(O)c(C=O)c1O"


def _canonical(smiles: str) -> str:
    assert Chem is not None
    molecule = Chem.MolFromSmiles(smiles)
    assert molecule is not None
    return str(Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False))


@unittest.skipIf(Chem is None, "RDKit is not available")
class COFidTests(unittest.TestCase):
    def test_read_cofid_from_cif_ignores_comment_suffix_tokens(self):
        cofid = f"3:aldehyde:{_canonical(TP)}.2:amine:{_canonical(PPD)}&&hcb&&imine"

        with tempfile.TemporaryDirectory() as temp_dir:
            cif_path = Path(temp_dir) / "stacked.cif"
            cif_path.write_text(
                cofid_comment_line(cofid, suffix="stacking=AA") + "\n" + "data_demo\n",
                encoding="utf-8",
            )

            parsed = read_cofid_from_cif(cif_path)

        self.assertEqual(parsed, cofid)

    def test_cofid_spec_reactive_group_table_covers_builtin_motif_kinds(self):
        spec_path = Path(__file__).resolve().parents[1] / "docs" / "COFid_Specification_v1.2.md"
        spec_groups = {
            line.split("|")[1].strip().strip("`")
            for line in spec_path.read_text(encoding="utf-8").splitlines()
            if line.startswith("| `")
        }

        self.assertTrue(set(default_motif_kind_registry().supported_kinds()).issubset(spec_groups))

    def test_cofid_build_request_requires_explicit_keto_aldehyde_group_for_beta_ketoenamine(self):
        cofid = f"3:keto_aldehyde:{_canonical(TP)}.2:amine:{_canonical(PPD)}&&hcb&&bken"

        request = cofid_to_build_request(cofid)

        self.assertEqual(request.template_id, "keto_enamine_bridge")
        self.assertEqual(request.target_dimensionality, "2D")
        self.assertEqual(tuple(monomer.motif_kind for monomer in request.monomers), ("keto_aldehyde", "amine"))

    def test_cofid_build_request_rejects_beta_ketoenamine_when_group_is_generic_aldehyde(self):
        cofid = f"3:aldehyde:{_canonical(TP)}.2:amine:{_canonical(PPD)}&&hcb&&bken"

        with self.assertRaisesRegex(ValueError, "do not match linkage"):
            cofid_to_build_request(cofid)

    def test_batch_summary_generates_spec_style_bken_cofid(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                allowed_reactions=("keto_enamine_bridge",),
                rdkit_num_conformers=2,
                single_node_topology_ids=("hcb",),
                write_cif=False,
            )
        )
        amine = BatchMonomerRecord(
            id="ppd",
            name="ppd",
            smiles=PPD,
            motif_kind="amine",
            expected_connectivity=2,
        )
        keto_aldehyde = BatchMonomerRecord(
            id="tp",
            name="tp",
            smiles=TP,
            motif_kind="keto_aldehyde",
            expected_connectivity=3,
        )

        summary, candidate = generator.generate_pair_candidate(amine, keto_aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertIsNotNone(candidate)
        self.assertEqual(summary.metadata["cofid"], f"3:keto_aldehyde:{_canonical(TP)}.2:amine:{_canonical(PPD)}&&hcb&&bken")

    def test_generate_cofid_with_cause_preserves_failure_cause(self):
        # A12 / T2-12: best-effort COFid generation must retain the
        # "TypeName: message" cause instead of returning a bare None.
        from cofkit.cofid import try_generate_cofid
        from cofkit.model import AssemblyState, Candidate

        # No net_plan topology in metadata: generate_cofid raises ValueError.
        candidate = Candidate(id="broken", score=None, state=AssemblyState(), events=(), metadata={})

        outcome = generate_cofid_with_cause(candidate, {})

        self.assertIsNone(outcome.cofid)
        self.assertIsNotNone(outcome.error)
        self.assertTrue(outcome.error.startswith("ValueError: "))
        self.assertIn("topology", outcome.error)
        # The legacy helper keeps its None-on-failure contract.
        self.assertIsNone(try_generate_cofid(candidate, {}))

    def test_batch_summary_records_cofid_error_when_generation_fails(self):
        # A12 / T2-12: when COFid generation fails for a built structure, the
        # manifest record carries the cause as metadata.cofid_error.
        import cofkit.batch as batch_module
        from cofkit.cofid import COFidGenerationOutcome

        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                rdkit_num_conformers=1,
                single_node_topology_ids=("hcb",),
                write_cif=False,
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
            id="tfb",
            name="tfb",
            smiles="C1=C(C=C(C=C1C=O)C=O)C=O",
            motif_kind="aldehyde",
            expected_connectivity=3,
        )
        real_outcome = generate_cofid_with_cause

        def failing_outcome(candidate, specs):
            outcome = real_outcome(candidate, specs)
            if outcome.error is not None:
                return outcome
            return COFidGenerationOutcome(cofid=None, error="RuntimeError: injected cofid failure")

        with unittest.mock.patch.object(batch_module, "generate_cofid_with_cause", failing_outcome):
            summary, candidate = generator.generate_pair_candidate(amine, aldehyde)

        self.assertEqual(summary.status, "ok")
        self.assertIsNotNone(candidate)
        self.assertNotIn("cofid", summary.metadata)
        self.assertEqual(summary.metadata["cofid_error"], "RuntimeError: injected cofid failure")


if __name__ == "__main__":
    unittest.main()
