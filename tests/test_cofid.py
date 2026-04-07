import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit import BatchGenerationConfig, BatchMonomerRecord, BatchStructureGenerator, default_motif_kind_registry
from cofkit.cofid import cofid_to_build_request

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
    def test_cofid_spec_reactive_group_table_covers_builtin_motif_kinds(self):
        spec_path = Path(__file__).resolve().parents[1] / "docs" / "COFid_Specification_v1.1.md"
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


if __name__ == "__main__":
    unittest.main()
