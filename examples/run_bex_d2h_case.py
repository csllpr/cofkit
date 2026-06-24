from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit import BatchGenerationConfig, BatchMonomerRecord, BatchStructureGenerator, CoarseValidationThresholds


BEX_D2H_ALDEHYDE = (
    "C1C=C(N(C2C=CC(C3C4C(=NON=4)C(C4C=CC(N(C5C=CC(C=O)=CC=5)C5C=CC(C=O)=CC=5)=CC=4)=CC=3)=CC=2)"
    "C2C=CC(C=O)=CC=2)C=CC=1C=O"
)
BEX_D2H_AMINE = (
    "Nc1ccc(-c2ccc3c(c2)C(=C2c4cc(-c5ccc(N)cc5)ccc4-c4ccc(-c5ccc(N)cc5)cc42)c2cc(-c4ccc(N)cc4)ccc2-3)cc1"
)


def main() -> None:
    output_dir = Path(__file__).resolve().parents[1] / "out" / "bex_d2h"
    generator = BatchStructureGenerator(
        BatchGenerationConfig(
            rdkit_num_conformers=4,
            retain_top_results=5,
            topology_ids=("bex",),
            hard_hard_max_bridge_distance=10.0,
            validation_thresholds=CoarseValidationThresholds(hard_hard_max_bridge_distance=10.0),
            write_cif=True,
        )
    )
    summary, candidate = generator.generate_pair_candidate(
        BatchMonomerRecord(
            id="bex_d2h_amine",
            name="bex_d2h_amine",
            smiles=BEX_D2H_AMINE,
            motif_kind="amine",
            expected_connectivity=4,
        ),
        BatchMonomerRecord(
            id="bex_d2h_aldehyde",
            name="bex_d2h_aldehyde",
            smiles=BEX_D2H_ALDEHYDE,
            motif_kind="aldehyde",
            expected_connectivity=4,
        ),
        out_dir=output_dir,
        write_cif=True,
    )
    print(f"status: {summary.status}")
    print(f"pair_mode: {summary.pair_mode}")
    print(f"topology: {summary.topology_id}")
    print(f"cif_path: {summary.cif_path}")
    print(f"validation: {summary.metadata.get('validation')}")
    if candidate is not None:
        print(f"placement_mode: {candidate.metadata['embedding']['placement_mode']}")
        print(f"assignment: {candidate.metadata['assignment']}")
        print(f"graph_summary: {candidate.metadata['graph_summary']}")
        bridge_distances = [
            round(metric["actual_distance"], 3)
            for metric in candidate.metadata["score_metadata"]["bridge_event_metrics"]
        ]
        print(f"bridge_distances_A: {bridge_distances}")
        print("note: this is a rigid-conformer decorated-bex placement; hard_invalid bridge metrics indicate that conformer relaxation is still needed.")


if __name__ == "__main__":
    main()
