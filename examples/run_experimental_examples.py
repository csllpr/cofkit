from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit import BatchGenerationConfig, BatchStructureGenerator, CIFWriter, build_rdkit_monomer, get_topology_hint


@dataclass(frozen=True)
class ExperimentalCase:
    case_id: str
    reaction_template: str
    target_dimensionality: str
    requested_topology: str
    monomer_a_name: str
    monomer_a_smiles: str
    monomer_a_kind: str
    monomer_b_name: str
    monomer_b_smiles: str
    monomer_b_kind: str
    notes: str = ""


CASES: tuple[ExperimentalCase, ...] = (
    ExperimentalCase(
        case_id="cof_300",
        reaction_template="imine_bridge",
        target_dimensionality="3D",
        requested_topology="dia-c5",
        monomer_a_name="TAPM",
        monomer_a_smiles="Nc1ccc(C(c2ccc(N)cc2)(c2ccc(N)cc2)c2ccc(N)cc2)cc1",
        monomer_a_kind="amine",
        monomer_b_name="terephthaldehyde",
        monomer_b_smiles="O=Cc1ccc(C=O)cc1",
        monomer_b_kind="aldehyde",
        notes="The bundled topology index exposes dia rather than dia-c5 as the supported 4+2 diamond-net target.",
    ),
    ExperimentalCase(
        case_id="tp_pa_1",
        reaction_template="keto_enamine_bridge",
        target_dimensionality="2D",
        requested_topology="hcb",
        monomer_a_name="Tp",
        monomer_a_smiles="O=Cc1c(O)c(C=O)c(O)c(C=O)c1O",
        monomer_a_kind="keto_aldehyde",
        monomer_b_name="Pa-1",
        monomer_b_smiles="Nc1ccc(N)cc1",
        monomer_b_kind="amine",
    ),
    ExperimentalCase(
        case_id="cof_42",
        reaction_template="hydrazone_bridge",
        target_dimensionality="2D",
        requested_topology="hcb",
        monomer_a_name="2,5-diethoxyterephthalohydrazide",
        monomer_a_smiles="CCOc1cc(C(=O)NN)cc(C(=O)NN)c1OCC",
        monomer_a_kind="hydrazide",
        monomer_b_name="1,3,5-triformylbenzene",
        monomer_b_smiles="O=Cc1cc(C=O)cc(C=O)c1",
        monomer_b_kind="aldehyde",
        notes="The earlier bnn/gra note was removed after rechecking the assignment; this example now targets the 3+2 hcb net.",
    ),
    ExperimentalCase(
        case_id="azine_tfb_hydrazine",
        reaction_template="azine_bridge",
        target_dimensionality="2D",
        requested_topology="hcb",
        monomer_a_name="hydrazine",
        monomer_a_smiles="NN",
        monomer_a_kind="hydrazine",
        monomer_b_name="1,3,5-triformylbenzene",
        monomer_b_smiles="O=Cc1cc(C=O)cc(C=O)c1",
        monomer_b_kind="aldehyde",
        notes="Experimental azine smoke case using the minimal ditopic hydrazine linker against a trigonal trialdehyde node on hcb.",
    ),
    ExperimentalCase(
        case_id="cof_5",
        reaction_template="boronate_ester_bridge",
        target_dimensionality="2D",
        requested_topology="hcb",
        monomer_a_name="BDBA",
        monomer_a_smiles="OB(O)c1ccc(B(O)O)cc1",
        monomer_a_kind="boronic_acid",
        monomer_b_name="HHTP",
        monomer_b_smiles="OC1=C(O)C=C2C(=C1)C1=CC(O)=C(O)C=C1C1=CC(O)=C(O)C=C21",
        monomer_b_kind="catechol",
        notes="The primary COF-5 paper discusses eclipsed layered stacking; the current builder only generates a single-layer 2D sheet without the layered stacking model.",
    ),
    ExperimentalCase(
        case_id="v_cof_1",
        reaction_template="vinylene_bridge",
        target_dimensionality="2D",
        requested_topology="hcb",
        monomer_a_name="TMT",
        monomer_a_smiles="Cc1nc(C)nc(C)n1",
        monomer_a_kind="activated_methylene",
        monomer_b_name="terephthaldehyde",
        monomer_b_smiles="O=Cc1ccc(C=O)cc1",
        monomer_b_kind="aldehyde",
        notes="This checks the new activated-methylene support for aza-heteroaryl methyl precursors in the current coarse vinylene builder.",
    ),
)


def _resolve_requested_topology(topology_id: str) -> tuple[str, tuple[str, ...]]:
    if topology_id == "dia-c5":
        return "dia", ("requested dia-c5, using bundled dia entry",)
    return topology_id, ()


def run_case(case: ExperimentalCase, *, output_dir: Path, num_conformers: int) -> dict[str, object]:
    case_dir = output_dir / case.case_id
    case_dir.mkdir(parents=True, exist_ok=True)
    topology_id, topology_notes = _resolve_requested_topology(case.requested_topology)
    record: dict[str, object] = {
        "case_id": case.case_id,
        "requested_topology": case.requested_topology,
        "resolved_topology": topology_id,
        "reaction_template": case.reaction_template,
        "target_dimensionality": case.target_dimensionality,
        "notes": list(topology_notes),
    }
    try:
        monomer_a = build_rdkit_monomer(
            case.monomer_a_name,
            case.monomer_a_name,
            case.monomer_a_smiles,
            case.monomer_a_kind,
            num_conformers=num_conformers,
        )
        monomer_b = build_rdkit_monomer(
            case.monomer_b_name,
            case.monomer_b_name,
            case.monomer_b_smiles,
            case.monomer_b_kind,
            num_conformers=num_conformers,
        )
        record["motif_counts"] = {
            case.monomer_a_name: len(monomer_a.motifs),
            case.monomer_b_name: len(monomer_b.motifs),
        }
    except Exception as exc:
        record["status"] = "monomer-build-failed"
        record["error"] = f"{type(exc).__name__}: {exc}"
        return record

    try:
        hint = get_topology_hint(topology_id)
        record["topology_metadata"] = {
            "dimensionality": hint.dimensionality,
            "two_monomer_node_node_modes": tuple(hint.metadata.get("two_monomer_node_node_modes", ())),
            "two_monomer_node_linker_modes": tuple(hint.metadata.get("two_monomer_node_linker_modes", ())),
            "two_monomer_reason": hint.metadata.get("two_monomer_reason"),
        }
    except Exception as exc:
        record["status"] = "topology-metadata-failed"
        record["error"] = f"{type(exc).__name__}: {exc}"
        return record

    generator = BatchStructureGenerator(
        BatchGenerationConfig(
            allowed_reactions=(case.reaction_template,),
            topology_ids=(topology_id,),
            write_cif=False,
        )
    )
    try:
        summaries, candidates, attempted_structures = generator.generate_monomer_pair_candidates(
            monomer_a,
            monomer_b,
            write_cif=False,
        )
    except Exception as exc:
        record["status"] = "generation-failed"
        record["error"] = f"{type(exc).__name__}: {exc}"
        return record

    record["attempted_structures"] = attempted_structures
    if not candidates:
        summary = summaries[0] if summaries else None
        record["status"] = "generation-failed"
        if summary is not None:
            error = summary.metadata.get("error")
            if error is None:
                error = summary.metadata.get("failed_topologies")
            record["error"] = str(error or summary.status)
            record["pair_mode"] = summary.pair_mode
        return record

    candidate = candidates[0]
    cif_path = case_dir / f"{case.case_id}__{topology_id}.cif"
    export = CIFWriter().write_candidate(cif_path, candidate, (monomer_a, monomer_b), data_name=case.case_id)
    embedding_metadata = candidate.metadata.get("embedding", {})
    graph_summary = candidate.metadata.get("graph_summary", {})
    record["status"] = "ok"
    record["cif_path"] = str(cif_path)
    record["score"] = candidate.score
    record["topology"] = candidate.metadata["net_plan"]["topology"]
    record["placement_mode"] = embedding_metadata.get("placement_mode", embedding_metadata.get("mode", "unknown"))
    record["n_reaction_events"] = graph_summary.get("n_reaction_events", len(candidate.events))
    record["n_monomer_instances"] = graph_summary.get("n_monomer_instances")
    record["cif_mode"] = export.mode
    record["reaction_realization"] = export.metadata.get("reaction_realization")
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description="Run literature-inspired single-pair COF examples.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("out/experimental_examples"),
        help="Directory to write CIFs and the run report into.",
    )
    parser.add_argument(
        "--num-conformers",
        type=int,
        default=4,
        help="RDKit conformers to sample per monomer.",
    )
    args = parser.parse_args()

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    records = [run_case(case, output_dir=output_dir, num_conformers=args.num_conformers) for case in CASES]
    (output_dir / "report.json").write_text(json.dumps(records, indent=2), encoding="utf-8")

    lines = ["# Experimental Example Report", ""]
    for case, record in zip(CASES, records):
        lines.append(f"## {case.case_id}")
        lines.append(f"- Requested topology: `{case.requested_topology}`")
        lines.append(f"- Resolved topology: `{record.get('resolved_topology')}`")
        lines.append(f"- Reaction template: `{case.reaction_template}`")
        lines.append(f"- Status: `{record.get('status')}`")
        if case.notes:
            lines.append(f"- Literature note: {case.notes}")
        if record.get("error"):
            lines.append(f"- Error: `{record['error']}`")
        if record.get("topology_metadata"):
            topology_metadata = record["topology_metadata"]
            lines.append(f"- Node-node modes: `{topology_metadata['two_monomer_node_node_modes']}`")
            lines.append(f"- Node-linker modes: `{topology_metadata['two_monomer_node_linker_modes']}`")
            lines.append(f"- Compatibility note: {topology_metadata['two_monomer_reason']}")
        if record.get("cif_path"):
            lines.append(f"- CIF: `{record['cif_path']}`")
            lines.append(f"- Placement mode: `{record['placement_mode']}`")
            lines.append(f"- Reaction events: `{record['n_reaction_events']}`")
            lines.append(f"- Monomer instances: `{record['n_monomer_instances']}`")
            lines.append(f"- Score: `{record['score']:.3f}`")
        lines.append("")
    (output_dir / "report.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
