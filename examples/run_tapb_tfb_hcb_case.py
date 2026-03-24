from __future__ import annotations

import configparser
import sys
from dataclasses import dataclass
from math import dist
from pathlib import Path
from typing import Iterable

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit import CIFWriter, COFEngine, COFProject, MonomerSpec, build_rdkit_monomer


@dataclass(frozen=True)
class WorkspaceMonomerInput:
    section: str
    name: str
    smiles: str


def main() -> None:
    workspace = Path(__file__).resolve().parents[1]
    smiles_path = workspace / "out" / "tapb_tfb_smiles.txt"
    cif_path = workspace / "out" / "tapb_tfb_hcb_rdkit_imine_product.cif"
    report_path = workspace / "out" / "tapb_tfb_hcb_rdkit_imine_product_report.md"

    amine_input, aldehyde_input = load_workspace_inputs(smiles_path)
    tapb = build_rdkit_monomer("tapb", amine_input.name, amine_input.smiles, "amine")
    tfb = build_rdkit_monomer("tfb", aldehyde_input.name, aldehyde_input.smiles, "aldehyde")

    project = COFProject(
        monomers=(tapb, tfb),
        allowed_reactions=("imine_bridge",),
        target_dimensionality="2D",
        target_topologies=("hcb",),
    )
    candidate = COFEngine().run(project).top(1)[0]
    cif_result = CIFWriter().write_candidate(cif_path, candidate, project.monomers, data_name="tapb_tfb_hcb_rdkit")

    report_path.write_text(render_report(smiles_path, amine_input, aldehyde_input, tapb, tfb, candidate, cif_result))


def load_workspace_inputs(path: Path) -> tuple[WorkspaceMonomerInput, WorkspaceMonomerInput]:
    parser = configparser.ConfigParser()
    parser.read(path)
    return (
        WorkspaceMonomerInput(
            section="amine",
            name=parser["amine"]["name"].strip(),
            smiles=parser["amine"]["smiles"].strip(),
        ),
        WorkspaceMonomerInput(
            section="aldehyde",
            name=parser["aldehyde"]["name"].strip(),
            smiles=parser["aldehyde"]["smiles"].strip(),
        ),
    )

def render_report(
    smiles_path: Path,
    amine_input: WorkspaceMonomerInput,
    aldehyde_input: WorkspaceMonomerInput,
    tapb: MonomerSpec,
    tfb: MonomerSpec,
    candidate,
    cif_result,
) -> str:
    pose_a = candidate.state.monomer_poses["m1"]
    pose_b = candidate.state.monomer_poses["m2"]
    center_distance = dist(pose_a.translation, pose_b.translation)
    bridge_metrics = candidate.metadata["score_metadata"]["bridge_event_metrics"]
    periodic_images = [event.participants[1].periodic_image for event in candidate.events]
    lines = [
        "# TAPB/TFB hcb RDKit imine-product report",
        "",
        f"- SMILES source: `{smiles_path}`",
        f"- Amine: {amine_input.name}",
        f"- Aldehyde: {aldehyde_input.name}",
        f"- Requested topology: hcb",
        f"- Selected topology: {candidate.metadata['net_plan']['topology']}",
        f"- Embedding mode: {candidate.metadata['embedding']['placement_mode']}",
        f"- CIF export mode: {cif_result.mode}",
        f"- Reacted events realized in CIF: {cif_result.metadata.get('reaction_realization', {}).get('applied_event_count', 0)}",
        "",
        "## Input SMILES",
        "",
        f"- amine: `{amine_input.smiles}`",
        f"- aldehyde: `{aldehyde_input.smiles}`",
        "",
        "## Geometry summary",
        "",
        f"- TAPB RDKit motifs: {len(tapb.motifs)} reactive amines across {tapb.metadata['n_atoms']} atoms",
        f"- TFB RDKit motifs: {len(tfb.motifs)} reactive aldehydes across {tfb.metadata['n_atoms']} atoms",
        f"- TAPB conformer energy ({tapb.metadata['forcefield']}): {float(tapb.metadata['selected_conformer_energy']):.6f}",
        f"- TFB conformer energy ({tfb.metadata['forcefield']}): {float(tfb.metadata['selected_conformer_energy']):.6f}",
        f"- TAPB motif detector: {tapb.metadata['motif_detection']}",
        f"- TFB motif detector: {tfb.metadata['motif_detection']}",
        f"- Center-to-center nearest-neighbor distance: {center_distance:.6f} A",
        f"- Target imine reactive-site separation: {candidate.metadata['embedding']['target_distance']:.6f} A",
        f"- Hexagonal reactive-site-based separation: {candidate.metadata['embedding']['reactive_site_distance']:.6f} A",
        f"- Bridge periodic images used for the TFB node: {periodic_images}",
        f"- Removed atom symbols during imine realization: {cif_result.metadata.get('reaction_realization', {}).get('removed_atom_symbols', {})}",
        "",
        "## Bridge metrics",
        "",
        *format_bridge_metrics(bridge_metrics),
        "",
        "## Files",
        "",
        "- CIF: `out/tapb_tfb_hcb_rdkit_imine_product.cif`",
        "- Report: `out/tapb_tfb_hcb_rdkit_imine_product_report.md`",
        f"- Candidate score: {candidate.score:.6f}",
        "",
        "## Chemical representation",
        "",
        "- Each imine event removes the reacting aldehyde oxygen and one hydrogen from the reacting amine nitrogen.",
        "- The CIF now includes explicit bond-loop records for the retained monomer connectivity plus the realized inter-monomer imine bridges, including periodic-image bridges in the hcb cell.",
        "- The retained carbon-bound aldehydic hydrogen becomes the imine carbon hydrogen in the exported product geometry.",
        "",
        "## Scientific limitations",
        "",
        "- The product is still rigid-body approximate: cofkit removes the leaving atoms and records the reacted connectivity, but it does not reoptimize local bond angles, refine CIF bond orders, or model the water byproduct geometry after condensation.",
        "- RDKit supplies chemically contextual motif detection and a low-energy conformer, but cofkit still treats each monomer as a rigid body during embedding and optimization.",
        "- The 2D hcb placement uses a planarized reactive-site model; out-of-plane torsions and stacking remain intentionally out of scope.",
        "- Stacking remains deliberately disabled and out of scope.",
    ]
    return "\n".join(lines) + "\n"


def format_bridge_metrics(metrics: Iterable[dict[str, object]]) -> list[str]:
    lines: list[str] = []
    for metric in metrics:
        lines.append(
            "- "
            f"{metric['event_id']}: image-aware reactive-site distance {float(metric['actual_distance']):.6f} A, "
            f"target {float(metric['target_distance']):.6f} A, "
            f"residual {float(metric['total_residual']):.6e}"
        )
    return lines


if __name__ == "__main__":
    main()
