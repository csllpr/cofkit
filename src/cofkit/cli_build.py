from __future__ import annotations

import argparse
import json
import shutil
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable

from .batch import BatchGenerationConfig, BatchStructureGenerator
from .chem.rdkit import build_rdkit_monomer
from .monomer_library import MonomerRoleResolver
from .reactions import ReactionLibrary


def add_build_group(subparsers) -> None:
    parser = subparsers.add_parser(
        "build",
        help="Build COF structures and build-time inputs.",
        description="Build COF structures, inspect build-capable templates, and generate detector-scanned monomer libraries.",
    )
    _set_help_default(parser)

    build_subparsers = parser.add_subparsers(dest="build_command")
    _add_single_pair_parser(build_subparsers)
    _add_batch_binary_bridge_parser(build_subparsers)
    _add_batch_all_binary_bridges_parser(build_subparsers)
    _add_default_library_parser(build_subparsers)
    _add_list_templates_parser(build_subparsers)


def _set_help_default(parser: argparse.ArgumentParser) -> None:
    def _show_help(_: argparse.Namespace) -> None:
        parser.print_help()

    parser.set_defaults(func=_show_help)


def _add_common_batch_generation_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--target-dimensionality",
        choices=("2D", "3D"),
        default="2D",
        help="Target dimensionality for topology selection. Defaults to 2D.",
    )
    parser.add_argument("--input-dir", required=False, help="Directory containing monomer libraries.")
    parser.add_argument("--output-dir", required=False, help="Directory for manifests, summaries, and CIF exports.")
    parser.add_argument("--max-pairs", type=int, default=None, help="Optional cap on attempted pairings.")
    parser.add_argument(
        "--write-cif",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Write CIF files for generated structures. Enabled by default; pass --no-write-cif to disable.",
    )
    parser.add_argument(
        "--max-cif-exports",
        type=int,
        default=None,
        help="Optional maximum number of CIF files to export.",
    )
    parser.add_argument(
        "--num-conformers",
        type=int,
        default=4,
        help="Number of RDKit conformers to sample per monomer during batch construction.",
    )
    parser.add_argument(
        "--all-topologies",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Enumerate all applicable topologies per monomer pair. Enabled by default.",
    )
    parser.add_argument(
        "--topology",
        action="append",
        default=[],
        help="Explicit topology id to evaluate. Repeat to request multiple indexed topologies.",
    )
    parser.add_argument(
        "--use-indexed-topology-defaults",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Augment legacy defaults with indexed two-monomer topology discovery for compatible pairs.",
    )
    parser.add_argument(
        "--auto-detect-libraries",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Infer monomer role and connectivity from SMILES instead of relying on role/count filename prefixes.",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=8,
        help="Worker budget for pair generation. Defaults to 8.",
    )


def _configure_generator(args: argparse.Namespace, *, template_id: str | None = None) -> BatchStructureGenerator:
    allowed_reactions = (template_id or getattr(args, "template_id", "imine_bridge"),)
    return BatchStructureGenerator(
        BatchGenerationConfig(
            allowed_reactions=allowed_reactions,
            target_dimensionality=args.target_dimensionality,
            topology_ids=tuple(args.topology),
            use_indexed_topology_defaults=args.use_indexed_topology_defaults,
            rdkit_num_conformers=args.num_conformers,
            enumerate_all_topologies=args.all_topologies,
            post_build_conversions=tuple(getattr(args, "annotate_post_build_conversion", ())),
            write_cif=args.write_cif,
            max_cif_exports=args.max_cif_exports,
            max_workers=args.max_workers,
        )
    )


def _add_single_pair_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "single-pair",
        help="Generate one COF candidate set from two monomers.",
        description="Generate one single-pair COF candidate set from two monomers using an explicit or autodetected binary-bridge template.",
    )
    parser.add_argument("--template-id", default="imine_bridge", help="Binary bridge template to use.")
    parser.add_argument("--first-smiles", required=True, help="SMILES for the first monomer.")
    parser.add_argument("--second-smiles", required=True, help="SMILES for the second monomer.")
    parser.add_argument("--first-id", default="monomer_a", help="Identifier for the first monomer.")
    parser.add_argument("--second-id", default="monomer_b", help="Identifier for the second monomer.")
    parser.add_argument("--first-name", default=None, help="Display name for the first monomer.")
    parser.add_argument("--second-name", default=None, help="Display name for the second monomer.")
    parser.add_argument("--first-motif-kind", default=None, help="Explicit motif kind for the first monomer.")
    parser.add_argument("--second-motif-kind", default=None, help="Explicit motif kind for the second monomer.")
    parser.add_argument(
        "--auto-detect-motifs",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Infer motif kinds from SMILES when explicit motif kinds are not provided. Enabled by default.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(Path(__file__).resolve().parents[2] / "out" / "single_pair_generation"),
        help="Directory for the single-pair summary and CIF outputs.",
    )
    parser.add_argument(
        "--write-cif",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Write CIF files for generated structures. Enabled by default; pass --no-write-cif to disable.",
    )
    parser.add_argument(
        "--num-conformers",
        type=int,
        default=4,
        help="Number of RDKit conformers to sample per monomer during monomer construction.",
    )
    parser.add_argument(
        "--all-topologies",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Enumerate all applicable topologies per monomer pair. Enabled by default.",
    )
    parser.add_argument(
        "--topology",
        action="append",
        default=[],
        help="Explicit topology id to evaluate. Repeat to request multiple indexed topologies.",
    )
    parser.add_argument(
        "--use-indexed-topology-defaults",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Augment legacy defaults with indexed two-monomer topology discovery for compatible pairs.",
    )
    parser.add_argument(
        "--target-dimensionality",
        choices=("2D", "3D"),
        default="2D",
        help="Target dimensionality for topology selection. Defaults to 2D.",
    )
    parser.add_argument(
        "--max-cif-exports",
        type=int,
        default=None,
        help="Optional maximum number of CIF files to export. Defaults to no limit.",
    )
    parser.add_argument("--max-workers", type=int, default=1, help="Worker budget. Defaults to 1 for single-pair mode.")
    parser.set_defaults(func=_run_single_pair)


def _run_single_pair(args: argparse.Namespace) -> None:
    library = ReactionLibrary.builtin()
    profile = library.linkage_profile(args.template_id)
    if profile is None or not profile.supports_binary_bridge_pair_generation:
        raise SystemExit(f"template {args.template_id!r} does not support current single-pair binary-bridge generation")

    allowed_motif_kinds = tuple(role.motif_kind for role in profile.binary_bridge_roles)
    resolver = MonomerRoleResolver.builtin()

    first_kind = _resolve_single_pair_motif_kind(
        smiles=args.first_smiles,
        explicit_kind=args.first_motif_kind,
        monomer_id=args.first_id,
        allowed_motif_kinds=allowed_motif_kinds,
        auto_detect=args.auto_detect_motifs,
        resolver=resolver,
        num_conformers=args.num_conformers,
    )
    second_kind = _resolve_single_pair_motif_kind(
        smiles=args.second_smiles,
        explicit_kind=args.second_motif_kind,
        monomer_id=args.second_id,
        allowed_motif_kinds=allowed_motif_kinds,
        auto_detect=args.auto_detect_motifs,
        resolver=resolver,
        num_conformers=args.num_conformers,
    )

    first = build_rdkit_monomer(
        args.first_id,
        args.first_name or args.first_id,
        args.first_smiles,
        first_kind,
        num_conformers=args.num_conformers,
    )
    second = build_rdkit_monomer(
        args.second_id,
        args.second_name or args.second_id,
        args.second_smiles,
        second_kind,
        num_conformers=args.num_conformers,
    )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    generator = _configure_generator(args, template_id=args.template_id)
    if args.all_topologies:
        summaries, candidates, attempted_structures = generator.generate_monomer_pair_candidates(
            first,
            second,
            pair_id=f"{args.first_id}__{args.second_id}",
            out_dir=output_dir / "cifs",
            write_cif=args.write_cif,
            cif_export_start_index=0,
        )
    else:
        summary, candidate = generator.generate_monomer_pair_candidate(
            first,
            second,
            pair_id=f"{args.first_id}__{args.second_id}",
            out_dir=output_dir / "cifs",
            write_cif=args.write_cif,
            cif_export_index=0,
        )
        summaries = (summary,)
        candidates = () if candidate is None else (candidate,)
        attempted_structures = 1 if summary.status == "ok" else 0

    report = {
        "template_id": args.template_id,
        "target_dimensionality": args.target_dimensionality,
        "first": {"id": args.first_id, "name": args.first_name or args.first_id, "motif_kind": first_kind, "motif_count": len(first.motifs)},
        "second": {"id": args.second_id, "name": args.second_name or args.second_id, "motif_kind": second_kind, "motif_count": len(second.motifs)},
        "attempted_structures": attempted_structures,
        "successful_structures": sum(1 for summary in summaries if summary.status == "ok"),
        "cifs_written": sum(1 for summary in summaries if summary.cif_path is not None),
        "results": [_summary_to_single_pair_result(summary) for summary in summaries],
    }
    (output_dir / "summary.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("template_id:", args.template_id)
    print("target_dimensionality:", args.target_dimensionality)
    print("first_monomer:", args.first_id, first_kind, len(first.motifs))
    print("second_monomer:", args.second_id, second_kind, len(second.motifs))
    print("attempted_structures:", attempted_structures)
    print("successful_structures:", report["successful_structures"])
    print("cifs_written:", report["cifs_written"])
    print("summary:", output_dir / "summary.json")
    for summary in summaries:
        print(
            "result:",
            summary.structure_id,
            summary.status,
            summary.topology_id or "none",
            f"{summary.score:.6f}" if summary.score is not None else "n/a",
            summary.cif_path or "-",
        )


def _resolve_single_pair_motif_kind(
    *,
    smiles: str,
    explicit_kind: str | None,
    monomer_id: str,
    allowed_motif_kinds: Iterable[str],
    auto_detect: bool,
    resolver: MonomerRoleResolver,
    num_conformers: int,
) -> str:
    if explicit_kind is not None:
        return explicit_kind
    if not auto_detect:
        raise SystemExit(f"monomer {monomer_id!r} needs an explicit motif kind or --auto-detect-motifs")
    record = resolver.infer_record(
        smiles,
        record_id=monomer_id,
        allowed_motif_kinds=allowed_motif_kinds,
        num_conformers=num_conformers,
    )
    return record.motif_kind


def _summary_to_single_pair_result(summary) -> dict[str, object]:
    return {
        "structure_id": summary.structure_id,
        "pair_id": summary.pair_id,
        "pair_mode": summary.pair_mode,
        "status": summary.status,
        "reactant_record_ids": dict(summary.reactant_record_ids),
        "reactant_connectivities": dict(summary.reactant_connectivities),
        "topology_id": summary.topology_id,
        "score": summary.score,
        "flags": list(summary.flags),
        "cif_path": summary.cif_path,
        "metadata": summary.metadata,
    }


def _print_batch_summary(summary, *, template_id: str | None = None) -> None:
    if template_id is not None:
        print("template_id:", template_id)
    print("input_dir:", summary.input_dir)
    print("output_dir:", summary.output_dir)
    print("attempted_pairs:", summary.attempted_pairs)
    print("successful_pairs:", summary.successful_pairs)
    print("attempted_structures:", summary.attempted_structures)
    print("successful_structures:", summary.successful_structures)
    print("built_monomers:", summary.built_monomers)
    print("failed_monomers:", summary.failed_monomers)
    print("mode_counts:", dict(summary.mode_counts))
    print("topology_counts:", dict(summary.topology_counts))
    print("cifs_written:", summary.cifs_written)
    print("manifest:", summary.manifest_path)
    for result in summary.top_results[:10]:
        print(
            "top_result:",
            result.structure_id,
            result.pair_mode,
            f"{result.score:.6f}" if result.score is not None else "n/a",
            result.topology_id or "none",
            dict(result.reactant_record_ids),
        )


def _add_batch_binary_bridge_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "batch-binary-bridge",
        help="Run batch generation for one binary-bridge template.",
        description="Run batch binary-bridge COF generation over monomer libraries.",
    )
    parser.add_argument("--template-id", default="imine_bridge", help="Binary bridge template to use.")
    parser.set_defaults(func=_run_batch_binary_bridge)
    _add_common_batch_generation_arguments(parser)
    parser.set_defaults(
        input_dir=str(Path(__file__).resolve().parents[2] / "examples" / "batch_test_monomers"),
        output_dir=str(Path(__file__).resolve().parents[2] / "out" / "binary_bridge_generation"),
    )


def _run_batch_binary_bridge(args: argparse.Namespace) -> None:
    generator = _configure_generator(args, template_id=args.template_id)
    summary = generator.run_binary_bridge_batch(
        args.input_dir,
        args.output_dir,
        template_id=args.template_id,
        max_pairs=args.max_pairs,
        write_cif=args.write_cif,
        auto_detect_libraries=args.auto_detect_libraries,
    )
    _print_batch_summary(summary, template_id=args.template_id)


def _add_batch_all_binary_bridges_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "batch-all-binary-bridges",
        help="Run all currently available binary-bridge templates over one library directory.",
        description="Run all currently available binary-bridge batch generators over one monomer-library directory.",
    )
    _add_common_batch_generation_arguments(parser)
    parser.add_argument(
        "--template-workers",
        type=int,
        default=1,
        help="How many linkage-template batches to run concurrently. Defaults to 1.",
    )
    parser.set_defaults(
        func=_run_batch_all_binary_bridges,
        input_dir=str(Path(__file__).resolve().parents[2] / "examples" / "default_monomers_library"),
        output_dir=str(Path(__file__).resolve().parents[2] / "out" / "available_binary_bridge_batches"),
        num_conformers=2,
    )


def _run_batch_all_binary_bridges(args: argparse.Namespace) -> None:
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    discovery = BatchStructureGenerator(BatchGenerationConfig(max_workers=max(1, args.max_workers)))
    available_templates = discovery.available_binary_bridge_template_ids(
        input_dir,
        auto_detect_libraries=args.auto_detect_libraries,
    )
    if not available_templates:
        raise SystemExit("no available binary-bridge templates were detected for the input directory")

    template_workers = min(max(1, args.template_workers), len(available_templates))
    summaries = {}
    with ThreadPoolExecutor(max_workers=template_workers) as executor:
        futures = {
            executor.submit(_run_template_batch, template_id, args, input_dir, output_dir): template_id
            for template_id in available_templates
        }
        for future in as_completed(futures):
            template_id, summary = future.result()
            summaries[template_id] = summary

    combined_summary = {
        "available_templates": available_templates,
        "template_workers": template_workers,
        "template_summaries": {template_id: _summary_to_json(summary) for template_id, summary in summaries.items()},
    }
    (output_dir / "combined_summary.json").write_text(json.dumps(combined_summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("input_dir:", input_dir)
    print("output_dir:", output_dir)
    print("available_templates:", available_templates)
    print("template_workers:", template_workers)
    for template_id in available_templates:
        _print_batch_summary(summaries[template_id], template_id=template_id)


def _run_template_batch(
    template_id: str,
    args: argparse.Namespace,
    input_dir: Path,
    output_dir: Path,
):
    generator = _configure_generator(args, template_id=template_id)
    summary = generator.run_binary_bridge_batch(
        input_dir,
        output_dir / template_id,
        template_id=template_id,
        max_pairs=args.max_pairs,
        write_cif=args.write_cif,
        auto_detect_libraries=args.auto_detect_libraries,
    )
    return template_id, summary


def _summary_to_json(summary) -> dict[str, object]:
    return {
        "input_dir": summary.input_dir,
        "output_dir": summary.output_dir,
        "attempted_pairs": summary.attempted_pairs,
        "successful_pairs": summary.successful_pairs,
        "attempted_structures": summary.attempted_structures,
        "successful_structures": summary.successful_structures,
        "cifs_written": summary.cifs_written,
        "built_monomers": summary.built_monomers,
        "failed_monomers": summary.failed_monomers,
        "build_failures": dict(summary.build_failures),
        "mode_counts": dict(summary.mode_counts),
        "topology_counts": dict(summary.topology_counts),
        "manifest_path": summary.manifest_path,
    }


def _add_default_library_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "default-library",
        help="Build a detector-scanned default monomer library from source SMILES files.",
    )
    parser.add_argument(
        "--input-dir",
        default=str(Path(__file__).resolve().parents[2] / "examples" / "batch_test_monomers"),
    )
    parser.add_argument(
        "--output-dir",
        default=str(Path(__file__).resolve().parents[2] / "examples" / "default_monomers_library"),
    )
    parser.add_argument("--num-conformers", type=int, default=2)
    parser.set_defaults(func=_run_default_library)


def _run_default_library(args: argparse.Namespace) -> None:
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    generator = BatchStructureGenerator(BatchGenerationConfig(rdkit_num_conformers=args.num_conformers))
    registered_records: list[dict[str, object]] = []
    failed_records: list[dict[str, object]] = []
    grouped_smiles: dict[str, list[str]] = {}
    grouped_entries: dict[str, list[dict[str, object]]] = {}

    for path in sorted(input_dir.glob("*.txt")):
        declared_connectivity = _source_declared_connectivity(path)
        declared_motif_kind = _source_declared_motif_kind(path, generator.reaction_library)
        seen = 0
        with path.open(encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.strip()
                if not line:
                    continue
                if line_number == 1 and line.lower() == "smiles":
                    continue
                seen += 1
                source_record_id = f"{path.stem}_{seen:04d}"
                base_row = {
                    "source_file": path.name,
                    "source_line": line_number,
                    "source_record_id": source_record_id,
                    "source_declared_motif_kind": declared_motif_kind,
                    "source_declared_connectivity": declared_connectivity,
                    "smiles": line,
                }
                try:
                    record = generator.infer_monomer_record(
                        line,
                        record_id=source_record_id,
                        source_path=str(path),
                        source_line=line_number,
                        library_stem=path.stem,
                    )
                except Exception as exc:
                    failed_records.append({**base_row, "status": "failed", "error": f"{type(exc).__name__}: {exc}"})
                    continue

                library_key = f"{generator._library_prefix_for_motif_kind(record.motif_kind)}_count_{record.expected_connectivity}"
                grouped_smiles.setdefault(library_key, []).append(record.smiles)
                entry = {
                    **base_row,
                    "status": "registered",
                    "detected_motif_kind": record.motif_kind,
                    "detected_connectivity": record.expected_connectivity,
                    "detected_motif_kinds": list(record.metadata.get("detected_motif_kinds", ())),
                    "detected_connectivities": dict(record.metadata.get("detected_connectivities", {})),
                    "library_key": library_key,
                    "source_role_match": declared_motif_kind == record.motif_kind if declared_motif_kind is not None else None,
                    "source_connectivity_match": declared_connectivity == record.expected_connectivity if declared_connectivity is not None else None,
                }
                grouped_entries.setdefault(library_key, []).append(entry)
                registered_records.append(entry)

    for library_key, smiles_rows in sorted(grouped_smiles.items()):
        path = output_dir / f"{library_key}.txt"
        path.write_text("smiles\n" + "\n".join(smiles_rows) + "\n", encoding="utf-8")

    registry_rows: list[dict[str, object]] = []
    for library_key, entries in sorted(grouped_entries.items()):
        for index, entry in enumerate(entries, start=1):
            registry_rows.append(
                {
                    **entry,
                    "registered_record_id": f"{library_key}_{index:04d}",
                    "registered_library_path": f"{library_key}.txt",
                }
            )

    _jsonl_write(output_dir / "registry.jsonl", registry_rows)
    _jsonl_write(output_dir / "failures.jsonl", failed_records)

    library_counts = Counter(entry["library_key"] for entry in registered_records)
    mismatch_count = sum(
        1
        for entry in registered_records
        if entry["source_role_match"] is False or entry["source_connectivity_match"] is False
    )
    summary_lines = [
        "# Default Monomers Library",
        "",
        "This directory is auto-generated using the current monomer role detector.",
        "",
        f"- source_dir: `{input_dir}`",
        f"- registered monomers: {len(registered_records)}",
        f"- failed or ambiguous detections: {len(failed_records)}",
        f"- source-label mismatches among registered monomers: {mismatch_count}",
        "",
        "## Registered Libraries",
        "",
    ]
    for library_key, count in sorted(library_counts.items()):
        summary_lines.append(f"- `{library_key}.txt`: {count}")
    summary_lines.extend(
        (
            "",
            "## Metadata Files",
            "",
            "- `registry.jsonl`: registered monomers with detector metadata and source provenance",
            "- `failures.jsonl`: failed or ambiguous autodetections",
        )
    )
    (output_dir / "README.md").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    print("output_dir:", output_dir)
    print("registered_monomers:", len(registered_records))
    print("failed_monomers:", len(failed_records))
    print("libraries:", dict(sorted(library_counts.items())))


def _jsonl_write(path: Path, rows: Iterable[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")


def _source_declared_connectivity(path: Path) -> int | None:
    stem = path.stem
    if "_count_" not in stem:
        return None
    try:
        return int(stem.rsplit("_", 1)[1])
    except ValueError:
        return None


def _source_declared_motif_kind(path: Path, reaction_library: ReactionLibrary) -> str | None:
    stem = path.stem
    if "_count_" not in stem:
        return None
    prefix = stem.split("_count_", 1)[0]
    for profile in reaction_library.linkage_profiles.values():
        for role in profile.binary_bridge_roles:
            library_prefix = role.library_prefix or f"{role.motif_kind}s"
            if library_prefix == prefix:
                return role.motif_kind
    if prefix.endswith("s"):
        return prefix[:-1]
    return prefix


def _add_list_templates_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "list-templates",
        help="List registered reaction templates and their execution capabilities.",
    )
    parser.add_argument("--json", action="store_true", help="Emit JSON instead of text.")
    parser.set_defaults(func=_run_list_templates)


def _run_list_templates(args: argparse.Namespace) -> None:
    library = ReactionLibrary.builtin()
    rows = []
    for template_id, template in sorted(library.templates.items()):
        profile = library.linkage_profile(template_id)
        rows.append(
            {
                "template_id": template_id,
                "arity": template.arity,
                "reactant_motif_kinds": template.reactant_motif_kinds,
                "workflow_family": None if profile is None else profile.workflow_family,
                "supports_pair_generation": False if profile is None else profile.supports_binary_bridge_pair_generation,
                "supports_atomistic_realization": False if profile is None else profile.supports_atomistic_realization,
                "supports_topology_guided_generation": False if profile is None else profile.supports_topology_guided_generation,
                "topology_assignment_mode": None if profile is None else profile.topology_assignment_mode,
                "roles": () if profile is None else tuple(role.motif_kind for role in profile.binary_bridge_roles),
            }
        )
    if args.json:
        print(json.dumps(rows, indent=2, sort_keys=True))
        return
    for row in rows:
        print(
            row["template_id"],
            "arity=" + str(row["arity"]),
            "family=" + str(row["workflow_family"]),
            "pair_generation=" + str(row["supports_pair_generation"]).lower(),
            "atomistic=" + str(row["supports_atomistic_realization"]).lower(),
            "topology_guided=" + str(row["supports_topology_guided_generation"]).lower(),
            "roles=" + ",".join(row["roles"]),
        )
