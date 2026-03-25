from __future__ import annotations

import argparse
from pathlib import Path

from .validation import CoarseValidationThresholds, classify_batch_output


def add_analyze_group(subparsers) -> None:
    parser = subparsers.add_parser(
        "analyze",
        help="Analyze and classify existing structures and output trees.",
        description="Analyze existing COF structures and batch outputs.",
    )
    _set_help_default(parser)

    analyze_subparsers = parser.add_subparsers(dest="analyze_command")
    _add_classify_output_parser(analyze_subparsers)


def _set_help_default(parser: argparse.ArgumentParser) -> None:
    def _show_help(_: argparse.Namespace) -> None:
        parser.print_help()

    parser.set_defaults(func=_show_help)


def _add_classify_output_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "classify-output",
        help="Classify a batch output into valid, warning, hard-invalid, and hard-hard-invalid sets.",
    )
    parser.add_argument("source_dir", help="Batch output directory containing manifest.jsonl and cifs/.")
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory for the classified output. Defaults to <source_dir>_coarse_validation.",
    )
    parser.add_argument("--link-mode", choices=("symlink", "hardlink", "copy"), default="symlink")
    parser.add_argument("--max-workers", type=int, default=None)
    parser.add_argument("--max-structures", type=int, default=None)
    parser.add_argument("--warning-max-bridge-distance-residual", type=float, default=0.75)
    parser.add_argument("--warning-mean-bridge-distance-residual", type=float, default=0.35)
    parser.add_argument("--warning-bad-bridge-distance-residual", type=float, default=0.50)
    parser.add_argument("--warning-max-bad-bridge-fraction", type=float, default=0.25)
    parser.add_argument("--hard-hard-max-bridge-distance", type=float, default=2.5)
    parser.add_argument("--hard-max-bridge-distance-residual", type=float, default=1.0)
    parser.add_argument("--hard-mean-bridge-distance-residual", type=float, default=0.60)
    parser.add_argument("--hard-min-bridge-distance-ratio", type=float, default=0.70)
    parser.add_argument("--hard-max-bridge-distance-ratio", type=float, default=1.60)
    parser.add_argument("--min-nonbonded-heavy-distance", type=float, default=1.05)
    parser.add_argument("--min-2d-cell-area", type=float, default=10.0)
    parser.add_argument("--min-3d-cell-volume", type=float, default=20.0)
    parser.add_argument("--no-skip-cif-checks-when-metadata-invalid", action="store_true")
    parser.set_defaults(func=_run_classify_output)


def _run_classify_output(args: argparse.Namespace) -> None:
    source_dir = Path(args.source_dir)
    output_dir = Path(args.output_dir) if args.output_dir else Path(f"{source_dir}_coarse_validation")
    thresholds = CoarseValidationThresholds(
        warning_max_bridge_distance_residual=args.warning_max_bridge_distance_residual,
        warning_mean_bridge_distance_residual=args.warning_mean_bridge_distance_residual,
        warning_bad_bridge_distance_residual=args.warning_bad_bridge_distance_residual,
        warning_max_bad_bridge_fraction=args.warning_max_bad_bridge_fraction,
        hard_hard_max_bridge_distance=args.hard_hard_max_bridge_distance,
        hard_max_bridge_distance_residual=args.hard_max_bridge_distance_residual,
        hard_mean_bridge_distance_residual=args.hard_mean_bridge_distance_residual,
        hard_min_bridge_distance_ratio=args.hard_min_bridge_distance_ratio,
        hard_max_bridge_distance_ratio=args.hard_max_bridge_distance_ratio,
        min_nonbonded_heavy_distance=args.min_nonbonded_heavy_distance,
        min_2d_cell_area=args.min_2d_cell_area,
        min_3d_cell_volume=args.min_3d_cell_volume,
        skip_cif_checks_when_metadata_invalid=not args.no_skip_cif_checks_when_metadata_invalid,
    )
    summary = classify_batch_output(
        source_dir,
        output_dir,
        thresholds=thresholds,
        link_mode=args.link_mode,
        max_workers=args.max_workers,
        max_structures=args.max_structures,
    )
    print("source_dir:", summary.source_dir)
    print("output_dir:", summary.output_dir)
    print("total_structures:", summary.total_structures)
    print("valid_structures:", summary.valid_structures)
    print("warning_structures:", summary.warning_structures)
    print("hard_hard_invalid_structures:", summary.hard_hard_invalid_structures)
    print("hard_invalid_structures:", summary.hard_invalid_structures)
    print("warning_reason_counts:", dict(summary.warning_reason_counts))
    print("hard_hard_invalid_reason_counts:", dict(summary.hard_hard_invalid_reason_counts))
    print("hard_invalid_reason_counts:", dict(summary.hard_invalid_reason_counts))
    print("classification_manifest:", summary.classification_manifest_path)
