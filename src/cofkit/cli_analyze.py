from __future__ import annotations

import argparse
import json
from pathlib import Path

from .validation import CoarseValidationThresholds, classify_batch_output
from .zeopp import ZeoppError, analyze_zeopp_pore_properties


def add_analyze_group(subparsers) -> None:
    parser = subparsers.add_parser(
        "analyze",
        help="Analyze and classify existing structures and output trees.",
        description="Analyze existing COF structures and batch outputs.",
    )
    _set_help_default(parser)

    analyze_subparsers = parser.add_subparsers(dest="analyze_command")
    _add_classify_output_parser(analyze_subparsers)
    _add_zeopp_parser(analyze_subparsers)


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


def _add_zeopp_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "zeopp",
        help="Run Zeo++ pore-property analysis on one CIF file.",
        description=(
            "Run Zeo++ on one CIF file and extract basic pore properties, point-probe "
            "surface/volume baselines, and optional accessibility-aware probe scans."
        ),
    )
    parser.add_argument("cif_path", help="Input CIF file to analyze with Zeo++.")
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory for raw Zeo++ outputs and the JSON report. Defaults to <input_stem>_zeopp beside the CIF.",
    )
    parser.add_argument(
        "--probe-radius",
        type=float,
        action="append",
        default=None,
        help=(
            "Optional probe radius in angstrom for accessibility-aware scans. Repeat to run more than one scan. "
            "When omitted, cofkit writes only the point-probe baseline."
        ),
    )
    parser.add_argument(
        "--channel-radius",
        type=float,
        default=None,
        help="Optional shared channel-radius override for every requested probe scan. Defaults to each probe radius.",
    )
    parser.add_argument(
        "--surface-samples-per-atom",
        type=int,
        default=250,
        help="Monte Carlo samples per atom for Zeo++ surface-area runs. Default: 250.",
    )
    parser.add_argument(
        "--volume-samples-total",
        type=int,
        default=5000,
        help="Monte Carlo total sample count for Zeo++ volume runs. Default: 5000.",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=float,
        default=300.0,
        help="Per-command timeout for Zeo++ subprocess calls. Default: 300.",
    )
    parser.add_argument(
        "--zeopp-path",
        default=None,
        help="Optional explicit path to the Zeo++ network binary. Defaults to COFKIT_ZEOPP_PATH.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print the full Zeo++ report as JSON instead of a short human-readable summary.",
    )
    parser.set_defaults(func=_run_zeopp)


def _run_zeopp(args: argparse.Namespace) -> None:
    try:
        result = analyze_zeopp_pore_properties(
            args.cif_path,
            output_dir=args.output_dir,
            probe_radii=tuple(args.probe_radius or ()),
            channel_radius=args.channel_radius,
            surface_samples_per_atom=args.surface_samples_per_atom,
            volume_samples_total=args.volume_samples_total,
            zeopp_path=args.zeopp_path,
            timeout_seconds=args.timeout_seconds,
        )
    except (ZeoppError, FileNotFoundError, ValueError) as exc:
        raise SystemExit(f"error: {exc}") from exc

    if args.json:
        print(json.dumps(result.to_dict(), indent=2))
        return

    properties = result.baseline.basic_pore_properties
    point_probe_channels = result.baseline.point_probe_channels
    point_probe_surface_area = result.baseline.point_probe_surface_area
    point_probe_volume = result.baseline.point_probe_volume
    successful_probe_scans = sum(scan.status == "ok" for scan in result.probe_scans)
    print("input_cif:", result.input_cif)
    print("zeopp_binary:", result.zeopp_binary)
    print("output_dir:", result.output_dir)
    print("largest_included_sphere:", properties.largest_included_sphere)
    print("largest_free_sphere:", properties.largest_free_sphere)
    print("largest_included_sphere_along_free_path:", properties.largest_included_sphere_along_free_path)
    print("axis_aligned_free_sphere:", dict(properties.axis_aligned_free_sphere))
    print(
        "axis_aligned_included_sphere_along_free_path:",
        dict(properties.axis_aligned_included_sphere_along_free_path),
    )
    print("point_probe_n_channels:", point_probe_channels.n_channels)
    print("point_probe_n_pockets:", point_probe_channels.n_pockets)
    print("point_probe_channel_dimensionality:", point_probe_channels.channel_dimensionality)
    print("point_probe_accessible_surface_area_a2:", point_probe_surface_area.accessible_surface_area_a2)
    print("point_probe_accessible_volume_a3:", point_probe_volume.accessible_volume_a3)
    print("probe_scans_requested:", len(result.probe_scans))
    print("probe_scans_successful:", successful_probe_scans)
    print("report_path:", result.report_path)
