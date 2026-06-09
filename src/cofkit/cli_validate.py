from __future__ import annotations

import argparse
import json

from .graspa import EqeqExecutionError, GraspaConfigurationError
from .lammps import (
    LammpsConfigurationError,
    LammpsError,
    LammpsExecutionError,
    LammpsInputError,
    LammpsParseError,
)
from .validate import (
    COFidValidationResult,
    validate_cif_against_cofid,
    validate_lammps_optimized_cif_against_cofid,
)


def add_validate_group(subparsers) -> None:
    parser = subparsers.add_parser(
        "validate",
        help="Validate CIF structures against supplied COFid strings.",
        description="Validate one CIF against one expected COFid.",
    )
    _set_help_default(parser)

    validate_subparsers = parser.add_subparsers(dest="validate_command")
    _add_simple_parser(validate_subparsers)
    _add_optimize_parser(validate_subparsers)


def _set_help_default(parser: argparse.ArgumentParser) -> None:
    def _show_help(_: argparse.Namespace) -> None:
        parser.print_help()

    parser.set_defaults(func=_show_help)


def _add_simple_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "simple",
        help="Decompose a CIF with distance-inferred bonds and compare it to a COFid.",
    )
    parser.add_argument("cofid", help="Expected COFid. Topology is used as a placeholder but not compared.")
    parser.add_argument("cif_path", help="Input CIF file to validate.")
    parser.add_argument("--json", action="store_true", help="Print the validation result as JSON.")
    parser.set_defaults(func=_run_simple)


def _add_optimize_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "optimize",
        help="Run default LAMMPS optimization, then decompose and compare the optimized CIF to a COFid.",
    )
    parser.add_argument("cofid", help="Expected COFid. Topology is used as a placeholder but not compared.")
    parser.add_argument("cif_path", help="Input CIF file to optimize and validate.")
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory for LAMMPS artifacts. Defaults to <input_stem>_lammps beside the CIF.",
    )
    parser.add_argument(
        "--lmp-path",
        default=None,
        help="Optional explicit path to the LAMMPS executable or alias. Defaults to COFKIT_LMP_PATH.",
    )
    parser.add_argument(
        "--eqeq-path",
        default=None,
        help="Optional explicit path to the EQeq executable used by the default LAMMPS charge stage.",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=float,
        default=300.0,
        help="Per-run timeout for the LAMMPS subprocess. Default: 300.",
    )
    parser.add_argument(
        "--eqeq-timeout-seconds",
        type=float,
        default=300.0,
        help="Per-run timeout for the default EQeq charge-assignment subprocess. Default: 300.",
    )
    parser.add_argument("--json", action="store_true", help="Print the validation result as JSON.")
    parser.set_defaults(func=_run_optimize)


def _run_simple(args: argparse.Namespace) -> None:
    try:
        result = validate_cif_against_cofid(args.cofid, args.cif_path)
    except (FileNotFoundError, ValueError) as exc:
        raise SystemExit(f"error: {exc}") from exc
    _emit_result(result, json_output=args.json)


def _run_optimize(args: argparse.Namespace) -> None:
    try:
        result = validate_lammps_optimized_cif_against_cofid(
            args.cofid,
            args.cif_path,
            output_dir=args.output_dir,
            lmp_path=args.lmp_path,
            eqeq_path=args.eqeq_path,
            timeout_seconds=args.timeout_seconds,
            eqeq_timeout_seconds=args.eqeq_timeout_seconds,
        )
    except (
        FileNotFoundError,
        ValueError,
        GraspaConfigurationError,
        EqeqExecutionError,
        LammpsConfigurationError,
        LammpsInputError,
        LammpsExecutionError,
        LammpsParseError,
        LammpsError,
    ) as exc:
        raise SystemExit(f"error: {exc}") from exc
    _emit_result(result, json_output=args.json)


def _emit_result(result: COFidValidationResult, *, json_output: bool) -> None:
    if json_output:
        print(json.dumps(result.to_dict(), indent=2, sort_keys=True))
    else:
        _print_human_result(result)
    if not result.ok:
        raise SystemExit(2 if result.status == "failed" else 1)


def _print_human_result(result: COFidValidationResult) -> None:
    decompose_metadata = dict(result.decomposition.metadata) if result.decomposition is not None else {}
    print("status:", result.status)
    print("input_cif:", result.input_cif)
    print("checked_cif:", result.checked_cif)
    if result.optimized_cif is not None:
        print("optimized_cif:", result.optimized_cif)
    print("expected_cofid:", result.input_cofid)
    print("recovered_cofid:", result.recovered_cofid or "")
    print("topology_compared:", str(result.topology_compared).lower())
    print("monomers_match:", str(result.monomers_match).lower())
    print("linkage_match:", str(result.linkage_match).lower())
    print("bond_source:", decompose_metadata.get("bond_source", ""))
    print("n_distance_inferred_bonds:", decompose_metadata.get("n_distance_inferred_bonds", ""))
    print("n_bond_orders_inferred:", decompose_metadata.get("n_bond_orders_inferred", ""))
    if result.reason:
        print("reason:", result.reason)


__all__ = ["add_validate_group"]
