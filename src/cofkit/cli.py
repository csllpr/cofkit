from __future__ import annotations

import argparse
import sys
from typing import Sequence

from ._version import __version__
from .cli_analyze import add_analyze_group
from .cli_build import add_build_group
from .cli_calculate import add_calculate_group

_LEGACY_COMMAND_ALIASES: dict[str, tuple[str, ...]] = {
    "single-pair": ("build", "single-pair"),
    "batch-binary-bridge": ("build", "batch-binary-bridge"),
    "batch-imine": ("build", "batch-binary-bridge", "--template-id", "imine_bridge"),
    "batch-all-binary-bridges": ("build", "batch-all-binary-bridges"),
    "classify-output": ("analyze", "classify-output"),
    "build-default-library": ("build", "default-library"),
    "list-templates": ("build", "list-templates"),
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cofkit",
        description="CLI for COF building, calculation, screening, analysis, and validation.",
    )
    _set_help_default(parser)
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command")
    add_build_group(subparsers)
    add_analyze_group(subparsers)
    add_calculate_group(subparsers)
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    normalized_argv = _normalize_argv(sys.argv[1:] if argv is None else argv)
    args = build_parser().parse_args(normalized_argv)
    args.func(args)


def _normalize_argv(argv: Sequence[str]) -> list[str]:
    normalized = list(argv)
    if not normalized:
        return normalized

    replacement = _LEGACY_COMMAND_ALIASES.get(normalized[0])
    if replacement is None:
        return normalized

    print(
        "warning: "
        f"`cofkit {normalized[0]}` is deprecated; use `cofkit {' '.join(replacement)}` instead.",
        file=sys.stderr,
    )
    return [*replacement, *normalized[1:]]


def _set_help_default(parser: argparse.ArgumentParser) -> None:
    def _show_help(_: argparse.Namespace) -> None:
        parser.print_help()

    parser.set_defaults(func=_show_help)


if __name__ == "__main__":
    main()
