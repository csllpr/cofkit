from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
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
    _load_dotenv()
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


def _load_dotenv(start_dir: Path | None = None) -> Path | None:
    dotenv_path = _find_dotenv(start_dir)
    if dotenv_path is None:
        return None

    for raw_line in dotenv_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("export "):
            line = line[len("export ") :].lstrip()
        if "=" not in line:
            continue
        key, raw_value = line.split("=", 1)
        key = key.strip()
        if not key or any(character.isspace() for character in key):
            continue
        os.environ.setdefault(key, _parse_dotenv_value(raw_value.strip()))

    return dotenv_path


def _find_dotenv(start_dir: Path | None = None) -> Path | None:
    current = (start_dir or Path.cwd()).expanduser().resolve()
    for directory in (current, *current.parents):
        candidate = directory / ".env"
        if candidate.is_file():
            return candidate
    return None


def _parse_dotenv_value(raw_value: str) -> str:
    if len(raw_value) >= 2 and raw_value[0] == raw_value[-1] and raw_value[0] in {"'", '"'}:
        return raw_value[1:-1]
    return raw_value


if __name__ == "__main__":
    main()
