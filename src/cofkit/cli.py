from __future__ import annotations

import argparse
import os
import signal
import sys
from pathlib import Path
from typing import Sequence

from ._version import __version__
from .cli_errors import COMMAND_FAILURE_EXIT_CODE
from .cli_analyze import add_analyze_group
from .cli_build import add_build_group
from .cli_calculate import add_calculate_group
from .cli_validate import add_validate_group

_LEGACY_COMMAND_ALIASES: dict[str, tuple[str, ...]] = {
    "single-pair": ("build", "single-pair"),
    "batch-binary-bridge": ("build", "batch-binary-bridge"),
    "batch-imine": ("build", "batch-binary-bridge", "--template-id", "imine_bridge"),
    "batch-all-binary-bridges": ("build", "batch-all-binary-bridges"),
    "classify-output": ("analyze", "classify-output"),
    "build-default-library": ("build", "default-library"),
    "list-templates": ("build", "list-templates"),
}


class _ArgumentParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # SUPPRESS prevents child parsers from clearing a parent's debug flag.
        self.add_argument("--debug", action="store_true", default=argparse.SUPPRESS,
                          help="Show the original traceback when a command raises an exception.")


def build_parser() -> argparse.ArgumentParser:
    parser = _ArgumentParser(
        prog="cofkit",
        description="CLI for COF building, calculation, screening, analysis, and validation.",
    )
    _set_help_default(parser)
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command")
    add_build_group(subparsers)
    add_analyze_group(subparsers)
    add_calculate_group(subparsers)
    add_validate_group(subparsers)
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    debug = "--debug" in raw_argv
    try:
        _load_dotenv()
        normalized_argv = _normalize_argv(raw_argv)
        args = build_parser().parse_args(normalized_argv)
        debug = getattr(args, "debug", False)
        args.func(args)
    except SystemExit as exc:
        # Existing handlers may wrap expected failures. Recover the original
        # exception for debugging and normalize string exits to stderr + code 1.
        if debug and exc.__cause__ is not None:
            raise exc.__cause__
        if isinstance(exc.code, str):
            message = exc.code.removeprefix("error: ")
            print(f"error: {message}", file=sys.stderr)
            raise SystemExit(COMMAND_FAILURE_EXIT_CODE) from None
        raise
    except KeyboardInterrupt:
        print("error: interrupted", file=sys.stderr)
        # Derived from the shell convention: 128 + the interrupt signal number.
        raise SystemExit(128 + signal.SIGINT) from None
    except Exception as exc:
        if debug:
            raise
        print(f"error: {type(exc).__name__}: {exc}", file=sys.stderr)
        if not isinstance(exc, (OSError, ValueError, RuntimeError)):
            print("hint: rerun with --debug for the original traceback", file=sys.stderr)
        raise SystemExit(COMMAND_FAILURE_EXIT_CODE) from None


def _normalize_argv(argv: Sequence[str]) -> list[str]:
    normalized = list(argv)
    global_flags = []
    while normalized and normalized[0] == "--debug":
        global_flags.append(normalized.pop(0))
    if not normalized:
        return global_flags

    replacement = _LEGACY_COMMAND_ALIASES.get(normalized[0])
    if replacement is None:
        return [*global_flags, *normalized]

    print(
        "warning: "
        f"`cofkit {normalized[0]}` is deprecated; use `cofkit {' '.join(replacement)}` instead.",
        file=sys.stderr,
    )
    return [*global_flags, *replacement, *normalized[1:]]


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
