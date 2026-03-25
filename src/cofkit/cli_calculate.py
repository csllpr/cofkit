from __future__ import annotations

import argparse


def add_calculate_group(subparsers) -> None:
    parser = subparsers.add_parser(
        "calculate",
        help="Run external calculation-tool workflows.",
        description="External calculation-tool integrations are planned but not implemented yet.",
    )
    parser.set_defaults(func=lambda _: _show_calculate_help(parser))


def _show_calculate_help(parser: argparse.ArgumentParser) -> None:
    parser.print_help()
    print()
    print("No calculation workflows are implemented yet. This namespace is reserved for future external-tool integrations.")
