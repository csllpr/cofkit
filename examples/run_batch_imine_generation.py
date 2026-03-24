from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit.cli import main as cofkit_cli_main


def main() -> None:
    cofkit_cli_main(["batch-imine", *sys.argv[1:]])


if __name__ == "__main__":
    main()
