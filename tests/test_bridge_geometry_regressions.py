"""Seven independent inputs from the bridge-geometry benchmark report.

The other six labels reuse case 35's COFid. Only input identities are stored;
no benchmark output or external simulation installation is required.
"""

import json
from pathlib import Path

import gemmi
import pytest

from cofkit.cli import main
from cofkit.validation import (
    CoarseStructureValidator,
    measure_inter_instance_bond_distances,
)

CASES = json.loads(
    (Path(__file__).parent / "fixtures/bridge_geometry/cofids.json").read_text()
)


@pytest.mark.parametrize("case", CASES, ids=lambda row: row["case"])
def test_reported_bridge_geometry(case, tmp_path, capsys):
    assert main(
        [
            "build",
            "single-pair",
            "--cofid",
            case["cofid"],
            "--output-dir",
            str(tmp_path),
        ]
    ) in (None, 0)
    summary = json.loads((tmp_path / "summary.json").read_text())
    assert len(summary["results"]) == 1
    record = summary["results"][0]
    assert record["status"] == "ok", record
    path = Path(record["cif_path"])
    metrics, reasons, _ = CoarseStructureValidator()._validate_cif(
        path, topology_id="hcb"
    )
    assert "heavy_atom_clash" not in reasons, metrics.get("heavy_atom_clash_min")
    assert "excluded_pair_severe_overlap" not in reasons, metrics.get(
        "excluded_pair_severe_overlap_min"
    )
    block = gemmi.cif.read_file(str(path)).sole_block()
    bonds = measure_inter_instance_bond_distances(block)
    assert bonds
    # These reaction families form C=N/C=C bonds, whose template target is
    # 1.3 A. This broad seed window rejects stretched/disconnected bridges.
    assert all(1.15 < bond.distance < 1.6 for bond in bonds), bonds
    assert record["metadata"]["validation"]["classification"] not in (
        "hard_invalid",
        "hard_hard_invalid",
    )
