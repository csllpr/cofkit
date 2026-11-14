"""Compare exhaustive/spatial distance bonding and complete decomposition.

Run: uv run python tests/benchmark_decompose_neighbors.py --output-dir out/decompose-neighbors
Inputs are generated precursor seeds, not experimental accuracy labels.
"""

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
from statistics import median
from time import perf_counter
from unittest.mock import patch

from cofkit import decompose
from cofkit.cofid import parse_cofid
from cofkit.decompose_cif import read_periodic_cif_atoms
from decompose_reference import exhaustive_bond_candidates
from soft_relax_cof_cases import COF_CASES, replicate_cif, write_cof_seed

# Heuristic measurement budget: three samples permit a median and alternate
# backend order to limit warm-cache bias without a long benchmark run.
TIMING_REPETITIONS = 3
SUPERCELLS = {"cof-5": ((2, 2, 1), (4, 4, 1)), "cof-300-single-net": ((2, 2, 2),)}


def compare(case, seed, directory):
    atoms = read_periodic_cif_atoms(seed)
    timings = {"exhaustive": [], "spatial": []}
    methods = {"exhaustive": exhaustive_bond_candidates, "spatial": decompose._distance_inferred_bond_candidates}
    bonds = {}
    for iteration in range(TIMING_REPETITIONS):
        for name in (tuple(methods) if iteration % 2 == 0 else tuple(reversed(methods))):
            started = perf_counter()
            bonds[name] = methods[name](atoms)
            timings[name].append(perf_counter() - started)
    assert bonds["spatial"] == bonds["exhaustive"], "ordered bond candidates differ"
    reports = {}
    full_timings = {}
    for name, method in methods.items():
        started = perf_counter()
        with patch.object(decompose, "_distance_inferred_bond_candidates", method):
            result = decompose.decompose_cif_to_cofid(
                seed, topology=parse_cofid(case["cofid"]).topology, bond_mode="distance",
            )
        full_timings[name] = perf_counter() - started
        reports[name] = result.to_dict()
        (directory / f"{name}.json").write_text(json.dumps(reports[name], indent=2))
    assert reports["spatial"] == reports["exhaustive"], "complete decomposition reports differ"
    (directory / "bonds.json").write_text(json.dumps([asdict(bond) for bond in bonds["spatial"]], indent=2))
    return {
        "status": "ok", "atoms": len(atoms), "bonds": len(bonds["spatial"]),
        "ordered_bonds_equal": True, "complete_reports_equal": True,
        "decomposition_status": reports["spatial"]["status"],
        "event_status": reports["spatial"]["metadata"].get("event_status"),
        "pair_seconds": {name: median(values) for name, values in timings.items()},
        "pair_samples": timings, "full_seconds": full_timings,
        "pair_speedup": median(timings["exhaustive"]) / median(timings["spatial"]),
        "full_speedup": full_timings["exhaustive"] / full_timings["spatial"],
        "seed_sha256": hashlib.sha256(seed.read_bytes()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("out/decompose-neighbors"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for case in COF_CASES:
        variants = [(case["id"], None)] + [
            (case["id"] + "__" + "x".join(map(str, repeat)), repeat)
            for repeat in SUPERCELLS.get(case["id"], ())
        ]
        base_seed = args.output_dir / case["id"] / "seed.cif"
        for identifier, repeat in variants:
            directory = args.output_dir / identifier
            directory.mkdir(parents=True, exist_ok=True)
            seed = directory / "seed.cif"
            try:
                if repeat is None:
                    write_cof_seed(case, seed)
                else:
                    replicate_cif(base_seed, seed, repeat)
                row = {"case": identifier, **compare(case, seed, directory)}
            except Exception as exc:
                row = {"case": identifier, "status": "error", "error": f"{type(exc).__name__}: {exc}"}
            rows.append(row)
            print(json.dumps(row), flush=True)
    provenance = {
        name: hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest()
        for name, module in (("decompose", decompose),)
    }
    provenance["exhaustive_reference"] = hashlib.sha256(Path(__file__).with_name("decompose_reference.py").read_bytes()).hexdigest()
    (args.output_dir / "results.json").write_text(json.dumps({
        "notes": "Generated seeds; known topology supplied. Equality includes unsuccessful verdicts; no experimental accuracy assertion.",
        "source_sha256": provenance, "cases": rows,
    }, indent=2))
    return int(any(row["status"] != "ok" for row in rows))


if __name__ == "__main__":
    raise SystemExit(main())
