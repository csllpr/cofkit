"""Side-by-side common-COF and supercell benchmark for spatial pair search.

Run with `uv run python tests/benchmark_soft_relax_neighbors.py --output-dir
out/soft-relax-neighbor-search/comparison`. Generated CIFs and reports stay
outside version control. Any case failure is recorded and the cohort continues.
"""

import argparse
import csv
from dataclasses import asdict
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
from statistics import median
import sys
from time import perf_counter
from unittest.mock import patch

import gemmi

from cofkit import soft_relax
from cofkit.cofid import cofid_to_build_request
from cofkit.validation import CoarseStructureValidator
from soft_relax_cof_cases import (
    COF_CASES, COF_TEST_CONFORMERS, CROWDED_LAYER_CASE_IDS, CROWDED_LAYER_GAP_ANGSTROM,
    normalized_report, replicate_cif, write_cof_seed, write_crowded_layers,
)
from soft_relax_reference import exhaustive_pair_list


# Heuristic — pending calibration: a median of three interleaved repetitions
# reduces timing noise without turning this validation into a long benchmark.
PAIR_TIMING_REPEATS = 3
SUPERCELLS = (
    ("cof-5", (2, 2, 1)),
    ("cof-5", (4, 4, 1)),
    ("cof-300-single-net", (2, 2, 2)),
)


def _pair_comparison(source, config):
    block = gemmi.cif.read_file(str(source)).sole_block()
    system, _ = soft_relax._parse_system(block, config)
    timings = {"exhaustive": [], "spatial": []}
    backends = {"exhaustive": exhaustive_pair_list, "spatial": soft_relax._build_pair_list}
    results = {}
    for repetition in range(PAIR_TIMING_REPEATS):
        order = ("exhaustive", "spatial") if repetition % 2 == 0 else ("spatial", "exhaustive")
        for label in order:
            start = perf_counter()
            pairs = backends[label](system, config)
            timings[label].append(perf_counter() - start)
            results[label] = pairs
        if results["exhaustive"] != results["spatial"]:
            raise AssertionError("ordered pair lists or periodic images differ")
    cart = [soft_relax._frac_to_cart(system.orth, f) for f in system.frac]
    energies_forces = {
        label: soft_relax._energy_and_forces(cart, [], pairs, config.repulsion_coefficients[-1], config)
        for label, pairs in results.items()
    }
    if energies_forces["exhaustive"] != energies_forces["spatial"]:
        raise AssertionError("pair energies or forces differ")
    radii = [config.hydrogen_radius if s == "H" else soft_relax._element_radius(s, config) for s in system.symbols]
    cutoff = 2 * max(radii, default=config.fallback_radius) * config.repulsion_scale + config.neighbor_margin
    n = len(system.frac)
    elapsed = {label: median(values) for label, values in timings.items()}
    return {
        "atoms": n, "pairs": len(results["spatial"]), "pair_lists_identical": True,
        "pair_energies_forces_identical": True,
        "exhaustive_candidates": n * (n + 1) // 2,
        "spatial_candidates": sum(1 for _ in soft_relax._nearby_atom_pairs(system, cutoff)),
        "pair_seconds": elapsed, "pair_samples_seconds": timings,
        "pair_speedup": elapsed["exhaustive"] / elapsed["spatial"],
    }


def _relax_comparison(source, directory, case, config):
    directory.mkdir(parents=True, exist_ok=True)
    outputs, reports, timings = {}, {}, {}
    for label in ("exhaustive", "spatial"):
        target = directory / f"{label}.cif"
        start = perf_counter()
        if label == "exhaustive":
            with patch.object(soft_relax, "_build_pair_list", exhaustive_pair_list):
                report = soft_relax.relax_cif_clashes(source, target, config)
        else:
            report = soft_relax.relax_cif_clashes(source, target, config)
        timings[label] = perf_counter() - start
        outputs[label] = target
        reports[label] = report
        (directory / f"{label}.json").write_text(json.dumps(asdict(report), indent=2) + "\n")
    if normalized_report(reports["exhaustive"]) != normalized_report(reports["spatial"]):
        raise AssertionError("relaxation diagnostics differ")
    if outputs["exhaustive"].read_bytes() != outputs["spatial"].read_bytes():
        raise AssertionError("relaxed CIF bytes differ")
    request = cofid_to_build_request(case["cofid"])
    validator = CoarseStructureValidator()
    validation = {
        label: validator._validate_cif(path, topology_id=request.topology, template_id=request.template_id)
        for label, path in outputs.items()
    }
    # Ignore the filename-only metric when comparing the final geometry scan.
    for metrics, _, _ in validation.values():
        metrics.pop("cif_path", None)
    if validation["exhaustive"] != validation["spatial"]:
        raise AssertionError("final geometry validation differs")
    report = reports["spatial"]
    return {
        "relax_seconds": timings, "relax_speedup": timings["exhaustive"] / timings["spatial"],
        "relax_reports_identical": True, "relaxed_cifs_identical": True,
        "final_validation_identical": True,
        "converged": report.converged,
        "clashes_before": report.clashes_before, "clashes_after": report.clashes_after,
        "max_bond_drift": report.max_bond_drift,
        "geometry_reasons": list(validation["spatial"][1]),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    config = soft_relax.SoftRelaxConfig()
    records = []
    for case in COF_CASES:
        directory = output / case["id"]
        directory.mkdir(exist_ok=True)
        record = {"id": case["id"], "name": case["name"], "cofid": case["cofid"], "notes": case["notes"], "status": "ok"}
        try:
            source = directory / "seed.cif"
            write_cof_seed(case, source)
            record["input_sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
            record.update(_pair_comparison(source, config))
            record.update(_relax_comparison(source, directory, case, config))
        except Exception as exc:
            record.update(status="failed", error=f"{type(exc).__name__}: {exc}")
            print(f"warning: benchmark {case['id']} failed ({record['error']}); continuing", file=sys.stderr)
        records.append(record)
        print(json.dumps(record), flush=True)
    for case_id in CROWDED_LAYER_CASE_IDS:
        case = next(case for case in COF_CASES if case["id"] == case_id)
        identifier = case_id + "__crowded_layers"
        directory = output / identifier
        directory.mkdir(exist_ok=True)
        record = {
            "id": identifier, "name": case["name"], "status": "ok",
            "artificial_layer_gap_angstrom": CROWDED_LAYER_GAP_ANGSTROM,
            "notes": "Intact layer copies deliberately placed too close; an adversarial relaxation input, not experimental geometry.",
        }
        try:
            source = directory / "seed.cif"
            write_crowded_layers(output / case_id / "seed.cif", source)
            record["input_sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
            record.update(_pair_comparison(source, config))
            record.update(_relax_comparison(source, directory, case, config))
            if record["clashes_before"] <= 0 or record["clashes_after"] >= record["clashes_before"]:
                raise AssertionError("crowded-layer fixture did not exercise clash relief")
        except Exception as exc:
            record.update(status="failed", error=f"{type(exc).__name__}: {exc}")
            print(f"warning: benchmark {identifier} failed ({record['error']}); continuing", file=sys.stderr)
        records.append(record)
        print(json.dumps(record), flush=True)
    for case_id, repeats in SUPERCELLS:
        case = next(case for case in COF_CASES if case["id"] == case_id)
        identifier = case_id + "__" + "x".join(str(n) for n in repeats)
        directory = output / identifier
        directory.mkdir(exist_ok=True)
        record = {"id": identifier, "name": case["name"], "repeats": repeats, "status": "ok"}
        try:
            source = directory / "seed.cif"
            replicate_cif(output / case_id / "seed.cif", source, repeats)
            record["input_sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
            record.update(_pair_comparison(source, config))
            record.update(_relax_comparison(source, directory, case, config))
        except Exception as exc:
            record.update(status="failed", error=f"{type(exc).__name__}: {exc}")
            print(f"warning: benchmark {identifier} failed ({record['error']}); continuing", file=sys.stderr)
        records.append(record)
        print(json.dumps(record), flush=True)
    summary = {
        "scope": "Generated common-COF precursor seeds and exact diagonal supercells; pair/force/full-relaxation equivalence, not physical accuracy. Crowded layers are intentional adversarial inputs.",
        "environment": {"python": sys.version, **{name: version(name) for name in ("cofkit", "gemmi", "rdkit")}},
        "configuration": asdict(config), "seed_conformers": COF_TEST_CONFORMERS,
        "spatial_source_sha256": hashlib.sha256(Path(soft_relax.__file__).read_bytes()).hexdigest(),
        "reference_sha256": hashlib.sha256(Path(__file__).with_name("soft_relax_reference.py").read_bytes()).hexdigest(),
        "records": records,
    }
    (output / "results.json").write_text(json.dumps(summary, indent=2) + "\n")
    fields = ("id", "atoms", "pairs", "status", "pair_speedup", "relax_speedup", "converged", "clashes_before", "clashes_after", "max_bond_drift", "error")
    with (output / "results.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)
    return int(any(record["status"] != "ok" for record in records))


if __name__ == "__main__":
    raise SystemExit(main())
