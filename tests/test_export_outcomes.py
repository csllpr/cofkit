"""Delivery outcomes reconcile artifacts without implying chemical validity."""

import json
from collections import Counter
from dataclasses import replace
from pathlib import Path

import pytest

from cofkit.batch import BatchGenerationConfig, BatchStructureGenerator
from cofkit.batch_models import BatchMonomerRecord
from cofkit.build_workflows.ring_forming import RingFormingStructureGenerator
from cofkit.chem.rdkit import build_rdkit_monomer
from cofkit.cif import CIFWriter
from cofkit.cli import main


AMINE = "Nc1cc(N)cc(N)c1"
ALDEHYDE = "O=Cc1cc(C=O)cc(C=O)c1"


@pytest.fixture
def records():
    return (
        BatchMonomerRecord("a", "a", AMINE, "amine", 3),
        BatchMonomerRecord("b", "b", ALDEHYDE, "aldehyde", 3),
    )


def _generator(**overrides):
    return BatchStructureGenerator(
        BatchGenerationConfig(
            **{
                "rdkit_num_conformers": 1,
                "shape_aware_conformer": False,
                "topology_ids": ("hcb",),
                "max_workers": 1,
                **overrides,
            }
        )
    )


@pytest.mark.parametrize("api", ["records", "monomers", "multiple"])
@pytest.mark.parametrize(
    "scenario", ["written", "blocked", "disabled", "quota", "failed", "no_directory"]
)
def test_binary_export_outcomes_match_delivery(
    tmp_path, monkeypatch, records, api, scenario
):
    options = {}
    if scenario == "blocked":
        options["hard_hard_max_bridge_distance"] = 0.1
    if scenario == "quota":
        options["max_cif_exports"] = 0
    generator = _generator(**options)
    if scenario == "failed":

        def fail(*args, **kwargs):
            raise RuntimeError("synthetic export failure")

        monkeypatch.setattr(generator.cif_writer, "write_candidate", fail)
    kwargs = {
        "out_dir": None if scenario == "no_directory" else tmp_path,
        "write_cif": scenario != "disabled",
    }
    if api == "records":
        summary, _ = generator.generate_pair_candidate(*records, **kwargs)
    else:
        monomers = tuple(generator.build_monomer(r).monomer for r in records)
        if api == "monomers":
            summary, _ = generator.generate_monomer_pair_candidate(*monomers, **kwargs)
        else:
            summaries, _, _ = generator.generate_monomer_pair_candidates(
                *monomers, **kwargs
            )
            assert len(summaries) == 1
            summary = summaries[0]
    expected = {
        "disabled": "not_requested",
        "quota": "limit_reached",
        "no_directory": "not_attempted",
    }.get(scenario, scenario)
    assert summary.cif_export.status == expected
    assert (summary.cif_path is not None) == (expected == "written")
    if summary.cif_path:
        assert Path(summary.cif_path).is_file()
    if scenario == "blocked":
        assert summary.status == "ok"  # legacy construction semantics survive
        assert summary.validation_classification == "hard_hard_invalid"
        assert summary.cif_export.reasons == (
            "bridge_distance_exceeds_cif_export_limit",
        )
    if scenario == "failed":
        assert summary.status == "export-failed"
        assert summary.cif_export.error == "RuntimeError: synthetic export failure"
    serialized = generator._summary_to_dict(summary)
    assert serialized["cif_export"]["status"] == expected


def test_disabled_export_takes_precedence_over_geometry_gate(tmp_path, records):
    summary, _ = _generator(hard_hard_max_bridge_distance=0.1).generate_pair_candidate(
        *records,
        out_dir=tmp_path,
        write_cif=False,
    )
    assert summary.cif_export.status == "not_requested"
    assert summary.metadata["cif_export_blocked"]
    assert not list(tmp_path.rglob("*.cif"))


@pytest.mark.parametrize("classification", ["hard_hard_invalid", "hard_invalid"])
def test_final_validation_distinguishes_blocked_and_diagnostic_exports(
    tmp_path, monkeypatch, records, classification
):
    generator = _generator()
    validate = generator.structure_validator.validate_manifest_record

    def altered_report(*args, **kwargs):
        report = validate(*args, **kwargs)
        return replace(
            report,
            classification=classification,
            blocks_cif_export=classification == "hard_hard_invalid",
            reasons=("synthetic_final_coordinate_rejection",),
        )

    monkeypatch.setattr(
        generator.structure_validator, "validate_manifest_record", altered_report
    )
    summary, _ = generator.generate_pair_candidate(
        *records, out_dir=tmp_path, write_cif=True
    )
    assert summary.status == "ok"
    assert summary.validation_classification == classification
    outcome = summary.cif_export
    assert outcome.reasons == ("synthetic_final_coordinate_rejection",)
    if classification == "hard_hard_invalid":
        assert outcome.status == "blocked"
        assert outcome.error is None
        assert summary.cif_path is None
        assert not list(tmp_path.rglob("*.cif"))
    else:
        assert outcome.status == "written"
        assert outcome.diagnostic
        assert Path(summary.cif_path).is_file()


def test_preconstruction_failure_does_not_claim_export(tmp_path, records, monkeypatch):
    generator = _generator()

    def fail(*args, **kwargs):
        raise ValueError("synthetic monomer failure")

    monkeypatch.setattr(generator, "smiles_monomer_builder", fail)
    summary, _ = generator.generate_pair_candidate(*records, out_dir=tmp_path)
    assert summary.status == "monomer-build-failed"
    assert summary.cif_export.status == "not_attempted"


@pytest.mark.parametrize("workers", [1, 2])
@pytest.mark.parametrize("quota", [None, 1])
def test_batch_export_counts_cover_every_manifest_row(tmp_path, workers, quota):
    library = tmp_path / "library"
    library.mkdir()
    (library / "amines_count_3.txt").write_text(f"smiles\n{AMINE}\n{AMINE}\n")
    (library / "aldehydes_count_3.txt").write_text(f"smiles\n{ALDEHYDE}\n{ALDEHYDE}\n")
    out = tmp_path / "build"
    summary = _generator(
        max_workers=workers, retain_top_results=1, max_cif_exports=quota
    ).run_binary_bridge_batch(library, out)
    rows = [
        json.loads(line) for line in (out / "manifest.jsonl").read_text().splitlines()
    ]
    assert len(rows) == 4
    assert len(summary.top_results) == 1
    assert summary.export_counts == Counter(r["cif_export"]["status"] for r in rows)
    expected_counts = (
        {"written": 4} if quota is None else {"written": 1, "limit_reached": 3}
    )
    assert summary.export_counts == expected_counts
    assert summary.cifs_written == (4 if quota is None else 1)
    assert summary.diagnostic_cifs == sum(r["cif_export"]["diagnostic"] for r in rows)
    assert (
        "CIF export outcomes (all manifest records)" in (out / "summary.md").read_text()
    )


@pytest.mark.parametrize("write", [False, True])
def test_single_pair_cli_reports_delivery_counts(tmp_path, write):
    main(
        [
            "build",
            "single-pair",
            "--cofid",
            f"3:amine:{AMINE}.3:aldehyde:{ALDEHYDE}&&hcb&&imine",
            "--num-conformers",
            "1",
            "--no-shape-aware-conformer",
            "--output-dir",
            str(tmp_path),
            "--write-cif" if write else "--no-write-cif",
        ]
    )
    report = json.loads((tmp_path / "summary.json").read_text())
    expected = "written" if write else "not_requested"
    assert report["export_counts"] == {expected: 1}
    assert report["results"][0]["cif_export"]["status"] == expected
    assert report["exported_structures"] == report["cifs_written"] == int(write)
    assert report["successful_structures"] == report["constructed_structures"] == 1


@pytest.fixture
def ring(monkeypatch):
    monomer = build_rdkit_monomer(
        "ring", "ring", "N#Cc1ccc(C#N)cc1", "nitrile", num_conformers=1
    )
    candidate = RingFormingStructureGenerator().generate(
        monomer, "triazine_trimerization"
    )
    monkeypatch.setattr(
        "cofkit.cli_build.build_rdkit_monomer", lambda *a, **kw: monomer
    )
    return candidate


def _run_ring(out, *extra):
    main(
        [
            "build",
            "ring-forming",
            "--cofid",
            "2:nitrile:N#Cc1ccc(C#N)cc1&&hcb&&triazine",
            "--output-dir",
            str(out),
            *extra,
        ]
    )
    return json.loads((out / "summary.json").read_text())


@pytest.mark.parametrize("write", [False, True])
def test_rejected_ring_delivery_is_explicit(tmp_path, monkeypatch, capsys, ring, write):
    ring = replace(
        ring,
        flags=(*ring.flags, "ring_geometry_rejected"),
        metadata={
            **ring.metadata,
            "ring_validation": {
                **ring.metadata["ring_validation"],
                "classification": "rejected",
                "arrangement_classification": "rejected",
                "reasons": ["displaced_participant"],
            },
        },
    )
    monkeypatch.setattr(
        RingFormingStructureGenerator, "generate_candidates", lambda *a, **kw: (ring,)
    )
    report = _run_ring(tmp_path, "--write-cif" if write else "--no-write-cif")
    row = report["results"][0]
    assert row["status"] == "ok"
    assert row["ring_validation"]["classification"] == "rejected"
    assert row["cif_export"]["status"] == ("written" if write else "not_requested")
    assert row["cif_export"]["diagnostic"] is write
    assert report["cifs_written"] == report["exported_structures"] == int(write)
    assert report["diagnostic_cifs"] == int(write)
    stderr = capsys.readouterr().err
    if write:
        assert "diagnostic CIF" in Path(row["cif_path"]).read_text()
        assert row["cif_export"]["reasons"] == ["displaced_participant"]
    else:
        assert "CIF export was not requested" in stderr
        assert "was written" not in stderr
        assert not list(tmp_path.glob("*.cif"))


def test_ring_export_failure_isolated_and_previous_artifact_preserved(
    tmp_path, monkeypatch, ring
):
    candidates = (replace(ring, id="first"), replace(ring, id="second"))
    monkeypatch.setattr(
        RingFormingStructureGenerator,
        "generate_candidates",
        lambda *a, **kw: candidates,
    )
    previous = tmp_path / "first.cif"
    previous.write_text("previous completed artifact")
    original = CIFWriter.write_candidate

    def write(self, path, *args, **kwargs):
        if Path(path).name == "first.cif":
            Path(path).write_text("partial candidate")
            raise OSError("synthetic disk error")
        return original(self, path, *args, **kwargs)

    monkeypatch.setattr(CIFWriter, "write_candidate", write)
    report = _run_ring(tmp_path, "--write-cif")
    first, second = report["results"]
    assert first["status"] == "export-failed"
    assert first["cif_path"] is None
    assert first["cif_export"]["error"] == "OSError: synthetic disk error"
    assert second["cif_export"]["status"] == "written"
    assert Path(second["cif_path"]).is_file()
    assert report["export_counts"] == {"failed": 1, "written": 1}
    assert previous.read_text() == "previous completed artifact"
    assert not list(tmp_path.glob(".cif-*"))
