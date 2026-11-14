"""A damaged record must not hide the rest of a classification run."""

import json

import pytest

from cofkit import validation


def make_input(tmp_path, records):
    source = tmp_path / "source"
    source.mkdir()
    (source / "manifest.jsonl").write_text("\n".join(records) + "\n", encoding="utf-8")
    return source


def row(identifier, **fields):
    return json.dumps({"status": "ok", "structure_id": identifier, **fields})


def read_rows(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def test_bad_json_shape_and_metadata_do_not_abort_later_records(tmp_path, capsys):
    source = make_input(tmp_path, [
        "", '{"truncated":', "[]", row("broken", metadata={"score_metadata": {"bridge_event_metrics": 1}}),
        row("after"),
    ])
    output = tmp_path / "classified"
    result = validation.classify_batch_output(source, output, max_workers=2)
    assert result.total_structures == 4
    assert result.failed_records == 3
    assert result.unvalidated_structures == 3
    assert result.hard_invalid_structures == 1  # the final record has no CIF
    assert [failure.line_number for failure in result.failures] == [2, 3, 4]
    assert [failure.stage for failure in result.failures] == ["manifest_read", "manifest_read", "validation"]
    failures = read_rows(output / "failures.jsonl")
    assert failures[0]["raw_line"] == '{"truncated":'
    assert failures[2]["structure_id"] == "broken"
    assert failures[2]["error"] == "TypeError: 'int' object is not iterable"
    records = read_rows(output / "classification_manifest.jsonl")
    assert {r["structure_id"] for r in records if r.get("status") == "ok"} == {"after"}
    failed = [r for r in records if r["status"] == "classification_failed"]
    assert all(r["validation"]["classification"] == "unvalidated" for r in failed)
    assert all(r["validation"]["is_valid"] is None for r in failed)
    assert "TypeError: 'int' object is not iterable" in (output / "summary.md").read_text()
    assert "classification line 4" in capsys.readouterr().err


@pytest.mark.parametrize("exception", [OSError("worker failure"), RuntimeError("unexpected")])
def test_worker_exception_keeps_all_later_results(tmp_path, monkeypatch, exception):
    source = make_input(tmp_path, [row(str(i)) for i in range(12)])
    worker = validation._validate_record_worker

    def fail_one(record, *args):
        if record["structure_id"] == "0":
            raise exception
        return worker(record, *args)

    monkeypatch.setattr(validation, "_validate_record_worker", fail_one)
    result = validation.classify_batch_output(source, tmp_path / "output", max_workers=2)
    assert result.total_structures == 12
    assert result.failed_records == result.unvalidated_structures == 1
    assert result.hard_invalid_structures == 11
    assert result.failures[0].error == f"{type(exception).__name__}: {exception}"


def test_artifact_delivery_exception_is_a_record_failure(tmp_path, monkeypatch):
    cif = tmp_path / "model.cif"
    cif.write_text("model")
    source = make_input(tmp_path, [row("bad-copy"), row("after")])

    def validate(record, *args):
        return validation.CoarseValidationReport(
            classification="valid", is_valid=True, passes_hard_validation=True,
            metrics={"cif_path": str(cif)} if record["structure_id"] == "bad-copy" else {},
        )

    def fail_copy(*args):
        raise PermissionError("destination is read-only")

    monkeypatch.setattr(validation, "_validate_record_worker", validate)
    monkeypatch.setattr(validation, "_materialize_link", fail_copy)
    result = validation.classify_batch_output(source, tmp_path / "output")
    assert result.valid_structures == result.unvalidated_structures == 1
    failure, = result.failures
    assert failure.stage == "artifact_delivery"
    assert failure.error == "PermissionError: destination is read-only"
    assert failure.record["completed_validation"]["classification"] == "valid"


def test_failure_counts_toward_limit_and_non_ok_rows_remain_skipped(tmp_path):
    source = make_input(tmp_path, [json.dumps({"status": "export-failed"}), "null", row("after")])
    result = validation.classify_batch_output(source, tmp_path / "output", max_structures=1)
    assert result.total_structures == result.failed_records == 1
    assert result.failures[0].line_number == 2


@pytest.mark.parametrize("kwargs", [{"max_workers": 0}, {"max_workers": -1}, {"max_structures": -1}, {"link_mode": "invalid"}])
def test_invalid_run_configuration_fails_before_output_creation(tmp_path, kwargs):
    source = make_input(tmp_path, [row("one")])
    output = tmp_path / "output"
    with pytest.raises(ValueError):
        validation.classify_batch_output(source, output, **kwargs)
    assert not output.exists()
