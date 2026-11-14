"""Reruns publish current views without overwriting unrelated files."""

import json
from pathlib import Path

import pytest

from cofkit import validation


_CIF = """data_model
_cell_length_a 20
_cell_length_b 20
_cell_length_c 20
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'P 1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
m_C1 C 0.2 0.2 0.2
"""


def setup_input(tmp_path):
    source = tmp_path / "input"
    source.mkdir()
    cif = source / "model.cif"
    cif.write_text(_CIF + "# first\n")
    manifest = source / "manifest.jsonl"
    manifest.write_text(json.dumps({"status": "ok", "structure_id": "model", "cif_path": str(cif)}) + "\n")
    return source, cif


@pytest.mark.parametrize("link_mode", ["copy", "hardlink", "symlink"])
def test_rerun_refreshes_files_and_removes_obsolete_classification_and_reason_views(tmp_path, link_mode):
    source, cif = setup_input(tmp_path)
    output = tmp_path / "output"
    first = validation.classify_batch_output(source, output, link_mode=link_mode)
    assert first.valid_structures == 1
    # Atomic source replacement also exercises stale hardlinks.
    replacement = source / "replacement.cif"
    replacement.write_text(_CIF + "# updated\n")
    replacement.replace(cif)
    second = validation.classify_batch_output(source, output, link_mode=link_mode)
    assert second.valid_structures == 1
    assert (output / "valid/cifs/model.cif").read_text().endswith("# updated\n")
    third = validation.classify_batch_output(source, output, link_mode=link_mode,
        thresholds=validation.CoarseValidationThresholds(min_3d_cell_volume=10000))
    assert third.hard_invalid_structures == 1
    assert not (output / "valid/cifs/model.cif").exists()
    assert (output / "hard_invalid/reasons/degenerate_cell/model.cif").is_file()
    empty = validation.classify_batch_output(source, output, max_structures=0)
    assert empty.total_structures == 0
    assert not list(output.rglob("*.cif"))
    assert (output / "classification_manifest.jsonl").read_text() == ""


def test_unrelated_files_survive_and_unmanaged_collisions_preserve_previous_run(tmp_path):
    source, _ = setup_input(tmp_path)
    output = tmp_path / "output"
    validation.classify_batch_output(source, output, max_structures=0)
    note = output / "valid/notes.txt"
    note.write_text("user notes")
    target = output / "valid/cifs/model.cif"
    target.parent.mkdir()
    target.write_text("unmanaged model")
    previous = (output / ".cofkit-classification.json").read_bytes()
    with pytest.raises(FileExistsError, match="unmanaged file"):
        validation.classify_batch_output(source, output)
    assert target.read_text() == "unmanaged model"
    assert note.read_text() == "user notes"
    assert (output / ".cofkit-classification.json").read_bytes() == previous
    target.unlink()
    validation.classify_batch_output(source, output)
    assert note.read_text() == "user notes"


def test_run_failure_preserves_existing_manifests_and_views(tmp_path, monkeypatch):
    source, _ = setup_input(tmp_path)
    output = tmp_path / "output"
    validation.classify_batch_output(source, output, link_mode="copy")
    before = {str(p.relative_to(output)): p.read_bytes() for p in output.rglob("*") if p.is_file()}

    def fail_summary(*args):
        raise OSError("report unavailable")

    monkeypatch.setattr(validation, "_write_classification_summary", fail_summary)
    with pytest.raises(OSError, match="report unavailable"):
        validation.classify_batch_output(source, output)
    after = {str(p.relative_to(output)): p.read_bytes() for p in output.rglob("*") if p.is_file()}
    assert after == before
    assert not list(tmp_path.glob(".cofkit-classification-*"))


def test_publication_failure_rolls_back_existing_files(tmp_path, monkeypatch):
    source, _ = setup_input(tmp_path)
    output = tmp_path / "output"
    validation.classify_batch_output(source, output, link_mode="copy")
    before = {str(p.relative_to(output)): p.read_bytes() for p in output.rglob("*") if p.is_file()}
    replace = Path.replace

    def fail_publish(path, destination):
        if ".cofkit-classification-stage-" in str(path) and str(path).endswith("valid/cifs/model.cif"):
            raise OSError("publication failed")
        return replace(path, destination)

    monkeypatch.setattr(Path, "replace", fail_publish)
    with pytest.raises(OSError, match="publication failed"):
        validation.classify_batch_output(source, output, link_mode="copy")
    after = {str(p.relative_to(output)): p.read_bytes() for p in output.rglob("*") if p.is_file()}
    assert after == before


def test_legacy_output_requires_fresh_directory_and_input_root_is_protected(tmp_path):
    source, _ = setup_input(tmp_path)
    output = tmp_path / "legacy"
    old = output / "valid/cifs/model.cif"
    old.parent.mkdir(parents=True)
    old.write_text("old")
    with pytest.raises(ValueError, match="fresh --output-dir"):
        validation.classify_batch_output(source, output)
    assert old.read_text() == "old"
    with pytest.raises(ValueError, match="must not contain"):
        validation.classify_batch_output(source, source)
    with pytest.raises(ValueError, match="must not contain"):
        validation.classify_batch_output(source, tmp_path)


def test_inventory_cannot_escape_output_and_symlink_parent_is_rejected(tmp_path):
    source, _ = setup_input(tmp_path)
    output = tmp_path / "output"
    validation.classify_batch_output(source, output, max_structures=0)
    marker = output / ".cofkit-classification.json"
    original = marker.read_text()
    marker.write_text(json.dumps({"schema_version": 1, "files": ["../input/model.cif"]}))
    with pytest.raises(ValueError, match="ownership path"):
        validation.classify_batch_output(source, output)
    marker.write_text(original)
    outside = tmp_path / "outside"
    outside.mkdir()
    (output / "valid/cifs").symlink_to(outside, target_is_directory=True)
    with pytest.raises(ValueError, match="parent must be a directory"):
        validation.classify_batch_output(source, output)
    assert not list(outside.iterdir())


def test_existing_lock_prevents_concurrent_publication(tmp_path):
    source, _ = setup_input(tmp_path)
    output = tmp_path / "output"
    output.mkdir()
    lock = output / ".cofkit-classification.lock"
    lock.touch()
    with pytest.raises(ValueError, match="in use"):
        validation.classify_batch_output(source, output)
    assert lock.exists()


def test_same_basename_cannot_silently_deliver_another_records_cif(tmp_path):
    source, cif = setup_input(tmp_path)
    other = source / "other/model.cif"
    other.parent.mkdir()
    other.write_text(_CIF + "# second\n")
    (source / "manifest.jsonl").write_text("".join(json.dumps({
        "status": "ok", "structure_id": str(i), "cif_path": str(path),
    }) + "\n" for i, path in enumerate((cif, other))))
    result = validation.classify_batch_output(source, tmp_path / "output", link_mode="copy")
    assert result.valid_structures == result.failed_records == 1
    assert "duplicate classification destination" in result.failures[0].error


def test_partial_copy_failure_does_not_publish_a_truncated_cif(tmp_path, monkeypatch):
    source, _ = setup_input(tmp_path)
    output = tmp_path / "output"

    def partial_copy(source, destination):
        destination.write_text("partial")
        raise OSError("copy interrupted")

    monkeypatch.setattr(validation.shutil, "copy2", partial_copy)
    result = validation.classify_batch_output(source, output, link_mode="copy")
    assert result.failed_records == result.unvalidated_structures == 1
    assert not list(output.rglob("*.cif"))
