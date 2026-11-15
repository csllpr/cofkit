"""Publish classifier-owned files and reconcile preceding snapshots safely."""

from __future__ import annotations

from contextlib import contextmanager
import json
from pathlib import Path
import shutil
import tempfile
from typing import Iterator, Sequence

from .calculation_io import atomic_write_text

_OWNERSHIP_FILE = ".cofkit-classification.json"
_LOCK_FILE = ".cofkit-classification.lock"


def _owned_paths(root: Path, reserved_names: Sequence[str]) -> set[Path]:
    marker = root / _OWNERSHIP_FILE
    if not marker.exists():
        if any((root / name).exists() or (root / name).is_symlink() for name in reserved_names):
            raise ValueError(
                f"classification output {root} contains unmanaged classification paths; "
                "choose a fresh --output-dir (older outputs have no ownership inventory)"
            )
        return set()
    if marker.is_symlink():
        raise ValueError(f"classification ownership file must not be a symlink: {marker}")
    payload = json.loads(marker.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or payload.get("schema_version") != 1 or not isinstance(payload.get("files"), list):
        raise ValueError(f"invalid classification ownership inventory: {marker}")
    paths = set()
    for name in payload["files"]:
        if not isinstance(name, str):
            raise ValueError(f"invalid classification ownership path: {name!r}")
        path = Path(name)
        if path.is_absolute() or not path.parts or any(part in {"..", "."} for part in path.parts) or path.parts[0] not in reserved_names:
            raise ValueError(f"invalid classification ownership path: {name!r}")
        paths.add(path)
    return paths


def _check_destination(root: Path, relative: Path) -> None:
    for parent in relative.parents:
        directory = root / parent
        if directory.is_symlink() or (directory.exists() and not directory.is_dir()):
            raise ValueError(f"classification output parent must be a directory, not a link/file: {directory}")
    destination = root / relative
    if destination.is_dir() and not destination.is_symlink():
        raise ValueError(f"classification file destination is a directory: {destination}")


def _publish(root: Path, stage: Path, owned: set[Path]) -> None:
    current = {p.relative_to(stage) for p in stage.rglob("*") if p.is_file() or p.is_symlink()}
    for relative in owned | current:
        _check_destination(root, relative)
        path = root / relative
        if relative in current and relative not in owned and (path.exists() or path.is_symlink()):
            raise FileExistsError(f"classification would overwrite an unmanaged file: {path}")
    # Keep previous files until every staged file and inventory can be installed.
    # Rename on the same filesystem and roll back if publication fails.
    backup = Path(tempfile.mkdtemp(prefix=".cofkit-classification-backup-", dir=root.parent))
    remove_backup = True
    try:
        moved = []
        installed = []
        try:
            for relative in sorted(owned):
                path = root / relative
                if path.exists() or path.is_symlink():
                    destination = backup / relative
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    path.replace(destination)
                    moved.append(relative)
            for relative in sorted(current):
                destination = root / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                (stage / relative).replace(destination)
                installed.append(relative)
            atomic_write_text(root / _OWNERSHIP_FILE, json.dumps({
                "schema_version": 1, "files": sorted(p.as_posix() for p in current),
            }, indent=2) + "\n")
        except BaseException as publication_error:
            # Do not let cleanup destroy the only remaining old copies if
            # restoration is interrupted or the filesystem keeps failing.
            remove_backup = False
            rollback_errors = []
            moved_paths = set(moved)
            for relative in reversed(installed):
                if relative in moved_paths:
                    # Restore owned files with replace directly, avoiding an
                    # unnecessary unlink and a window with no public copy.
                    continue
                try:
                    (root / relative).unlink(missing_ok=True)
                except OSError as exc:
                    rollback_errors.append(f"remove {relative}: {type(exc).__name__}: {exc}")
            for relative in moved:
                try:
                    (backup / relative).replace(root / relative)
                except OSError as exc:
                    rollback_errors.append(f"restore {relative}: {type(exc).__name__}: {exc}")
            if rollback_errors:
                raise OSError(
                    f"classification publication failed ({type(publication_error).__name__}: {publication_error}); "
                    f"rollback was incomplete ({'; '.join(rollback_errors)}). "
                    f"Recovery files retained at {backup.resolve()}"
                ) from publication_error
            remove_backup = True
            raise
    finally:
        if remove_backup:
            shutil.rmtree(backup)
    # Remove empty old views; retain unrelated files, including user notes.
    for parent in sorted({p for relative in owned for p in relative.parents if p != Path('.')}, key=lambda p: len(p.parts), reverse=True):
        try:
            (root / parent).rmdir()
        except OSError:
            pass


@contextmanager
def classification_run(root: Path, *, source: Path, reserved_names: Sequence[str]) -> Iterator[Path]:
    if source.resolve().is_relative_to(root.resolve()):
        raise ValueError("classification output must not contain the input directory")
    root.mkdir(parents=True, exist_ok=True)
    lock = root / _LOCK_FILE
    try:
        handle = lock.open("x", encoding="utf-8")
    except FileExistsError as exc:
        raise ValueError(f"classification output is in use: {root}; remove {lock.name} only if its previous run was interrupted") from exc
    try:
        with handle:
            handle.write("classification in progress\n")
        owned = _owned_paths(root, reserved_names)
        with tempfile.TemporaryDirectory(prefix=".cofkit-classification-stage-", dir=root.parent) as stage_name:
            stage = Path(stage_name)
            yield stage
            _publish(root, stage, owned)
    finally:
        lock.unlink(missing_ok=True)
