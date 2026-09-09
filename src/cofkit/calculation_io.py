"""Calculation attempt identity, atomic publication, and numeric input checks."""

from __future__ import annotations

import hashlib
import json
import math
import os
import platform
from dataclasses import fields, is_dataclass
from functools import lru_cache
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import get_args, get_origin, get_type_hints
from uuid import uuid4

from ._version import __version__

# LAMMPS RanMars rejects seeds above this value (random_mars.cpp). Using
# the strictest supported range also makes derived MC seeds portable.
MAX_ENGINE_SEED = 900_000_000


def validate_numbers(settings: object) -> None:
    """Reject nonfinite values and nonintegral integer settings before rendering."""
    hints = get_type_hints(type(settings))
    for field in fields(settings):
        value = getattr(settings, field.name)
        if value is None:
            continue
        annotation = hints[field.name]
        alternatives = (
            (annotation,)
            if get_origin(annotation) in (tuple, list)
            else get_args(annotation) or (annotation,)
        )
        if bool in alternatives and not isinstance(value, bool):
            raise ValueError(f"{field.name} must be a boolean.")
        if int in alternatives and (
            not isinstance(value, int) or isinstance(value, bool)
        ):
            raise ValueError(f"{field.name} must be an integer.")
        if (
            float in alternatives
            and str not in alternatives
            and (not isinstance(value, (float, int)) or isinstance(value, bool))
        ):
            raise ValueError(f"{field.name} must be numeric.")
        values = value if isinstance(value, (tuple, list)) else (value,)
        for item in values:
            if is_dataclass(item):
                validate_numbers(item)
            elif isinstance(item, (float, int)) and not math.isfinite(item):
                raise ValueError(f"{field.name} must contain only finite numbers.")


def file_identity(path: Path) -> dict[str, object]:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"path": str(path.resolve()), "sha256": digest.hexdigest()}


@lru_cache(maxsize=1)
def implementation_identity() -> dict[str, object]:
    """Fingerprint the installed adapter code as well as its release label."""
    root = Path(__file__).parent
    digest = hashlib.sha256()
    for path in sorted(root.rglob("*.py")):
        digest.update(str(path.relative_to(root)).encode())
        digest.update(path.read_bytes())
    dependencies = {}
    for package in ("gemmi", "rdkit", "openbabel", "pymatgen"):
        try:
            dependencies[package] = version(package)
        except PackageNotFoundError:
            dependencies[package] = None
    return {
        "version": __version__,
        "source_sha256": digest.hexdigest(),
        "python": platform.python_version(),
        "dependencies": dependencies,
    }


def atomic_write_text(path: Path, text: str, *, encoding: str = "utf-8") -> None:
    temporary = path.with_name(f".{path.name}.{uuid4().hex}.tmp")
    try:
        with temporary.open("w", encoding=encoding) as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def begin_attempt(
    root: Path, *, input_path: Path, settings: object, binary: Path | None
) -> Path:
    """Claim a fresh directory; preserve every existing directory and its contents.

    A new output root is itself the first attempt. Subsequent invocations create
    children; resume is intentionally not inferred from previous files.
    """
    try:
        root.mkdir(parents=True, exist_ok=False)
        attempt = root
    except FileExistsError:
        attempt = root / f"attempt-{uuid4().hex}"
        attempt.mkdir(parents=True, exist_ok=False)
    atomic_write_text(
        attempt / "attempt.json",
        json.dumps(
            {
                "schema_version": 1,
                "status": "incomplete",
                "input": file_identity(input_path),
                "executable": file_identity(binary) if binary else None,
                "settings": settings.to_dict(),
                "cofkit": implementation_identity(),
            },
            indent=2,
            allow_nan=False,
        ),
    )
    return attempt


def record_execution(directory: Path, command: list[str]) -> None:
    """Record the actual staged inputs, including parameter-file contents."""
    inputs = [
        file_identity(path)
        for path in sorted(directory.rglob("*"))
        if path.is_file()
        and path.suffix not in {".log"}
        and path.name not in {"execution.json", "attempt.json"}
    ]
    atomic_write_text(
        directory / "execution.json",
        json.dumps(
            {
                "command": command,
                "executable": file_identity(Path(command[0])),
                "inputs": inputs,
            },
            indent=2,
            allow_nan=False,
        ),
    )


def finish_attempt(directory: Path) -> None:
    manifest = directory / "attempt.json"
    payload = json.loads(manifest.read_text())
    payload["status"] = "completed"
    atomic_write_text(manifest, json.dumps(payload, indent=2, allow_nan=False))


def derive_seed(seed: int, *labels: object) -> int:
    """Stable child seed in the range accepted by LAMMPS RanMars and MC engines."""
    digest = hashlib.sha256(json.dumps([seed, *labels]).encode()).digest()
    return int.from_bytes(digest[:8], "big") % MAX_ENGINE_SEED + 1
