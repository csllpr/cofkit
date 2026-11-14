"""Run portable wheel checks outside the checkout using the current interpreter."""

import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


root = Path(__file__).resolve().parents[2]
env = dict(os.environ)
env.pop("PYTHONPATH", None)
with tempfile.TemporaryDirectory(prefix="cofkit-wheel-") as directory:
    work = Path(directory)
    shutil.copytree(root / "tests", work / "tests", ignore=shutil.ignore_patterns("__pycache__"))
    (work / "docs").mkdir()
    shutil.copy2(root / "docs/COFid_Specification_v1.2.md", work / "docs")
    subprocess.run([
        sys.executable, "-c",
        'import cofkit; from pathlib import Path; assert "site-packages" in Path(cofkit.__file__).parts, cofkit.__file__; print(cofkit.__file__)',
    ], cwd=work, env=env, check=True)
    subprocess.run([
        sys.executable, "-m", "pytest", "-q",
        "tests/test_distribution_smoke.py", "tests/test_core.py", "tests/test_cli_errors.py",
    ], cwd=work, env=env, check=True)
