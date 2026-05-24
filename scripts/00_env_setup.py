"""Step 00 — environment check and dependency install.

Verifies that the project root is writable, internet is reachable (Kaggle
needs this toggled on explicitly), and the lightweight Python deps that
aren't preinstalled in the Kaggle base image are available.

The heavy `cellxgene` package is installed lazily, only when the FBS
decoder is actually needed (see `t1d_coa.hpap`).
"""

from __future__ import annotations

import importlib
import os
import random
import shutil
import socket
import subprocess
import sys

import numpy as np

from t1d_coa.config import PROJECT_ROOT, ensure_dirs


SEED = 42
LIGHTWEIGHT_DEPS = {
    "httpx": "httpx>=0.27",
    "tenacity": "tenacity>=8.5",
    "filelock": "filelock>=3.15",
}
DISK_FLOOR_GB = 5.0


def check_internet(host: str = "www.ebi.ac.uk", port: int = 443) -> bool:
    try:
        socket.create_connection((host, port), timeout=5)
        return True
    except OSError:
        return False


def install_missing(specs: dict[str, str]) -> None:
    missing = []
    for mod, spec in specs.items():
        try:
            importlib.import_module(mod)
        except ImportError:
            missing.append(spec)
    if missing:
        print(f"Installing: {missing}")
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", "-q", *missing]
        )


def main() -> None:
    if not check_internet():
        raise RuntimeError(
            "No internet. On Kaggle: notebook Settings → Internet → ON."
        )
    print("✓ Internet reachable")

    install_missing(LIGHTWEIGHT_DEPS)
    print("✓ Dependencies ready")

    random.seed(SEED)
    np.random.seed(SEED)
    os.environ["PYTHONHASHSEED"] = str(SEED)

    ensure_dirs()
    free_gb = shutil.disk_usage(PROJECT_ROOT.parent).free / 1e9
    print(f"✓ Project root: {PROJECT_ROOT}")
    print(f"✓ Free space: {free_gb:.1f} GB")
    if free_gb < DISK_FLOOR_GB:
        print(f"⚠ Less than {DISK_FLOOR_GB} GB free; HPAP fetch may fail.")


if __name__ == "__main__":
    main()
