"""00_env_setup: Set up working directories and project paths.

Extracted from notebook cell 0. Part of the t1d-celltype-of-action pipeline.
Run modules in numeric order — each consumes outputs from earlier steps.
"""

# Before running: Notebook settings → Internet: ON
# Optional: Notebook settings → Accelerator: GPU T4 x2 (for week 3+)

import sys, os, subprocess, importlib
from pathlib import Path

# Kaggle paths
WORKING = Path("/kaggle/working")
INPUT = Path("/kaggle/input")
PROJECT = WORKING / "t1d_mech"
RAW = PROJECT / "data" / "raw"
INTERIM = PROJECT / "data" / "interim"
PROCESSED = PROJECT / "data" / "processed"
RESULTS = PROJECT / "results"
for d in (RAW, INTERIM, PROCESSED, RESULTS):
    d.mkdir(parents=True, exist_ok=True)

# Verify internet (Kaggle disables it by default)
def check_internet():
    import socket
    try:
        socket.create_connection(("www.ebi.ac.uk", 443), timeout=5)
        return True
    except OSError:
        return False

if not check_internet():
    raise RuntimeError(
        "No internet. Open notebook Settings (right sidebar) → "
        "Internet → toggle ON, then re-run."
    )
print("✓ Internet on")

# Install only what isn't preinstalled. Kaggle base image already has
# pandas, numpy, scipy, scikit-learn, matplotlib, requests, tqdm.
needed = {
    "httpx": "httpx>=0.27",
    "tenacity": "tenacity>=8.5",
    "filelock": "filelock>=3.15",
}
to_install = []
for mod, spec in needed.items():
    try:
        importlib.import_module(mod)
    except ImportError:
        to_install.append(spec)

if to_install:
    print(f"Installing: {to_install}")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", *to_install])
print("✓ Deps ready")

# Fix random seed
import random, numpy as np
SEED = 42
random.seed(SEED); np.random.seed(SEED); os.environ["PYTHONHASHSEED"] = str(SEED)

# Disk-space awareness — fail loud before downloading 15 GB
import shutil
free_gb = shutil.disk_usage(WORKING).free / 1e9
print(f"✓ /kaggle/working free space: {free_gb:.1f} GB")
if free_gb < 5:
    print("⚠ Less than 5 GB free; HPAP scRNA-seq pull (cell 4) will fail.")

print(f"✓ Project root: {PROJECT}")
