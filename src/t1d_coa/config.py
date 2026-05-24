"""Paths and constants used across the pipeline.

A single import surface so no module hard-codes a Kaggle path or an HPAP URL.
Override PROJECT_ROOT by setting the T1D_COA_ROOT environment variable.
"""

from __future__ import annotations

import os
from pathlib import Path

# Project root. Defaults to a Kaggle-style layout for reproducibility with the
# original notebook, but any local checkout works once T1D_COA_ROOT is set.
PROJECT_ROOT = Path(os.environ.get("T1D_COA_ROOT", "/kaggle/working/t1d_mech"))

RAW = PROJECT_ROOT / "data" / "raw"
INTERIM = PROJECT_ROOT / "data" / "interim"
PROCESSED = PROJECT_ROOT / "data" / "processed"
RESULTS = PROJECT_ROOT / "results"
TOOLS = PROJECT_ROOT / "tools"
REFERENCE = PROJECT_ROOT / "reference"

EXPR_CACHE = RAW / "cellxgene_expr"
OBS_CACHE = RAW / "cellxgene_obs"

PROVENANCE = PROJECT_ROOT / "data" / "PROVENANCE.md"

# HPAP CellxGene REST endpoint (Faryabi lab, U Penn).
HPAP_BASE = "https://faryabi16.pmacs.upenn.edu/view/T1D_T2D_public.h5ad"
HPAP_API = f"{HPAP_BASE}/api/v0.2"

# Note: the upstream cert's SAN list omits this hostname; httpx with verify=True
# refuses, browsers tolerate. We disable verification for this host only and
# rely on downstream integrity checks (cell counts, known marker behavior).
SSL_VERIFY = False

# Expected atlas dimensions; asserted to catch silent schema drift.
EXPECTED_N_CELLS = 222_077

# Tau cell-type specificity threshold and expression floor.
TAU_THRESHOLD = 0.70
EXPR_FLOOR_LOG1P = 0.5

# Permutation test reproducibility.
PERM_SEED = 42
N_PERMUTATIONS = 10_000


def ensure_dirs() -> None:
    """Create the standard directory layout if any pieces are missing."""
    for d in (RAW, INTERIM, PROCESSED, RESULTS, EXPR_CACHE, OBS_CACHE):
        d.mkdir(parents=True, exist_ok=True)
