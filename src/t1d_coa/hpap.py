"""Thin client for the HPAP CellxGene REST API.

Encapsulates the two awkward parts of talking to this endpoint:

1. SSL verification must be disabled (cert hostname mismatch upstream).
2. Expression payloads come back as raw FlatBuffers, not Arrow IPC; the
   `cellxgene` package's internal `decode_matrix_fbs` is the only stable
   decoder published for the schema.
"""

from __future__ import annotations

import json
import subprocess
import sys
import warnings
from pathlib import Path
from typing import Any

import httpx
import numpy as np
import pandas as pd

from .config import EXPR_CACHE, HPAP_API, OBS_CACHE, SSL_VERIFY
from .provenance import log_artifact

warnings.filterwarnings("ignore", message="Unverified HTTPS request")


def _client() -> httpx.Client:
    return httpx.Client(timeout=300.0, verify=SSL_VERIFY, follow_redirects=True)


def _decoder():
    """Return the cellxgene FBS matrix decoder, installing it on demand."""
    try:
        from server.common.fbs.matrix import decode_matrix_fbs  # type: ignore
    except ImportError:
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", "-q",
             "--no-deps", "cellxgene", "flatbuffers"]
        )
        from server.common.fbs.matrix import decode_matrix_fbs  # type: ignore
    return decode_matrix_fbs


def _extract_json_array(blob: bytes, column_name: str) -> list[Any]:
    """Pull the data array for `column_name` out of an obs FBS payload."""
    marker = f'["{column_name}"]'.encode()
    idx = blob.find(marker)
    if idx < 0:
        raise ValueError(f"Column marker {marker!r} not found")
    start = blob.find(b'["', idx + len(marker))
    if start < 0:
        raise ValueError("No data array found after column marker")

    i, depth, in_str, esc = start, 0, False, False
    while i < len(blob):
        c = blob[i : i + 1]
        if esc:
            esc = False
        elif c == b"\\":
            esc = True
        elif c == b'"':
            in_str = not in_str
        elif not in_str:
            if c == b"[":
                depth += 1
            elif c == b"]":
                depth -= 1
                if depth == 0:
                    break
        i += 1
    return json.loads(blob[start : i + 1].decode("utf-8"))


def fetch_obs_column(name: str, client: httpx.Client | None = None) -> pd.Series:
    """Fetch one obs annotation column; cache the raw blob."""
    OBS_CACHE.mkdir(parents=True, exist_ok=True)
    cache = OBS_CACHE / f"obs_{name}.bin"

    if cache.exists() and cache.stat().st_size > 1000:
        blob = cache.read_bytes()
    else:
        c = client or _client()
        r = c.get(f"{HPAP_API}/annotations/obs",
                  params={"annotation-name": name})
        r.raise_for_status()
        blob = r.content
        cache.write_bytes(blob)
        log_artifact(
            cache,
            source=f"GET {HPAP_API}/annotations/obs?annotation-name={name}",
            notes=f"HPAP obs '{name}' (SSL verify off; see provenance note)",
        )

    return pd.Series(_extract_json_array(blob, name), name=name)


def fetch_gene_list(interim_cache: Path,
                    client: httpx.Client | None = None) -> list[str]:
    """Pull var.name_0 — the gene symbol vector in matrix row order."""
    if interim_cache.exists():
        return json.loads(interim_cache.read_text())

    c = client or _client()
    r = c.get(f"{HPAP_API}/annotations/var",
              params={"annotation-name": "name_0"})
    r.raise_for_status()
    genes = _extract_json_array(r.content, "name_0")

    interim_cache.parent.mkdir(parents=True, exist_ok=True)
    interim_cache.write_text(json.dumps(genes))
    log_artifact(
        interim_cache,
        source=f"GET {HPAP_API}/annotations/var?annotation-name=name_0",
        notes=f"{len(genes):,} gene symbols, var row order",
    )
    return genes


def fetch_gene_expression(symbol: str, var_idx: int,
                          client: httpx.Client | None = None) -> np.ndarray:
    """Fetch one gene's per-cell expression vector; cache the FBS blob."""
    EXPR_CACHE.mkdir(parents=True, exist_ok=True)
    cache = EXPR_CACHE / f"{symbol}_idx{var_idx}.fbs"

    if cache.exists() and cache.stat().st_size > 1000:
        blob = cache.read_bytes()
    else:
        c = client or _client()
        r = c.put(
            f"{HPAP_API}/data/var",
            json={"filter": {"var": {"index": [var_idx]}}},
            headers={"Accept": "application/octet-stream"},
        )
        r.raise_for_status()
        blob = r.content
        cache.write_bytes(blob)
        log_artifact(
            cache,
            source=f"PUT {HPAP_API}/data/var var.index=[{var_idx}]",
            notes=f"HPAP {symbol} expression FBS (SSL verify off)",
        )

    decode_matrix_fbs = _decoder()
    return decode_matrix_fbs(blob).iloc[:, 0].to_numpy()
