"""04_hpap_load: Load HPAP scRNA-seq atlas via CellxGene REST API (obs only).

Extracted from notebook cell 22. Part of the t1d-celltype-of-action pipeline.
Run modules in numeric order — each consumes outputs from earlier steps.
"""

# Penn's CellxGene host (faryabi16.pmacs.upenn.edu) serves a cert that
# doesn't include this hostname in its SAN list. Browsers tolerate this
# variably; Python's strict default refuses. We explicitly disable
# verification for THIS HOST ONLY. Documented in PROVENANCE.md.

import os, json, logging, hashlib, ssl
from pathlib import Path
from datetime import datetime, timezone
import warnings

import httpx
import pandas as pd
from tqdm.auto import tqdm

# Suppress just the urllib3/httpx InsecureRequestWarning noise for this host
warnings.filterwarnings("ignore", message="Unverified HTTPS request")

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S", force=True)
log = logging.getLogger("hpap")

PROJECT = Path("/kaggle/working/t1d_mech")
RAW = PROJECT / "data" / "raw"
PROCESSED = PROJECT / "data" / "processed"
for d in (RAW, PROCESSED):
    d.mkdir(parents=True, exist_ok=True)
PROVENANCE = PROJECT / "data" / "PROVENANCE.md"

CXG_BASE = "https://faryabi16.pmacs.upenn.edu/view/T1D_T2D_public.h5ad"
CXG_API = f"{CXG_BASE}/api/v0.2"

# --- Document the SSL choice explicitly -----------------------------------
SSL_NOTE = (
    "NOTE: requests to faryabi16.pmacs.upenn.edu use verify=False because the "
    "server's TLS cert hostname doesn't match. The endpoint is a public "
    "scientific data API (no credentials sent). Data integrity is validated "
    "downstream by checking expected cell counts (~222k) and known cell-type "
    "labels."
)
log.warning(SSL_NOTE)

# Single client with SSL verification off, reused for all requests
_client = httpx.Client(timeout=120.0, verify=False, follow_redirects=True)

# --- Provenance helper -----------------------------------------------------
def _sha(p, chunk=1 << 20):
    h = hashlib.sha256()
    with p.open("rb") as f:
        while (d := f.read(chunk)): h.update(d)
    return h.hexdigest()

def log_artifact(path, source, notes=""):
    if not PROVENANCE.exists():
        PROVENANCE.write_text(
            "# Data Provenance\n\n| Path | Source | Date (UTC) | sha256 | Bytes | Notes |\n"
            "|------|--------|------------|--------|-------|-------|\n"
        )
    sha = _sha(path); size = path.stat().st_size
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    rel = path.relative_to(PROJECT)
    with PROVENANCE.open("a") as f:
        f.write(f"| `{rel}` | {source.replace('|','\\|')} | {now} | "
                f"`{sha[:16]}…` | {size:,} | {notes.replace('|','\\|')} |\n")

# --- FlatBuffer JSON-array extractor (verified working on your obs blob) ---
def _extract_json_array(blob: bytes, column_name: str) -> list:
    marker = f'["{column_name}"]'.encode()
    idx = blob.find(marker)
    if idx < 0:
        raise ValueError(f"Column marker {marker!r} not found in response")
    start = blob.find(b'["', idx + len(marker))
    if start < 0:
        raise ValueError("No data array found after column marker")
    i = start; depth = 0; in_str = False; esc = False
    while i < len(blob):
        c = blob[i:i+1]
        if esc:    esc = False
        elif c == b'\\': esc = True
        elif c == b'"':  in_str = not in_str
        elif not in_str:
            if c == b'[': depth += 1
            elif c == b']':
                depth -= 1
                if depth == 0: break
        i += 1
    return json.loads(blob[start:i+1].decode("utf-8"))

# --- Sanity check first ----------------------------------------------------
log.info("Connectivity check...")
r = _client.get(CXG_API + "/config")
log.info("  /config -> %d, %d bytes, ct=%s",
         r.status_code, len(r.content), r.headers.get("content-type"))
if r.status_code != 200:
    raise RuntimeError(f"API not reachable: status {r.status_code}")

# --- Fetch obs columns with proper error reporting -------------------------
def fetch_obs_column(name: str) -> pd.Series:
    cache_path = RAW / "cellxgene_obs" / f"obs_{name}.bin"
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    if cache_path.exists() and cache_path.stat().st_size > 1000:
        log.info("  cached: %s (%.1f MB)", name, cache_path.stat().st_size / 1e6)
        blob = cache_path.read_bytes()
    else:
        log.info("  fetching: %s", name)
        r = _client.get(f"{CXG_API}/annotations/obs",
                        params={"annotation-name": name})
        r.raise_for_status()
        blob = r.content
        cache_path.write_bytes(blob)
        log_artifact(cache_path,
                     source=f"GET {CXG_API}/annotations/obs?annotation-name={name}",
                     notes=f"HPAP CellxGene obs '{name}' (SSL verify disabled, see PROVENANCE note)")

    return pd.Series(_extract_json_array(blob, name), name=name)

OBS_COLUMNS = ["cell_type", "disease_state", "donor_id", "age", "sex", "race", "assay"]
obs_data = {}
for col in OBS_COLUMNS:
    try:
        obs_data[col] = fetch_obs_column(col)
        log.info("  ✓ %s: %d values, %d unique",
                 col, len(obs_data[col]), obs_data[col].nunique())
    except httpx.HTTPStatusError as e:
        log.warning("  ✗ %s: HTTP %d", col, e.response.status_code)
    except Exception as e:
        log.warning("  ✗ %s: %s: %s", col, type(e).__name__, e)

if not obs_data:
    raise RuntimeError("All obs fetches failed — see warnings above")

obs_df = pd.DataFrame(obs_data)
obs_path = PROCESSED / "hpap_cellxgene_obs.tsv"
obs_df.to_csv(obs_path, sep="\t", index=False)
log_artifact(obs_path, source="Combined CellxGene obs columns",
             notes=f"{len(obs_df):,} cells × {obs_df.shape[1]} annotations")

# --- Headline numbers ------------------------------------------------------
print(f"\n=== HPAP atlas via CellxGene API ===")
print(f"  Total cells: {len(obs_df):,}")
print(f"  Annotations available: {list(obs_df.columns)}")

if "cell_type" in obs_df:
    print(f"\n  Cell type breakdown:")
    for ct, n in obs_df["cell_type"].value_counts().items():
        print(f"    {ct:35s} {n:>8,}")

if "disease_state" in obs_df:
    print(f"\n  Disease state breakdown:")
    for ds, n in obs_df["disease_state"].value_counts().items():
        print(f"    {ds:35s} {n:>8,}")

if "donor_id" in obs_df:
    print(f"\n  Unique donors: {obs_df['donor_id'].nunique()}")

if "cell_type" in obs_df and "disease_state" in obs_df:
    print(f"\n  Cell type × disease state crosstab:")
    print(pd.crosstab(obs_df["cell_type"], obs_df["disease_state"], margins=True).to_string())

print(f"\nSaved: {obs_path}")
