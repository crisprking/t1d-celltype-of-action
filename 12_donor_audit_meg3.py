"""05_hpap_expression: Decode per-gene expression vectors from CellxGene FBS blobs.

Extracted from notebook cell 30. Part of the t1d-celltype-of-action pipeline.
Run modules in numeric order — each consumes outputs from earlier steps.
"""

# Verified contract:
#   verb    = PUT
#   path    = {CXG_API}/data/var
#   body    = {"filter": {"var": {"index": [i, j, k, ...]}}}
#   accept  = application/octet-stream
#   decoder = server.common.fbs.matrix.decode_matrix_fbs(bytes) -> pd.DataFrame
#
# Smoke test passed: beta cells lead INS mean expression by ~11x over delta
# (the next highest), exactly as expected for the canonical beta marker.
#
# Output: gene × cell_type × disease_state long-form TSV with mean,
# median, percent-expressing, and cell counts.
# ============================================================

import os, sys, json, hashlib, logging, warnings, subprocess, importlib
from pathlib import Path
from datetime import datetime, timezone

import httpx
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

warnings.filterwarnings("ignore", message="Unverified HTTPS request")
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S", force=True)
log = logging.getLogger("hpap")

PROJECT = Path("/kaggle/working/t1d_mech")
RAW = PROJECT / "data" / "raw"
INTERIM = PROJECT / "data" / "interim"
PROCESSED = PROJECT / "data" / "processed"
PROVENANCE = PROJECT / "data" / "PROVENANCE.md"
for d in (RAW, INTERIM, PROCESSED):
    d.mkdir(parents=True, exist_ok=True)

CXG_BASE = "https://faryabi16.pmacs.upenn.edu/view/T1D_T2D_public.h5ad"
CXG_API = f"{CXG_BASE}/api/v0.2"

# --- Ensure cellxgene's FBS decoder is importable (installed last cell) ---
try:
    from server.common.fbs.matrix import decode_matrix_fbs
except ImportError:
    log.info("Installing cellxgene (FBS decoder)...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q",
                           "--no-deps", "cellxgene", "flatbuffers"])
    from server.common.fbs.matrix import decode_matrix_fbs
log.info("FBS decoder ready: %s", decode_matrix_fbs)

# --- Provenance helper ----------------------------------------------------
def _sha(p, chunk=1 << 20):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for d in iter(lambda: f.read(chunk), b""):
            h.update(d)
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

# --- HTTP client (SSL verify off for this host; cert SAN mismatch) --------
client = httpx.Client(timeout=300.0, verify=False, follow_redirects=True)

# --- 1. Load obs (cached) + recover gene → index map ----------------------
obs_path = PROCESSED / "hpap_cellxgene_obs.tsv"
obs = pd.read_csv(obs_path, sep="\t")
log.info("Loaded obs: %d cells × %d annotations", *obs.shape)
assert len(obs) == 222077, "obs cell count != 222,077 — re-run the obs cell"

def fetch_gene_list():
    """Return the list of gene symbols (var.name_0) in row order."""
    cache = INTERIM / "hpap_gene_list.json"
    if cache.exists():
        return json.loads(cache.read_text())
    log.info("Fetching var/name_0 from CellxGene...")
    r = client.get(f"{CXG_API}/annotations/var",
                   params={"annotation-name": "name_0"})
    r.raise_for_status()
    blob = r.content
    marker = b'["name_0"]'
    i = blob.find(marker)
    start = blob.find(b'["', i + len(marker))
    j, depth, in_str, esc = start, 0, False, False
    while j < len(blob):
        c = blob[j:j+1]
        if esc:    esc = False
        elif c == b'\\': esc = True
        elif c == b'"':  in_str = not in_str
        elif not in_str:
            if c == b'[': depth += 1
            elif c == b']':
                depth -= 1
                if depth == 0: break
        j += 1
    genes = json.loads(blob[start:j+1].decode("utf-8"))
    cache.write_text(json.dumps(genes))
    log_artifact(cache, source=f"GET {CXG_API}/annotations/var?annotation-name=name_0",
                 notes=f"{len(genes):,} gene symbols, var row order")
    return genes

genes_all = fetch_gene_list()
gene_to_idx = {g: i for i, g in enumerate(genes_all)}
log.info("Var dim: %d genes", len(genes_all))

# --- 2. Pick the pilot 20 genes from our 145 T1D loci ---------------------
loci = pd.read_csv(PROCESSED / "t1d_independent_loci_annotated.tsv", sep="\t")

def _expand(s):
    if pd.isna(s) or not str(s).strip(): return []
    s = str(s).replace(" - ", ", ")
    return [g.strip() for g in s.split(",") if g.strip()]

loci["genes_list"] = loci["MAPPED_GENE"].apply(_expand)
gene_rows = []
for _, r in loci.iterrows():
    for g in r["genes_list"]:
        gene_rows.append({"gene": g, "lead_snp": r["SNP"], "P": r["P"],
                          "CHR": r["CHR"], "n_studies": r.get("n_studies", 1)})
gene_df = pd.DataFrame(gene_rows)
gene_best = gene_df.sort_values("P").drop_duplicates("gene", keep="first")

top_15 = gene_best.head(15)
spread_idx = [i for i in (30, 60, 100, 140, 180) if i < len(gene_best)]
pilot = (pd.concat([top_15, gene_best.iloc[spread_idx]])
           .drop_duplicates("gene")
           .head(20)
           .reset_index(drop=True))

# Only keep genes that exist in the HPAP var
pilot["in_hpap"] = pilot["gene"].isin(gene_to_idx)
missing = pilot[~pilot["in_hpap"]]
if len(missing):
    log.warning("Pilot genes missing from HPAP var: %s", missing["gene"].tolist())
pilot = pilot[pilot["in_hpap"]].drop(columns=["in_hpap"]).reset_index(drop=True)
pilot["var_idx"] = pilot["gene"].map(gene_to_idx)

print("\nPilot genes (kept):")
print(pilot[["gene", "lead_snp", "CHR", "P", "var_idx"]].to_string(index=False))

# --- 3. Fetch expression (batched, cached) --------------------------------
EXPR_CACHE_DIR = RAW / "cellxgene_expr"
EXPR_CACHE_DIR.mkdir(parents=True, exist_ok=True)

def fetch_one_gene(symbol: str, idx: int) -> np.ndarray:
    """Fetch one gene's expression vector; cache the raw FBS blob to disk."""
    cache = EXPR_CACHE_DIR / f"{symbol}_idx{idx}.fbs"
    if cache.exists() and cache.stat().st_size > 1000:
        blob = cache.read_bytes()
    else:
        r = client.put(
            f"{CXG_API}/data/var",
            json={"filter": {"var": {"index": [idx]}}},
            headers={"Accept": "application/octet-stream"},
        )
        r.raise_for_status()
        blob = r.content
        cache.write_bytes(blob)
        log_artifact(cache, source=f"PUT {CXG_API}/data/var var.index=[{idx}]",
                     notes=f"HPAP {symbol} expression FBS blob (SSL verify off)")
    df = decode_matrix_fbs(blob)
    vec = df.iloc[:, 0].to_numpy()
    if len(vec) != len(obs):
        raise RuntimeError(
            f"{symbol}: returned {len(vec)} values, expected {len(obs)}"
        )
    return vec

log.info("Fetching %d gene expression vectors...", len(pilot))
expr = {}
for _, row in tqdm(pilot.iterrows(), total=len(pilot), desc="genes"):
    expr[row["gene"]] = fetch_one_gene(row["gene"], int(row["var_idx"]))

expr_df = pd.DataFrame(expr)  # cells × genes (in obs row order)
log.info("Expression matrix: %d cells × %d genes", *expr_df.shape)

# --- 4. Aggregate to gene × cell_type × disease_state ---------------------
log.info("Aggregating to gene × cell_type × disease_state...")
records = []
group_keys = ["cell_type", "disease_state"]
groups = obs.groupby(group_keys, observed=True).indices  # dict: key -> row positions

for (cell_type, disease_state), row_pos in tqdm(groups.items(), desc="groups"):
    n_cells = len(row_pos)
    sub = expr_df.iloc[row_pos]
    for gene in expr_df.columns:
        vals = sub[gene].to_numpy()
        records.append({
            "gene": gene,
            "cell_type": cell_type,
            "disease_state": disease_state,
            "n_cells": n_cells,
            "mean": float(vals.mean()),
            "median": float(np.median(vals)),
            "pct_expressing": float((vals > 0).mean() * 100),
            "n_nonzero": int((vals > 0).sum()),
        })

result = pd.DataFrame(records)
result = result.sort_values(
    ["gene", "cell_type", "disease_state"]
).reset_index(drop=True)

out_path = PROCESSED / "hpap_pilot20_expression_by_celltype_disease.tsv"
result.to_csv(out_path, sep="\t", index=False)
log_artifact(out_path,
             source="Aggregated from per-gene CellxGene blobs (see preceding rows)",
             notes=f"{len(result):,} rows: {pilot.gene.nunique()} genes × cell_type × disease_state")

# --- 5. Headline summary --------------------------------------------------
print(f"\n=== Pilot 20 T1D-locus genes × cell_type × disease_state ===")
print(f"  Rows: {len(result):,}")
print(f"  File: {out_path}")
print(f"\nCanonical sanity (INS mean by cell_type, all disease states pooled):")
sanity = (result[result["gene"] == "INS"]
          .groupby("cell_type")
          .apply(lambda g: (g["mean"] * g["n_cells"]).sum() / g["n_cells"].sum(),
                 include_groups=False)
          .sort_values(ascending=False))
print(sanity.to_string())
print("\nTop 5 (gene, cell_type, disease_state) by mean expression:")
print(result.nlargest(5, "mean")[
    ["gene", "cell_type", "disease_state", "n_cells", "mean", "pct_expressing"]
].to_string(index=False))
