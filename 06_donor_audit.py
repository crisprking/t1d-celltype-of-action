"""07_full_panel_fetch: Fetch + aggregate expression for all 150 testable T1D-locus genes.

Extracted from notebook cell 33. Part of the t1d-celltype-of-action pipeline.
Run modules in numeric order — each consumes outputs from earlier steps.
"""

# Reuses the per-gene PUT + decode_matrix_fbs path from the pilot cell.
# Idempotent: hits cached .fbs blobs on second run, so re-runs are
# essentially free.
#
# Expected wall-clock for a cold run:
#   ~250-400 unique genes × ~0.5-1.5s per PUT ≈ 3-8 minutes
# Disk: each blob is ~880KB, so total ~250-350 MB. Fits within Kaggle's
# 20 GB free space comfortably.
# ============================================================

import os, sys, json, hashlib, logging, warnings, subprocess
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

CXG_BASE = "https://faryabi16.pmacs.upenn.edu/view/T1D_T2D_public.h5ad"
CXG_API = f"{CXG_BASE}/api/v0.2"

# --- Decoder + provenance (re-import in case kernel restarted) -----------
try:
    from server.common.fbs.matrix import decode_matrix_fbs
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q",
                           "--no-deps", "cellxgene", "flatbuffers"])
    from server.common.fbs.matrix import decode_matrix_fbs

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

client = httpx.Client(timeout=300.0, verify=False, follow_redirects=True)
EXPR_CACHE = RAW / "cellxgene_expr"
EXPR_CACHE.mkdir(parents=True, exist_ok=True)

# --- Load obs + gene list ------------------------------------------------
obs = pd.read_csv(PROCESSED / "hpap_cellxgene_obs.tsv", sep="\t")
assert len(obs) == 222077
gene_list_path = INTERIM / "hpap_gene_list.json"
genes_all = json.loads(gene_list_path.read_text())
gene_to_idx = {g: i for i, g in enumerate(genes_all)}
log.info("Var dim: %d genes; obs: %d cells", len(genes_all), len(obs))

# --- 1. Expand all 145 loci → unique gene list ---------------------------
loci = pd.read_csv(PROCESSED / "t1d_independent_loci_annotated.tsv", sep="\t")

def _expand(s):
    if pd.isna(s) or not str(s).strip(): return []
    return [g.strip() for g in str(s).replace(" - ", ", ").split(",") if g.strip()]

loci["genes_list"] = loci["MAPPED_GENE"].apply(_expand)
gene_rows = []
for _, r in loci.iterrows():
    for g in r["genes_list"]:
        gene_rows.append({"gene": g, "lead_snp": r["SNP"], "P": r["P"],
                          "CHR": r["CHR"], "n_studies": r.get("n_studies", 1)})
gene_df = pd.DataFrame(gene_rows)

# Best (smallest) p-value per gene if same gene comes up via multiple loci
gene_best = gene_df.sort_values("P").drop_duplicates("gene", keep="first")
log.info("Unique genes across 145 loci: %d", len(gene_best))

# Filter to those present in HPAP var
gene_best["in_hpap"] = gene_best["gene"].isin(gene_to_idx)
missing = gene_best.loc[~gene_best["in_hpap"], "gene"].tolist()
log.info("In HPAP var: %d; missing: %d", gene_best["in_hpap"].sum(), len(missing))
if missing:
    (INTERIM / "missing_genes_from_hpap.txt").write_text("\n".join(missing))
    log.info("  missing gene list saved to interim/missing_genes_from_hpap.txt")

target = gene_best[gene_best["in_hpap"]].copy()
target["var_idx"] = target["gene"].map(gene_to_idx)
target = target.reset_index(drop=True)

# --- 2. Fetch loop (idempotent; per-gene cache) --------------------------
def fetch_one_gene(symbol: str, idx: int) -> np.ndarray:
    cache = EXPR_CACHE / f"{symbol}_idx{idx}.fbs"
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
                     notes=f"HPAP {symbol} expression FBS (SSL verify off)")
    df = decode_matrix_fbs(blob)
    return df.iloc[:, 0].to_numpy()

# Report cache hit-rate up front so user sees what's about to happen
n_cached = sum(1 for _, row in target.iterrows()
               if (EXPR_CACHE / f"{row['gene']}_idx{row['var_idx']}.fbs").exists())
log.info("Cache: %d/%d genes already fetched", n_cached, len(target))

expr = {}
failures = []
for _, row in tqdm(target.iterrows(), total=len(target), desc="genes"):
    try:
        expr[row["gene"]] = fetch_one_gene(row["gene"], int(row["var_idx"]))
    except Exception as e:
        failures.append((row["gene"], str(e)))
        log.warning("  failed: %s — %s", row["gene"], e)

if failures:
    log.warning("%d gene fetches failed (will be omitted from output)", len(failures))
    (INTERIM / "failed_gene_fetches.txt").write_text(
        "\n".join(f"{g}\t{e}" for g, e in failures)
    )

# --- 3. Aggregate to gene × cell_type × disease_state -------------------
log.info("Aggregating %d genes × cell_type × disease_state...", len(expr))
group_keys = ["cell_type", "disease_state"]
groups = obs.groupby(group_keys, observed=True).indices
gene_names = list(expr.keys())
# Stack expression columns once for speed
expr_mat = np.column_stack([expr[g] for g in gene_names])  # cells × genes

records = []
for (cell_type, disease_state), row_pos in tqdm(groups.items(), desc="groups"):
    n_cells = len(row_pos)
    sub = expr_mat[row_pos, :]
    means = sub.mean(axis=0)
    medians = np.median(sub, axis=0)
    n_nz = (sub > 0).sum(axis=0)
    for k, gene in enumerate(gene_names):
        records.append({
            "gene": gene,
            "cell_type": cell_type,
            "disease_state": disease_state,
            "n_cells": int(n_cells),
            "mean": float(means[k]),
            "median": float(medians[k]),
            "pct_expressing": float(n_nz[k] / n_cells * 100),
            "n_nonzero": int(n_nz[k]),
        })
result = pd.DataFrame(records).sort_values(
    ["gene", "cell_type", "disease_state"]
).reset_index(drop=True)

out_path = PROCESSED / "hpap_T1D145_expression_by_celltype_disease.tsv"
result.to_csv(out_path, sep="\t", index=False)
log_artifact(out_path,
             source="Aggregated from per-gene CellxGene blobs (T1D 145-locus gene set)",
             notes=f"{len(result):,} rows × {result.gene.nunique()} genes")

# Also save the matched gene→locus mapping (which lead SNP / chrom / P each gene came from)
locus_map = target[["gene", "lead_snp", "CHR", "P", "n_studies", "var_idx"]].copy()
locus_map.to_csv(PROCESSED / "T1D145_gene_locus_map.tsv", sep="\t", index=False)

print(f"\n=== Full T1D-locus gene panel ===")
print(f"  Loci:    145")
print(f"  Unique genes from loci: {len(gene_best):,}")
print(f"  In HPAP var: {len(target):,}")
print(f"  Fetched OK:  {len(expr):,}")
if failures: print(f"  Failed:      {len(failures)}")
print(f"  Output rows: {len(result):,}  ->  {out_path}")

# --- 4. Quick scientific peek: which genes show the strongest beta-cell
#        T1D-vs-Control shift?
print("\n=== Beta cells: T1D vs Control log2((T1D+1)/(Ctrl+1)) ===")
beta = result[result["cell_type"] == "beta cell"]
piv = beta.pivot(index="gene", columns="disease_state", values="mean")
if {"T1D", "Control"}.issubset(piv.columns):
    piv["log2_T1D_Ctrl"] = np.log2((piv["T1D"] + 1) / (piv["Control"] + 1))
    print("\nTop 15 UP in T1D beta cells:")
    print(piv.nlargest(15, "log2_T1D_Ctrl")[["Control", "T1D", "log2_T1D_Ctrl"]]
              .to_string())
    print("\nTop 15 DOWN in T1D beta cells:")
    print(piv.nsmallest(15, "log2_T1D_Ctrl")[["Control", "T1D", "log2_T1D_Ctrl"]]
              .to_string())

# --- 5. Cell-type specificity index: where is each gene most expressed?
print("\n=== Per-gene dominant cell type (pooled across disease states) ===")
def pooled_mean(g):
    return (g["mean"] * g["n_cells"]).sum() / g["n_cells"].sum()

pooled = (result.groupby(["gene", "cell_type"])
                .apply(pooled_mean, include_groups=False)
                .rename("mean").reset_index())
dominant = (pooled.sort_values("mean", ascending=False)
                  .drop_duplicates("gene", keep="first"))
print("\nCell-type-of-action distribution (top compartment per gene):")
print(dominant["cell_type"].value_counts().to_string())
dominant.to_csv(PROCESSED / "T1D145_gene_dominant_celltype.tsv", sep="\t", index=False)
