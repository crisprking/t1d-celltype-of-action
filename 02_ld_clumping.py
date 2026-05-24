"""10_cross_validation: Full-panel cross-validation of tau-based cell-type calls (the headline result).

Extracted from notebook cell 38. Part of the t1d-celltype-of-action pipeline.
Run modules in numeric order — each consumes outputs from earlier steps.
"""

# Question: if our tau analysis (Cell A) correctly identified the
# cell-type-of-action for T1D-locus genes, then AAB-vs-Control effect
# sizes should be largest in the called cell type, not random ones.
#
# This is a stronger test than running stats on a few hand-picked
# genes. We pool information across all 150 genes:
#   - For each gene, get its Cell A cell-type call (leukocyte, beta,
#     duct, acinar, etc., or "broad" / "below detection")
#   - For each gene × cell-type combination, compute AAB-vs-Control
#     effect size (donor-pseudobulk Δlog1p)
#   - Ask: do the AAB effects in the CALLED cell type systematically
#     exceed those in OTHER cell types?
#
# Design choices:
#   - AAB vs Control as the primary contrast (we established in Cell C
#     that AAB is where active insulitis lives in this dataset)
#   - All cell types tested per gene (not just the called one) so we
#     have a within-gene "background" of effect sizes
#   - Permutation test for the overall hypothesis to avoid relying on
#     per-gene FDR (which is underpowered here)
# ============================================================

import json, warnings, logging, subprocess, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from tqdm.auto import tqdm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", message="Unverified HTTPS request")

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S", force=True)
log = logging.getLogger("xval")

PROJECT = Path("/kaggle/working/t1d_mech")
RAW = PROJECT / "data" / "raw"
INTERIM = PROJECT / "data" / "interim"
PROCESSED = PROJECT / "data" / "processed"
RESULTS = PROJECT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)
EXPR_CACHE_DIR = RAW / "cellxgene_expr"

# --- 1. Load obs, specificity calls, gene index --------------------------
obs = pd.read_csv(PROCESSED / "hpap_cellxgene_obs.tsv", sep="\t")
spec = pd.read_csv(PROCESSED / "T1D145_celltype_specificity.tsv", sep="\t")
gene_list = json.loads((INTERIM / "hpap_gene_list.json").read_text())
gene_to_idx = {g: i for i, g in enumerate(gene_list)}

# Use ALL 150 genes that were successfully fetched, regardless of tau
cache_files = {p.stem.split("_idx")[0]: p
               for p in EXPR_CACHE_DIR.glob("*_idx*.fbs")}
panel_genes = [g for g in spec["gene"] if g in cache_files]
log.info("Panel size: %d genes (cache available)", len(panel_genes))

# --- 2. Decode all 150 genes once, into a sparse-ish per-gene dict ------
try:
    from server.common.fbs.matrix import decode_matrix_fbs
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q",
                           "--no-deps", "cellxgene", "flatbuffers"])
    from server.common.fbs.matrix import decode_matrix_fbs

log.info("Decoding %d FBS blobs (one-time cost)...", len(panel_genes))
expr_all = {}
for g in tqdm(panel_genes, desc="decoding"):
    blob = cache_files[g].read_bytes()
    df = decode_matrix_fbs(blob)
    expr_all[g] = df.iloc[:, 0].to_numpy()
log.info("All gene vectors loaded: %d genes × 222,077 cells", len(expr_all))

# --- 3. Donor-pseudobulk effect-size machinery ---------------------------
def donor_effect(expr_vec, mask_subset, obs_subset, state_a, state_b):
    """
    Δ log1p donor-mean expression (state_a - state_b) within a cell-type subset.
    Returns NaN if either group has < 2 donors with cells.
    """
    if mask_subset.sum() == 0:
        return np.nan
    sub_obs = obs_subset.loc[mask_subset, ["donor_id", "disease_state"]].copy()
    sub_obs["expr"] = expr_vec[mask_subset.to_numpy()]
    sub_obs = sub_obs[sub_obs["disease_state"].isin([state_a, state_b])]
    if len(sub_obs) == 0:
        return np.nan
    sub_obs["log_expr"] = np.log1p(sub_obs["expr"])
    donor_means = (sub_obs.groupby(["donor_id", "disease_state"], observed=True)["log_expr"]
                          .mean().reset_index())
    a = donor_means.loc[donor_means["disease_state"] == state_a, "log_expr"].to_numpy()
    b = donor_means.loc[donor_means["disease_state"] == state_b, "log_expr"].to_numpy()
    if len(a) < 2 or len(b) < 2:
        return np.nan
    return float(a.mean() - b.mean())

# --- 4. For each gene, compute Δlog1p in every cell type × AAB vs Control --
cell_types = sorted(obs["cell_type"].unique())
log.info("Computing effect sizes: %d genes × %d cell types", len(panel_genes), len(cell_types))

# Pre-compute cell-type masks once
ct_masks = {ct: (obs["cell_type"] == ct) for ct in cell_types}

effects = []  # rows: gene, cell_type, effect
for g in tqdm(panel_genes, desc="effects"):
    vec = expr_all[g]
    for ct in cell_types:
        eff = donor_effect(vec, ct_masks[ct], obs, "AAB", "Control")
        effects.append({"gene": g, "cell_type": ct, "AAB_vs_Control": eff})
eff_df = pd.DataFrame(effects)

# Join with the Cell A call
eff_df = eff_df.merge(
    spec[["gene", "tau", "top_celltype", "above_floor", "call"]],
    on="gene", how="left",
)

# --- 5. Stratify by whether this is the gene's called cell type ----------
eff_df["is_called_ct"] = eff_df["cell_type"] == eff_df["top_celltype"]
# Apply same filters as Cell A: tau ≥ 0.70 AND above expression floor
eff_df["confident_call"] = (eff_df["tau"] >= 0.70) & eff_df["above_floor"]

# Save the full effect table
out_path = PROCESSED / "T1D145_full_AAB_effects_by_celltype.tsv"
eff_df.to_csv(out_path, sep="\t", index=False)
print(f"Saved: {out_path}")

# --- 6. The cross-validation test ---------------------------------------
# For genes with confident Cell A calls, is the AAB effect in the called
# cell type systematically larger than in non-called cell types?

confident = eff_df[eff_df["confident_call"]].copy()
n_confident_genes = confident["gene"].nunique()
log.info("Genes with confident Cell A calls: %d", n_confident_genes)

# Per-gene "called-vs-other" delta: effect in called CT minus mean effect across non-called CTs
def signed_called_advantage(group):
    if group["is_called_ct"].sum() != 1:
        return np.nan
    called_eff = group.loc[group["is_called_ct"], "AAB_vs_Control"].iloc[0]
    other_effs = group.loc[~group["is_called_ct"], "AAB_vs_Control"].dropna()
    if pd.isna(called_eff) or len(other_effs) == 0:
        return np.nan
    return called_eff - other_effs.mean()

per_gene_advantage = (confident.groupby("gene")
                              .apply(signed_called_advantage, include_groups=False)
                              .dropna())
log.info("Per-gene 'called CT advantage' (n=%d genes):", len(per_gene_advantage))
print(per_gene_advantage.sort_values(ascending=False).to_string())

# Hypothesis: mean advantage > 0 (called CT systematically higher AAB effect)
# Permutation test: shuffle which cell type is "called" within each gene,
# recompute mean advantage; ask how often shuffled exceeds observed.
observed = per_gene_advantage.mean()
log.info("\nObserved mean called-CT advantage: %+.4f", observed)

rng = np.random.default_rng(42)
n_perm = 10_000
perm_means = np.empty(n_perm)
genes_in_test = per_gene_advantage.index.tolist()
gene_groups = {g: eff_df[eff_df["gene"] == g].dropna(subset=["AAB_vs_Control"]).copy()
               for g in genes_in_test}
for i in range(n_perm):
    fake_advantages = []
    for g in genes_in_test:
        group = gene_groups[g]
        if len(group) < 2:
            continue
        # Randomly pick one cell type as "called"
        idx = rng.integers(0, len(group))
        called = group["AAB_vs_Control"].iloc[idx]
        others = group["AAB_vs_Control"].drop(group.index[idx])
        if len(others) == 0:
            continue
        fake_advantages.append(called - others.mean())
    perm_means[i] = np.mean(fake_advantages) if fake_advantages else np.nan

valid_perms = perm_means[np.isfinite(perm_means)]
p_one_sided = (np.sum(valid_perms >= observed) + 1) / (len(valid_perms) + 1)
p_two_sided = (np.sum(np.abs(valid_perms) >= abs(observed)) + 1) / (len(valid_perms) + 1)

print(f"\nPermutation test (n={n_perm:,}):")
print(f"  Observed mean called-CT advantage: {observed:+.4f}")
print(f"  Null distribution: mean={valid_perms.mean():+.4f}, "
      f"SD={valid_perms.std():.4f}")
print(f"  One-sided p (called > random): {p_one_sided:.4f}")
print(f"  Two-sided p:                   {p_two_sided:.4f}")

# --- 7. Stratified by the called cell type (which compartments work?) ---
print("\n" + "=" * 72)
print("Stratified by called cell type: AAB effect in called CT vs others")
print("=" * 72)
stratified = []
for ct in confident["top_celltype"].dropna().unique():
    ct_genes = confident[confident["top_celltype"] == ct]["gene"].unique()
    if len(ct_genes) < 2:
        continue  # need ≥2 genes per CT for a test
    advantages_this_ct = [per_gene_advantage[g] for g in ct_genes
                          if g in per_gene_advantage.index]
    if len(advantages_this_ct) < 2:
        continue
    stratified.append({
        "called_celltype": ct,
        "n_genes": len(ct_genes),
        "n_genes_tested": len(advantages_this_ct),
        "mean_advantage": float(np.mean(advantages_this_ct)),
        "median_advantage": float(np.median(advantages_this_ct)),
    })
strat_df = pd.DataFrame(stratified).sort_values("mean_advantage", ascending=False)
print(strat_df.round(4).to_string(index=False))

# --- 8. Visualization: per-gene effects in their called CT vs alternatives -
fig, axes = plt.subplots(1, 2, figsize=(14, 7))

# Panel A: per-gene called-CT advantage histogram
ax = axes[0]
ax.hist(per_gene_advantage.values, bins=20, color="#3a7ca5",
        edgecolor="white", alpha=0.85)
ax.axvline(0, color="black", linewidth=0.8)
ax.axvline(observed, color="red", linewidth=1.5, linestyle="--",
           label=f"observed mean = {observed:+.3f}")
ax.set_xlabel("called-CT AAB effect − mean(other-CT AAB effects), per gene")
ax.set_ylabel("number of T1D-locus genes")
ax.set_title(f"Cross-validation of cell-type-of-action calls\n"
             f"({n_confident_genes} genes with tau ≥ 0.70, AAB vs Control)")
ax.legend()
ax.spines[["top", "right"]].set_visible(False)

# Panel B: scatter of called-CT effect vs max-other-CT effect, per gene
ax = axes[1]
plot_rows = []
for g in genes_in_test:
    group = eff_df[eff_df["gene"] == g]
    called = group.loc[group["is_called_ct"], "AAB_vs_Control"]
    if called.empty or called.isna().all():
        continue
    others = group.loc[~group["is_called_ct"], "AAB_vs_Control"].dropna()
    if others.empty:
        continue
    plot_rows.append({
        "gene": g,
        "called_eff": float(called.iloc[0]),
        "max_other_eff": float(others.abs().max() * np.sign(others.iloc[others.abs().argmax()])),
        "called_ct": group["top_celltype"].iloc[0],
    })
plot_df = pd.DataFrame(plot_rows)

# Color by called cell type
compartment_colors = {
    "beta cell": "#2ca02c",
    "alpha cell": "#9467bd",
    "delta cell": "#8c564b",
    "leukocyte": "#d62728",
    "duct epithelial cell": "#7f7f7f",
    "endothelial cell": "#17becf",
    "mesenchymal cell": "#bcbd22",
    "acinar cell": "#ff7f0e",
    "PP cell": "#e377c2",
    "epsilon cell": "#1f77b4",
}
for ct in plot_df["called_ct"].unique():
    sub = plot_df[plot_df["called_ct"] == ct]
    ax.scatter(sub["max_other_eff"], sub["called_eff"],
               s=60, color=compartment_colors.get(ct, "#999"),
               edgecolor="black", linewidth=0.5, alpha=0.85,
               label=f"{ct} (n={len(sub)})")

# y=x reference
lim = max(plot_df["called_eff"].abs().max(),
          plot_df["max_other_eff"].abs().max()) * 1.15
ax.plot([-lim, lim], [-lim, lim], color="gray", linestyle="--",
        linewidth=0.5, alpha=0.7, zorder=1)
ax.axhline(0, color="black", linewidth=0.3)
ax.axvline(0, color="black", linewidth=0.3)
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
ax.set_xlabel("largest AAB effect in any OTHER cell type")
ax.set_ylabel("AAB effect in the called cell type")
ax.set_title("Per-gene AAB effect: called cell type vs alternatives\n"
             "(above y=x line = called CT is stronger)")
ax.legend(loc="lower right", fontsize=8, frameon=True, ncol=1)
ax.spines[["top", "right"]].set_visible(False)

# Label outlier genes
for _, r in plot_df.iterrows():
    if abs(r["called_eff"]) > 0.3 or abs(r["max_other_eff"]) > 0.3:
        ax.annotate(r["gene"], (r["max_other_eff"], r["called_eff"]),
                    xytext=(3, 3), textcoords="offset points", fontsize=7)

plt.tight_layout()
plot_path = RESULTS / "T1D145_celltype_call_crossvalidation.png"
plt.savefig(plot_path, dpi=150, bbox_inches="tight")
plt.show()
print(f"\nSaved: {plot_path}")

# --- 9. Headline summary for the writeup --------------------------------
print("\n" + "=" * 72)
print("WRITEUP-READY SUMMARY")
print("=" * 72)
n_called_higher = (per_gene_advantage > 0).sum()
n_total = len(per_gene_advantage)
print(f"\nFor {n_total} T1D-locus genes with confident cell-type-of-action calls:")
print(f"  - {n_called_higher} ({n_called_higher/n_total:.0%}) have larger AAB effect "
      f"in the called cell type than the mean across other cell types")
print(f"  - Mean per-gene called-CT advantage: {observed:+.4f} log1p")
print(f"  - Permutation p (vs random cell-type assignment): "
      f"one-sided {p_one_sided:.4f}, two-sided {p_two_sided:.4f}")
print(f"\nInterpretation:")
if p_one_sided < 0.05:
    print("  Cell-type-of-action calls from tau (Cell A) cross-validate: AAB-stage")
    print("  expression changes systematically concentrate in the called cell types.")
elif observed > 0 and p_one_sided < 0.15:
    print("  Trend toward cross-validation but underpowered; effect direction is")
    print("  consistent with cell-type calls being biologically meaningful.")
else:
    print("  Tau-based cell-type calls do NOT predict where AAB effects show up.")
    print("  This could mean: (1) AAB effects are diffuse / not cell-type-localized,")
    print("  (2) the calls capture baseline specificity but not disease activation,")
    print("  or (3) the AAB cohort is too small to detect localized effects.")

# Show the top 10 individual hits
print(f"\nTop 10 genes by called-CT AAB advantage:")
top10 = per_gene_advantage.sort_values(ascending=False).head(10)
for g in top10.index:
    ct = spec.loc[spec["gene"] == g, "top_celltype"].iloc[0]
    print(f"  {g:15s}  called CT: {ct:25s}  advantage = {top10[g]:+.4f}")
