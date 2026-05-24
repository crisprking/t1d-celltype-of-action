"""08_tau_specificity: Compute tau cell-type specificity + cell-type-of-action calls.

Extracted from notebook cell 35. Part of the t1d-celltype-of-action pipeline.
Run modules in numeric order — each consumes outputs from earlier steps.
"""

# Replaces the previous "top cell type by raw mean" call, which rewarded
# housekeeping noise and was inflated by ambient mRNA contamination in
# duct cells. Tau is the standard tissue-specificity index from Yanai
# et al. 2005, redefined here on cell types rather than tissues.
#
# tau = sum_i (1 - x_i_hat) / (N - 1)
#   where x_i_hat = x_i / max(x), N = number of cell types
#
# tau = 1.0 -> expressed in exactly one cell type
# tau = 0.0 -> uniformly expressed across all cell types
#
# We pool across disease states for the specificity call, because we
# want the constitutive cell-type-of-action; disease-state effects are
# handled separately in cell B.

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

PROJECT = Path("/kaggle/working/t1d_mech")
PROCESSED = PROJECT / "data" / "processed"
RESULTS = PROJECT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

# --- 1. Load the full 145-locus aggregated table -------------------------
agg = pd.read_csv(PROCESSED / "hpap_T1D145_expression_by_celltype_disease.tsv",
                  sep="\t")
print(f"Loaded: {len(agg):,} rows, {agg['gene'].nunique()} genes")

# --- 2. Collapse to gene × cell_type (weighted by n_cells) ---------------
# We weight by n_cells so disease states with few cells (e.g. T1D beta
# cells, n=715) don't dominate the per-cell-type mean.
def weighted_mean(g):
    return (g["mean"] * g["n_cells"]).sum() / g["n_cells"].sum()

mat = (agg.groupby(["gene", "cell_type"])
          .apply(weighted_mean, include_groups=False)
          .unstack("cell_type")
          .fillna(0))
print(f"Gene × cell_type matrix: {mat.shape}")

# Log1p before computing tau. Raw counts have huge dynamic range (INS
# is in thousands, most loci are < 1), so tau on raw values is dominated
# by which cell type has the largest absolute count, not the most
# specific expression pattern.
mat_log = np.log1p(mat)

# --- 3. Compute tau per gene ---------------------------------------------
def tau(row):
    """Tau specificity index. NaN if all-zero (no expression anywhere)."""
    x = row.to_numpy(dtype=float)
    xmax = x.max()
    if xmax == 0:
        return np.nan
    x_hat = x / xmax
    n = len(x)
    return (1.0 - x_hat).sum() / (n - 1)

tau_series = mat_log.apply(tau, axis=1)

# Best cell type = the one with the highest log1p mean (only meaningful
# if tau is reasonably high; we'll annotate this below)
top_celltype = mat_log.idxmax(axis=1)
top_value = mat_log.max(axis=1)

spec = pd.DataFrame({
    "gene": tau_series.index,
    "tau": tau_series.values,
    "top_celltype": top_celltype.values,
    "top_celltype_log1p_mean": top_value.values,
    "any_expression": mat_log.max(axis=1).values > 0,
}).sort_values("tau", ascending=False).reset_index(drop=True)

# --- 4. Bin tau into interpretable buckets -------------------------------
def tau_bucket(t):
    if pd.isna(t):                  return "no expression"
    if t >= 0.85:                   return "highly specific (≥0.85)"
    if t >= 0.70:                   return "specific (0.70–0.85)"
    if t >= 0.50:                   return "moderate (0.50–0.70)"
    return "broad / ubiquitous (<0.50)"

spec["bucket"] = spec["tau"].apply(tau_bucket)

# --- 5. Apply an expression floor ----------------------------------------
# A gene with max log1p of 0.01 across all cell types is essentially
# noise. Specificity is only meaningful if the gene is actually expressed
# somewhere. Floor: top_celltype log1p mean > 0.5 (i.e. raw mean ~ 0.65+).
EXPR_FLOOR = 0.5
spec["above_floor"] = spec["top_celltype_log1p_mean"] >= EXPR_FLOOR
spec["call"] = np.where(
    spec["above_floor"] & (spec["tau"] >= 0.70),
    spec["top_celltype"],
    np.where(spec["above_floor"], "broad / multi-compartment", "below detection"),
)

print(f"\n=== Tau distribution ===")
print(spec["bucket"].value_counts().to_string())

print(f"\n=== Cell-type-of-action calls (tau ≥ 0.70, above expression floor) ===")
print(spec["call"].value_counts().to_string())

# Save the full table
out_path = PROCESSED / "T1D145_celltype_specificity.tsv"
spec.to_csv(out_path, sep="\t", index=False)
print(f"\nSaved: {out_path}")

# --- 6. The headline table for the writeup -------------------------------
print(f"\n=== T1D loci with confident cell-type-of-action (top 30) ===")
confident = spec[(spec["tau"] >= 0.70) & (spec["above_floor"])].head(30)
print(confident[["gene", "tau", "top_celltype", "top_celltype_log1p_mean"]].to_string(index=False))

# --- 7. Plots -------------------------------------------------------------
# Plot 1: tau distribution histogram
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))

ax = axes[0]
plot_data = spec[spec["above_floor"]]["tau"].dropna()
ax.hist(plot_data, bins=25, color="#3a7ca5", edgecolor="white")
ax.axvline(0.70, color="red", linestyle="--", linewidth=1.2,
           label="tau = 0.70 (specificity threshold)")
ax.set_xlabel("tau (cell-type specificity)")
ax.set_ylabel("number of T1D-locus genes")
ax.set_title(f"Specificity of T1D-locus genes\n(n={len(plot_data)} above expression floor)")
ax.legend()
ax.spines[["top", "right"]].set_visible(False)

# Plot 2: cell-type-of-action call distribution
ax = axes[1]
call_counts = (spec[spec["above_floor"]]["call"]
               .value_counts()
               .sort_values(ascending=True))
colors = ["#888888" if c in ("broad / multi-compartment", "below detection")
          else "#d62728" if c == "leukocyte"
          else "#2ca02c" if c == "beta cell"
          else "#1f77b4" for c in call_counts.index]
ax.barh(call_counts.index, call_counts.values, color=colors, edgecolor="white")
ax.set_xlabel("number of T1D-locus genes")
ax.set_title("Cell-type-of-action (tau ≥ 0.70)")
ax.spines[["top", "right"]].set_visible(False)

plt.tight_layout()
plot_path = RESULTS / "T1D145_specificity_distribution.png"
plt.savefig(plot_path, dpi=150, bbox_inches="tight")
plt.show()
print(f"\nSaved: {plot_path}")

# Plot 3: focused look at high-confidence calls by compartment
fig, ax = plt.subplots(figsize=(10, max(6, len(confident) * 0.28)))
top_for_plot = confident.head(40).iloc[::-1]  # reverse so highest tau at top
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
bar_colors = [compartment_colors.get(ct, "#999999") for ct in top_for_plot["top_celltype"]]
ax.barh(top_for_plot["gene"], top_for_plot["tau"], color=bar_colors, edgecolor="white")
ax.set_xlabel("tau (cell-type specificity)")
ax.set_xlim(0.65, 1.0)
ax.set_title("T1D-locus genes with confident cell-type-of-action")
ax.spines[["top", "right"]].set_visible(False)

# Add cell-type legend
import matplotlib.patches as mpatches
present = top_for_plot["top_celltype"].unique()
legend_handles = [mpatches.Patch(color=compartment_colors.get(ct, "#999999"), label=ct)
                  for ct in present]
ax.legend(handles=legend_handles, bbox_to_anchor=(1.02, 1), loc="upper left",
          frameon=False, fontsize=9)

plt.tight_layout()
plot_path2 = RESULTS / "T1D145_topgenes_by_compartment.png"
plt.savefig(plot_path2, dpi=150, bbox_inches="tight")
plt.show()
print(f"Saved: {plot_path2}")

# --- 8. Compartment-stratified summary for the article ------------------
print(f"\n=== For the writeup: T1D loci by cell-type-of-action ===")
for ct in ["beta cell", "alpha cell", "delta cell", "leukocyte",
           "duct epithelial cell", "endothelial cell",
           "mesenchymal cell", "acinar cell", "epsilon cell"]:
    genes_in_ct = confident[confident["top_celltype"] == ct]["gene"].tolist()
    if genes_in_ct:
        print(f"  {ct} (n={len(genes_in_ct)}): {', '.join(genes_in_ct[:20])}"
              + ("..." if len(genes_in_ct) > 20 else ""))
