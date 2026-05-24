"""09_leukocyte_stage: Leukocyte-compartment expression of immune-arm T1D loci by disease stage.

Extracted from notebook cell 36. Part of the t1d-celltype-of-action pipeline.
Run modules in numeric order — each consumes outputs from earlier steps.
"""

#         — testing whether the signal is insulitis-specific
# ============================================================
# Hypothesis: leukocyte-restricted T1D-locus genes (HLA-II, RUNX3,
# SH2B3, SLC15A3 — all called by Cell A) should be more highly
# expressed in T1D leukocytes than control leukocytes, reflecting
# autoimmune infiltrate (insulitis).
#
# Key design:
#   - Two contrasts:
#       T1D vs Control  -> "is there an immune signal in T1D?"
#       T1D vs T2D      -> "is it autoimmunity-specific, or just
#                          generic disease inflammation?"
#   - Positive controls: CD3D (T cells) and CD68 (macrophages).
#     If insulitis is the right story, CD3D should be elevated in
#     T1D leukocytes specifically (T cell infiltrate).
#   - Same donor-pseudobulk + LMM dual-test strategy as Cell B.
#   - Sample size warning: 229 T1D leukocytes, 6-9 donors. Small.
# ============================================================

import json, warnings, logging, subprocess, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
import httpx
from tqdm.auto import tqdm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", message=".*Maximum likelihood.*")
warnings.filterwarnings("ignore", message=".*did not converge.*")
warnings.filterwarnings("ignore", message="Unverified HTTPS request")

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S", force=True)
log = logging.getLogger("leuko")

PROJECT = Path("/kaggle/working/t1d_mech")
RAW = PROJECT / "data" / "raw"
INTERIM = PROJECT / "data" / "interim"
PROCESSED = PROJECT / "data" / "processed"
RESULTS = PROJECT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)
EXPR_CACHE_DIR = RAW / "cellxgene_expr"

CXG_BASE = "https://faryabi16.pmacs.upenn.edu/view/T1D_T2D_public.h5ad"
CXG_API = f"{CXG_BASE}/api/v0.2"

# --- 1. Gene panel: 7 leukocyte calls from Cell A + immune marker controls
LEUKOCYTE_T1D_GENES = ["RUNX3", "HLA-DQA1", "SLC15A3", "HLA-DRB5",
                       "HLA-DQB1", "HLA-DRB1", "SH2B3"]
POSITIVE_CONTROL_MARKERS = {
    "CD3D":  "T cell receptor complex (T cell marker)",
    "CD3E":  "T cell receptor complex (T cell marker)",
    "CD8A":  "CD8+ cytotoxic T cells",
    "CD4":   "CD4+ helper T cells",
    "CD68":  "macrophage marker",
    "CD19":  "B cell marker",
    "NKG7":  "NK cell / cytotoxic granule",
    "IL2RA": "T1D risk locus + activated T cell marker",
    "PTPN22":"T1D risk locus + T/B cell signaling",
    "CTLA4": "T1D risk locus + T cell checkpoint",
}
ALL_GENES = LEUKOCYTE_T1D_GENES + list(POSITIVE_CONTROL_MARKERS.keys())
print(f"Testing {len(ALL_GENES)} genes:")
print(f"  T1D leukocyte calls: {LEUKOCYTE_T1D_GENES}")
print(f"  Immune marker positive controls: {list(POSITIVE_CONTROL_MARKERS.keys())}")

# --- 2. Set up: obs, decoder, gene index ---------------------------------
obs = pd.read_csv(PROCESSED / "hpap_cellxgene_obs.tsv", sep="\t")
assert len(obs) == 222077

try:
    from server.common.fbs.matrix import decode_matrix_fbs
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q",
                           "--no-deps", "cellxgene", "flatbuffers"])
    from server.common.fbs.matrix import decode_matrix_fbs

gene_list = json.loads((INTERIM / "hpap_gene_list.json").read_text())
gene_to_idx = {g: i for i, g in enumerate(gene_list)}
missing = [g for g in ALL_GENES if g not in gene_to_idx]
if missing:
    log.warning("Genes missing from HPAP var (dropped): %s", missing)
ALL_GENES = [g for g in ALL_GENES if g in gene_to_idx]
log.info("Gene panel size after var filter: %d", len(ALL_GENES))

# --- 3. Fetch any genes not yet cached -----------------------------------
client = httpx.Client(timeout=300.0, verify=False, follow_redirects=True)
cache_files = {p.stem.split("_idx")[0]: p
               for p in EXPR_CACHE_DIR.glob("*_idx*.fbs")}

for g in tqdm(ALL_GENES, desc="ensure cache"):
    if g in cache_files:
        continue
    idx = gene_to_idx[g]
    r = client.put(
        f"{CXG_API}/data/var",
        json={"filter": {"var": {"index": [idx]}}},
        headers={"Accept": "application/octet-stream"},
    )
    r.raise_for_status()
    out = EXPR_CACHE_DIR / f"{g}_idx{idx}.fbs"
    out.write_bytes(r.content)
    cache_files[g] = out

# --- 4. Subset to leukocytes; decode all genes for those cells ----------
leuko_mask = (obs["cell_type"] == "leukocyte").to_numpy()
leuko_obs = obs[leuko_mask].reset_index(drop=True)
log.info("Leukocyte cells: %d", len(leuko_obs))
log.info("By disease state:\n%s",
         leuko_obs["disease_state"].value_counts().to_string())
log.info("By donor × disease state (leukocyte counts):")
donor_ct = pd.crosstab(leuko_obs["donor_id"], leuko_obs["disease_state"])
# Only print donors with any leukocytes
donor_ct = donor_ct[donor_ct.sum(axis=1) > 0]
print(donor_ct.to_string())

log.info("Decoding %d genes for leukocyte cells...", len(ALL_GENES))
expr_leuko = {}
for g in tqdm(ALL_GENES, desc="decoding"):
    blob = cache_files[g].read_bytes()
    df = decode_matrix_fbs(blob)
    vec = df.iloc[:, 0].to_numpy()
    expr_leuko[g] = vec[leuko_mask]
expr_leuko_df = pd.DataFrame(expr_leuko)
log.info("Leukocyte expression matrix: %d cells × %d genes", *expr_leuko_df.shape)

# Stack with obs
leuko_full = pd.concat(
    [leuko_obs[["donor_id", "disease_state", "sex", "assay"]].reset_index(drop=True),
     expr_leuko_df.reset_index(drop=True)],
    axis=1,
)

# --- 5. Donor pseudobulk + LMM dual-test machinery -----------------------
def donor_pseudobulk_test(df, gene, state_a, state_b):
    sub = df[df["disease_state"].isin([state_a, state_b])][
        ["donor_id", "disease_state", gene]].copy()
    sub["log_expr"] = np.log1p(sub[gene])
    donor_means = (sub.groupby(["donor_id", "disease_state"], observed=True)["log_expr"]
                      .mean().reset_index())
    a = donor_means.loc[donor_means["disease_state"] == state_a, "log_expr"].to_numpy()
    b = donor_means.loc[donor_means["disease_state"] == state_b, "log_expr"].to_numpy()
    if len(a) < 2 or len(b) < 2:
        return {"pb_effect": np.nan, "pb_p": np.nan,
                "n_donors_a": len(a), "n_donors_b": len(b),
                "mean_a": float(a.mean()) if len(a) else np.nan,
                "mean_b": float(b.mean()) if len(b) else np.nan}
    t = stats.ttest_ind(a, b, equal_var=False)
    return {
        "pb_effect": float(a.mean() - b.mean()),
        "pb_p": float(t.pvalue),
        "n_donors_a": len(a), "n_donors_b": len(b),
        "mean_a": float(a.mean()), "mean_b": float(b.mean()),
    }

def mixed_model_test(df, gene, state_a, state_b):
    sub = df[df["disease_state"].isin([state_a, state_b])][
        ["donor_id", "disease_state", gene]].copy()
    sub["log_expr"] = np.log1p(sub[gene])
    sub["donor_id"] = sub["donor_id"].astype(str)
    sub["disease_state"] = pd.Categorical(
        sub["disease_state"], categories=[state_b, state_a], ordered=False)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fit = smf.mixedlm("log_expr ~ disease_state",
                              data=sub, groups=sub["donor_id"]
                              ).fit(method="lbfgs", reml=True)
        coef_name = f"disease_state[T.{state_a}]"
        beta = float(fit.params.get(coef_name, np.nan))
        se = float(fit.bse.get(coef_name, np.nan))
        # Sanity flag: SE > 10x |beta| or > 1.0 means model is unstable
        unstable = (np.isnan(se) or se > 1.0 or
                    (np.isfinite(beta) and abs(beta) > 0 and se / max(abs(beta), 1e-9) > 10))
        return {"lmm_beta": beta, "lmm_se": se,
                "lmm_p": float(fit.pvalues.get(coef_name, np.nan)),
                "lmm_converged": bool(fit.converged),
                "lmm_unstable": unstable}
    except Exception as e:
        return {"lmm_beta": np.nan, "lmm_se": np.nan, "lmm_p": np.nan,
                "lmm_converged": False, "lmm_unstable": True,
                "lmm_error": str(e)[:80]}

# --- 6. Run all contrasts ------------------------------------------------
contrasts = [
    ("T1D", "Control"),    # main hypothesis: insulitis vs healthy
    ("T1D", "T2D"),        # specificity check: autoimmune vs generic disease
    ("AAB", "Control"),    # pre-clinical: insulitis present histologically
    ("T2D", "Control"),    # nuisance check: does generic disease elevate immune genes?
]

results = []
for gene in tqdm(ALL_GENES, desc="testing"):
    for state_a, state_b in contrasts:
        is_t1d_locus_call = gene in LEUKOCYTE_T1D_GENES
        is_positive_control = gene in POSITIVE_CONTROL_MARKERS
        row = {
            "gene": gene,
            "category": ("T1D locus (Cell A leukocyte call)" if is_t1d_locus_call
                        else "Immune marker control"),
            "contrast": f"{state_a} vs {state_b}",
        }
        row.update(donor_pseudobulk_test(leuko_full, gene, state_a, state_b))
        row.update(mixed_model_test(leuko_full, gene, state_a, state_b))
        results.append(row)

res = pd.DataFrame(results)

# BH-FDR within each contrast, separately within each category
res["pb_fdr"] = np.nan
res["lmm_fdr"] = np.nan
for c in res["contrast"].unique():
    for cat in res["category"].unique():
        mask = (res["contrast"] == c) & (res["category"] == cat)
        for col_p, col_q in [("pb_p", "pb_fdr"), ("lmm_p", "lmm_fdr")]:
            p = res.loc[mask, col_p].to_numpy()
            valid = np.isfinite(p)
            if valid.any():
                adj = np.full_like(p, np.nan, dtype=float)
                adj[valid] = multipletests(p[valid], method="fdr_bh")[1]
                res.loc[mask, col_q] = adj

# Save
out_path = PROCESSED / "leukocyte_immune_tests.tsv"
res.to_csv(out_path, sep="\t", index=False)
print(f"\nSaved: {out_path}")

# --- 7. Reporting --------------------------------------------------------
def report(contrast, panel_genes, panel_name):
    sub = (res[(res["contrast"] == contrast) & (res["gene"].isin(panel_genes))]
           .sort_values("pb_p"))
    if not len(sub):
        return
    print(f"\n--- {contrast} | {panel_name} ---")
    cols = ["gene", "pb_effect", "pb_p", "pb_fdr",
            "lmm_beta", "lmm_se", "lmm_p", "lmm_fdr",
            "n_donors_a", "n_donors_b", "lmm_unstable"]
    disp = sub[cols].rename(columns={
        "pb_effect": "PB_Δlog1p", "pb_p": "PB_p", "pb_fdr": "PB_FDR",
        "lmm_beta": "LMM_β", "lmm_se": "LMM_SE", "lmm_p": "LMM_p", "lmm_fdr": "LMM_FDR",
        "n_donors_a": "n_a", "n_donors_b": "n_b", "lmm_unstable": "unstable"})
    print(disp.round(4).to_string(index=False))

print("\n" + "=" * 72)
print("RESULTS: T1D vs Control leukocytes")
print("=" * 72)
report("T1D vs Control", LEUKOCYTE_T1D_GENES, "T1D loci (Cell A calls)")
report("T1D vs Control", list(POSITIVE_CONTROL_MARKERS), "Immune marker controls")

print("\n" + "=" * 72)
print("RESULTS: T1D vs T2D leukocytes (autoimmunity-specific?)")
print("=" * 72)
report("T1D vs T2D", LEUKOCYTE_T1D_GENES, "T1D loci (Cell A calls)")
report("T1D vs T2D", list(POSITIVE_CONTROL_MARKERS), "Immune marker controls")

print("\n" + "=" * 72)
print("RESULTS: AAB vs Control leukocytes (pre-clinical insulitis?)")
print("=" * 72)
report("AAB vs Control", LEUKOCYTE_T1D_GENES, "T1D loci (Cell A calls)")

print("\n" + "=" * 72)
print("RESULTS: T2D vs Control leukocytes (nuisance / generic disease)")
print("=" * 72)
report("T2D vs Control", LEUKOCYTE_T1D_GENES, "T1D loci (Cell A calls)")

# --- 8. The headline figure: effect-size dot plot -----------------------
# For each gene × contrast: plot pseudobulk effect size with 95% CI from
# a bootstrap over donors (more robust than LMM SE here).
def bootstrap_ci(df, gene, state_a, state_b, n_boot=2000, seed=42):
    sub = df[df["disease_state"].isin([state_a, state_b])][
        ["donor_id", "disease_state", gene]].copy()
    sub["log_expr"] = np.log1p(sub[gene])
    donor_means = (sub.groupby(["donor_id", "disease_state"], observed=True)["log_expr"]
                      .mean().reset_index())
    a = donor_means.loc[donor_means["disease_state"] == state_a, "log_expr"].to_numpy()
    b = donor_means.loc[donor_means["disease_state"] == state_b, "log_expr"].to_numpy()
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(seed)
    boot_effects = np.empty(n_boot)
    for i in range(n_boot):
        ai = rng.choice(a, size=len(a), replace=True)
        bi = rng.choice(b, size=len(b), replace=True)
        boot_effects[i] = ai.mean() - bi.mean()
    return float(a.mean() - b.mean()), float(np.percentile(boot_effects, 2.5)), \
           float(np.percentile(boot_effects, 97.5))

# Compute bootstrap CIs for the headline plot
plot_rows = []
for gene in ALL_GENES:
    for state_a, state_b in [("T1D", "Control"), ("T1D", "T2D"), ("AAB", "Control")]:
        eff, lo, hi = bootstrap_ci(leuko_full, gene, state_a, state_b)
        plot_rows.append({
            "gene": gene,
            "contrast": f"{state_a} vs {state_b}",
            "effect": eff, "ci_lo": lo, "ci_hi": hi,
            "category": ("T1D locus" if gene in LEUKOCYTE_T1D_GENES
                        else "Immune marker"),
        })
plot_df = pd.DataFrame(plot_rows)

# Order genes: T1D loci first (alphabetical), then immune markers (alphabetical)
gene_order = (sorted([g for g in ALL_GENES if g in LEUKOCYTE_T1D_GENES])
              + sorted([g for g in ALL_GENES if g not in LEUKOCYTE_T1D_GENES]))

fig, axes = plt.subplots(1, 3, figsize=(15, 8), sharey=True)
contrast_titles = {
    "T1D vs Control": "T1D vs Control\n(does autoimmunity show up?)",
    "T1D vs T2D":     "T1D vs T2D\n(is it autoimmune-specific?)",
    "AAB vs Control": "AAB vs Control\n(pre-clinical signal?)",
}
for ax, contrast in zip(axes, ["T1D vs Control", "T1D vs T2D", "AAB vs Control"]):
    sub = plot_df[plot_df["contrast"] == contrast].set_index("gene").reindex(gene_order)
    y = np.arange(len(sub))
    colors = ["#d62728" if cat == "T1D locus" else "#666666"
              for cat in sub["category"]]
    # Error bars
    err_lo = sub["effect"] - sub["ci_lo"]
    err_hi = sub["ci_hi"] - sub["effect"]
    ax.errorbar(sub["effect"], y, xerr=[err_lo, err_hi],
                fmt="o", color="black", markerfacecolor="white",
                markersize=0, capsize=3, linewidth=1, zorder=2)
    ax.scatter(sub["effect"], y, c=colors, s=70, zorder=3,
               edgecolor="black", linewidth=0.5)
    ax.axvline(0, color="gray", linestyle="--", linewidth=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(sub.index)
    ax.set_xlabel("Δ log1p mean (donor-pseudobulk)")
    ax.set_title(contrast_titles[contrast])
    ax.spines[["top", "right"]].set_visible(False)
    # Divider between T1D loci and immune markers
    boundary = sum(1 for g in gene_order if g in LEUKOCYTE_T1D_GENES) - 0.5
    ax.axhline(boundary, color="black", linewidth=0.6, alpha=0.3)

# Legend
red_p = mpatches.Patch(color="#d62728", label="T1D locus (Cell A leukocyte call)")
gray_p = mpatches.Patch(color="#666666", label="Immune marker positive control")
axes[0].legend(handles=[red_p, gray_p], loc="lower left", fontsize=9, frameon=True)
plt.suptitle("Leukocyte-compartment expression of immune-arm T1D loci, by disease state",
             fontsize=13)
plt.tight_layout()
plot_path = RESULTS / "leukocyte_immune_t1d_loci_by_disease.png"
plt.savefig(plot_path, dpi=150, bbox_inches="tight")
plt.show()
print(f"\nSaved: {plot_path}")

# --- 9. Summary for the writeup ------------------------------------------
print("\n" + "=" * 72)
print("WRITEUP-READY SUMMARY")
print("=" * 72)

def summarize_contrast(contrast):
    sub_loci = plot_df[(plot_df["contrast"] == contrast) &
                       (plot_df["category"] == "T1D locus")].dropna(subset=["effect"])
    sub_ctrl = plot_df[(plot_df["contrast"] == contrast) &
                       (plot_df["category"] == "Immune marker")].dropna(subset=["effect"])
    n_up_loci = (sub_loci["effect"] > 0).sum()
    n_up_ci_excl_zero = ((sub_loci["effect"] > 0) & (sub_loci["ci_lo"] > 0)).sum()
    n_up_ctrl = (sub_ctrl["effect"] > 0).sum()
    n_up_ctrl_ci_excl_zero = ((sub_ctrl["effect"] > 0) & (sub_ctrl["ci_lo"] > 0)).sum()
    print(f"\n{contrast}:")
    print(f"  T1D locus genes (n={len(sub_loci)}):")
    print(f"    upregulated direction: {n_up_loci}/{len(sub_loci)}")
    print(f"    bootstrap 95% CI excludes zero (up): {n_up_ci_excl_zero}/{len(sub_loci)}")
    print(f"  Immune marker controls (n={len(sub_ctrl)}):")
    print(f"    upregulated direction: {n_up_ctrl}/{len(sub_ctrl)}")
    print(f"    bootstrap 95% CI excludes zero (up): {n_up_ctrl_ci_excl_zero}/{len(sub_ctrl)}")

for c in ["T1D vs Control", "T1D vs T2D", "AAB vs Control", "T2D vs Control"]:
    summarize_contrast(c)

# Optional: a one-shot sign test on the T1D locus genes
# Null: 50/50 chance of being up vs down if no insulitis signal
print("\n" + "=" * 72)
print("BINOMIAL SIGN TEST: are T1D locus genes preferentially UP in T1D leukocytes?")
print("=" * 72)
for c in ["T1D vs Control", "T1D vs T2D", "AAB vs Control"]:
    sub = plot_df[(plot_df["contrast"] == c) &
                  (plot_df["category"] == "T1D locus")].dropna(subset=["effect"])
    n_up = int((sub["effect"] > 0).sum())
    n_tot = len(sub)
    p = stats.binomtest(n_up, n_tot, p=0.5, alternative="greater").pvalue
    print(f"  {c}: {n_up}/{n_tot} up, one-sided binomial p = {p:.3g}")
