"""Step 09 — leukocyte-compartment expression of immune-arm T1D loci.

Tests whether the 7 leukocyte calls from step 08 (HLA class II + RUNX3 +
SH2B3 + SLC15A3) plus a panel of canonical immune markers show
disease-stage shifts consistent with insulitis. Three contrasts:
T1D vs Control, T1D vs T2D (autoimmune vs generic inflammation), and
AAB vs Control (pre-clinical).
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.multitest import multipletests
from tqdm.auto import tqdm

from t1d_coa.config import EXPECTED_N_CELLS, INTERIM, PROCESSED, RESULTS
from t1d_coa.hpap import fetch_gene_expression, fetch_gene_list


LEUKO_T1D_GENES = ("RUNX3", "HLA-DQA1", "SLC15A3", "HLA-DRB5",
                   "HLA-DQB1", "HLA-DRB1", "SH2B3")
IMMUNE_MARKERS = {
    "CD3D": "T cell receptor complex",
    "CD3E": "T cell receptor complex",
    "CD8A": "CD8+ cytotoxic T cells",
    "CD4": "CD4+ helper T cells",
    "CD68": "macrophage marker",
    "CD19": "B cell marker",
    "NKG7": "NK cell granule",
    "IL2RA": "activated T cells",
    "PTPN22": "T/B cell signaling",
    "CTLA4": "T cell checkpoint",
}
CONTRASTS = [("T1D", "Control"), ("T1D", "T2D"),
             ("AAB", "Control"), ("T2D", "Control")]


def pseudobulk_test(df: pd.DataFrame, gene: str, a: str, b: str) -> dict:
    sub = df[df["disease_state"].isin([a, b])][
        ["donor_id", "disease_state", gene]].copy()
    sub["log_expr"] = np.log1p(sub[gene])
    means = (sub.groupby(["donor_id", "disease_state"], observed=True)["log_expr"]
             .mean().reset_index())
    ga = means.loc[means["disease_state"] == a, "log_expr"].to_numpy()
    gb = means.loc[means["disease_state"] == b, "log_expr"].to_numpy()
    if len(ga) < 2 or len(gb) < 2:
        return {"pb_effect": np.nan, "pb_p": np.nan,
                "n_donors_a": len(ga), "n_donors_b": len(gb)}
    return {
        "pb_effect": float(ga.mean() - gb.mean()),
        "pb_p": float(stats.ttest_ind(ga, gb, equal_var=False).pvalue),
        "n_donors_a": len(ga), "n_donors_b": len(gb),
    }


def lmm_test(df: pd.DataFrame, gene: str, a: str, b: str) -> dict:
    sub = df[df["disease_state"].isin([a, b])][
        ["donor_id", "disease_state", gene]].copy()
    sub["log_expr"] = np.log1p(sub[gene])
    sub["donor_id"] = sub["donor_id"].astype(str)
    sub["disease_state"] = pd.Categorical(sub["disease_state"],
                                          categories=[b, a], ordered=False)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fit = smf.mixedlm("log_expr ~ disease_state", data=sub,
                              groups=sub["donor_id"]).fit(method="lbfgs", reml=True)
        coef = f"disease_state[T.{a}]"
        beta, se = float(fit.params.get(coef, np.nan)), float(fit.bse.get(coef, np.nan))
        unstable = (np.isnan(se) or se > 1.0 or
                    (abs(beta) > 0 and se / max(abs(beta), 1e-9) > 10))
        return {"lmm_beta": beta, "lmm_se": se,
                "lmm_p": float(fit.pvalues.get(coef, np.nan)),
                "lmm_unstable": unstable}
    except Exception:
        return {"lmm_beta": np.nan, "lmm_se": np.nan,
                "lmm_p": np.nan, "lmm_unstable": True}


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)

    obs = pd.read_csv(PROCESSED / "hpap_cellxgene_obs.tsv", sep="\t")
    assert len(obs) == EXPECTED_N_CELLS

    genes_all = fetch_gene_list(INTERIM / "hpap_gene_list.json")
    gene_to_idx = {g: i for i, g in enumerate(genes_all)}

    panel = [g for g in (*LEUKO_T1D_GENES, *IMMUNE_MARKERS) if g in gene_to_idx]
    print(f"Panel: {len(panel)} genes")

    leuko_mask = (obs["cell_type"] == "leukocyte").to_numpy()
    leuko_obs = obs[leuko_mask].reset_index(drop=True)
    print(f"Leukocyte cells: {len(leuko_obs)}")

    expr_leuko = {}
    for g in tqdm(panel, desc="decoding"):
        vec = fetch_gene_expression(g, gene_to_idx[g])
        expr_leuko[g] = vec[leuko_mask]

    leuko_full = pd.concat(
        [leuko_obs[["donor_id", "disease_state"]].reset_index(drop=True),
         pd.DataFrame(expr_leuko).reset_index(drop=True)],
        axis=1,
    )

    results = []
    for gene in tqdm(panel, desc="testing"):
        for a, b in CONTRASTS:
            row = {
                "gene": gene,
                "category": ("T1D locus" if gene in LEUKO_T1D_GENES
                             else "Immune marker"),
                "contrast": f"{a} vs {b}",
                **pseudobulk_test(leuko_full, gene, a, b),
                **lmm_test(leuko_full, gene, a, b),
            }
            results.append(row)
    res = pd.DataFrame(results)

    # BH-FDR within each contrast × category.
    for col_p, col_q in (("pb_p", "pb_fdr"), ("lmm_p", "lmm_fdr")):
        res[col_q] = np.nan
        for c in res["contrast"].unique():
            for cat in res["category"].unique():
                m = (res["contrast"] == c) & (res["category"] == cat)
                p = res.loc[m, col_p].to_numpy()
                valid = np.isfinite(p)
                if valid.any():
                    adj = np.full_like(p, np.nan, dtype=float)
                    adj[valid] = multipletests(p[valid], method="fdr_bh")[1]
                    res.loc[m, col_q] = adj

    out = PROCESSED / "leukocyte_immune_tests.tsv"
    res.to_csv(out, sep="\t", index=False)
    print(f"\nSaved: {out}")

    # Sign-test summary.
    print("\n=== Sign test: T1D-locus genes UP in T1D leukocytes? ===")
    for c in ("T1D vs Control", "T1D vs T2D", "AAB vs Control"):
        sub = res[(res["contrast"] == c) & (res["category"] == "T1D locus")].dropna(
            subset=["pb_effect"])
        n_up = int((sub["pb_effect"] > 0).sum())
        pval = stats.binomtest(n_up, len(sub), p=0.5,
                               alternative="greater").pvalue
        print(f"  {c}: {n_up}/{len(sub)} up, one-sided binomial p = {pval:.3g}")


if __name__ == "__main__":
    main()
