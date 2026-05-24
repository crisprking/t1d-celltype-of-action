"""Step 10 — cross-validate tau cell-type calls against AAB disease effects.

The headline test of the project. For each gene with a confident tau
call, ask: is the AAB-vs-Control effect in the called compartment
larger than the mean across other compartments? Permutation null:
shuffle which cell type is "called" within each gene, recompute the
mean advantage, repeat 10,000 times.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from t1d_coa.config import (
    EXPR_CACHE, INTERIM, N_PERMUTATIONS, PERM_SEED, PROCESSED, RESULTS,
    TAU_THRESHOLD,
)
from t1d_coa.hpap import _decoder
from t1d_coa.plotting import COMPARTMENT_COLORS, style_axes
from t1d_coa.stats import donor_pseudobulk_effect


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)

    obs = pd.read_csv(PROCESSED / "hpap_cellxgene_obs.tsv", sep="\t")
    spec = pd.read_csv(PROCESSED / "T1D145_celltype_specificity.tsv", sep="\t")

    cache_files = {p.stem.split("_idx")[0]: p
                   for p in EXPR_CACHE.glob("*_idx*.fbs")}
    panel = [g for g in spec["gene"] if g in cache_files]
    print(f"Panel size: {len(panel)} genes")

    decode = _decoder()
    expr = {g: decode(cache_files[g].read_bytes()).iloc[:, 0].to_numpy()
            for g in tqdm(panel, desc="decoding")}

    cell_types = sorted(obs["cell_type"].unique())
    ct_masks = {ct: (obs["cell_type"] == ct) for ct in cell_types}

    effects = []
    for g in tqdm(panel, desc="effects"):
        vec = expr[g]
        for ct in cell_types:
            effects.append({
                "gene": g, "cell_type": ct,
                "AAB_vs_Control": donor_pseudobulk_effect(
                    vec, ct_masks[ct], obs, "AAB", "Control"
                ),
            })
    eff = (pd.DataFrame(effects)
           .merge(spec[["gene", "tau", "top_celltype", "above_floor", "call"]],
                  on="gene", how="left"))
    eff["is_called_ct"] = eff["cell_type"] == eff["top_celltype"]
    eff["confident_call"] = (eff["tau"] >= TAU_THRESHOLD) & eff["above_floor"]

    out_eff = PROCESSED / "T1D145_full_AAB_effects_by_celltype.tsv"
    eff.to_csv(out_eff, sep="\t", index=False)
    print(f"Saved: {out_eff}")

    confident = eff[eff["confident_call"]]

    def called_advantage(group: pd.DataFrame) -> float:
        if group["is_called_ct"].sum() != 1:
            return np.nan
        called = group.loc[group["is_called_ct"], "AAB_vs_Control"].iloc[0]
        others = group.loc[~group["is_called_ct"], "AAB_vs_Control"].dropna()
        if pd.isna(called) or others.empty:
            return np.nan
        return called - others.mean()

    per_gene = (confident.groupby("gene")
                .apply(called_advantage, include_groups=False)
                .dropna())
    observed = per_gene.mean()
    print(f"\nGenes with confident calls: {len(per_gene)}")
    print(f"Observed mean called-CT advantage: {observed:+.4f}")

    # Permutation: shuffle which CT is "called" inside each gene.
    rng = np.random.default_rng(PERM_SEED)
    gene_groups = {g: eff[eff["gene"] == g].dropna(subset=["AAB_vs_Control"])
                   for g in per_gene.index}
    perm = np.empty(N_PERMUTATIONS)
    for i in range(N_PERMUTATIONS):
        fakes = []
        for g, gp in gene_groups.items():
            if len(gp) < 2:
                continue
            idx = rng.integers(0, len(gp))
            called = gp["AAB_vs_Control"].iloc[idx]
            others = gp["AAB_vs_Control"].drop(gp.index[idx])
            if others.empty:
                continue
            fakes.append(called - others.mean())
        perm[i] = np.mean(fakes) if fakes else np.nan
    perm = perm[np.isfinite(perm)]
    p_one = (np.sum(perm >= observed) + 1) / (len(perm) + 1)
    print(f"Permutation null: mean={perm.mean():+.4f}, SD={perm.std():.4f}")
    print(f"One-sided p (called > random): {p_one:.4f}")

    # ----- figure -----
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    ax = axes[0]
    ax.hist(per_gene.values, bins=20, color="#3a7ca5", edgecolor="white", alpha=0.85)
    ax.axvline(0, color="black", lw=0.8)
    ax.axvline(observed, color="red", lw=1.5, ls="--",
               label=f"observed = {observed:+.3f}")
    ax.set_xlabel("called-CT effect − mean(other-CT effects), per gene")
    ax.set_ylabel("number of T1D-locus genes")
    ax.set_title(f"Cross-validation of cell-type-of-action calls\n"
                 f"({len(per_gene)} genes, tau ≥ {TAU_THRESHOLD})")
    ax.legend()
    style_axes(ax)

    ax = axes[1]
    rows = []
    for g in per_gene.index:
        gp = eff[eff["gene"] == g]
        called = gp.loc[gp["is_called_ct"], "AAB_vs_Control"]
        others = gp.loc[~gp["is_called_ct"], "AAB_vs_Control"].dropna()
        if called.empty or others.empty:
            continue
        max_other_signed = float(others.abs().max() *
                                 np.sign(others.iloc[others.abs().argmax()]))
        rows.append({
            "gene": g,
            "called_eff": float(called.iloc[0]),
            "max_other": max_other_signed,
            "called_ct": gp["top_celltype"].iloc[0],
        })
    plot_df = pd.DataFrame(rows)
    for ct in plot_df["called_ct"].unique():
        sub = plot_df[plot_df["called_ct"] == ct]
        ax.scatter(sub["max_other"], sub["called_eff"],
                   s=60, color=COMPARTMENT_COLORS.get(ct, "#999"),
                   edgecolor="black", lw=0.5, alpha=0.85,
                   label=f"{ct} (n={len(sub)})")
    lim = max(plot_df["called_eff"].abs().max(),
              plot_df["max_other"].abs().max()) * 1.15
    ax.plot([-lim, lim], [-lim, lim], color="gray", ls="--", lw=0.5, alpha=0.7)
    ax.axhline(0, color="black", lw=0.3)
    ax.axvline(0, color="black", lw=0.3)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("largest AAB effect in any OTHER cell type")
    ax.set_ylabel("AAB effect in the called cell type")
    ax.set_title("Per-gene effect: called CT vs alternatives")
    ax.legend(loc="lower right", fontsize=8)
    style_axes(ax)
    for _, r in plot_df.iterrows():
        if abs(r["called_eff"]) > 0.3 or abs(r["max_other"]) > 0.3:
            ax.annotate(r["gene"], (r["max_other"], r["called_eff"]),
                        xytext=(3, 3), textcoords="offset points", fontsize=7)

    plt.tight_layout()
    plot_path = RESULTS / "T1D145_celltype_call_crossvalidation.png"
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {plot_path}")


if __name__ == "__main__":
    main()
