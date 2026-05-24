"""Step 08 — tau cell-type specificity + cell-type-of-action calls.

Tau (Yanai et al., 2005) is computed on log1p means weighted by cell
count, so disease states with few cells don't dominate. A gene gets
a "confident" call if tau ≥ 0.70 AND its top compartment is above the
expression floor; otherwise it's "broad" or "below detection".
"""

from __future__ import annotations

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from t1d_coa.config import (
    EXPR_FLOOR_LOG1P, PROCESSED, RESULTS, TAU_THRESHOLD,
)
from t1d_coa.plotting import COMPARTMENT_COLORS, style_axes
from t1d_coa.stats import tau


def _bucket(t: float) -> str:
    if pd.isna(t):
        return "no expression"
    if t >= 0.85:
        return "highly specific (≥0.85)"
    if t >= 0.70:
        return "specific (0.70–0.85)"
    if t >= 0.50:
        return "moderate (0.50–0.70)"
    return "broad / ubiquitous (<0.50)"


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    agg = pd.read_csv(
        PROCESSED / "hpap_T1D145_expression_by_celltype_disease.tsv", sep="\t"
    )
    print(f"Loaded: {len(agg):,} rows, {agg['gene'].nunique()} genes")

    # Collapse to gene × cell_type, weighted by n_cells.
    def weighted_mean(g: pd.DataFrame) -> float:
        return (g["mean"] * g["n_cells"]).sum() / g["n_cells"].sum()

    mat = (
        agg.groupby(["gene", "cell_type"])
        .apply(weighted_mean, include_groups=False)
        .unstack("cell_type")
        .fillna(0)
    )
    mat_log = np.log1p(mat)

    spec = pd.DataFrame({
        "gene": mat_log.index,
        "tau": mat_log.apply(tau, axis=1).values,
        "top_celltype": mat_log.idxmax(axis=1).values,
        "top_celltype_log1p_mean": mat_log.max(axis=1).values,
        "any_expression": mat_log.max(axis=1).values > 0,
    }).sort_values("tau", ascending=False).reset_index(drop=True)

    spec["bucket"] = spec["tau"].apply(_bucket)
    spec["above_floor"] = spec["top_celltype_log1p_mean"] >= EXPR_FLOOR_LOG1P
    spec["call"] = np.where(
        spec["above_floor"] & (spec["tau"] >= TAU_THRESHOLD),
        spec["top_celltype"],
        np.where(spec["above_floor"], "broad / multi-compartment",
                 "below detection"),
    )

    print("\n=== Tau distribution ===")
    print(spec["bucket"].value_counts().to_string())
    print(f"\n=== Cell-type-of-action calls (tau ≥ {TAU_THRESHOLD}, above floor) ===")
    print(spec["call"].value_counts().to_string())

    out = PROCESSED / "T1D145_celltype_specificity.tsv"
    spec.to_csv(out, sep="\t", index=False)
    print(f"\nSaved: {out}")

    # ----- figures -----
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))

    ax = axes[0]
    plot_data = spec[spec["above_floor"]]["tau"].dropna()
    ax.hist(plot_data, bins=25, color="#3a7ca5", edgecolor="white")
    ax.axvline(TAU_THRESHOLD, color="red", ls="--", lw=1.2,
               label=f"tau = {TAU_THRESHOLD}")
    ax.set_xlabel("tau (cell-type specificity)")
    ax.set_ylabel("number of T1D-locus genes")
    ax.set_title(f"Specificity of T1D-locus genes (n={len(plot_data)} above floor)")
    ax.legend()
    style_axes(ax)

    ax = axes[1]
    counts = (spec[spec["above_floor"]]["call"].value_counts()
              .sort_values(ascending=True))
    colors = [
        "#888" if c in ("broad / multi-compartment", "below detection")
        else COMPARTMENT_COLORS.get(c, "#1f77b4")
        for c in counts.index
    ]
    ax.barh(counts.index, counts.values, color=colors, edgecolor="white")
    ax.set_xlabel("number of T1D-locus genes")
    ax.set_title(f"Cell-type-of-action (tau ≥ {TAU_THRESHOLD})")
    style_axes(ax)

    plt.tight_layout()
    plot_path = RESULTS / "T1D145_specificity_distribution.png"
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {plot_path}")

    # Confident calls, ranked.
    confident = spec[(spec["tau"] >= TAU_THRESHOLD) & spec["above_floor"]].head(40)
    fig, ax = plt.subplots(figsize=(10, max(6, len(confident) * 0.28)))
    rev = confident.iloc[::-1]
    bar_colors = [COMPARTMENT_COLORS.get(ct, "#999") for ct in rev["top_celltype"]]
    ax.barh(rev["gene"], rev["tau"], color=bar_colors, edgecolor="white")
    ax.set_xlabel("tau (cell-type specificity)")
    ax.set_xlim(0.65, 1.0)
    ax.set_title("T1D-locus genes with confident cell-type-of-action")
    style_axes(ax)
    handles = [mpatches.Patch(color=COMPARTMENT_COLORS.get(ct, "#999"), label=ct)
               for ct in rev["top_celltype"].unique()]
    ax.legend(handles=handles, bbox_to_anchor=(1.02, 1), loc="upper left",
              frameon=False, fontsize=9)
    plt.tight_layout()
    plot_path2 = RESULTS / "T1D145_topgenes_by_compartment.png"
    plt.savefig(plot_path2, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {plot_path2}")


if __name__ == "__main__":
    main()
