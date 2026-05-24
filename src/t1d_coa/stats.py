"""Statistics used across the cell-type-of-action analysis.

Three primitives, each one defined once and reused:

* `tau` — Yanai 2005 tissue-specificity index, applied to cell-type means.
* `donor_pseudobulk_effect` — Δ log1p donor-mean expression between two states,
  the unit on which all disease-stage comparisons are built.
* `permutation_pvalue` — empirical p-value with a documented null definition.
"""

from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd


def tau(row: Iterable[float]) -> float:
    """Cell-type specificity index. NaN if every entry is zero.

    tau = sum_i (1 - x_i / x_max) / (N - 1).  Bounded [0, 1]; 1 means
    expression in exactly one cell type, 0 means uniform expression.
    """
    x = np.asarray(list(row), dtype=float)
    xmax = x.max()
    if xmax == 0:
        return float("nan")
    return float((1.0 - x / xmax).sum() / (len(x) - 1))


def donor_pseudobulk_effect(
    expr: np.ndarray,
    mask: pd.Series,
    obs: pd.DataFrame,
    state_a: str,
    state_b: str,
    min_donors: int = 2,
) -> float:
    """Δ log1p donor-mean expression (state_a minus state_b) inside a mask.

    Donor means are computed first, then the difference of group means —
    so a single cell-rich donor cannot dominate the statistic.
    """
    if mask.sum() == 0:
        return float("nan")
    sub = obs.loc[mask, ["donor_id", "disease_state"]].copy()
    sub["expr"] = expr[mask.to_numpy()]
    sub = sub[sub["disease_state"].isin([state_a, state_b])]
    if sub.empty:
        return float("nan")
    sub["log_expr"] = np.log1p(sub["expr"])
    donor_means = (
        sub.groupby(["donor_id", "disease_state"], observed=True)["log_expr"]
        .mean()
        .reset_index()
    )
    a = donor_means.loc[donor_means["disease_state"] == state_a, "log_expr"].to_numpy()
    b = donor_means.loc[donor_means["disease_state"] == state_b, "log_expr"].to_numpy()
    if len(a) < min_donors or len(b) < min_donors:
        return float("nan")
    return float(a.mean() - b.mean())


def permutation_pvalue(
    observed: float,
    null_samples: np.ndarray,
    alternative: str = "greater",
) -> float:
    """Empirical p-value with the standard +1/+1 small-sample correction.

    alternative ∈ {"greater", "less", "two_sided"}.
    """
    null_samples = null_samples[np.isfinite(null_samples)]
    n = len(null_samples)
    if n == 0:
        return float("nan")
    if alternative == "greater":
        k = int(np.sum(null_samples >= observed))
    elif alternative == "less":
        k = int(np.sum(null_samples <= observed))
    elif alternative == "two_sided":
        k = int(np.sum(np.abs(null_samples) >= abs(observed)))
    else:
        raise ValueError(f"unknown alternative: {alternative!r}")
    return (k + 1) / (n + 1)
