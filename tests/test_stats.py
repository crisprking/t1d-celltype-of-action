"""Unit tests for t1d_coa.stats."""

import math

import numpy as np
import pandas as pd

from t1d_coa.stats import donor_pseudobulk_effect, permutation_pvalue, tau


# ---- tau ----------------------------------------------------------------

def test_tau_single_cell_type_is_one():
    assert tau([0, 0, 5, 0]) == 1.0


def test_tau_uniform_is_zero():
    assert tau([3, 3, 3, 3]) == 0.0


def test_tau_all_zero_is_nan():
    assert math.isnan(tau([0, 0, 0]))


def test_tau_intermediate():
    # Two equal peaks, two zeros → tau = (0 + 0 + 1 + 1) / 3 = 2/3
    assert math.isclose(tau([5, 5, 0, 0]), 2 / 3, rel_tol=1e-9)


# ---- permutation_pvalue --------------------------------------------------

def test_permutation_p_greater_when_observed_at_max():
    null = np.array([0.0, 0.1, 0.2, 0.3, 0.4])
    # Observed equals or exceeds 1 of 5 samples (just itself at 0.4).
    p = permutation_pvalue(0.4, null, "greater")
    assert math.isclose(p, (1 + 1) / (5 + 1), rel_tol=1e-9)


def test_permutation_p_below_min():
    null = np.array([1.0, 2.0, 3.0])
    # Observed exceeds zero null samples.
    p = permutation_pvalue(0.5, null, "greater")
    assert math.isclose(p, 1 / 4)


def test_permutation_two_sided():
    null = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
    p = permutation_pvalue(1.5, null, "two_sided")
    # |null| >= 1.5: only -2 and 2 → 2 samples.
    assert math.isclose(p, (2 + 1) / (5 + 1))


# ---- donor_pseudobulk_effect --------------------------------------------

def test_donor_pseudobulk_effect_clean_signal():
    # Two donors per state; per-cell expression equal within a donor.
    obs = pd.DataFrame({
        "donor_id":      ["A1", "A1", "A2", "B1", "B1", "B2"],
        "disease_state": ["A",  "A",  "A",  "B",  "B",  "B"],
        "cell_type":     ["x"] * 6,
    })
    expr = np.array([np.expm1(2), np.expm1(2), np.expm1(2),
                     np.expm1(1), np.expm1(1), np.expm1(1)])
    mask = pd.Series([True] * 6)
    effect = donor_pseudobulk_effect(expr, mask, obs, "A", "B")
    assert math.isclose(effect, 1.0, rel_tol=1e-6)


def test_donor_pseudobulk_effect_insufficient_donors():
    obs = pd.DataFrame({
        "donor_id":      ["A1", "B1"],
        "disease_state": ["A",  "B"],
        "cell_type":     ["x", "x"],
    })
    effect = donor_pseudobulk_effect(np.array([1.0, 1.0]),
                                     pd.Series([True, True]),
                                     obs, "A", "B")
    assert math.isnan(effect)
