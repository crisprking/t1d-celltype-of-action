"""Step 12 — per-donor INS and MEG3 audit (HPAP107 + MEG3 robustness).

Two quality-control questions:

1. Does HPAP107 (multi-AAB+, 4,606 beta cells) dominate the AAB signal
   the way HPAP084 dominates the T1D signal?
2. Does MEG3's elevation in T1D beta cells hold across multiple donors,
   or is it one donor's transcriptome?
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from t1d_coa.config import EXPR_CACHE, PROCESSED
from t1d_coa.donors import AAB_CLASS, T1D_CLASS
from t1d_coa.hpap import _decoder


def _load_gene(symbol: str) -> np.ndarray:
    files = list(EXPR_CACHE.glob(f"{symbol}_idx*.fbs"))
    if not files:
        raise FileNotFoundError(f"No cached blob for {symbol}")
    return _decoder()(files[0].read_bytes()).iloc[:, 0].to_numpy()


def _per_donor(df: pd.DataFrame, class_map: dict[str, str],
               label: str) -> pd.DataFrame:
    sub = df[df["donor_id"].isin(class_map)].copy()
    sub["class"] = sub["donor_id"].map(class_map)
    out = (
        sub.groupby(["donor_id", "class"], observed=True)
        .agg(
            n_beta=("INS", "size"),
            INS_log1p_mean=("INS", lambda x: np.log1p(x).mean()),
            MEG3_log1p_mean=("MEG3", lambda x: np.log1p(x).mean()),
        )
        .reset_index()
        .sort_values("n_beta", ascending=False)
    )
    print(f"\n=== {label} beta cells, per donor ===")
    print(out.round(2).to_string(index=False))
    return out


def _dominance(table: pd.DataFrame, label: str) -> None:
    total = table["n_beta"].sum()
    top = table.iloc[0]
    print(f"{label} dominance: {top['donor_id']} ({top['class']}) → "
          f"{top['n_beta']}/{total} = {top['n_beta'] / total:.1%}")


def _coherence_table(table: pd.DataFrame, label: str,
                     ctrl_meg3: float, ctrl_ins: float) -> None:
    print(f"\n=== {label} per-donor (control log1p(MEG3) = {ctrl_meg3:.2f}, "
          f"INS = {ctrl_ins:.2f}) ===")
    print(f"  donor      n_beta   MEG3   ΔMEG3   INS    ΔINS   class")
    for _, r in table.iterrows():
        dm = r["MEG3_log1p_mean"] - ctrl_meg3
        di = r["INS_log1p_mean"] - ctrl_ins
        a_m = "↑" if dm > 0.1 else ("↓" if dm < -0.1 else "·")
        a_i = "↑" if di > 0.1 else ("↓" if di < -0.1 else "·")
        print(f"  {r['donor_id']}  {int(r['n_beta']):5d}    "
              f"{r['MEG3_log1p_mean']:.2f}  {dm:+.2f} {a_m}  "
              f"{r['INS_log1p_mean']:.2f}  {di:+.2f} {a_i}  ({r['class']})")


def main() -> None:
    obs = pd.read_csv(PROCESSED / "hpap_cellxgene_obs.tsv", sep="\t")
    ins = _load_gene("INS")
    meg3 = _load_gene("MEG3")

    beta = obs[obs["cell_type"] == "beta cell"].copy()
    beta["INS"] = ins[beta.index]
    beta["MEG3"] = meg3[beta.index]

    t1d_tbl = _per_donor(beta, T1D_CLASS, "T1D")
    aab_tbl = _per_donor(beta, AAB_CLASS, "AAB")

    print()
    _dominance(t1d_tbl, "T1D")
    _dominance(aab_tbl, "AAB")

    ctrl = beta[beta["disease_state"] == "Control"]
    ctrl_meg3 = np.log1p(ctrl["MEG3"]).mean()
    ctrl_ins = np.log1p(ctrl["INS"]).mean()

    _coherence_table(t1d_tbl, "T1D", ctrl_meg3, ctrl_ins)
    _coherence_table(aab_tbl, "AAB", ctrl_meg3, ctrl_ins)


if __name__ == "__main__":
    main()
