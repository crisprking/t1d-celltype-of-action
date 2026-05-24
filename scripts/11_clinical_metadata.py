"""Step 11 — integrate HPAP donor clinical metadata.

Joins curated donor metadata (Faryabi 2023 Supp Table S2) onto our cell
counts and quantifies how much of the "T1D beta cell" signal is driven
by HPAP084 — a 12-year-old with C-peptide 2.20 ng/mL whose pancreas
contributes 67% of all surviving T1D beta cells in the cohort.
"""

from __future__ import annotations

import pandas as pd

from t1d_coa.config import PROCESSED
from t1d_coa.donors import AAB_META, T1D_META


def main() -> None:
    obs = pd.read_csv(PROCESSED / "hpap_cellxgene_obs.tsv", sep="\t")
    beta_per_donor = (
        obs[obs["cell_type"] == "beta cell"]
        .groupby("donor_id").size().rename("n_beta_cells")
    )
    total_per_donor = obs.groupby("donor_id").size().rename("n_cells_total")

    t1d = (T1D_META
           .merge(beta_per_donor, left_on="donor_id", right_index=True, how="left")
           .merge(total_per_donor, left_on="donor_id", right_index=True, how="left"))
    t1d["n_beta_cells"] = t1d["n_beta_cells"].fillna(0).astype(int)

    aab = (AAB_META
           .merge(beta_per_donor, left_on="donor_id", right_index=True, how="left")
           .merge(total_per_donor, left_on="donor_id", right_index=True, how="left"))

    print("=== T1D donor clinical metadata ===")
    print(t1d[["donor_id", "age", "sex", "disease_duration", "HbA1c",
               "c_peptide_ngml", "n_beta_cells", "n_cells_total",
               "clinical_class"]].to_string(index=False))

    print("\n=== AAB donor metadata ===")
    print(aab[["donor_id", "age", "sex", "HbA1c", "GAD_titer",
               "aab_class", "n_beta_cells", "n_cells_total"]].to_string(index=False))

    t1d.to_csv(PROCESSED / "T1D_donor_metadata.tsv", sep="\t", index=False)
    aab.to_csv(PROCESSED / "AAB_donor_metadata.tsv", sep="\t", index=False)

    total = t1d["n_beta_cells"].sum()
    h084 = t1d.loc[t1d["donor_id"] == "HPAP084", "n_beta_cells"].iloc[0]
    print(f"\n=== HPAP084 dominance audit ===")
    print(f"Total T1D beta cells:   {total}")
    print(f"From HPAP084 alone:     {h084} ({h084 / total:.1%})")

    by_class = t1d.groupby("clinical_class", observed=True)["n_beta_cells"].agg(
        ["sum", "count"])
    print("\nBeta cells by clinical class:")
    print(by_class.to_string())


if __name__ == "__main__":
    main()
