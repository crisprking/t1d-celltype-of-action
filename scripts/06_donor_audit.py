"""Step 06 — donor distribution audit for every cell-type × disease cell.

Before any disease-stage statistics run, check that each
(cell_type, disease_state) bucket draws from ≥ 2 donors. Anywhere it
doesn't, we can't generalize — flag it and let downstream steps decide
how to handle it.

Beta cells are the bottleneck (post-T1D ablation): we report per-donor
counts in that compartment so the dominance of single high-cell donors
is visible from the audit alone.
"""

from __future__ import annotations

import pandas as pd

from t1d_coa.config import PROCESSED


def main() -> None:
    obs = pd.read_csv(PROCESSED / "hpap_cellxgene_obs.tsv", sep="\t")
    print(f"obs: {len(obs):,} cells × {obs.shape[1]} annotations")

    # Sanity: no donor should span multiple disease states.
    multi = (obs.groupby("donor_id")["disease_state"].nunique() > 1).sum()
    assert multi == 0, f"{multi} donor(s) labeled with multiple disease states"

    donors = (
        obs.groupby(["cell_type", "disease_state"])["donor_id"]
        .nunique()
        .unstack(fill_value=0)
    )
    cells = (
        obs.groupby(["cell_type", "disease_state"])
        .size()
        .unstack(fill_value=0)
    )

    print("\n=== Donors per (cell_type, disease_state) ===")
    print(donors.to_string())
    print("\n=== Cells per (cell_type, disease_state) ===")
    print(cells.to_string())

    flat = donors.stack().rename("n_donors").reset_index()
    weak = flat[flat["n_donors"] <= 1].sort_values("n_donors")
    print("\n=== Combinations with ≤1 donor (NOT generalizable) ===")
    print(weak.to_string(index=False) if len(weak) else
          "  None — every combo has ≥2 donors")

    beta = (
        obs[obs["cell_type"] == "beta cell"]
        .groupby(["disease_state", "donor_id"]).size()
        .rename("n_beta_cells").reset_index()
        .sort_values(["disease_state", "n_beta_cells"], ascending=[True, False])
    )
    print("\n=== Beta cells per donor (the T1D-stats bottleneck) ===")
    print(beta.to_string(index=False))

    donors.to_csv(PROCESSED / "audit_ndonors_celltype_disease.tsv", sep="\t")
    cells.to_csv(PROCESSED / "audit_ncells_celltype_disease.tsv", sep="\t")
    beta.to_csv(PROCESSED / "audit_beta_donors_per_disease.tsv",
                sep="\t", index=False)
    print(f"\nSaved audit tables to {PROCESSED}/audit_*.tsv")


if __name__ == "__main__":
    main()
