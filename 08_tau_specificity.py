"""06_donor_audit: Audit donor distribution per (cell_type, disease_state).

Extracted from notebook cell 31. Part of the t1d-celltype-of-action pipeline.
Run modules in numeric order — each consumes outputs from earlier steps.
"""

# Before any disease-state statistics: how many donors per
# (disease_state, cell_type), and is any cell × disease combination
# coming from a single donor (which would mean we can't generalize)?
#
# Outputs:
#   - donor × disease_state table (1 row per donor)
#   - n_donors per (cell_type, disease_state)
#   - n_cells per (donor_id, cell_type, disease_state)
#   - flags any beta-cell-disease combination with <2 donors
# ============================================================

from pathlib import Path
import pandas as pd
import numpy as np

PROJECT = Path("/kaggle/working/t1d_mech")
PROCESSED = PROJECT / "data" / "processed"

obs = pd.read_csv(PROCESSED / "hpap_cellxgene_obs.tsv", sep="\t")
print(f"obs: {len(obs):,} cells × {obs.shape[1]} annotations")

# --- 1. Donors per disease_state (each donor is one row) -------------------
print("\n=== 1. Donors by disease_state ===")
donor_disease = (obs.drop_duplicates("donor_id")
                    .groupby("disease_state").size()
                    .rename("n_donors").reset_index())
print(donor_disease.to_string(index=False))
print(f"  Total donors: {obs['donor_id'].nunique()}")

# Sanity: no donor should appear under two disease states
multi_state = (obs.groupby("donor_id")["disease_state"].nunique() > 1).sum()
assert multi_state == 0, f"{multi_state} donor(s) labeled with multiple disease states!"

# --- 2. n_donors per (cell_type, disease_state) ----------------------------
print("\n=== 2. n_donors per (cell_type, disease_state) ===")
donors_per_cell_disease = (
    obs.groupby(["cell_type", "disease_state"])["donor_id"].nunique().unstack(fill_value=0)
)
print(donors_per_cell_disease.to_string())

# --- 3. n_cells per (cell_type, disease_state) -----------------------------
print("\n=== 3. n_cells per (cell_type, disease_state) ===")
cells_per = (
    obs.groupby(["cell_type", "disease_state"]).size().unstack(fill_value=0)
)
print(cells_per.to_string())

# --- 4. Single-donor flags (where we have ≤1 donor contributing) -----------
print("\n=== 4. Combinations with ≤1 donor (NOT generalizable) ===")
flat = donors_per_cell_disease.stack().rename("n_donors").reset_index()
weak = flat[flat["n_donors"] <= 1].sort_values("n_donors")
if len(weak):
    print(weak.to_string(index=False))
else:
    print("  None — every (cell_type, disease_state) combo has ≥2 donors")

# --- 5. Beta-cell donor breakdown across disease states (the key cell) -----
print("\n=== 5. Beta cells per donor (the bottleneck for T1D stats) ===")
beta = obs[obs["cell_type"] == "beta cell"]
beta_breakdown = (
    beta.groupby(["disease_state", "donor_id"]).size().rename("n_beta_cells")
        .reset_index()
        .sort_values(["disease_state", "n_beta_cells"], ascending=[True, False])
)
print(beta_breakdown.to_string(index=False))
print(f"\n  Summary per disease_state:")
print(
    beta_breakdown.groupby("disease_state")["n_beta_cells"]
        .agg(n_donors="count", total_cells="sum",
             median_per_donor="median", min_per_donor="min", max_per_donor="max")
        .to_string()
)

# --- 6. Save audit tables for downstream stats -----------------------------
donors_per_cell_disease.to_csv(PROCESSED / "audit_ndonors_celltype_disease.tsv", sep="\t")
cells_per.to_csv(PROCESSED / "audit_ncells_celltype_disease.tsv", sep="\t")
beta_breakdown.to_csv(PROCESSED / "audit_beta_donors_per_disease.tsv", sep="\t", index=False)
print(f"\nSaved audit tables to {PROCESSED}/audit_*.tsv")
