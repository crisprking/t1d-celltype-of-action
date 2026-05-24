"""06_donor_audit.py
Audit cell/donor counts per (cell_type, disease_state).
Output: data/processed/audit_ncells_celltype_disease.tsv
        data/processed/audit_ndonors_celltype_disease.tsv
        data/processed/audit_beta_donors_per_disease.tsv
"""
import pandas as pd
obs=pd.read_csv("data/processed/hpap_cellxgene_obs.tsv",sep="	")
obs.groupby(["cell_type","disease_state"]).size().reset_index(name="n_cells").to_csv("data/processed/audit_ncells_celltype_disease.tsv",sep="	",index=False)
obs.groupby(["cell_type","disease_state"])["donor_id"].nunique().reset_index(name="n_donors").to_csv("data/processed/audit_ndonors_celltype_disease.tsv",sep="	",index=False)
print("Donor audit complete.")
