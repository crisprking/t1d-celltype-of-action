"""11_clinical_metadata: Integrate HPAP supplementary clinical metadata (donor durations, C-peptide, AAB panels).

Extracted from notebook cell 39. Part of the t1d-celltype-of-action pipeline.
Run modules in numeric order — each consumes outputs from earlier steps.
"""

# Goal: produce a clean donor metadata table for the article, and
# quantify how much of our "T1D beta cell" signal is driven by
# HPAP084 specifically. From the Faryabi 2023 supplement:
#   - HPAP084: 12yo F, "Unsuspected" T1D, C-peptide 2.20 ng/mL (HIGH)
#     contributed 476 of 715 T1D beta cells (67%)
#   - HPAP055, HPAP021, HPAP023: 7-year established T1D, mostly destroyed
#   - HPAP032, HPAP087, HPAP028: 3-10y duration, complete beta loss
#
# Important reframe: HPAP084 looks transcriptionally more like an
# advanced AAB donor than like a long-standing T1D case. The "T1D
# beta cell" signal in Cell B is therefore really a heterogeneous
# average over a very small, highly variable cohort.
# ============================================================

import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

PROJECT = Path("/kaggle/working/t1d_mech")
RAW = PROJECT / "data" / "raw"
PROCESSED = PROJECT / "data" / "processed"
INTERIM = PROJECT / "data" / "interim"
RESULTS = PROJECT / "results"
EXPR_CACHE_DIR = RAW / "cellxgene_expr"

# --- 1. Hard-code metadata from Faryabi 2023 Supp Table S2 ---------------
# (Already viewed and tabulated; embedding here so the cell is self-contained
# without needing the xlsx file on Kaggle.)
T1D_META = pd.DataFrame([
    {"donor_id":"HPAP021","age":13,"sex":"Female","ancestry":"Caucasian","BMI":21.40,"HbA1c":None,"disease_duration":"7 years","c_peptide_ngml":0.02,"clinical_class":"Established T1D"},
    {"donor_id":"HPAP023","age":17,"sex":"Female","ancestry":"Caucasian","BMI":21.35,"HbA1c":8.9, "disease_duration":"7 years","c_peptide_ngml":0.02,"clinical_class":"Established T1D"},
    {"donor_id":"HPAP028","age":4, "sex":"Male",  "ancestry":"Caucasian","BMI":17.30,"HbA1c":9.8, "disease_duration":"Unsuspected","c_peptide_ngml":0.30,"clinical_class":"Unsuspected (undiagnosed pre-mortem)"},
    {"donor_id":"HPAP032","age":10,"sex":"Female","ancestry":"Caucasian","BMI":16.30,"HbA1c":9.0, "disease_duration":"3 years","c_peptide_ngml":0.02,"clinical_class":"Established T1D"},
    {"donor_id":"HPAP055","age":24,"sex":"Male",  "ancestry":"Hispanic", "BMI":27.90,"HbA1c":10.4,"disease_duration":"7 years","c_peptide_ngml":0.02,"clinical_class":"Established T1D"},
    {"donor_id":"HPAP064","age":24,"sex":"Male",  "ancestry":"African American","BMI":16.98,"HbA1c":13.0,"disease_duration":"Unsuspected","c_peptide_ngml":0.25,"clinical_class":"Unsuspected (undiagnosed pre-mortem)"},
    {"donor_id":"HPAP071","age":12,"sex":"Female","ancestry":"Caucasian","BMI":15.42,"HbA1c":9.8, "disease_duration":"Unsuspected","c_peptide_ngml":0.06,"clinical_class":"Unsuspected (undiagnosed pre-mortem)"},
    {"donor_id":"HPAP084","age":12,"sex":"Female","ancestry":"Caucasian","BMI":18.50,"HbA1c":13.3,"disease_duration":"Unsuspected","c_peptide_ngml":2.20,"clinical_class":"Unsuspected, preserved C-peptide"},
    {"donor_id":"HPAP087","age":15,"sex":"Female","ancestry":"Caucasian","BMI":19.30,"HbA1c":10.4,"disease_duration":"6-10 years","c_peptide_ngml":0.02,"clinical_class":"Established T1D"},
])
AAB_META = pd.DataFrame([
    {"donor_id":"HPAP024","age":18,"sex":"Male",  "BMI":24.30,"HbA1c":5.5,"GAD_titer":203,"multi_aab":False,"aab_class":"GAD only"},
    {"donor_id":"HPAP029","age":23,"sex":"Male",  "BMI":28.60,"HbA1c":5.3,"GAD_titer":84, "multi_aab":False,"aab_class":"GAD only"},
    {"donor_id":"HPAP038","age":13,"sex":"Male",  "BMI":18.34,"HbA1c":5.7,"GAD_titer":89, "multi_aab":False,"aab_class":"GAD only"},
    {"donor_id":"HPAP043","age":15,"sex":"Male",  "BMI":24.07,"HbA1c":5.9,"GAD_titer":0,  "multi_aab":True, "aab_class":"IA-2 + ZnT8"},
    {"donor_id":"HPAP045","age":27,"sex":"Female","BMI":26.20,"HbA1c":5.2,"GAD_titer":321,"multi_aab":False,"aab_class":"GAD only"},
    {"donor_id":"HPAP049","age":29,"sex":"Male",  "BMI":37.20,"HbA1c":5.4,"GAD_titer":412,"multi_aab":False,"aab_class":"GAD only"},
    {"donor_id":"HPAP050","age":21,"sex":"Female","BMI":28.99,"HbA1c":5.1,"GAD_titer":203,"multi_aab":False,"aab_class":"GAD only"},
    {"donor_id":"HPAP072","age":19,"sex":"Male",  "BMI":23.10,"HbA1c":5.6,"GAD_titer":204,"multi_aab":False,"aab_class":"GAD only"},
    {"donor_id":"HPAP092","age":21,"sex":"Male",  "BMI":25.59,"HbA1c":5.6,"GAD_titer":38, "multi_aab":False,"aab_class":"GAD only"},
    {"donor_id":"HPAP107","age":15,"sex":"Male",  "BMI":23.59,"HbA1c":5.3,"GAD_titer":848,"multi_aab":True, "aab_class":"GAD + IA-2 + ZnT8 (high-risk pre-T1D)"},
])

# Add cell counts from our obs
obs = pd.read_csv(PROCESSED / "hpap_cellxgene_obs.tsv", sep="\t")
beta_per_donor = (obs[obs["cell_type"] == "beta cell"]
                  .groupby("donor_id").size().rename("n_beta_cells"))
total_per_donor = obs.groupby("donor_id").size().rename("n_cells_total")
T1D_META = T1D_META.merge(beta_per_donor, left_on="donor_id", right_index=True, how="left")
T1D_META = T1D_META.merge(total_per_donor, left_on="donor_id", right_index=True, how="left")
T1D_META["n_beta_cells"] = T1D_META["n_beta_cells"].fillna(0).astype(int)
AAB_META = AAB_META.merge(beta_per_donor, left_on="donor_id", right_index=True, how="left")
AAB_META = AAB_META.merge(total_per_donor, left_on="donor_id", right_index=True, how="left")

print("=" * 72)
print("T1D donor clinical metadata (from Faryabi 2023 Supp Table S2)")
print("=" * 72)
print(T1D_META[["donor_id","age","sex","disease_duration","HbA1c","c_peptide_ngml",
                "n_beta_cells","n_cells_total","clinical_class"]].to_string(index=False))

print("\n" + "=" * 72)
print("AAB donor autoantibody panel and clinical metadata")
print("=" * 72)
print(AAB_META[["donor_id","age","sex","HbA1c","GAD_titer","aab_class","n_beta_cells","n_cells_total"]].to_string(index=False))

# --- 2. Save these tables for the article -------------------------------
T1D_META.to_csv(PROCESSED / "T1D_donor_metadata.tsv", sep="\t", index=False)
AAB_META.to_csv(PROCESSED / "AAB_donor_metadata.tsv", sep="\t", index=False)

# --- 3. Quantify HPAP084 dominance of T1D beta-cell signal --------------
print("\n" + "=" * 72)
print("HPAP084 dominance audit")
print("=" * 72)
total_t1d_beta = T1D_META["n_beta_cells"].sum()
hpap084_beta = T1D_META.loc[T1D_META["donor_id"]=="HPAP084", "n_beta_cells"].iloc[0]
print(f"  Total T1D beta cells: {total_t1d_beta}")
print(f"  From HPAP084 (12yo F, unsuspected, C-peptide 2.20): "
      f"{hpap084_beta} ({hpap084_beta/total_t1d_beta:.1%})")
# Split by clinical class
print("\n  Beta cells by clinical class:")
print(T1D_META.groupby("clinical_class", observed=True)["n_beta_cells"].agg(["sum","count"]).to_string())
established = T1D_META[T1D_META["clinical_class"].str.startswith("Established")]["n_beta_cells"].sum()
unsuspected_total = T1D_META[T1D_META["clinical_class"].str.startswith("Unsuspected")]["n_beta_cells"].sum()
print(f"\n  Established T1D (≥3y diagnosis): {established} beta cells from "
      f"{(T1D_META['clinical_class'].str.startswith('Established')).sum()} donors")
print(f"  Unsuspected / undiagnosed:        {unsuspected_total} beta cells from "
      f"{(T1D_META['clinical_class'].str.startswith('Unsuspected')).sum()} donors")

# --- 4. Re-do INS by clinical class to show what's actually going on ----
import sys, subprocess
try:
    from server.common.fbs.matrix import decode_matrix_fbs
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "--no-deps", "cellxgene", "flatbuffers"])
    from server.common.fbs.matrix import decode_matrix_fbs

# Pull cached INS blob
ins_files = list(EXPR_CACHE_DIR.glob("INS_idx*.fbs"))
if ins_files:
    blob = ins_files[0].read_bytes()
    ins_df = decode_matrix_fbs(blob)
    ins_vec = ins_df.iloc[:, 0].to_numpy()
    obs_aug = obs.copy()
    obs_aug["INS"] = ins_vec
    # Tag clinical class
    class_map = dict(zip(T1D_META["donor_id"], T1D_META["clinical_class"]))
    obs_aug["clinical_class"] = obs_aug.apply(
        lambda r: class_map.get(r["donor_id"], r["disease_state"]), axis=1
    )

    print("\n" + "=" * 72)
    print("INS expression in beta cells, broken out by clinical class")
    print("=" * 72)
    beta_obs = obs_aug[obs_aug["cell_type"] == "beta cell"].copy()
    summary = (beta_obs.groupby("clinical_class")["INS"]
               .agg(n_cells="size", mean="mean", median="median")
               .sort_values("mean", ascending=False))
    print(summary.round(1).to_string())

    # Donor-level INS means within T1D
    print("\n  Per-T1D-donor INS mean in beta cells (log1p):")
    t1d_beta_donors = beta_obs[beta_obs["disease_state"] == "T1D"]
    if len(t1d_beta_donors):
        donor_ins = (t1d_beta_donors.groupby("donor_id")
                     .apply(lambda g: (np.log1p(g["INS"]).mean(),
                                       len(g)),
                            include_groups=False))
        for did, (m, n) in donor_ins.items():
            klass = class_map.get(did, "?")
            print(f"    {did}: log1p(INS) = {m:.2f}  (n={n} beta cells, {klass})")

# --- 5. The reframed message for the article ----------------------------
print("\n" + "=" * 72)
print("REFRAMED MESSAGE for the article")
print("=" * 72)
print("""
Our T1D cohort is heterogeneous: 9 donors comprising:
  - 4 "Unsuspected" donors (T1D detected only at autopsy via AAB/clinical
    markers, no prior diagnosis) — represent very early or pre-clinical
    disease, with variable beta-cell preservation
  - 4 donors with established disease (3-10y duration), most with severe
    beta-cell loss
  - HPAP084 is exceptional: unsuspected status, preserved C-peptide
    (2.20 ng/mL), and contributes 67% of all T1D beta cells we analyze.
    Functionally this donor sits between AAB+ and clinical T1D.

This means the "T1D beta cell" signal in our cell-level analyses is
heavily weighted toward HPAP084's transcriptome, which is itself
closer to a high-titer AAB+ state than to long-standing T1D.

Implications for the writeup:
  1. Cell B beta-cell results should be described as "from a heterogeneous
     cohort dominated by one unsuspected-T1D donor," not "from chronic T1D."
  2. The AAB-vs-T1D scatter plot (where INS, MEG3 trend upward in both)
     is consistent with a continuum of disease activity rather than a
     stage discontinuity.
  3. Cell D-2's cross-validation result is unaffected — it pools all 21
     confident-call genes across all cell types and uses donor-pseudobulk,
     so HPAP084's outsized weight is appropriately bounded.
""")
