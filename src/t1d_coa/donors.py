"""Hand-curated clinical metadata for HPAP T1D and AAB donors.

From Faryabi et al. 2023 Supplementary Table S2. Encoded here so any
script can `from t1d_coa.donors import T1D_META, AAB_META` instead of
re-typing the table; one source of truth, easy to diff if HPAP updates.
"""

from __future__ import annotations

import pandas as pd

T1D_META = pd.DataFrame(
    [
        # donor_id, age, sex, ancestry, BMI, HbA1c, duration, c_peptide, clinical_class
        ("HPAP021", 13, "Female", "Caucasian", 21.40, None, "7 years", 0.02, "Established T1D"),
        ("HPAP023", 17, "Female", "Caucasian", 21.35, 8.9, "7 years", 0.02, "Established T1D"),
        ("HPAP028", 4, "Male", "Caucasian", 17.30, 9.8, "Unsuspected", 0.30, "Unsuspected (undiagnosed pre-mortem)"),
        ("HPAP032", 10, "Female", "Caucasian", 16.30, 9.0, "3 years", 0.02, "Established T1D"),
        ("HPAP055", 24, "Male", "Hispanic", 27.90, 10.4, "7 years", 0.02, "Established T1D"),
        ("HPAP064", 24, "Male", "African American", 16.98, 13.0, "Unsuspected", 0.25, "Unsuspected (undiagnosed pre-mortem)"),
        ("HPAP071", 12, "Female", "Caucasian", 15.42, 9.8, "Unsuspected", 0.06, "Unsuspected (undiagnosed pre-mortem)"),
        ("HPAP084", 12, "Female", "Caucasian", 18.50, 13.3, "Unsuspected", 2.20, "Unsuspected, preserved C-peptide"),
        ("HPAP087", 15, "Female", "Caucasian", 19.30, 10.4, "6-10 years", 0.02, "Established T1D"),
    ],
    columns=["donor_id", "age", "sex", "ancestry", "BMI", "HbA1c",
             "disease_duration", "c_peptide_ngml", "clinical_class"],
)

AAB_META = pd.DataFrame(
    [
        # donor_id, age, sex, BMI, HbA1c, GAD_titer, multi_aab, aab_class
        ("HPAP024", 18, "Male", 24.30, 5.5, 203, False, "GAD only"),
        ("HPAP029", 23, "Male", 28.60, 5.3, 84, False, "GAD only"),
        ("HPAP038", 13, "Male", 18.34, 5.7, 89, False, "GAD only"),
        ("HPAP043", 15, "Male", 24.07, 5.9, 0, True, "IA-2 + ZnT8"),
        ("HPAP045", 27, "Female", 26.20, 5.2, 321, False, "GAD only"),
        ("HPAP049", 29, "Male", 37.20, 5.4, 412, False, "GAD only"),
        ("HPAP050", 21, "Female", 28.99, 5.1, 203, False, "GAD only"),
        ("HPAP072", 19, "Male", 23.10, 5.6, 204, False, "GAD only"),
        ("HPAP092", 21, "Male", 25.59, 5.6, 38, False, "GAD only"),
        ("HPAP107", 15, "Male", 23.59, 5.3, 848, True,
         "GAD + IA-2 + ZnT8 (high-risk pre-T1D)"),
    ],
    columns=["donor_id", "age", "sex", "BMI", "HbA1c", "GAD_titer",
             "multi_aab", "aab_class"],
)

# Per-donor clinical-class lookup used in donor audits.
T1D_CLASS = {
    "HPAP021": "Established 7y", "HPAP023": "Established 7y",
    "HPAP032": "Established 3y", "HPAP055": "Established 7y",
    "HPAP087": "Established 6-10y",
    "HPAP028": "Unsuspected", "HPAP064": "Unsuspected",
    "HPAP071": "Unsuspected",
    "HPAP084": "Unsuspected (preserved C-pep)",
}

AAB_CLASS = {
    "HPAP024": "GAD only", "HPAP029": "GAD only", "HPAP038": "GAD only",
    "HPAP045": "GAD only", "HPAP049": "GAD only", "HPAP050": "GAD only",
    "HPAP072": "GAD only", "HPAP092": "GAD only",
    "HPAP043": "Multi-AAB (IA-2+ZnT8)",
    "HPAP107": "Multi-AAB (GAD+IA-2+ZnT8) HIGH-RISK",
}
