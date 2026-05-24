"""Step 01 — filter GWAS Catalog bulk associations to pure T1D.

The Catalog mixes pure T1D associations with composite traits ("T1D and
celiac disease", "T1D nephropathy", …). Pure T1D is the target. The
filter is conservative: exact MAPPED_TRAIT name match, then exclude any
row whose trait string mentions a known composite or complication.
"""

from __future__ import annotations

import pandas as pd

from t1d_coa.config import RAW
from t1d_coa.provenance import log_artifact


T1D_NAMES = {"type 1 diabetes mellitus", "type 1 diabetes"}
EXCLUDE_SUBSTR = (
    "nephropathy", "retinopathy", "neuropathy", "celiac", "coeliac",
    "and ", "ketoacidosis", "complication",
)
GW_SIG = 5e-8


def main() -> None:
    src = RAW / "gwas-catalog-download-associations-alt-full.tsv"
    full = pd.read_csv(src, sep="\t", low_memory=False, dtype=str)
    print(f"Loaded: {len(full):,} rows × {full.shape[1]} columns")

    mt = full["MAPPED_TRAIT"].fillna("").str.lower().str.strip()
    in_t1d = mt.isin(T1D_NAMES)
    excluded = mt.apply(lambda s: any(x in s for x in EXCLUDE_SUBSTR))
    t1d = full.loc[in_t1d & ~excluded].copy()
    print(f"Pure T1D associations: {len(t1d):,}")

    t1d_path = RAW / "gwas_catalog_t1d_associations.tsv"
    t1d.to_csv(t1d_path, sep="\t", index=False)
    log_artifact(
        t1d_path,
        source=f"filtered {src.name} on MAPPED_TRAIT in T1D_NAMES",
        notes=f"{len(t1d)} pure-T1D associations; excluded composites",
    )

    p_num = pd.to_numeric(t1d["P-VALUE"], errors="coerce")
    gw = t1d.loc[p_num <= GW_SIG].copy()
    gw_path = RAW / "gwas_catalog_t1d_gwsig.tsv"
    gw.to_csv(gw_path, sep="\t", index=False)
    log_artifact(
        gw_path,
        source=t1d_path.name,
        notes=f"{len(gw)} GW-significant T1D associations (p ≤ {GW_SIG:g})",
    )

    print(f"\nGW-sig:        {len(gw):,}")
    print(f"Unique SNPs:   {gw['SNPS'].nunique():,}")
    print(f"Chromosomes:   {gw['CHR_ID'].dropna().nunique()}")
    print(f"\nSaved:")
    print(f"  {t1d_path}")
    print(f"  {gw_path}")


if __name__ == "__main__":
    main()
