"""01_gwas_t1d_filter: Filter GWAS Catalog associations to T1D (MONDO:0005147).

Extracted from notebook cell 14. Part of the t1d-celltype-of-action pipeline.
Run modules in numeric order — each consumes outputs from earlier steps.
"""

import pandas as pd
from pathlib import Path

RAW = Path("/kaggle/working/t1d_mech/data/raw")
src = RAW / "gwas-catalog-download-associations-alt-full.tsv"

full = pd.read_csv(src, sep="\t", low_memory=False, dtype=str)
log.info("Loaded: %d rows × %d columns", len(full), full.shape[1])

# Filter on MAPPED_TRAIT name with explicit T1D variants. We deliberately
# exclude T2D, gestational, and related-but-distinct conditions.
T1D_NAMES = {
    "type 1 diabetes mellitus",
    "type 1 diabetes",
}
# Substrings to exclude even if they contain "type 1 diabetes" (e.g.
# "type 1 diabetes nephropathy" or "type 1 diabetes and celiac disease"
# are composite traits that aren't pure T1D susceptibility).
EXCLUDE_SUBSTR = [
    "nephropathy", "retinopathy", "neuropathy", "celiac", "coeliac",
    "and ", "ketoacidosis", "complication",
]

mt = full["MAPPED_TRAIT"].fillna("").str.lower().str.strip()
in_t1d = mt.isin(T1D_NAMES)
excluded = mt.apply(lambda s: any(x in s for x in EXCLUDE_SUBSTR))
t1d = full.loc[in_t1d & ~excluded].copy()

log.info("Pure T1D associations: %d", len(t1d))

# What ontology URIs map to T1D in this release? (informational)
print("\nOntology URIs found for T1D:")
uri_counts = t1d["MAPPED_TRAIT_URI"].value_counts()
print(uri_counts.head(10).to_string())

# Show what we excluded — useful to verify we're not over-filtering
print(f"\nExcluded composite/complication traits (sample):")
composite = full.loc[mt.str.contains("type 1 diabetes", na=False) & excluded, "MAPPED_TRAIT"]
print(composite.value_counts().head(10).to_string())

# Save
t1d_path = RAW / "gwas_catalog_t1d_associations.tsv"
t1d.to_csv(t1d_path, sep="\t", index=False)
log_artifact(t1d_path, source=f"filtered {src.name} on MAPPED_TRAIT == 'type 1 diabetes'",
             notes=f"{len(t1d)} pure-T1D associations; excluded composites")

# Headline numbers
p_num = pd.to_numeric(t1d["P-VALUE"], errors="coerce")
gw = t1d.loc[p_num <= 5e-8].copy()

print(f"\n=== T1D headline numbers ===")
print(f"  Pure T1D associations:     {len(t1d):,}")
print(f"  Unique studies:            {t1d['STUDY ACCESSION'].nunique():,}")
print(f"  Unique SNPs (raw lead):    {t1d['SNPS'].nunique():,}")
print(f"  Genome-wide sig (p≤5e-8):  {len(gw):,}")
print(f"  Unique SNPs at GW sig:     {gw['SNPS'].nunique():,}")
print(f"  Chromosomes covered:       {t1d['CHR_ID'].dropna().nunique()}")
print(f"  Median p-value:            {p_num.median():.2e}")

# Most-replicated SNPs across studies — useful sanity check
print(f"\nTop 10 most-replicated lead SNPs in T1D:")
print(
    t1d.groupby("SNPS")["STUDY ACCESSION"].nunique()
    .sort_values(ascending=False).head(10).to_string()
)

# Save the genome-wide-significant subset separately — this is what you'll
# actually use downstream for locus inventory.
gw_path = RAW / "gwas_catalog_t1d_gwsig.tsv"
gw.to_csv(gw_path, sep="\t", index=False)
log_artifact(gw_path, source=t1d_path.name,
             notes=f"{len(gw)} GW-significant T1D associations (p ≤ 5e-8)")
print(f"\nSaved:")
print(f"  All T1D:         {t1d_path}")
print(f"  GW-sig subset:   {gw_path}")
