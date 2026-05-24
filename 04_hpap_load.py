"""03_locus_annotation: Annotate independent loci with mapped gene names.

Extracted from notebook cell 17. Part of the t1d-celltype-of-action pipeline.
Run modules in numeric order — each consumes outputs from earlier steps.
"""

import pandas as pd
from pathlib import Path

RAW = Path("/kaggle/working/t1d_mech/data/raw")
PROCESSED = Path("/kaggle/working/t1d_mech/data/processed")

loci = pd.read_csv(PROCESSED / "t1d_independent_loci.tsv", sep="\t")
t1d_all = pd.read_csv(RAW / "gwas_catalog_t1d_associations.tsv",
                      sep="\t", low_memory=False, dtype=str)

# For each clump lead SNP, pull the mapped gene and reported gene from the
# GWAS Catalog rows (collapse across studies — take the most common annotation)
ann = (
    t1d_all.dropna(subset=["SNPS"])
           .groupby("SNPS")
           .agg(MAPPED_GENE=("MAPPED_GENE",
                             lambda s: ", ".join(sorted(set(x for x in s if pd.notna(x) and x))[:3])),
                REPORTED_GENE=("REPORTED GENE(S)",
                               lambda s: ", ".join(sorted(set(x for x in s if pd.notna(x) and x))[:3])),
                MAPPED_TRAIT=("MAPPED_TRAIT",
                              lambda s: s.mode().iloc[0] if len(s.mode()) else ""),
                n_studies=("STUDY ACCESSION", "nunique"))
           .reset_index()
           .rename(columns={"SNPS": "SNP"})
)

loci_annotated = loci.merge(ann, on="SNP", how="left")
loci_annotated = loci_annotated.sort_values("P").reset_index(drop=True)

out = PROCESSED / "t1d_independent_loci_annotated.tsv"
loci_annotated.to_csv(out, sep="\t", index=False)
log_artifact(out, source="PLINK clumps annotated with GWAS Catalog gene mappings",
             notes=f"{len(loci_annotated)} loci with MAPPED_GENE")

# Quick look at the annotated top loci
print(loci_annotated[["CHR", "SNP", "BP", "P", "MAPPED_GENE", "n_studies", "n_absorbed"]].head(20).to_string(index=False))
print(f"\nSaved: {out}")
