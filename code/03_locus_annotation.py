"""03_locus_annotation.py
Annotate clumped loci with GWAS Catalog gene mappings.
Output: data/processed/t1d_independent_loci_annotated.tsv
"""
import pandas as pd
loci=pd.read_csv("data/interim/t1d_clumped.clumped",sep="\s+")
t1d=pd.read_csv("data/raw/gwas_catalog_t1d_associations.tsv",sep="	",dtype=str)
ann=t1d.groupby("SNPS").agg(MAPPED_GENE=("MAPPED_GENE","first"),n_studies=("STUDY ACCESSION","nunique")).reset_index().rename(columns={"SNPS":"SNP"})
loci.merge(ann,on="SNP",how="left").to_csv("data/processed/t1d_independent_loci_annotated.tsv",sep="	",index=False)
