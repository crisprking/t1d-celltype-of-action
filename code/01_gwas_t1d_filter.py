"""01_gwas_t1d_filter.py
Filter GWAS Catalog bulk TSV to pure T1D (MONDO_0005147).
Input:  data/raw/gwas-catalog-download-associations-alt-full.tsv
Output: data/raw/gwas_catalog_t1d_associations.tsv
        data/raw/gwas_catalog_t1d_gwsig.tsv
"""
import pathlib,pandas as pd
T1D_URI="http://purl.obolibrary.org/obo/MONDO_0005147"
bulk=pathlib.Path("data/raw/gwas-catalog-download-associations-alt-full.tsv")
df=pd.read_csv(bulk,sep="	",low_memory=False,dtype=str)
t1d=df[df["MAPPED_TRAIT_URI"].str.contains(T1D_URI,na=False)]
t1d=t1d[t1d["MAPPED_TRAIT"].str.lower().str.count(",")==0]
t1d.to_csv("data/raw/gwas_catalog_t1d_associations.tsv",sep="	",index=False)
gwsig=t1d[pd.to_numeric(t1d["PVALUE_MLOG"],errors="coerce")>=7.3]
gwsig.to_csv("data/raw/gwas_catalog_t1d_gwsig.tsv",sep="	",index=False)
print(f"T1D:{len(t1d):,}  GW-sig:{len(gwsig):,}")
