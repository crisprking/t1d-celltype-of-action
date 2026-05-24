"""02_ld_clumping.py
LD-clump T1D GW-sig SNPs with PLINK + 1000G EUR.
Output: data/interim/t1d_clumped.clumped  (145 independent loci)
"""
import subprocess,pathlib,pandas as pd
df=pd.read_csv("data/raw/gwas_catalog_t1d_gwsig.tsv",sep="	",dtype=str)
df=df.rename(columns={"SNPS":"SNP","P-VALUE":"P"})[["SNP","P"]].dropna()
df.to_csv("data/interim/t1d_clump_input.tsv",sep="	",index=False)
subprocess.run(["tools/plink","--bfile","reference/g1000_eur",
    "--clump","data/interim/t1d_clump_input.tsv",
    "--clump-snp-field","SNP","--clump-field","P",
    "--clump-p1","5e-8","--clump-p2","1e-5",
    "--clump-r2","0.1","--clump-kb","1000",
    "--out","data/interim/t1d_clumped"],check=True)
