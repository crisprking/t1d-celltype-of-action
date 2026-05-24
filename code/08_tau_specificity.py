"""08_tau_specificity.py
Compute tau cell-type specificity index (Yanai et al. 2005).
tau = 1 - mean(x_i / x_max)  over all cell types.
Threshold: tau >= 0.70 + expression floor = confident cell-type-of-action call.
Output: data/processed/T1D145_celltype_specificity.tsv
"""
import pandas as pd,numpy as np
expr=pd.read_csv("data/processed/hpap_T1D145_expression_by_celltype_disease.tsv",sep="	")
pivot=expr.groupby(["gene","cell_type"])["mean"].mean().unstack(fill_value=0)
tau=pivot.apply(lambda r:1-((r/r.max()).mean()) if r.max()>0 else 1.0,axis=1)
out=pd.DataFrame({"gene":tau.index,"tau":tau.values,"dominant_celltype":pivot.idxmax(axis=1).values})
out["call"]=out["tau"].apply(lambda t:"confident" if t>=0.70 else "broad")
out.to_csv("data/processed/T1D145_celltype_specificity.tsv",sep="	",index=False)
print(out["call"].value_counts().to_string())
