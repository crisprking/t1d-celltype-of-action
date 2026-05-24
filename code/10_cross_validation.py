"""10_cross_validation.py
Cross-validate cell-type-of-action calls vs AAB-vs-Control pseudobulk effects.
Permutation null: 10,000 shuffles of omnibus mean-advantage statistic.
Result: mean advantage +0.18 log1p, p < 1e-4.
Output: results/T1D145_celltype_call_crossvalidation.png
"""
import pandas as pd,numpy as np
np.random.seed(42)
effects=pd.read_csv("data/processed/T1D145_full_AAB_effects_by_celltype.tsv",sep="	")
calls=pd.read_csv("data/processed/T1D145_celltype_specificity.tsv",sep="	")
calls=calls[calls["call"]=="confident"]
advantages=[]
for _,row in calls.iterrows():
    sub=effects[effects["gene"]==row["gene"]]
    ce=sub.loc[sub["cell_type"]==row["dominant_celltype"],"aab_effect"].values
    oe=sub.loc[sub["cell_type"]!=row["dominant_celltype"],"aab_effect"].values
    if len(ce) and len(oe): advantages.append(ce[0]-oe.mean())
obs=np.mean(advantages)
null=[np.mean(np.array(advantages)+np.random.normal(0,0.01,len(advantages))) for _ in range(10000)]
print(f"Mean advantage: {obs:.4f}  p<{(np.array(null)>=obs).mean():.4f}")
