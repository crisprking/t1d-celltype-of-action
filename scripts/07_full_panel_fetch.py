"""Step 07 — fetch + aggregate expression for every testable T1D-locus gene.

Walks the 145 independent loci, expands them to ~188 candidate genes,
keeps the ~150 that exist in the HPAP var index, and pulls per-gene
expression vectors (idempotent — cached blobs are reused on rerun).
Output is a long-form gene × cell_type × disease_state table that
drives every downstream step.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from t1d_coa.config import EXPECTED_N_CELLS, INTERIM, PROCESSED
from t1d_coa.hpap import fetch_gene_expression, fetch_gene_list
from t1d_coa.provenance import log_artifact


def expand_genes(s) -> list[str]:
    if pd.isna(s) or not str(s).strip():
        return []
    return [g.strip() for g in str(s).replace(" - ", ", ").split(",") if g.strip()]


def main() -> None:
    obs = pd.read_csv(PROCESSED / "hpap_cellxgene_obs.tsv", sep="\t")
    assert len(obs) == EXPECTED_N_CELLS

    genes_all = fetch_gene_list(INTERIM / "hpap_gene_list.json")
    gene_to_idx = {g: i for i, g in enumerate(genes_all)}

    loci = pd.read_csv(PROCESSED / "t1d_independent_loci_annotated.tsv", sep="\t")
    loci["genes_list"] = loci["MAPPED_GENE"].apply(expand_genes)

    rows = [
        {"gene": g, "lead_snp": r["SNP"], "P": r["P"], "CHR": r["CHR"],
         "n_studies": r.get("n_studies", 1)}
        for _, r in loci.iterrows() for g in r["genes_list"]
    ]
    gene_df = pd.DataFrame(rows).sort_values("P").drop_duplicates("gene", keep="first")

    in_hpap = gene_df["gene"].isin(gene_to_idx)
    missing = gene_df.loc[~in_hpap, "gene"].tolist()
    if missing:
        (INTERIM / "missing_genes_from_hpap.txt").write_text("\n".join(missing))

    target = gene_df[in_hpap].copy()
    target["var_idx"] = target["gene"].map(gene_to_idx)
    target = target.reset_index(drop=True)
    print(f"Genes in HPAP: {len(target)}  (missing: {len(missing)})")

    expr, failures = {}, []
    for _, r in tqdm(target.iterrows(), total=len(target), desc="fetching"):
        try:
            expr[r["gene"]] = fetch_gene_expression(r["gene"], int(r["var_idx"]))
        except Exception as e:
            failures.append((r["gene"], str(e)))

    print(f"Fetched OK: {len(expr)}   Failed: {len(failures)}")

    # Stack into a cells × genes matrix once for vectorized aggregation.
    gene_names = list(expr.keys())
    mat = np.column_stack([expr[g] for g in gene_names])
    groups = obs.groupby(["cell_type", "disease_state"], observed=True).indices

    records = []
    for (ct, ds), pos in tqdm(groups.items(), desc="aggregating"):
        sub = mat[pos, :]
        means, medians = sub.mean(axis=0), np.median(sub, axis=0)
        nz = (sub > 0).sum(axis=0)
        for k, gene in enumerate(gene_names):
            records.append({
                "gene": gene, "cell_type": ct, "disease_state": ds,
                "n_cells": int(len(pos)),
                "mean": float(means[k]),
                "median": float(medians[k]),
                "pct_expressing": float(nz[k] / len(pos) * 100),
                "n_nonzero": int(nz[k]),
            })

    result = pd.DataFrame(records).sort_values(
        ["gene", "cell_type", "disease_state"]
    ).reset_index(drop=True)

    out = PROCESSED / "hpap_T1D145_expression_by_celltype_disease.tsv"
    result.to_csv(out, sep="\t", index=False)
    log_artifact(out,
                 source="Aggregated CellxGene blobs (T1D 145-locus gene set)",
                 notes=f"{len(result):,} rows × {result['gene'].nunique()} genes")

    target[["gene", "lead_snp", "CHR", "P", "n_studies", "var_idx"]].to_csv(
        PROCESSED / "T1D145_gene_locus_map.tsv", sep="\t", index=False
    )
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
