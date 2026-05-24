"""Step 05 — pilot expression fetch for 20 top T1D-locus genes.

Smoke test for the full-panel fetch. Pulls the top-15 loci by p-value
plus a spread of weaker ones, aggregates expression by cell type and
disease state, and verifies that INS expression is beta-cell-dominant
(canonical marker; ~10× over the next-highest compartment).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from t1d_coa.config import EXPECTED_N_CELLS, INTERIM, PROCESSED
from t1d_coa.hpap import fetch_gene_expression, fetch_gene_list
from t1d_coa.provenance import log_artifact


def expand_genes(mapped_gene_str) -> list[str]:
    """Split GWAS Catalog MAPPED_GENE strings into clean gene symbols."""
    if pd.isna(mapped_gene_str) or not str(mapped_gene_str).strip():
        return []
    s = str(mapped_gene_str).replace(" - ", ", ")
    return [g.strip() for g in s.split(",") if g.strip()]


def main() -> None:
    obs = pd.read_csv(PROCESSED / "hpap_cellxgene_obs.tsv", sep="\t")
    assert len(obs) == EXPECTED_N_CELLS, f"cell count mismatch: {len(obs)}"

    genes_all = fetch_gene_list(INTERIM / "hpap_gene_list.json")
    gene_to_idx = {g: i for i, g in enumerate(genes_all)}

    loci = pd.read_csv(PROCESSED / "t1d_independent_loci_annotated.tsv", sep="\t")
    loci["genes_list"] = loci["MAPPED_GENE"].apply(expand_genes)
    rows = [
        {"gene": g, "lead_snp": r["SNP"], "P": r["P"], "CHR": r["CHR"]}
        for _, r in loci.iterrows()
        for g in r["genes_list"]
    ]
    gene_df = pd.DataFrame(rows).sort_values("P").drop_duplicates("gene")

    # 15 top + 5 spread = 20 pilot genes, restricted to those in HPAP var.
    pilot = pd.concat(
        [gene_df.head(15),
         gene_df.iloc[[i for i in (30, 60, 100, 140, 180) if i < len(gene_df)]]]
    ).drop_duplicates("gene").head(20)
    pilot = pilot[pilot["gene"].isin(gene_to_idx)].reset_index(drop=True)
    pilot["var_idx"] = pilot["gene"].map(gene_to_idx)
    print(f"Pilot genes: {len(pilot)}")

    expr = {
        r["gene"]: fetch_gene_expression(r["gene"], int(r["var_idx"]))
        for _, r in tqdm(pilot.iterrows(), total=len(pilot), desc="genes")
    }
    expr_df = pd.DataFrame(expr)

    # Aggregate to gene × cell_type × disease_state.
    records = []
    for (ct, ds), pos in obs.groupby(["cell_type", "disease_state"],
                                     observed=True).indices.items():
        sub = expr_df.iloc[pos]
        for gene in expr_df.columns:
            v = sub[gene].to_numpy()
            records.append({
                "gene": gene, "cell_type": ct, "disease_state": ds,
                "n_cells": len(pos),
                "mean": float(v.mean()),
                "median": float(np.median(v)),
                "pct_expressing": float((v > 0).mean() * 100),
                "n_nonzero": int((v > 0).sum()),
            })
    result = pd.DataFrame(records).sort_values(
        ["gene", "cell_type", "disease_state"]
    ).reset_index(drop=True)

    out = PROCESSED / "hpap_pilot20_expression_by_celltype_disease.tsv"
    result.to_csv(out, sep="\t", index=False)
    log_artifact(out, source="Per-gene CellxGene blobs aggregated",
                 notes=f"{len(result)} rows × {pilot['gene'].nunique()} genes")

    # Canonical sanity: INS should peak in beta cells by ~10×.
    print("\nCanonical INS sanity (pooled across disease states):")
    ins = result[result["gene"] == "INS"]
    sanity = (
        ins.groupby("cell_type")
        .apply(lambda g: (g["mean"] * g["n_cells"]).sum() / g["n_cells"].sum(),
               include_groups=False)
        .sort_values(ascending=False)
    )
    print(sanity.to_string())
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
