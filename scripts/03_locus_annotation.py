"""Step 03 — annotate clumped loci with GWAS Catalog gene mappings.

For each lead SNP, collapse all Catalog rows for that SNP into a single
record: mapped gene(s), reported gene(s), most-common mapped trait, and
the number of independent studies. Output drives every downstream
expression query.
"""

from __future__ import annotations

import pandas as pd

from t1d_coa.config import PROCESSED, RAW
from t1d_coa.provenance import log_artifact


def _top_set(series: pd.Series, n: int = 3) -> str:
    return ", ".join(sorted({x for x in series if pd.notna(x) and x})[:n])


def main() -> None:
    loci = pd.read_csv(PROCESSED / "t1d_independent_loci.tsv", sep="\t")
    t1d_all = pd.read_csv(RAW / "gwas_catalog_t1d_associations.tsv",
                          sep="\t", low_memory=False, dtype=str)

    ann = (
        t1d_all.dropna(subset=["SNPS"])
        .groupby("SNPS")
        .agg(
            MAPPED_GENE=("MAPPED_GENE", _top_set),
            REPORTED_GENE=("REPORTED GENE(S)", _top_set),
            MAPPED_TRAIT=("MAPPED_TRAIT",
                          lambda s: s.mode().iloc[0] if len(s.mode()) else ""),
            n_studies=("STUDY ACCESSION", "nunique"),
        )
        .reset_index()
        .rename(columns={"SNPS": "SNP"})
    )

    annotated = loci.merge(ann, on="SNP", how="left").sort_values("P").reset_index(drop=True)
    annotated["n_absorbed"] = (
        annotated["SP2"].fillna("").apply(
            lambda s: 0 if s in ("", "NONE") else len(s.split(","))
        )
    )

    out = PROCESSED / "t1d_independent_loci_annotated.tsv"
    annotated.to_csv(out, sep="\t", index=False)
    log_artifact(
        out,
        source="PLINK clumps annotated with GWAS Catalog gene mappings",
        notes=f"{len(annotated)} loci with MAPPED_GENE",
    )

    cols = ["CHR", "SNP", "BP", "P", "MAPPED_GENE", "n_studies", "n_absorbed"]
    print(annotated[cols].head(20).to_string(index=False))
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
