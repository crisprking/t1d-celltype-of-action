"""Step 04 — load HPAP cell-level annotations from the CellxGene REST API.

Pulls the seven obs columns we need (cell type, disease state, donor,
demographics, assay) for all 222,077 cells. Each column is fetched as
a FlatBuffers blob, cached on disk, and joined into a single TSV.
"""

from __future__ import annotations

import pandas as pd

from t1d_coa.config import EXPECTED_N_CELLS, PROCESSED
from t1d_coa.hpap import _client, fetch_obs_column
from t1d_coa.provenance import log_artifact


OBS_COLUMNS = ("cell_type", "disease_state", "donor_id", "age",
               "sex", "race", "assay")


def main() -> None:
    client = _client()
    obs_data = {}
    for col in OBS_COLUMNS:
        try:
            obs_data[col] = fetch_obs_column(col, client=client)
            print(f"  ✓ {col:15s} {len(obs_data[col]):>7,} values, "
                  f"{obs_data[col].nunique()} unique")
        except Exception as e:
            print(f"  ✗ {col:15s} {type(e).__name__}: {e}")

    if not obs_data:
        raise RuntimeError("All obs fetches failed; see errors above")

    obs = pd.DataFrame(obs_data)
    if len(obs) != EXPECTED_N_CELLS:
        print(f"⚠ Expected {EXPECTED_N_CELLS:,} cells, got {len(obs):,}")

    out = PROCESSED / "hpap_cellxgene_obs.tsv"
    obs.to_csv(out, sep="\t", index=False)
    log_artifact(out, source="Combined CellxGene obs columns",
                 notes=f"{len(obs):,} cells × {obs.shape[1]} annotations")

    print(f"\n=== HPAP atlas via CellxGene API ===")
    print(f"  Cells: {len(obs):,}   Donors: {obs['donor_id'].nunique()}")
    print(f"\n  Disease states:")
    for ds, n in obs["disease_state"].value_counts().items():
        print(f"    {ds:15s} {n:>8,}")
    print(f"\n  Cell types:")
    for ct, n in obs["cell_type"].value_counts().items():
        print(f"    {ct:25s} {n:>8,}")
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
