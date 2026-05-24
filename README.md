# t1d-celltype-of-action

**Mapping 145 type 1 diabetes GWAS loci to the pancreatic cell types they likely act in, using cell-type specificity in the HPAP scRNA-seq atlas — cross-validated against pre-clinical (autoantibody-positive) transcriptional changes.**

Companion to the Substack article: *Where do type 1 diabetes risk genes actually live?* (`docs/substack_article.md`).

## TL;DR

- 145 independent T1D risk loci (GWAS Catalog v1.0.3.1) → 188 candidate genes → 150 testable in the HPAP atlas (222,077 cells, 67 donors).
- Tau cell-type specificity (Yanai 2005) assigns a cell-type-of-action to each gene from healthy expression.
- **21 of 150 genes received a confident call:** 7 leukocyte, 6 duct epithelial, 5 acinar, 2 mesenchymal, 1 endothelial. Zero beta-cell — a tau exclusivity artifact, not biology.
- **Headline:** the calls cross-validate. 20/21 confidently-called genes show a larger AAB-vs-Control effect in the called compartment than the mean across other compartments. Mean advantage +0.18 log1p; permutation p < 10⁻⁴.
- **Caveats kept honest:** AAB contributes ~14% of the matrix used for tau (cross-validation not fully out-of-sample); permutation null is uniform over cell types, not marginal-preserving; effect is carried by ~5 strong genes (CEL, CTRB1, HLA-DRB1, HLA-DQA1, PLPP1).

## Repository layout

```
.
├── src/t1d_coa/             reusable library — config, HPAP client, stats, donors
│   ├── config.py            paths, thresholds, HPAP endpoint
│   ├── hpap.py              CellxGene REST + FBS decoder
│   ├── stats.py             tau, donor-pseudobulk effect, permutation p-value
│   ├── donors.py            curated HPAP donor metadata
│   ├── plotting.py          shared compartment colors & axes style
│   └── provenance.py        append-only sha256 artifact log
├── scripts/                 13 pipeline steps, runnable in order
│   ├── 00_env_setup.py
│   ├── 01_gwas_t1d_filter.py
│   ├── 02_ld_clumping.py
│   ├── 03_locus_annotation.py
│   ├── 04_hpap_load.py
│   ├── 05_hpap_expression.py
│   ├── 06_donor_audit.py
│   ├── 07_full_panel_fetch.py
│   ├── 08_tau_specificity.py
│   ├── 09_leukocyte_stage.py
│   ├── 10_cross_validation.py
│   ├── 11_clinical_metadata.py
│   └── 12_donor_audit_meg3.py
├── tests/                   pytest unit tests for stats + config
├── docs/substack_article.md
├── pyproject.toml
├── requirements.txt
└── .github/workflows/ci.yml
```

## Quick start

```bash
# Install
pip install -e .

# Optional: HPAP FBS decoder (only needed for full pipeline runs)
pip install --no-deps cellxgene flatbuffers

# Set project root (defaults to /kaggle/working/t1d_mech if unset)
export T1D_COA_ROOT=$PWD/workdir

# Run the pipeline
python scripts/00_env_setup.py
python scripts/01_gwas_t1d_filter.py
# … through 12_donor_audit_meg3.py
```

Each script is self-contained, single-purpose, and runs in 30 seconds to a few minutes depending on cache state. Reruns are cheap — every HPAP fetch caches the raw FlatBuffers blob to disk and skips the network round-trip on the next pass.

## Design choices

**Why a library + scripts split?** The original notebook had the same provenance helper, FBS decoder, and HPAP path constants pasted into five different cells. Moving that to `src/t1d_coa/` means each pipeline step reads like the analysis it performs, not boilerplate.

**Why tau and not a fancier specificity metric?** Tau is the standard tissue-specificity index and has a built-in failure mode that's easy to interpret: it under-calls genes that are biologically specific but have non-zero ambient expression elsewhere (the reason zero beta-cell calls land at this threshold). Other indices like the Tau-style Gini variants would recover some of those at the cost of a less interpretable threshold.

**Why donor-pseudobulk and not cell-level mixed models for the headline?** Cell-level LMMs were tested in step 09; they're underpowered for the AAB-vs-Control contrast where the within-donor variance dominates. Donor-pseudobulk Δ log1p is the unit on which both the headline and the permutation null are computed.

## Data provenance

Every artifact the pipeline writes is logged to `data/PROVENANCE.md` with a sha256 hash and the exact source URL or command that produced it. This is appended automatically — no need to maintain it by hand.

- GWAS Catalog: release v1.0.3.1, filtered to pure T1D (excludes composite traits like T1D-nephropathy).
- LD reference: MAGMA pre-built 1000 Genomes Phase 3 EUR PLINK panel.
- HPAP atlas: 222,077 cells from 67 donors via the Faryabi lab CellxGene REST endpoint. Underlying dataset: Faryabi et al., *bioRxiv* 2023.01.03.522578.

## Known limitations

1. Tau pools across all disease states, so the cross-validation against AAB-vs-Control is not strictly out-of-sample. A Control-only tau re-run is the first item on the followup list.
2. The permutation null picks uniformly random cell types; preserving the marginal distribution of called cell types would be a stricter test.
3. The headline mean is carried by ~5 strong genes. "20 of 21 genes" should be read as the omnibus passes, not as every individual gene cross-validates.
4. T1D beta-cell statistics are underpowered (6/9 donors have any surviving beta cells; HPAP084 contributes 67%).
5. The HPAP "leukocyte" bucket is resident macrophages plus a smaller mix of immune cells — not the infiltrating cytotoxic T cells that drive insulitis.

## License

MIT — see `LICENSE`. Underlying data is third-party (GWAS Catalog, HPAP/PANC-DB) under their respective terms.
