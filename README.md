# t1d-celltype-of-action

**Mapping the 145 type 1 diabetes GWAS loci to the pancreatic cell types they likely act in, using cell-type specificity in the HPAP scRNA-seq atlas — and cross-validating those calls against pre-clinical (autoantibody-positive) transcriptional changes.**

Companion to the Substack article: *[Where do type 1 diabetes risk genes actually live?](docs/substack_article.md)*

---

## TL;DR

- Pulled the 145 independent T1D risk loci from the GWAS Catalog (release v1.0.3.1, May 2026) → 188 candidate genes → 150 testable in the HPAP atlas (222,077 cells, 67 donors).
- Used the tau specificity index (Yanai et al., 2005) to assign each gene a cell-type-of-action based on healthy expression specificity.
- **21 of 150 genes received a confident cell-type call:** 7 leukocyte (HLA class II + RUNX3, SH2B3, SLC15A3), 6 duct epithelial (read with caution — see below), 5 acinar (CEL, CTRB1, CTRB2, PGM1, SLC25A37), 2 mesenchymal, 1 endothelial. **Zero confidently beta cell** — a limitation of tau's exclusivity requirement, not an absence of beta-cell biology in the locus set.
- **The headline result:** the cell-type-of-action calls cross-validate against pre-clinical (AAB-vs-Control) transcriptional changes. Across the 21 confidently-called genes, 20 had a called-cell-type AAB effect larger than the average effect across the other compartments; the omnibus mean called-cell-type advantage of +0.18 log1p exceeded every shuffled mean across 10,000 permutations (p < 10⁻⁴).
- **Honest caveats:** the cross-validation isn't fully out-of-sample (AAB cells contributed ~14% of the matrix used to compute tau), the permutation null is weaker than ideal (uniform over cell types rather than preserving the marginal of called cell types), and the effect is driven disproportionately by ~5 strong genes (CEL, CTRB1, HLA-DRB1, HLA-DQA1, PLPP1) with the rest weakly consistent.

## Headline figures

### Tau distribution and cell-type-of-action calls

![tau distribution and calls](figures/fig03_tau_distribution_and_calls.png)

Of the 150 testable genes, 34 sat above the expression floor. Of those, 21 had tau ≥ 0.70 (a confident single-compartment call). The other 13 were broadly expressed across compartments. The remaining 116 of the original 150 were below detection in any islet cell type.

### The 21 confident calls, ranked by tau

![confident calls by compartment](figures/fig04_confident_calls_ranked.png)

The leukocyte cluster at the top is the antigen-presentation machinery (HLA class II) plus three lymphocyte-program transcription factors. The acinar cluster is pancreatic digestive enzymes — interpretable in light of the nPOD consortium's finding that T1D pancreata show acinar atrophy. The duct cluster is mixed: CFTR is real biology, but GLIS3, CTSH, and NOTCH2 are independently known beta-cell genes whose duct-call here probably reflects scRNA-seq ambient-mRNA contamination plus tau's preference for the most abundant non-endocrine epithelial compartment.

### Cross-validation against AAB-stage effects

![cross-validation result](figures/fig05_cross_validation.png)

Left panel: histogram of per-gene called-cell-type advantage (called-CT AAB effect minus mean of other CTs) across the 21 confidently-called genes. Right panel: per-gene scatter of called-CT effect vs the largest effect in any other compartment. Genes above the y=x line are cases where the called compartment shows the strongest disease-stage effect; points below the line are cases where some other compartment beats it. The mean across all 21 genes is +0.18 log1p (red dashed line in the histogram), shifting the whole distribution to the right of zero.

Stratified by called cell type, the acinar genes are doing most of the work (mean advantage +0.33 across 5 genes), with leukocyte genes a clear second (+0.17, 7 genes). Duct and mesenchymal calls contribute weakly. The headline result is carried by a minority of strong genes rather than uniformly across the panel — see the article for the longer hedge.

---

## Repository layout

```
.
├── README.md                       this file
├── LICENSE                         MIT
├── .gitignore
├── notebook.ipynb                  cleaned Kaggle notebook (20 cells, final pipeline only)
├── push_to_github.sh               interactive push helper
├── code/                           13 modular pipeline scripts (00–12)
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
├── tools/
│   ├── make_figures.py             regenerate all 9 article figures
│   └── make_data.py                regenerate the processed data tables
├── data/
│   └── processed/
│       ├── T1D145_celltype_specificity.tsv
│       ├── T1D145_celltype_call_crossvalidation.tsv
│       └── donor_metadata_t1d_aab.tsv
├── docs/
│   └── substack_article.md         the article, ready to paste into Substack
└── figures/
    ├── fig01_pilot_heatmap_pooled.png
    ├── fig02_pilot_heatmap_by_disease.png
    ├── fig03_tau_distribution_and_calls.png
    ├── fig04_confident_calls_ranked.png
    ├── fig05_cross_validation.png
    ├── fig06_leukocyte_by_state.png
    ├── fig07_beta_t1d_vs_ctrl_bar.png
    ├── fig08_beta_aab_vs_t1d_scatter.png
    └── fig09_meg3_per_donor.png
```

## Reproducing the figures and tables

The figures and processed-data tables in `figures/` and `data/processed/` can be regenerated end-to-end with no external data dependencies:

```bash
pip install pandas numpy matplotlib
python tools/make_figures.py
python tools/make_data.py
```

The summary statistics that drive these scripts are hard-coded from the original notebook outputs, so anyone can rebuild the article-quality figures locally without re-fetching the HPAP atlas.

## Reproducing the full analysis from raw data

The hard part isn't the analysis, it's loading the data. HPAP's CellxGene instance returns expression payloads as raw FlatBuffers (not Arrow IPC), and the schema isn't documented publicly. The `05_hpap_expression.py` module shows the path that worked: install the `cellxgene` package solely for its `server.common.fbs.matrix.decode_matrix_fbs` helper.

```bash
# Minimum dependencies
pip install pandas numpy scipy statsmodels pyarrow httpx tenacity tqdm matplotlib Pillow
pip install --no-deps cellxgene flatbuffers   # only for the FBS decoder

# Run the pipeline
cd code/
python 00_env_setup.py
python 01_gwas_t1d_filter.py
# ... through 12_donor_audit_meg3.py
```

The pipeline expects ~21 GB of disk space (mostly for the GWAS Catalog bulk TSV and the cached HPAP expression blobs) and runs end-to-end in 30–45 minutes on a Kaggle CPU instance.

## Pushing your edits to GitHub

`push_to_github.sh` is interactive — run it and it'll prompt for username, repository name, and a Personal Access Token. It verifies the token, confirms the repo exists, warns before force-pushing over existing content, and scrubs the token from `.git/config` once the push succeeds. No file ever contains the token.

```bash
bash push_to_github.sh
```

## Data provenance

- **GWAS Catalog:** release v1.0.3.1 (May 2026), filtered to MONDO:0005147 (pure T1D, excluding composite autoimmune-disease studies and T1D nephropathy).
- **LD reference:** 1000 Genomes Phase 3 EUR panel via the MAGMA bundle, used for PLINK clumping (`--clump-r2 0.1 --clump-kb 1000 --clump-p1 5e-8 --clump-p2 1e-5`).
- **HPAP atlas:** 222,077 cells from 67 donors, accessed through the Faryabi lab's CellxGene REST API. Underlying dataset: Faryabi et al., *bioRxiv* 2023.01.03.522578.
- **Clinical metadata:** HPAP Supplementary Table S2 from the same preprint (donor demographics, autoantibody panels, HbA1c, C-peptide, disease duration).

## Limitations the article doesn't fully cover

1. **Tau calls aren't fully out-of-sample for the cross-validation.** The tau matrix pools across all four disease states with cell-count weighting. AAB contributes ~14% of cells; the cross-validation against AAB-vs-Control effects is therefore not strictly out-of-sample. A clean re-run on Control-only cells is the first item on the article's "what I'd want to see" list.
2. **The permutation null is weaker than ideal.** It picks a uniformly random cell type as "called" from each gene's available cell types, but the empirical distribution of tau calls concentrates in the same compartments where AAB effects are largest. A null that preserved the marginal distribution of called cell types would be a stronger test.
3. **The effect is driven by ~5 strong genes.** CEL (+0.71), CTRB1 (+0.51), HLA-DRB1 (+0.50), HLA-DQA1 (+0.35), PLPP1 (+0.23) carry the mean; the remaining 16 genes cluster between 0 and +0.15. "20 of 21 genes" should be read as the omnibus passes, not as every individual gene cross-validates.
4. **T1D beta-cell statistics are underpowered.** 6 of 9 T1D donors have any surviving beta cells at all; HPAP084 (an "Unsuspected" pre-clinical case) contributes 67% of the pool. No FDR-corrected significant hits in beta cells.
5. **HPAP "leukocyte" isn't infiltrating T cells.** It's resident macrophages plus a smaller mix of other immune cells. The T cell-mediated autoimmune attack that defines T1D is happening in lymph nodes and at the islet edge during a window of disease activity that's mostly closed by the time HPAP donors die.
6. **Duct compartment calls are a mix of biology and artifact.** See the article's longer hedge.

## Citation

If this is useful in your own work, cite the Substack article and this repository together. The methodological contribution — using tau-derived cell-type-of-action calls with a built-in cross-validation against disease-stage data — is what would carry over to other GWAS-and-scRNA pairings.

## License

MIT for the code. The data is third-party (GWAS Catalog, HPAP/PANC-DB) under their respective terms.
