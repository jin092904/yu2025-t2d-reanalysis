# Reanalysis of Yu et al. (2025), *Sci. Rep.* 15:44904: type 2 diabetes and obesity in adipose biomarker data

Code and results for the Matters Arising comment "Type 2 diabetes is confounded with obesity in the data used to propose AIM2 and FHOD3 as adipose biomarkers" (H. Lee, submitted to *Scientific Reports*), arising from Yu, Y.-S., Lee, D. S., Lim, J. H. & Lee, Y. J. Identifying potential biomarkers for type 2 diabetes in the adipose tissue of older adults via multiple machine learning algorithms. *Sci. Rep.* **15**, 44904 (2025). https://doi.org/10.1038/s41598-025-29141-9

## What is reanalysed

- The published list of 210 differentially expressed genes (DEGs; Supplementary Tables S1–S2 of Yu et al.) under four ComBat-seq/DESeq2 settings
- Feature selection by LASSO, SVM-RFE and random forest in the training–test split and over 100 splits of two kinds, with single-gene test-set AUCs
- DEG selection from the training samples only (leakage)
- Within-cohort comparisons of lean, obese non-T2D and obese T2D participants in GSE141432 and GSE166047, without batch correction
- An independent Python re-implementation of each step, run on the same splits and cross-validation folds

## Contents

- `src/r/`, `src/py/`: the R implementation (primary) and the Python re-implementation
- `data/raw/MANIFEST.tsv`: name, size, MD5 checksum, URL and download time of every input file (GEO files and the original article's supplementary files); the files themselves are downloaded by `src/py/00_download.py` and are not redistributed here
- `data/processed/samples.tsv`: all 91 GEO samples of the three series with their groups; the 31 analysed samples are the rows with a `paper_group` value and `in_ncbi_matrix` = `True`
- `data/processed/`: the 31-sample count matrix (columns taken from the NCBI-generated count matrices), ComBat-seq outputs, machine-learning inputs and split assignments
- `results/tables/`: every number reported in the comment; `02_paper_vs_R_group_pc.tsv` includes the published values of Supplementary Tables S1–S2 of Yu et al. (CC BY 4.0) for comparison
- `manuscript/figures/`: Figures 1–2 and Supplementary Figure S1 (vector PDF, 600-dpi PNG and JPEG)

## Figures and tables

- Figure 1 (FHOD3 and AIM2 within cohorts): `05_within_cohort.R` → `05_within_cohort_aim2_fhod3_expr.tsv`, `05_within_cohort_aim2_fhod3_tests_R.tsv`
- Figure 2a (DEGs by setting): `02_combat_deseq2.R`, `02_compare.py` → `02_deg_by_spec.tsv`
- Figure 2b (genes selected by all three methods over 100 splits): `03_ml.R stability` with both split types → `03_seed_stability_vst_table1_{sample,pergroup}_no.tsv`
- Figure 2c (leakage): `05_leakage.R` → `05_leakage_100seeds.tsv`
- Supplementary Figure S1 (published DEGs by direction within cohorts): `05_within_cohort.R`, `05_figures.py` → `05_within_cohort_deg210_summary.tsv`
- Comparison with Supplementary Tables S1–S2: `02_compare.py` → `02_paper_vs_R_group_pc.tsv`
- R–Python agreement: `02_deg_by_spec.tsv`, `03_table2_comparison.tsv`, `04_roc_R_vs_py.tsv`, `05_within_cohort_aim2_fhod3_tests_py.tsv`
- All figures are drawn by `src/r/07_manuscript_figures.R` from these tables.

## How to run (from the repository root)

```bash
python src/py/00_download.py                  # GEO files and the original article's files (into paper/), checksums into MANIFEST.tsv
python src/py/01_data.py                      # sample table, 31-sample count matrix
Rscript src/r/02_combat_deseq2.R              # pre-specified: ComBat-seq batch only, no gene filter
Rscript src/r/02_combat_deseq2.R pc           # protein-coding genes only
Rscript src/r/02_combat_deseq2.R group        # disease-label covariate in ComBat-seq
Rscript src/r/02_combat_deseq2.R group pc     # both (reproduces Supplementary Tables S1-S2)
python src/py/02_combat_deseq2.py
python src/py/02_combat_deseq2.py group pc
python src/py/02_compare.py                   # DEG counts by setting, comparison with the published list
Rscript src/r/03_ml_inputs.R                  # ML inputs for the 210 published DEGs
Rscript src/r/03_ml.R                         # pre-specified split (set.seed(100); simple random), LASSO, SVM-RFE, random forest
Rscript src/r/03_ml.R split=pergroup          # split within each group (65% of each group)
Rscript src/r/03_ml.R search                  # LASSO lambda under 32 combinations of input, order, split and re-seeding
Rscript src/r/03_ml.R stability               # seeds 1-100, simple random splits (8 cores)
Rscript src/r/03_ml.R stability split=pergroup  # seeds 1-100, within-group splits (8 cores)
Rscript src/r/03_lambda_path_check.R
python src/py/03_ml.py                        # Python re-implementation on the R split and folds
python src/py/03_compare.py
Rscript src/r/04_roc.R                        # AIM2 and FHOD3 AUCs, test-set AUCs over 100 splits
python src/py/04_roc.py
Rscript src/r/05_within_cohort.R              # within-cohort comparisons
Rscript src/r/05_leakage.R                    # DEG selection from training samples only (8 cores)
python src/py/05_figures.py                   # Python check of the within-cohort tests, direction summary
Rscript src/r/07_pergroup_split.R             # summary of within-group splits
python src/py/07_split_composition.py         # simple random splits grouped by test composition
Rscript src/r/07_gene_auc_100splits.R         # number of splits with a test AUC of 1, per gene
Rscript src/r/07_manuscript_figures.R         # Figures 1-2 and Supplementary Figure S1 (PDF, PNG)
python src/py/07_figure_jpeg.py               # JPEG copies of the figures
```

- Run R scripts from the repository root (they use relative paths). The 8-core steps use `parallel::mclapply`, which does not run in parallel on Windows.
- `02_compare.py`, `03_compare.py`, `04_roc.py` and `05_figures.py` also write working figures with Korean labels to `results/figures/`; these are not included.
- Every random step uses a fixed seed.

## Environment

Python 3.10 with `requirements.txt`; R 4.6.1 with the packages and versions listed at the top of `requirements.txt`. Figures use the cairo graphics device.

## Use of AI tools

Claude (Anthropic) assisted with writing and running the analysis code and with drafting and revising the manuscript, and Gemini (Google) with literature and data searches and with language editing and translation. All numbers in the manuscript were checked against `results/tables/` with `src/py/07_check_numbers.py` (which reads the manuscript text, not included here), and the author reviewed every analysis setting and takes full responsibility for the content.

## Notes

- Code comments and console messages are in Korean because the work began as a university course project.

## License

Code: MIT (see `LICENSE`).
