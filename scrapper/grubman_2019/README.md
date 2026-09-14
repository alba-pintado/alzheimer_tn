# Astrocytes from Grubman et al. 2019 — `scrapper` pipeline

One component of the MSc thesis *Transcriptional noise in Alzheimer's disease
progression* (CBMSO). This folder covers a single analysis: the full
single-nucleus RNA-seq workflow applied to cortical astrocytes from one dataset.
Other datasets and other cell types live in sibling folders.

## Dataset

Grubman et al. 2019, *Nature Neuroscience* (PMID 31768052). Single-nucleus
RNA-seq of human entorhinal cortex, 6 libraries (3 Alzheimer's disease, 3
control). Input to this script is the astrocyte subset after cell-type
assignment, stored as a `SingleCellExperiment`.

## Pipeline

`scrapper_grubman.R` follows the `scrapper` book
(https://bioconductor.org/books/release/scrapbook/), chapter by chapter:

| Block | Chapter | What it does |
|---|---|---|
| 1 | 1 | Quality control, per-block median ± 3 MAD thresholds |
| 2 | 2 | Size factors per library, log-normalisation |
| 3 | 3 | Highly variable genes from the mean-variance residual |
| 4 | 4 | PCA, plus variance attribution per component |
| 5-7 | 5-7 | Visualisation, SNN clustering and markers, uncorrected |
| 8 | 8 | MNN batch correction, then the same three steps again |
| 8.5 | 8.5 | Pseudobulk differential expression and abundance (`edgeR`) |
| 9 | 9 | Protein multi-omics — not applicable, RNA only |

## Experimental design, and why it constrains the conclusions

Every library in Grubman is either 100 % AD or 100 % control. Batch and
diagnosis are therefore the same variable, and the effective sample size is 3
vs 3, not the ~2000 nuclei that the matrix suggests. Two consequences are
handled explicitly in the code:

- MNN correction removes part of the disease effect along with the batch
  effect. The corrected space is used only to visualise and to cluster;
  expression is always measured on the normalised counts.
- Library cannot enter the differential expression model as a covariate — the
  design matrix loses rank. The script computes the rank to show this rather
  than asserting it.

Grubman validates the pipeline. Disease conclusions require cohorts where batch
and condition are not confounded (Mathys, Gazestani).

## Running it

```bash
conda env create -f ../../envs/alzheimer_snRNA.yml
conda activate alzheimer_snRNA
Rscript scrapper_grubman.R
```

Set the `project` variable at the top of the script to your local path first.

Figures and tables are written to `results/01_exploration/Grubman_31768052/`.
Data, results and PDFs are not tracked in this repository.

## Requirements

R ≥ 4.3 with `SingleCellExperiment`, `scrapper`, `scater`, `edgeR`, `ggplot2`,
`gridExtra` and `Matrix`.