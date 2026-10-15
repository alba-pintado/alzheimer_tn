# Transcriptional noise in Alzheimer's disease progression

Code for the MSc thesis *Transcriptional noise in Alzheimer's disease
progression* (CBMSO). Alba Pintado, MSc in Bioinformatics and Data Science
(ISCIII).

The aim is to measure transcriptional noise (cell-to-cell variability in gene
expression) in single-nucleus RNA-seq of the human brain and compare it between
Alzheimer's disease (AD) and control donors.

## Contents

| Folder | Dataset | What it does |
|---|---|---|
| `scrapper/grubman_2019` | Grubman et al. 2019, *Nat Neurosci* | Analysis with the scrapper book |
| `scrapper/lau_2020` | Lau et al. 2020, *PNAS* | Preprocessing and analysis with the scrapper book (chapters 1–9) |

Each folder has its own README. The Lau scripts are provided in English and in
Spanish; the code is identical, only comments and plot labels change.

## Requirements

R ≥ 4.3 with `SingleCellExperiment`, `scrapper`, `celda`, `scDblFinder`,
`scater`, `edgeR`, `org.Hs.eg.db`, `AnnotationDbi` and `BiocParallel`
(Bioconductor), and `mclust` and `ggplot2` (CRAN).

## Data

Data and results are not tracked in this repository. Set the `proyecto`
variable at the top of each script to the local project path.

## References

- Lun A. *Single-cell RNA-seq analysis with scrapper*. https://bioconductor.org/books/release/scrapbook/
- Amezquita RA et al. (2020). Orchestrating single-cell analysis with Bioconductor. *Nat Methods* 17:137–145.