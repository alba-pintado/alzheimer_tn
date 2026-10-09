# Lau et al. 2020 — snRNA-seq preprocessing

Part of the MSc thesis *Transcriptional noise in Alzheimer's disease
progression* (CBMSO). This folder contains the preprocessing of one dataset.
The transcriptional noise analysis lives in a separate script.

Script: `scrapper_lau.R`

## Dataset

Lau et al. 2020, *PNAS* (doi:10.1073/pnas.2008762117). Single-nucleus RNA-seq
of human prefrontal cortex from 21 donors (12 Alzheimer's disease, 9 control),
with 179,392 nuclei and seven annotated cell types.

## Steps

| Step | What it does | Tool |
|---|---|---|
| 1 | Quality control, median ± 3 MAD thresholds per donor | `scrapper` |
| 2 | Ambient RNA correction and filtering (contamination > 0.5) | `decontX` (`celda`) |
| 3 | Doublet detection on the original counts, per donor | `scDblFinder` |
| 4 | Normalisation on the decontX-corrected counts | `scrapper` |
| 5 | Removal of genes detected in fewer than 3 nuclei | — |
| 6 | Summary table and saving of the final object | — |

decontX is run donor by donor to keep memory use low; this gives the same
result as running it with one batch per donor.

## Results

| Step | Nuclei |
|---|---|
| Original | 179,392 |
| After QC | 163,710 |
| After contamination filter | 158,231 |
| After doublet removal | 142,997 |

After gene filtering, 29,467 genes remain. The proportion of nuclei removed is
similar in AD and control at every step.

## Running it

Open the script in RStudio, set the `project` variable at the top to your
local path, and run it (Ctrl+A, Ctrl+Enter). The full run takes over an hour
and needs about 20 GB of RAM at its peak.

Outputs:

- `results/01_preprocessing/Lau_32989152/01_contamination_histogram.png`
- `results/01_preprocessing/Lau_32989152/02_filtering_summary.csv`
- `results/objects/lau_sce_preprocessed.rds`

Data and results are not tracked in this repository.

## Requirements

R ≥ 4.3 with `SingleCellExperiment`, `scrapper`, `celda`, `scDblFinder` and
`BiocParallel` (all from Bioconductor).

## References

- Lun A. *Single-cell RNA-seq analysis with scrapper*. https://bioconductor.org/books/release/scrapbook/
- Amezquita RA et al. (2020). Orchestrating single-cell analysis with Bioconductor. *Nat Methods* 17:137–145.
- Yang S et al. (2020). Decontamination of ambient RNA in single-cell RNA-seq with DecontX. *Genome Biol* 21:57.
- Germain PL et al. (2021). Doublet identification in single-cell sequencing data using scDblFinder. *F1000Research* 10:979.
- Lau SF et al. (2020). Single-nucleus transcriptome analysis reveals dysregulation of angiogenic endothelial cells and neuroprotective glia in Alzheimer's disease. *PNAS* 117:25800–25809.