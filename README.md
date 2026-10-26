# Transcriptional noise in Alzheimer's disease progression

Code for the MSc thesis *Transcriptional noise in Alzheimer's disease
progression* (CBMSO). Alba Pintado, MSc in Bioinformatics and Data Science
(ISCIII).

The aim is to measure transcriptional noise (cell-to-cell variability in gene
expression within a cell type) in single-nucleus RNA-seq of the human brain and
compare it between Alzheimer's disease (AD) and control donors.

## Datasets

| Dataset | Tissue | Donors | Nuclei (raw → analysed) | Folder |
|---|---|---|---|---|
| Grubman et al. 2019, *Nat Neurosci* | Entorhinal cortex | 12 (6 AD, 6 control) | 13,214 (2,171 astrocytes) | `scrapper/grubman_2019` |
| Lau et al. 2020, *PNAS* | Prefrontal cortex | 21 (12 AD, 9 control) | 179,392 → 137,303 | `scrapper/lau_2020`, `scanpy/lau_2020` |

**Grubman et al. 2019.** Analysis of the astrocytes following the scrapper book.
Disease and batch are confounded in this design, so it was used to build and
test the pipeline.

**Lau et al. 2020.** The main dataset. The workflow has three parts, run in
this order:

1. **R (`scrapper/lau_2020`).** Preprocessing (QC with median ± 2.5 MAD per
   donor, ambient RNA with decontX, doublets with scDblFinder, normalisation,
   gene filtering and a check that every nucleus keeps its metadata) and
   analysis (chapters 3–9 of the scrapper book: HVGs, PCA blocked by donor, MNN,
   UMAP, clustering, markers, pseudobulk AD vs control and differential
   abundance). Nuclei without a cell type are removed. All 21 donors are kept.
2. **Export to Python (`scanpy/lau_2020`).** The final R object is written to
   `.h5ad` in two steps to fit in memory. The `.h5ad` is the input for the
   Python analysis and for DECIBEL.
3. **Transcriptional noise (DECIBEL).** *Work in progress:* the Python analysis
   and the noise scripts will be added once they have been run and reviewed.

See the README in each folder for details.

## Current status (27-10-2026)

- Final object: 137,303 nuclei × 29,447 genes, 21 donors (12 AD, 9 control).
- Cell type explains 90–97 % of the variance of PC1–PC5; donor ≤ 7 %.
- Pseudobulk: one gene in astrocytes (PPDPF) and one in microglia (PTPRG) at
  FDR < 0.05; endothelial cells tend to increase in AD (FDR 0.11).
- Next: re-run the Python analysis on the final `.h5ad` and run DECIBEL on all
  cell types (target: end of November 2026).

## Requirements

- R ≥ 4.3 with `SingleCellExperiment`, `scrapper`, `celda`, `scDblFinder`,
  `scater`, `edgeR`, `org.Hs.eg.db`, `AnnotationDbi`, `BiocParallel`, `rhdf5`
  and `dbscan`, plus `mclust` and `ggplot2` (CRAN).
- Python ≥ 3.10 with `scanpy`, `anndata`, `h5py`, `harmonypy`, `leidenalg`,
  `pydeseq2` and `statsmodels` (see `envs/`). DECIBEL is not on pip: clone it
  into `decibel/decibel_repo` (`git clone https://gitlab.com/olgaibanez/decibel.git decibel/decibel_repo`).

## Data

Data, results, objects (`.rds`, `.h5ad`, `.qs`), logs and the cloned DECIBEL
repository are not tracked (`.gitignore`). Set the `proyecto` variable at the
top of each R script to the local project path; Python scripts use the current
directory, so run everything from the project root, for example:

```bash
Rscript scrapper/lau_2020/01_lau_preprocesado.R > log_preprocesado.txt 2>&1
```

## Key references

**Datasets**

- Grubman A et al. (2019). A single-cell atlas of entorhinal cortex from individuals with Alzheimer's disease reveals cell-type-specific gene expression regulation. *Nat Neurosci* 22:2087–2097. [doi:10.1038/s41593-019-0539-4](https://doi.org/10.1038/s41593-019-0539-4)
- Lau SF et al. (2020). Single-nucleus transcriptome analysis reveals dysregulation of angiogenic endothelial cells and neuroprotective glia in Alzheimer's disease. *PNAS* 117:25800–25809. [doi:10.1073/pnas.2008762117](https://doi.org/10.1073/pnas.2008762117)

**Alzheimer's disease at single-cell resolution**

- Mathys H et al. (2019). Single-cell transcriptomic analysis of Alzheimer's disease. *Nature* 570:332–337. [doi:10.1038/s41586-019-1195-2](https://doi.org/10.1038/s41586-019-1195-2)
- Gazestani V et al. (2023). Early Alzheimer's disease pathology in human cortex involves transient cell states. *Cell* 186:4438–4453. [doi:10.1016/j.cell.2023.08.005](https://doi.org/10.1016/j.cell.2023.08.005)
- Mathys H et al. (2024). Single-cell multiregion dissection of Alzheimer's disease. *Nature* 632:858–868. [doi:10.1038/s41586-024-07606-7](https://doi.org/10.1038/s41586-024-07606-7)

**Transcriptional noise**

- Ibañez-Solé O et al. (2022). Lack of evidence for increased transcriptional noise in aged tissues. *eLife* 11:e80380. [doi:10.7554/eLife.80380](https://doi.org/10.7554/eLife.80380)
- Bartz J et al. (2023). Progress in discovering transcriptional noise in aging. *Int J Mol Sci* 24:3701. [doi:10.3390/ijms24043701](https://doi.org/10.3390/ijms24043701)
- Rosales-Alvarez RE et al. (2023). VarID2 quantifies gene expression noise dynamics and unveils functional heterogeneity of ageing hematopoietic stem cells. *Genome Biol* 24:148. [doi:10.1186/s13059-023-02974-1](https://doi.org/10.1186/s13059-023-02974-1)
- Amit G et al. (2022). Global coordination level in single-cell transcriptomic data. *Sci Rep* 12:7547. [doi:10.1038/s41598-022-11507-y](https://doi.org/10.1038/s41598-022-11507-y)

**Methods**

- Lun A. *Single-cell RNA-seq analysis with scrapper*. https://bioconductor.org/books/release/scrapbook/
- Amezquita RA et al. (2020). Orchestrating single-cell analysis with Bioconductor. *Nat Methods* 17:137–145. [doi:10.1038/s41592-019-0654-x](https://doi.org/10.1038/s41592-019-0654-x)
- Wolf FA, Angerer P, Theis FJ (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol* 19:15. [doi:10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0)
