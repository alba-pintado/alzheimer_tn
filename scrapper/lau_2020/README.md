# Lau et al. 2020

Lau et al. 2020, *PNAS* (doi:10.1073/pnas.2008762117). Single-nucleus RNA-seq
of human prefrontal cortex from 21 donors (12 AD, 9 control), 179,392 nuclei.
All 21 donors are kept; the final object has 137,303 annotated nuclei.

| Script | English | Spanish |
|---|---|---|
| 1. Preprocessing | `01_lau_preprocessing.R` | `01_lau_preprocesado.R` |
| 2. Analysis | `02_lau_analysis.R` | `02_lau_analisis.R` |
| 3. UMAP diagnostic (exploratory) | — | `03_lau_diagnostico_umap.R` |
| 4. Harmony test (exploratory) | — | `04_lau_harmony.R` |

The analysis script reads the object saved by the preprocessing script, so they
must be run in that order. Scripts 3 and 4 read the analysis object and do not
modify it: they are checks of the UMAP and a test of an alternative batch
correction.

## 1. Preprocessing

| Step | What it does | Tool |
|---|---|---|
| 1 | Quality control, median ± 2.5 MAD thresholds per donor (`NUM_MADS`) | `scrapper` |
| 2 | Ambient RNA correction, run donor by donor; nuclei with contamination > 0.5 removed | `decontX` (`celda`) |
| 2.4 | Loss per donor in the contamination filter (see below) | — |
| 3 | Doublet detection on the original counts, per donor | `scDblFinder` |
| 4 | Normalisation on the decontX-corrected counts | `scrapper` |
| 5 | Removal of genes detected in fewer than 3 nuclei | — |
| 5.1 | Check of the nucleus order (see below) | — |
| 6 | Summary table and saving of the final object | — |

**AD6.** An earlier version excluded AD6 because it lost 24 % of its nuclei in
the contamination filter. In the current version it loses 9 %, like AD4
(9 %), so it no longer stands out and all 21 donors are kept. The loss per
donor is still printed to check that no donor separates from the rest.

**Stricter QC.** Thresholds use 2.5 MADs instead of the default 3. This removes
about 4,400 more nuclei, in similar proportions in AD (10.7 %) and control
(11.8 %), removes a small low-quality endothelial cluster and does not change
the structure of the data or the AD vs control results.

**Check of the nucleus order.** decontX and scDblFinder process each donor
separately and the results are joined back, so the script checks that every
column is still the same nucleus. It stops with an error if any of these fails:

1. nucleus names are unique and all exist in the original object;
2. the metadata of each nucleus (donor, condition, cell type…) is identical to the original;
3. `counts`, `decontXcounts` and `logcounts` have the nuclei in the same order;
4. in every nucleus, the decontX-corrected counts do not exceed the original counts;
5. the counts of a random sample of nuclei are identical to those in the original object.

It also saves the number of nuclei per donor at each step.

| Step | Nuclei |
|---|---|
| Original | 179,392 |
| After QC (2.5 MAD) | 159,295 |
| After contamination filter | 154,154 |
| After doublet removal | 139,430 |
| Without unannotated nuclei (analysis) | 137,303 |

29,447 genes remain. The proportion of nuclei removed is similar in AD and
control at every step, and all order checks pass.

## 2. Analysis

| Step | What it does | scrapper book |
|---|---|---|
| 0 | Removal of nuclei without cell type (NA in the annotation) | — |
| 1 | 3,000 highly variable genes, blocking by donor | ch. 3 |
| 2 | PCA (25 components) blocked by donor, and R² of each PC against cell type and donor | ch. 4, 8 |
| 3 | MNN correction by donor, used only for UMAP and clustering | ch. 8 |
| 4 | UMAP, with nuclei drawn in random order | ch. 5 |
| 5 | SNN clustering, purity, ARI and donor composition per cluster | ch. 6 |
| 6 | Markers per cluster on uncorrected logcounts | ch. 7 |
| 7 | Pseudobulk AD vs control per cell type (edgeR) and differential abundance | ch. 8 |
| 8 | Multi-omics: checked, not applicable (snRNA-seq only) | ch. 9 |

Because every donor is either AD or control, correcting by donor also removes
part of the disease effect. The MNN space is therefore used only for
visualisation and clustering; markers and AD vs control use uncorrected data.

Main results of the current version: cell type explains 90–97 % of the variance
of PC1–PC5 and donor ≤ 7 %; 30 clusters (ARI 0.30 against the authors'
annotation), none dominated by one donor; pseudobulk finds one gene in
astrocytes (PPDPF) and one in microglia (PTPRG) at FDR < 0.05; endothelial cells
tend to increase in AD (logFC 1.54, FDR 0.11).

## 3. UMAP diagnostic (exploratory)

`03_lau_diagnostico_umap.R` finds the islands of the MNN UMAP with DBSCAN
(`dbscan`, eps = 0.3, minPts = 20) and, for each island, reports size, cell
type, donors, AD fraction and quality metrics (detected genes, mitochondrial %,
decontX contamination, scDblFinder score). It also recomputes the UMAP on the
uncorrected PCA.

- Current version: 38 islands with MNN and 20 without correction. Blocking the
  PCA by donor and stricter QC do not change the number of MNN islands, so the
  fragmentation comes from MNN, which splits neuronal subtypes, not from
  nucleus quality.
- Islands contain nuclei from most donors and none is dominated by one donor.
- A group of about 1,350 excitatory nuclei with few genes and mitochondrial
  genes among its markers remains; it can be filtered before the noise
  analysis.

## 4. Harmony test (exploratory)

`04_lau_harmony.R` replaces MNN with Harmony (by donor) on the same 25-PC PCA
and recomputes R², UMAP and SNN clustering. It was run on an earlier version
of the object (3 MADs, without AD6) as a comparison; MNN is kept in the main
analysis.

| | MNN | Harmony |
|---|---|---|
| UMAP islands | 42 | 17 |
| SNN clusters | 33 | 25 |
| ARI vs authors' annotation | 0.30 | 0.44 |
| PC1 R² vs cell type | 0.977 | 0.974 |
| PC1 R² vs donor | 0.066 | 0.067 |

Markers are computed per cluster, so they would change with Harmony; pseudobulk
and differential abundance use uncorrected counts and the authors' cell types,
so they do not depend on this choice.

## Running it

From a terminal rather than an editor (the preprocessing takes over an hour and
needs about 20 GB of RAM at its peak):

```bash
Rscript 01_lau_preprocesado.R > log_preprocesado.txt 2>&1 && \
Rscript 02_lau_analisis.R > log_analisis.txt 2>&1
Rscript 03_lau_diagnostico_umap.R > log_diagnostico.txt 2>&1
Rscript 04_lau_harmony.R > log_harmony.txt 2>&1
```

Outputs:

- `results/01_preprocesado/Lau_32989152/`: contamination histogram, filtering summary and nuclei per donor.
- `results/02_analisis/Lau_32989152/`: mean–variance plot, scree plot, UMAPs, markers, pseudobulk results and differential abundance.
- `results/02_analisis/Lau_32989152/diagnostico_umap/` and `.../harmony/`: island table, UMAPs and the MNN vs Harmony comparison.
- `results/objects/lau_sce_preprocesado.rds` and `results/objects/lau_sce_analisis.rds`.

## References

- Lau SF et al. (2020). Single-nucleus transcriptome analysis reveals dysregulation of angiogenic endothelial cells and neuroprotective glia in Alzheimer's disease. *PNAS* 117:25800–25809.
- Yang S et al. (2020). Decontamination of ambient RNA in single-cell RNA-seq with DecontX. *Genome Biol* 21:57.
- Germain PL et al. (2021). Doublet identification in single-cell sequencing data using scDblFinder. *F1000Research* 10:979.
- Haghverdi L et al. (2018). Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors. *Nat Biotechnol* 36:421–427.
- Crowell HL et al. (2020). muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data. *Nat Commun* 11:6077.
- Squair JW et al. (2021). Confronting false discoveries in single-cell differential expression. *Nat Commun* 12:5692.
- Korsunsky I et al. (2019). Fast, sensitive and accurate integration of single-cell data with Harmony. *Nat Methods* 16:1289–1296.
- Hahsler M, Piekenbrock M, Doran D (2019). dbscan: Fast density-based clustering with R. *J Stat Softw* 91(1):1–30.
