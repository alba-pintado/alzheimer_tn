# Decibel

Measuring transcriptional noise in single-cell data. Check our manuscript [Lack of evidence for increased transcriptional noise in aged tissues](https://www.biorxiv.org/content/10.1101/2022.05.18.492432v1).

## Official documentation

https://decibel.readthedocs.io/en/latest/index.html

## Getting started

Clone repository:
```python
git clone https://gitlab.com/olgaibanez/decibel.git
```
Our implementation of the methods for transcriptional noise quantification from single-cell RNAseq data can be found in ```module```. 
```python
cd decibel/module
```

## Requirements
Official documentation for Decibel can be found [here](https://decibel.readthedocs.io/en/latest/). Decibel is designed to be used within the [scanpy](https://scanpy.readthedocs.io/en/stable/) framework. Check documentation for [scallop](https://scallop.readthedocs.io/en/latest/) (our method for transcriptional noise quantification), and [triku](https://triku/readthedocs.io/en/latest/) (our feature selection method for single-cell RNAseq data). 

```pandas```
```numpy``` 
```scanpy```
```scanpy.external```
```tqdm```
```scallop```
```triku```


## Basic usage

Import required modules:

``` python
import pandas as pd
import numpy as np
import scanpy as sc
import scanpy.external as sce
import tqdm as tqdm
import scallop as sl
import triku as tk
import sys
sys.path.append('path_tp_decibel/')
import decibel as dcb

```
Load and process data (normalization, log-transformation, QC filtering on cells and genes): 

``` python
adata = sc.read('path_to_file/filename.h5ad')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=100)
```
Run PCA, feature selection (triku), batch effect correction (harmony) and dimensionality reduction (UMAP). 

``` python
sc.pp.pca(adata)
sc.pp.neighbors(adata)
tk.tl.triku(adata)
sc.pp.pca(adata)
sce.pp.harmony_integrate(adata, 'batch')
sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.umap(adata)
```

Compute transcriptional noise as in Enge et al, (2017):
``` python
dcb.enge_transcriptional_noise(adata, 'batch')

```
## Functions

**enge_transcriptional_noise**

It computes the transcriptional noise as the biological variation over the technical variation, as in Enge et al (2017).
This can only be computed in datasets with ERCC spike-ins. The biological variation is computed as the correlation distance 
between each cell and the average gene expression of all the cells of the same cell type and the same batch. The technical variation is 
computed as the correlation distance between each cell and the mean ERCC spike-in expression of all the cell of the same cell type and batch. 
The transcriptional noise is computed as the biological variation over the technical variation. 

---

**distance_to_celltype_mean**

It computes the distance between each cell and the mean expression of its cell type in the 
same batch (donor/mouse). It computes three distances: euclidean, correlation and manhattan.

---

**enge_euclidean_dist**

It computes the Euclidean distance to the average expression across cell types using a set of invariant genes. 
The invariant genes are selected as follows: 
1. Create equally sized bins of genes according to their mean expression
2. Discard the two most extreme bins (lowest and highest mean expression)
3. Select the 10% with the lowest coefficient of variation within each of the remaining bins

---

**pairwise_euclidean_sample**
    
It computes the euclidean distance between n randomly selected cells of the same cell type.

---

**hernando_herraez**

It computes the correlation distance of each cell to the cell type median using the 500 most variably expressed genes.

---

**distance_to_celltype_mean_invariant**

It computes the distance between each cell and the mean expression of its cell type in the same batch (donor/mouse), using a set of invariant genes as in Enge (2017). It computes three distances: euclidean, correlation and manhattan.

---

**gcl**

It computes the global coordination level (GCL), following the original GCL.m script provided by the authors (https://github.com/guy531/gcl).

---

**gcl_per_cell_type_and_batch**

It computes the GCL for each cell type and batch in ```adata.obs['batch'].```

---

**rerun_preprocessing**

It re-runs preprocessing steps: filter lowly expressed genes, compute HVGs, run batch-effect corrected PCA (harmony), neighbors.

---

**scallop_pipeline**
    
It computes transcriptional noise as 1 - membership score (averaged over aa range of resolution values). 
It runs the whole Scallop pipeline: 
1. Create separate annData object per condition (young/old, smoker/non-smoker) and cell type. 
2. Re-run preprocessing on annData
3. Run Scallop over range of resolution values
4. Compute average membership score across resolutions
5. Transcriptional noise = 1 - mean membership score

## How to cite

Lack of evidence for increased transcriptional noise in aged tissues
Olga Ibáñez-Solé, Alex M. Ascensión, Marcos J. Araúzo-Bravo, Ander Izeta
bioRxiv 2022.05.18.492432; doi: https://doi.org/10.1101/2022.05.18.492432 
    
## Support
olgaibanezsole@gmail.com

