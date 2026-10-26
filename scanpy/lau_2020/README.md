# Lau et al. 2020 — export to Python

Export of the processed R object to AnnData (`.h5ad`), with the same nuclei and
genes as R (137,303 nuclei × 29,447 genes, 21 donors: 12 AD and 9 control). The `.h5ad`
is the input for the transcriptional noise analysis with DECIBEL, which is
written in Python.

| Script | Language | What it does |
|---|---|---|
| `00_exportar_h5ad.R` | R | Writes the counts and logcounts matrices piece by piece to HDF5 (`rhdf5`), plus nucleus metadata, genes and R's PCA/UMAP as CSV |
| `00b_construir_h5ad.py` | Python | Rebuilds the matrices and saves `lau_sce_analisis.h5ad` |

Scripts must be run in that order, from the project root.

## Why the export has two steps

`zellkonverter::writeH5AD` copies the whole object into an embedded Python
before writing, so two copies of the matrices sit in memory at once and the
process was killed for lack of RAM. Here R writes each sparse matrix
(`x`, `i`, `p` slots of the `dgCMatrix`) to HDF5 one at a time, removing it from
the object before the next, and Python assembles the AnnData afterwards. Each
step holds a single copy in memory.

The resulting `.h5ad` contains:

- `X`: logcounts (nuclei × genes, sparse, float32)
- `layers["counts"]`: raw counts
- `obs`: `donor`, `condition`, `celltype`, `cluster` (R clusters), `subsets_mito_percent`, `sum`, `detected`, `decontX_contamination`
- `var`: `simbolo` (gene symbol; Ensembl ID when there is none)
- `obsm`: `X_pca_R`, `X_umap_R` (coordinates from the R analysis)

`00b_construir_h5ad.py` checks that names are unique, that there are 21 donors.

## Environments

- **R step:** the system R (4.6.1, with `scrapper`), not the conda R. Run
  `conda deactivate` first and check with `Rscript -e 'R.version.string'`.
  Needs `rhdf5` (`BiocManager::install("rhdf5")`).
- **Python steps:** a conda environment with `scanpy`, `anndata`, `h5py`,
  for example the one in `envs/env_python_scanpy.yml`:

```bash
conda env create -f envs/env_python_scanpy.yml
conda activate alzheimer_python
python -c "import scanpy, anndata, h5py; print('ok')"
```

## Running it

Set the `proyecto` path at the top of `00_exportar_h5ad.R` first. The Python
scripts use the current directory, so run everything from the project root.

```bash
# 1. R: export to HDF5 (system R)
conda deactivate
nohup Rscript scanpy/lau_2020/00_exportar_h5ad.R > log_exportar.txt 2>&1 &

# 2. Python: build the .h5ad
conda activate alzheimer_python
python scanpy/lau_2020/00b_construir_h5ad.py > log_construir.txt 2>&1
```

Step 1 needs
about 10 GB of RAM at its peak; close other programs while it runs.

## Outputs

- `results/objects/lau_2020/lau_export.h5` and `lau_export_*.csv`: intermediate files from step 1; they can be deleted once the `.h5ad` exists.
- `results/objects/lau_2020/lau_sce_analisis.h5ad`: the object for Python and DECIBEL.

## References

- Virshup I et al. (2024). anndata: Access and store annotated data matrices. *J Open Source Softw* 9:4371.
- Wolf FA, Angerer P, Theis FJ (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol* 19:15.