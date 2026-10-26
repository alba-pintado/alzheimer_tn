# Lau et al. 2020: construir el .h5ad a partir de lo exportado por 00_exportar_h5ad.R
#
# Ejecutar desde la raíz del proyecto, con el entorno de Python activado:
#   python scanpy/lau_2020/00b_construir_h5ad.py > log_construir.txt 2>&1

import os
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

base = os.path.join(os.getcwd(), "results/objects/lau_2020")

f = h5py.File(os.path.join(base, "lau_export.h5"), "r")
n_genes, n_nucleos = (int(v) for v in f["shape"][:])
print("genes:", n_genes, "| núcleos:", n_nucleos)

def leer(nombre):
    """Matriz dispersa de R (genes x núcleos, por columnas) -> núcleos x genes."""
    g = f[nombre]
    m = sp.csc_matrix((g["x"][:], g["i"][:], g["p"][:]), shape=(n_genes, n_nucleos))
    return m.T.tocsr().astype(np.float32)

logcounts = leer("logcounts")
counts = leer("counts")
f.close()

obs = pd.read_csv(os.path.join(base, "lau_export_obs.csv"), index_col=0)
var = pd.read_csv(os.path.join(base, "lau_export_var.csv"), index_col=0)
assert obs.shape[0] == n_nucleos and var.shape[0] == n_genes

adata = ad.AnnData(X=logcounts, obs=obs, var=var, layers={"counts": counts})
adata.obsm["X_pca_R"] = pd.read_csv(os.path.join(base, "lau_export_pca_R.csv"), index_col=0).values
adata.obsm["X_umap_R"] = pd.read_csv(os.path.join(base, "lau_export_umap_R.csv"), index_col=0).values

# Comprobaciones: mismas dimensiones que en R y 21 donantes (12 AD, 9 control).
print(adata)
assert adata.obs_names.is_unique and adata.var_names.is_unique
assert adata.obs["donor"].nunique() == 21
print(pd.crosstab(adata.obs["donor"], adata.obs["condition"]))

salida = os.path.join(base, "lau_sce_analisis.h5ad")
adata.write_h5ad(salida, compression="gzip")
print("Guardado:", salida)