# Lau et al. 2020: exportar el objeto de R a HDF5 sin zellkonverter

# zellkonverter copia el objeto entero a Python antes de escribir el .h5ad, así
# que en memoria hay dos copias a la vez y el sistema mata el proceso (killed).
# Este script escribe las matrices dispersas por partes (x, i, p) con rhdf5 y la
# metadata en CSV; 00b_construir_h5ad.py monta después el .h5ad en Python.

# Ejecutar desde la raíz del proyecto:
#   nohup Rscript scanpy/lau_2020/00_exportar_h5ad.R > log_exportar.txt 2>&1 &

library(SingleCellExperiment)
library(Matrix)
library(rhdf5)                     # BiocManager::install("rhdf5")


proyecto <- "~/Documentos/proyectos/alzheimer_tn"
dir_obj  <- file.path(proyecto, "results/objects/lau_2020")
salida   <- file.path(dir_obj, "lau_export.h5")

sce <- readRDS(file.path(dir_obj, "lau_sce_analisis.rds"))
dim(sce)                           # 29456 genes x 140727 núcleos
assay(sce, "decontXcounts") <- NULL   # Python no la usa
gc()

if (file.exists(salida)) file.remove(salida)
h5createFile(salida)
h5write(dim(sce), salida, "shape")     # genes, núcleos

escribir <- function(v, nombre, modo) {
  h5createDataset(salida, nombre, dims = length(v), storage.mode = modo,
                  chunk = min(1e6, length(v)), level = 4)
  h5write(v, salida, nombre)
}

# Cada matriz se guarda en formato columna disperso (dgCMatrix: genes x núcleos)
# y se borra del objeto antes de pasar a la siguiente, para no duplicar memoria.
for (a in c("counts", "logcounts")) {
  print(a)
  m <- as(assay(sce, a), "CsparseMatrix")
  h5createGroup(salida, a)
  escribir(m@x, paste0(a, "/x"), "double")
  escribir(m@i, paste0(a, "/i"), "integer")
  escribir(m@p, paste0(a, "/p"), "integer")
  rm(m); assay(sce, a) <- NULL; gc()
}
h5closeAll()

# Metadata de los núcleos y de los genes, y coordenadas de R (PCA y UMAP con MNN).
columnas <- c("donor", "condition", "celltype", "cluster", "subsets_mito_percent",
              "sum", "detected", "decontX_contamination")
obs <- as.data.frame(colData(sce)[, intersect(columnas, colnames(colData(sce)))])
write.csv(obs, file.path(dir_obj, "lau_export_obs.csv"))
write.csv(data.frame(simbolo = rowData(sce)$simbolo, row.names = rownames(sce)),
          file.path(dir_obj, "lau_export_var.csv"))
write.csv(reducedDim(sce, "PCA"),  file.path(dir_obj, "lau_export_pca_R.csv"))
write.csv(reducedDim(sce, "UMAP"), file.path(dir_obj, "lau_export_umap_R.csv"))

print(paste("Guardado:", salida))
sessionInfo()