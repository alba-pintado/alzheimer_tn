# Lau et al. 2020: UMAP y clustering con Harmony en lugar de MNN
#
# TFM "Transcriptional Noise in Alzheimer's Disease Progression" (CBMSO)
#
# Comparación con el análisis principal (02_lau_analisis.R), que corrige
# por donante con MNN. El diagnóstico del UMAP mostró que MNN fragmenta los
# subtipos neuronales (42 islas con MNN, 17 sin corregir). Aquí se prueba
# Harmony, una corrección más suave, sobre el MISMO PCA. No modifica el objeto
# del análisis principal: marcadores, pseudobulk y abundancia no cambian,
# porque se calculan sobre datos sin corregir.
#
#   (1) Harmony sobre el PCA, por donante
#   (2) R² de las componentes frente a tipo celular y donante
#   (3) UMAP y clustering sobre el espacio Harmony
#   (4) Pureza, ARI, donantes por cluster e islas del UMAP (como el diagnóstico)
#   (5) Comparación con los clusters de MNN
#
# Fuente: Korsunsky et al. 2019, Nat Methods 16:1289-1296,
# doi:10.1038/s41592-019-0619-0; https://cran.r-project.org/package=harmony

library(SingleCellExperiment)
library(scrapper)
library(scater)
library(ggplot2)
library(harmony)                   # install.packages("harmony")
library(mclust)                    # ARI
library(dbscan)                    # islas del UMAP

set.seed(1234)
proyecto <- "~/Documentos/proyectos/alzheimer_tn"   # carpeta del proyecto
outdir <- file.path(proyecto, "results/02_analisis/Lau_32989152/harmony")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

sce <- readRDS(file.path(proyecto, "results/objects/lau_2020/lau_sce_analisis.rds"))
dim(sce)
anotada <- !is.na(sce$celltype)



# (1) HARMONY SOBRE EL PCA
#     Mismo PCA de 25 componentes del análisis principal (núcleos x PCs).

pca <- reducedDim(sce, "PCA")
harm <- RunHarmony(pca, data.frame(donor = sce$donor), vars_use = "donor")
reducedDim(sce, "HARMONY") <- harm



# (2) R² FRENTE A TIPO CELULAR Y DONANTE (como en el análisis principal)

r2 <- data.frame(pc = 1:6, tipo = NA, donante = NA)
for (k in 1:6) {
  pc_k <- harm[anotada, k]
  r2$tipo[k]    <- summary(lm(pc_k ~ sce$celltype[anotada]))$r.squared
  r2$donante[k] <- summary(lm(pc_k ~ sce$donor[anotada]))$r.squared
}
round(r2, 3)
write.csv(round(r2, 3), file.path(outdir, "01_r2_harmony.csv"), row.names = FALSE)



# (3) UMAP Y CLUSTERING SOBRE HARMONY
#     scrapper espera las componentes en filas (PCs x núcleos), por eso t().

reducedDim(sce, "UMAP_HARMONY") <- runUmap(t(harm), num.threads = 2)

snn <- buildSnnGraph(t(harm), num.threads = 2)
sce$cluster_harmony <- factor(clusterGraph(snn)$membership)
rm(snn); gc()
table(sce$cluster_harmony)

orden <- sample(ncol(sce))
sce_plot <- sce[, orden]
ggsave(file.path(outdir, "02_umap_harmony_celltype.png"),
       plotReducedDim(sce_plot, "UMAP_HARMONY", colour_by = "celltype", point_size = 0.3),
       width = 8, height = 7)
ggsave(file.path(outdir, "03_umap_harmony_donante.png"),
       plotReducedDim(sce_plot, "UMAP_HARMONY", colour_by = "donor", point_size = 0.3),
       width = 8, height = 7)
ggsave(file.path(outdir, "04_umap_harmony_paneles.png"),
       plotReducedDim(sce_plot, "UMAP_HARMONY", colour_by = "celltype",
                      other_fields = "condition", point_size = 0.3) +
         facet_wrap(~ condition),
       width = 14, height = 7)
ggsave(file.path(outdir, "05_umap_harmony_cluster.png"),
       plotReducedDim(sce_plot, "UMAP_HARMONY", colour_by = "cluster_harmony",
                      text_by = "cluster_harmony", point_size = 0.3),
       width = 8, height = 7)



# (4) VALIDACIÓN: PUREZA, ARI, DONANTES E ISLAS

tab <- table(sce$cluster_harmony[anotada], sce$celltype[anotada])
pureza <- apply(prop.table(tab, 1), 1, max)
print(paste("Clusters:", nrow(tab), "| con pureza >= 0,99:", sum(pureza >= 0.99)))
round(sort(pureza)[1:5], 3)

print(paste("ARI frente a la anotación:",
            round(adjustedRandIndex(sce$cluster_harmony[anotada], sce$celltype[anotada]), 3)))

max_donante <- apply(prop.table(table(sce$cluster_harmony, sce$donor), 1), 1, max)
print("Clusters con más del 70 % de un donante:")
max_donante[max_donante > 0.7]

islas <- dbscan(reducedDim(sce, "UMAP_HARMONY"), eps = 0.3, minPts = 20)
print(paste("Islas del UMAP con Harmony:", length(unique(islas$cluster[islas$cluster > 0])),
            "(con MNN: 42; sin corregir: 17)"))



# (5) COMPARACIÓN CON LOS CLUSTERS DE MNN

print(paste("ARI entre clusters de Harmony y de MNN:",
            round(adjustedRandIndex(sce$cluster_harmony, sce$cluster), 3)))
write.csv(as.data.frame.matrix(table(harmony = sce$cluster_harmony, mnn = sce$cluster)),
          file.path(outdir, "06_harmony_vs_mnn.csv"))

sessionInfo()