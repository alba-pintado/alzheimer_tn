# Lau et al. 2020: analisis con scrapper sobre el objeto preprocesado

# TFM "Transcriptional Noise in Alzheimer's Disease Progression" (CBMSO)
# Master en Bioinformatica y Ciencia de Datos (ISCIII)

# Segundo script. Parte de lau_sce_preprocesado.rds (salida de
# scrapper_lau_preprocesado.R: QC, decontX, dobletes, normalizacion y genes >= 3
# nucleos) y aplica el resto de capitulos del libro de scrapper:

#   (1) Genes variables (HVGs)            scrapbook cap. 3
#   (2) PCA                               scrapbook cap. 4
#   (3) Correccion de lote con MNN        scrapbook cap. 8
#   (4) UMAP                              scrapbook cap. 5
#   (5) Clustering                        scrapbook cap. 6
#   (6) Marcadores                        scrapbook cap. 7
#   (7) AD frente a control: pseudobulk   scrapbook cap. 8
#   (8) Multiomica                        scrapbook cap. 9
#   (9) Guardado

# Ejecutar desde una terminal (no desde el editor), por la memoria:
#   Rscript scrapper_lau_analisis.R > log_analisis.txt 2>&1

# Fuentes generales:
#   Libro de scrapper: https://bioconductor.org/books/release/scrapbook/
#   OSCA (Amezquita et al. 2020, Nat Methods, doi:10.1038/s41592-019-0654-x):
#     https://bioconductor.org/books/release/OSCA/
#   Datos: Lau et al. 2020, PNAS, doi:10.1073/pnas.2008762117



# (0) PREPARACION

library(SingleCellExperiment)  # https://bioconductor.org/packages/SingleCellExperiment
library(scrapper)              # https://bioconductor.org/packages/scrapper
library(scater)                # solo para las figuras (McCarthy et al. 2017, doi:10.1093/bioinformatics/btw777)
library(ggplot2)
library(AnnotationDbi)         # https://bioconductor.org/packages/AnnotationDbi
library(org.Hs.eg.db)          # Ensembl -> simbolo: https://bioconductor.org/packages/org.Hs.eg.db
library(mclust)                # solo para el ARI: https://cran.r-project.org/package=mclust
library(edgeR)                 # pseudobulk: https://bioconductor.org/packages/edgeR

set.seed(1234)
nthreads <- 2

proyecto <- "/home/vant/Documentos/proyectos/alzheimer_tn"   # cambiar por la ruta local del proyecto
outdir   <- file.path(proyecto, "results/02_analisis/Lau_32989152")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

sce <- readRDS(file.path(proyecto, "results/objects/lau_sce_preprocesado.rds"))
dim(sce)                         # genes x nucleos (20 donantes)
assayNames(sce)                  # "counts" "decontXcounts" "logcounts"

# (0.1) Simbolos de los genes (vienen en Ensembl), para leer los marcadores.
#       Fuente: ?mapIds (AnnotationDbi).
simbolo <- mapIds(org.Hs.eg.db, keys = rownames(sce), keytype = "ENSEMBL",
                  column = "SYMBOL", multiVals = "first")
# Si un gen no tiene simbolo, se deja su ID de Ensembl (antes salia "NA").
simbolo <- ifelse(is.na(simbolo), rownames(sce), simbolo)
names(simbolo) <- rownames(sce)
rowData(sce)$simbolo <- simbolo

# (0.2) Columnas de diseno.
sce$celltype <- sce$celltype_compendium
anotada <- !is.na(sce$celltype)
table(sce$celltype, useNA = "ifany")



# (1) GENES VARIABLES (scrapbook cap. 3)
#     Se ajusta la varianza esperada segun la media y se elige el residuo.
#     block = donante para que la diferencia entre donantes no cuente como
#     senal. N = 3000, el mismo criterio de la cobertura de marcadores
#     canonicos que se decidio en la tutoria del 23-09-2026.
#     Fuente: ?modelGeneVariances, ?chooseHighlyVariableGenes; OSCA "Feature selection".

var_genes <- modelGeneVariances(logcounts(sce), block = sce$donor,
                                use.min.width = TRUE, num.threads = nthreads)$statistics
hvgs <- chooseHighlyVariableGenes(var_genes$residuals, top = 3000)
length(hvgs)

png(file.path(outdir, "01_varianza_media.png"), 1400, 1200, res = 150)
plot(var_genes$means, var_genes$variances, pch = 16, cex = 0.4, col = "grey",
     xlab = "Media (log)", ylab = "Varianza")
points(var_genes$means[hvgs], var_genes$variances[hvgs], pch = 16, cex = 0.4)
dev.off()



# (2) PCA (scrapbook cap. 4)
#     25 componentes sobre los 3000 HVGs.
#     Fuente: runPca; OSCA "Dimensionality reduction".

pca <- runPca(logcounts(sce)[hvgs, ], number = 25, num.threads = nthreads)
var_pct <- pca$variance.explained / pca$total.variance * 100
round(head(var_pct, 6), 2)
reducedDim(sce, "PCA") <- t(pca$components)

png(file.path(outdir, "02_scree_plot.png"), 1400, 1200, res = 150)
plot(var_pct, log = "y", type = "b", xlab = "PC", ylab = "Varianza explicada (%)")
dev.off()

# (2.1) R^2 de las 6 primeras PCs frente a tipo celular y frente a donante,
#       solo en los nucleos anotados. Si el tipo celular explica mas que el
#       donante, la biologia manda sobre el lote.
r2_pca <- data.frame(pc = 1:6, tipo = NA, donante = NA)
for (k in 1:6) {
  pc_k <- pca$components[k, anotada]
  r2_pca$tipo[k]    <- summary(lm(pc_k ~ sce$celltype[anotada]))$r.squared
  r2_pca$donante[k] <- summary(lm(pc_k ~ sce$donor[anotada]))$r.squared
}
round(r2_pca, 3)



# (3) CORRECCION DE LOTE CON MNN (scrapbook cap. 8)
#     Alinea los donantes con pares de nucleos vecinos mutuos.
#     AVISO: cada donante es solo AD o solo control, asi que corregir por
#     donante tambien quita parte del efecto de la enfermedad. El espacio
#     corregido se usa SOLO para el UMAP y el clustering. Marcadores,
#     AD frente a control y ruido se calculan sobre logcounts sin corregir.
#     Fuente: Haghverdi et al. 2018, Nat Biotechnol 36:421-427,
#     doi:10.1038/nbt.4091; ?correctMnn; OSCA "Integrating datasets".

mnn <- correctMnn(pca$components, block = sce$donor, num.threads = nthreads)
reducedDim(sce, "MNN") <- t(mnn$corrected)

# (3.1) Mismo R^2 que en (2.1), ahora sobre el espacio corregido.
r2_mnn <- data.frame(pc = 1:6, tipo = NA, donante = NA)
for (k in 1:6) {
  pc_k <- mnn$corrected[k, anotada]
  r2_mnn$tipo[k]    <- summary(lm(pc_k ~ sce$celltype[anotada]))$r.squared
  r2_mnn$donante[k] <- summary(lm(pc_k ~ sce$donor[anotada]))$r.squared
}
round(r2_mnn, 3)



# (4) UMAP (scrapbook cap. 5)
#     Solo para ver: las distancias entre islas no son cuantitativas.
#     Fuente: ?runUmap; McInnes et al. 2018, arXiv:1802.03426.

reducedDim(sce, "UMAP") <- runUmap(mnn$corrected, num.threads = nthreads)

# (4.1) Los puntos se dibujan en el orden de las columnas: como los nucleos
#       estan agrupados por donante, el ultimo donante (y los controles)
#       tapaban al resto. Se dibujan en orden aleatorio; el objeto no se toca.
orden <- sample(ncol(sce))
sce_plot <- sce[, orden]

ggsave(file.path(outdir, "03_umap_por_celltype.png"),
       plotReducedDim(sce_plot, "UMAP", colour_by = "celltype", point_size = 0.3),
       width = 8, height = 7)
ggsave(file.path(outdir, "04_umap_por_donante.png"),
       plotReducedDim(sce_plot, "UMAP", colour_by = "donor", point_size = 0.3),
       width = 8, height = 7)
ggsave(file.path(outdir, "05_umap_por_diagnostico.png"),
       plotReducedDim(sce_plot, "UMAP", colour_by = "condition", point_size = 0.3),
       width = 8, height = 7)

# (4.2) Diagnostico en dos paneles (AD | Ctrl): asi se ve que ambas
#       condiciones ocupan todas las islas, sin que una tape a la otra.
ggsave(file.path(outdir, "05b_umap_por_diagnostico_paneles.png"),
       plotReducedDim(sce_plot, "UMAP", colour_by = "celltype",
                      other_fields = "condition", point_size = 0.3) +
         facet_wrap(~ condition),
       width = 14, height = 7)



# (5) CLUSTERING (scrapbook cap. 6)
#     Grafo de vecinos compartidos (SNN) y deteccion de comunidades sobre el
#     espacio corregido, sin mirar ninguna etiqueta.
#     Fuente: ?buildSnnGraph, ?clusterGraph; OSCA "Clustering".

snn <- buildSnnGraph(mnn$corrected, num.threads = nthreads)
sce$cluster <- factor(clusterGraph(snn)$membership)
rm(snn)
gc()
table(sce$cluster)

ggsave(file.path(outdir, "06_umap_por_cluster.png"),
       plotReducedDim(sce[, orden], "UMAP", colour_by = "cluster", text_by = "cluster"),
       width = 8, height = 7)

# (5.1) Pureza: fraccion del tipo celular mayoritario en cada cluster.
tab <- table(cluster = sce$cluster[anotada], tipo = sce$celltype[anotada])
pureza <- apply(prop.table(tab, 1), 1, max)
tipo_cluster <- colnames(tab)[apply(tab, 1, which.max)]
names(tipo_cluster) <- rownames(tab)
sum(pureza > 0.7)
round(sort(pureza)[1:5], 3)      # los 5 clusters menos puros

# (5.2) ARI frente a la anotacion de los autores (Hubert y Arabie 1985,
#       doi:10.1007/BF01908075). Es bajo si un tipo se parte en varios clusters.
adjustedRandIndex(sce$cluster[anotada], sce$celltype[anotada])

# (5.3) Ningun cluster deberia venir casi entero de un donante (> 70 %;
#       Mathys et al. 2019, Nature, doi:10.1038/s41586-019-1195-2).
max_donante <- apply(prop.table(table(sce$cluster, sce$donor), 1), 1, max)
max_donante[max_donante > 0.7]

# (5.4) Proporcion mitocondrial por cluster: los 5 clusters con mas ARN
#       mitocondrial, candidatos a nucleos de baja calidad.
round(sort(tapply(sce$subsets_mito_percent, sce$cluster, median), decreasing = TRUE)[1:5], 2)



# (6) MARCADORES (scrapbook cap. 7)
#     Sobre logcounts sin corregir, comparando dentro de cada donante
#     (block). Se ordena por AUC.
#     Fuente: ?scoreMarkers, ?reportGroupMarkerStatistics; OSCA "Marker detection".

markers <- scoreMarkers(logcounts(sce), groups = sce$cluster,
                        block = sce$donor, num.threads = nthreads)

resumen_marcadores <- data.frame(cluster = levels(sce$cluster), tipo = NA,
                                 pureza = NA, top10 = NA)
for (i in seq_along(levels(sce$cluster))) {
  k <- levels(sce$cluster)[i]
  res_k <- reportGroupMarkerStatistics(markers, k)
  top_k <- head(rownames(res_k)[order(res_k$auc.mean, decreasing = TRUE)], 10)
  resumen_marcadores$tipo[i]   <- tipo_cluster[k]
  resumen_marcadores$pureza[i] <- round(pureza[k], 3)
  resumen_marcadores$top10[i]  <- paste(simbolo[top_k], collapse = ", ")
}
resumen_marcadores
write.csv(resumen_marcadores, file.path(outdir, "07_marcadores_top10.csv"), row.names = FALSE)



# (7) AD FRENTE A CONTROL: PSEUDOBULK (scrapbook cap. 8)
#     Se suman las cuentas crudas (counts) de cada donante dentro de cada
#     tipo celular: la replica es el donante, no el nucleo. Se usan las
#     crudas porque edgeR modela cuentas enteras. Control es la referencia:
#     logFC > 0 = mas alto en AD. AD6 ya se excluyo en el
#     preprocesado (paso 2.4).
#     Fuentes: Crowell et al. 2020, Nat Commun 11:6077, doi:10.1038/s41467-020-19894-4;
#     Squair et al. 2021, Nat Commun 12:5692, doi:10.1038/s41467-021-25960-2;
#     edgeR (Chen, Lun y Smyth 2016, doi:10.12688/f1000research.8987.2).

sce$condition <- relevel(factor(sce$condition), ref = "Ctrl")
tipos <- sort(unique(sce$celltype[anotada]))
resumen_DE <- data.frame()

for (tipo in tipos) {
  print(tipo)
  sel <- which(sce$celltype == tipo)
  donantes <- unique(sce$donor[sel])

  # (7.1) Matriz pseudobulk: una columna por donante.
  sumas <- matrix(0, nrow(sce), length(donantes),
                  dimnames = list(rownames(sce), donantes))
  for (d in donantes) {
    sumas[, d] <- rowSums(counts(sce)[, sel[sce$donor[sel] == d], drop = FALSE])
  }
  info <- data.frame(donante = donantes,
                     diag = sce$condition[match(donantes, sce$donor)])

  # (7.2) Modelo de edgeR (quasi-likelihood).
  y <- DGEList(sumas, samples = info)
  y <- y[filterByExpr(y, group = y$samples$diag), , keep.lib.sizes = FALSE]
  y <- normLibSizes(y)
  diseno <- model.matrix(~ diag, y$samples)
  y <- estimateDisp(y, diseno)
  ajuste <- glmQLFit(y, diseno, robust = TRUE)
  res <- topTags(glmQLFTest(ajuste, coef = 2), n = Inf)$table
  res$simbolo <- simbolo[rownames(res)]
  write.csv(res, file.path(outdir, paste0("08_DE_", tipo, ".csv")))

  resumen_DE <- rbind(resumen_DE, data.frame(
    tipo = tipo, donantes = ncol(y), genes = nrow(y),
    significativos = sum(res$FDR < 0.05)))
}
resumen_DE
write.csv(resumen_DE, file.path(outdir, "09_resumen_DE.csv"), row.names = FALSE)

# (7.3) Abundancia diferencial: numero de nucleos de cada tipo en cada
#       donante, con el mismo esquema de edgeR (OSCA, "Differential abundance").
abund <- unclass(table(sce$celltype[anotada], sce$donor[anotada]))
info_ab <- data.frame(donante = colnames(abund),
                      diag = sce$condition[match(colnames(abund), sce$donor)])
y_ab <- DGEList(abund, samples = info_ab)
d_ab <- model.matrix(~ diag, y_ab$samples)
y_ab <- estimateDisp(y_ab, d_ab, trend = "none")
ajuste_ab <- glmQLFit(y_ab, d_ab, robust = TRUE, abundance.trend = FALSE)
abundancia <- topTags(glmQLFTest(ajuste_ab, coef = 2))$table
abundancia
write.csv(abundancia, file.path(outdir, "10_abundancia_diferencial.csv"))



# (8) MULTIOMICA (scrapbook cap. 9)
#     El capitulo trata datos con proteina (CITE-seq) guardada como
#     experimento alternativo. Lau es solo snRNA-seq: se comprueba y no aplica.

altExpNames(sce)                 # character(0)



# (9) GUARDADO

saveRDS(sce, file.path(proyecto, "results/objects/lau_sce_analisis.rds"))
sessionInfo()