# Lau et al. 2020: scrapper analysis of the preprocessed object

# MSc thesis "Transcriptional Noise in Alzheimer's Disease Progression" (CBMSO)
# MSc in Bioinformatics and Data Science (ISCIII)

# Second script. Starts from lau_sce_preprocesado.rds (output of
# scrapper_lau_preprocessing.R: QC, decontX, doublets, normalisation and
# genes >= 3 nuclei) and applies the remaining chapters of the scrapper book:
#   (1) Highly variable genes (HVGs)      scrapbook ch. 3
#   (2) PCA                               scrapbook ch. 4
#   (3) Batch correction with MNN         scrapbook ch. 8
#   (4) UMAP                              scrapbook ch. 5
#   (5) Clustering                        scrapbook ch. 6
#   (6) Markers                           scrapbook ch. 7
#   (7) AD versus control: pseudobulk     scrapbook ch. 8
#   (8) Multi-omics                       scrapbook ch. 9
#   (9) Saving

# Run from a terminal (not from the editor), because of memory use:
#   Rscript scrapper_lau_analysis.R > log_analysis.txt 2>&1

# General sources:
#   scrapper book: https://bioconductor.org/books/release/scrapbook/
#   OSCA (Amezquita et al. 2020, Nat Methods, doi:10.1038/s41592-019-0654-x):
#     https://bioconductor.org/books/release/OSCA/
#   Data: Lau et al. 2020, PNAS, doi:10.1073/pnas.2008762117



# (0) SETUP

library(SingleCellExperiment)  # https://bioconductor.org/packages/SingleCellExperiment
library(scrapper)              # https://bioconductor.org/packages/scrapper
library(scater)                # figures only (McCarthy et al. 2017, doi:10.1093/bioinformatics/btw777)
library(ggplot2)
library(AnnotationDbi)         # https://bioconductor.org/packages/AnnotationDbi
library(org.Hs.eg.db)          # Ensembl -> symbol: https://bioconductor.org/packages/org.Hs.eg.db
library(mclust)                # ARI only: https://cran.r-project.org/package=mclust
library(edgeR)                 # pseudobulk: https://bioconductor.org/packages/edgeR

set.seed(1234)
nthreads <- 2

proyecto <- "/path/to/alzheimer_tn"   # change to the local project path
outdir   <- file.path(proyecto, "results/02_analisis/Lau_32989152")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

sce <- readRDS(file.path(proyecto, "results/objects/lau_sce_preprocesado.rds"))
dim(sce)                         # genes x nuclei (20 donors)
assayNames(sce)                  # "counts" "decontXcounts" "logcounts"

# (0.1) Gene symbols (genes are Ensembl IDs), to read the markers.
#       Source: ?mapIds (AnnotationDbi).
simbolo <- mapIds(org.Hs.eg.db, keys = rownames(sce), keytype = "ENSEMBL",
                  column = "SYMBOL", multiVals = "first")
# Genes without a symbol keep their Ensembl ID (instead of "NA").
simbolo <- ifelse(is.na(simbolo), rownames(sce), simbolo)
names(simbolo) <- rownames(sce)
rowData(sce)$simbolo <- simbolo

# (0.2) Design columns.
sce$celltype <- sce$celltype_compendium
anotada <- !is.na(sce$celltype)
table(sce$celltype, useNA = "ifany")



# (1) HIGHLY VARIABLE GENES (scrapbook ch. 3)
#     The expected variance is fitted as a function of the mean and genes
#     are chosen by the residual. block = donor so that differences between
#     donors do not count as signal. N = 3000, the same criterion as the
#     coverage of canonical markers agreed in the 23-09-2026 meeting.
#     Source: ?modelGeneVariances, ?chooseHighlyVariableGenes; OSCA "Feature selection".

var_genes <- modelGeneVariances(logcounts(sce), block = sce$donor,
                                use.min.width = TRUE, num.threads = nthreads)$statistics
hvgs <- chooseHighlyVariableGenes(var_genes$residuals, top = 3000)
length(hvgs)

png(file.path(outdir, "01_varianza_media.png"), 1400, 1200, res = 150)
plot(var_genes$means, var_genes$variances, pch = 16, cex = 0.4, col = "grey",
     xlab = "Mean (log)", ylab = "Variance")
points(var_genes$means[hvgs], var_genes$variances[hvgs], pch = 16, cex = 0.4)
dev.off()



# (2) PCA (scrapbook ch. 4)
#     25 components on the 3000 HVGs.
#     Source: ?runPca; OSCA "Dimensionality reduction".

pca <- runPca(logcounts(sce)[hvgs, ], number = 25, num.threads = nthreads)
var_pct <- pca$variance.explained / pca$total.variance * 100
round(head(var_pct, 6), 2)
reducedDim(sce, "PCA") <- t(pca$components)

png(file.path(outdir, "02_scree_plot.png"), 1400, 1200, res = 150)
plot(var_pct, log = "y", type = "b", xlab = "PC", ylab = "Variance explained (%)")
dev.off()

# (2.1) R^2 of the first 6 PCs against cell type and against donor, in the
#       annotated nuclei only. If cell type explains more than donor,
#       biology dominates over batch.
r2_pca <- data.frame(pc = 1:6, tipo = NA, donante = NA)
for (k in 1:6) {
  pc_k <- pca$components[k, anotada]
  r2_pca$tipo[k]    <- summary(lm(pc_k ~ sce$celltype[anotada]))$r.squared
  r2_pca$donante[k] <- summary(lm(pc_k ~ sce$donor[anotada]))$r.squared
}
round(r2_pca, 3)



# (3) BATCH CORRECTION WITH MNN (scrapbook ch. 8)
#     Aligns the donors using pairs of mutual nearest neighbour nuclei.
#     WARNING: each donor is either AD or control only, so correcting by
#     donor also removes part of the disease effect. The corrected space is
#     used ONLY for UMAP and clustering. Markers, AD versus control and
#     noise are computed on uncorrected logcounts.
#     Source: Haghverdi et al. 2018, Nat Biotechnol 36:421-427,
#     doi:10.1038/nbt.4091; ?correctMnn; OSCA "Integrating datasets".

mnn <- correctMnn(pca$components, block = sce$donor, num.threads = nthreads)
reducedDim(sce, "MNN") <- t(mnn$corrected)

# (3.1) Same R^2 as in (2.1), now on the corrected space.
r2_mnn <- data.frame(pc = 1:6, tipo = NA, donante = NA)
for (k in 1:6) {
  pc_k <- mnn$corrected[k, anotada]
  r2_mnn$tipo[k]    <- summary(lm(pc_k ~ sce$celltype[anotada]))$r.squared
  r2_mnn$donante[k] <- summary(lm(pc_k ~ sce$donor[anotada]))$r.squared
}
round(r2_mnn, 3)



# (4) UMAP (scrapbook ch. 5)
#     For visualisation only: distances between islands are not quantitative.
#     Source: ?runUmap; McInnes et al. 2018, arXiv:1802.03426.

reducedDim(sce, "UMAP") <- runUmap(mnn$corrected, num.threads = nthreads)

# (4.1) Points are drawn in column order: since nuclei are grouped by
#       donor, the last donor (and the controls) hid the rest. They are
#       drawn in random order; the object itself is not changed.
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

# (4.2) Diagnosis in two panels (AD | Ctrl): shows that both conditions
#       occupy every island, without one hiding the other.
ggsave(file.path(outdir, "05b_umap_por_diagnostico_paneles.png"),
       plotReducedDim(sce_plot, "UMAP", colour_by = "celltype",
                      other_fields = "condition", point_size = 0.3) +
         facet_wrap(~ condition),
       width = 14, height = 7)



# (5) CLUSTERING (scrapbook ch. 6)
#     Shared nearest neighbour (SNN) graph and community detection on the
#     corrected space, without looking at any label.
#     Source: ?buildSnnGraph, ?clusterGraph; OSCA "Clustering".

snn <- buildSnnGraph(mnn$corrected, num.threads = nthreads)
sce$cluster <- factor(clusterGraph(snn)$membership)
rm(snn)
gc()
table(sce$cluster)

ggsave(file.path(outdir, "06_umap_por_cluster.png"),
       plotReducedDim(sce[, orden], "UMAP", colour_by = "cluster", text_by = "cluster"),
       width = 8, height = 7)

# (5.1) Purity: fraction of the majority cell type in each cluster.
tab <- table(cluster = sce$cluster[anotada], tipo = sce$celltype[anotada])
pureza <- apply(prop.table(tab, 1), 1, max)
tipo_cluster <- colnames(tab)[apply(tab, 1, which.max)]
names(tipo_cluster) <- rownames(tab)
sum(pureza > 0.7)
round(sort(pureza)[1:5], 3)      # the 5 least pure clusters

# (5.2) ARI against the authors' annotation (Hubert and Arabie 1985,
#       doi:10.1007/BF01908075). It is low when a type splits into several clusters.
adjustedRandIndex(sce$cluster[anotada], sce$celltype[anotada])

# (5.3) No cluster should come almost entirely from one donor (> 70 %;
#       Mathys et al. 2019, Nature, doi:10.1038/s41586-019-1195-2).
max_donante <- apply(prop.table(table(sce$cluster, sce$donor), 1), 1, max)
max_donante[max_donante > 0.7]

# (5.4) Mitochondrial proportion per cluster: the 5 clusters with the most
#       mitochondrial RNA, candidates for low-quality nuclei.
round(sort(tapply(sce$subsets_mito_percent, sce$cluster, median), decreasing = TRUE)[1:5], 2)



# (6) MARKERS (scrapbook ch. 7)
#     On uncorrected logcounts, comparing within each donor (block).
#     Ranked by AUC.
#     Source: ?scoreMarkers, ?reportGroupMarkerStatistics; OSCA "Marker detection".

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



# (7) AD VERSUS CONTROL: PSEUDOBULK (scrapbook ch. 8)
#     Raw counts are summed per donor within each cell type: the replicate
#     is the donor, not the nucleus. Raw counts are used because edgeR
#     models integer counts. Control is the reference: logFC > 0 = higher
#     in AD. AD6 was already excluded in
#     preprocessing (step 2.4).
#     Sources: Crowell et al. 2020, Nat Commun 11:6077, doi:10.1038/s41467-020-19894-4;
#     Squair et al. 2021, Nat Commun 12:5692, doi:10.1038/s41467-021-25960-2;
#     edgeR (Chen, Lun and Smyth 2016, doi:10.12688/f1000research.8987.2).

sce$condition <- relevel(factor(sce$condition), ref = "Ctrl")
tipos <- sort(unique(sce$celltype[anotada]))
resumen_DE <- data.frame()

for (tipo in tipos) {
  print(tipo)
  sel <- which(sce$celltype == tipo)
  donantes <- unique(sce$donor[sel])

  # (7.1) Pseudobulk matrix: one column per donor.
  sumas <- matrix(0, nrow(sce), length(donantes),
                  dimnames = list(rownames(sce), donantes))
  for (d in donantes) {
    sumas[, d] <- rowSums(counts(sce)[, sel[sce$donor[sel] == d], drop = FALSE])
  }
  info <- data.frame(donante = donantes,
                     diag = sce$condition[match(donantes, sce$donor)])

  # (7.2) edgeR model (quasi-likelihood).
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

# (7.3) Differential abundance: number of nuclei of each type in each
#       donor, with the same edgeR scheme (OSCA, "Differential abundance").
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



# (8) MULTI-OMICS (scrapbook ch. 9)
#     The chapter covers data with protein (CITE-seq) stored as an
#     alternative experiment. Lau is snRNA-seq only: this is checked and
#     does not apply.

altExpNames(sce)                 # character(0)



# (9) SAVING

saveRDS(sce, file.path(proyecto, "results/objects/lau_sce_analisis.rds"))
sessionInfo()