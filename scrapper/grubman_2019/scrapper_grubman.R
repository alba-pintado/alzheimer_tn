# Astrocytes from Grubman 2019 with scrapper — full pipeline (chapters 1-9)
# MSc thesis "Transcriptional noise in Alzheimer's disease progression" — CBMSO
# Alba Pintado

# The scrapper book, applied to the astrocytes of Grubman:
# https://bioconductor.org/books/release/scrapbook/
# Chapters 1-7 are the first pass, WITHOUT batch correction; chapter 8 corrects
# the batch effect (MNN) and redoes visualisation, clustering and markers.
# Chapter 9 (multi-omics) does not apply: this dataset is RNA only.

# Each block carries a short note on what it does and which section of the book
# it comes from. Figures are written to results/01_exploration/Grubman_31768052.


# (0) Setup

library(SingleCellExperiment)
library(scater)
library(scrapper)
library(ggplot2)
library(gridExtra)
library(Matrix)

set.seed(1234)
nthreads <- 2

project  <- "/home/vant/Documentos/proyectos/alzheimer_tn"
rds.path <- file.path(project, "data/organized_data/human/brain/snRNA",
                      "Grubman_31768052", "Astro_data_arranged_idUpdated_final.rds")
outdir   <- file.path(project, "results/01_exploration/Grubman_31768052")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
fig <- function(x) file.path(outdir, x)

sce <- readRDS(rds.path)

# (0.1) Design variables. In Grubman every library is 100% AD or 100% control,
#       so batch and diagnosis are the same variable (true n = 3 vs 3).
COL_DONOR <- "individuals.ch1"
COL_DIAG  <- "oupSample.batchCond"
donor     <- factor(colData(sce)[[COL_DONOR]])
diagnosis <- factor(colData(sce)[[COL_DIAG]])   # not 'diag': that masks base::diag
counts    <- assay(sce, "counts")

table(donor, diagnosis)   # each library falls entirely in one group
range(table(donor))       # imbalance: 96 .. 853 nuclei per library


# (1) Quality control
# Source: scrapbook ch. 1. Median +/- 3 MAD thresholds per block (1.3.1, 1.5).

is.mito <- grep("^MT-", rownames(sce))
qc <- computeRnaQcMetrics(counts, subsets = list(MT = is.mito), num.threads = nthreads)

png(fig("01_profundidad_por_diagnostico.png"), 1400, 1200, res = 150)
boxplot(qc$detected ~ diagnosis, ylab = "Genes per nucleus", main = "Depth by diagnosis")
dev.off()
png(fig("02_profundidad_por_donante.png"), 1400, 1200, res = 150)
boxplot(qc$detected ~ donor, las = 2, ylab = "Genes per nucleus", main = "Depth by donor")
dev.off()

tapply(qc$detected, diagnosis, median)   # AD 857.5 vs ct 1004

# (1.1) Per-block thresholds and filtering. Check that no threshold degenerates (1.5).
thresholds <- suggestRnaQcThresholds(qc, block = donor)
qc.keep    <- filterRnaQcMetrics(thresholds, qc, block = donor)
sum(qc.keep); sum(!qc.keep)              # 1991 kept, 180 discarded

# (1.2) Which metric discards each nucleus (in nuclei almost everything goes on MT).
b <- as.character(donor)
c(UMIs  = sum(qc$sum      < thresholds$sum[b]),
  genes = sum(qc$detected < thresholds$detected[b]),
  MT    = sum(qc$subsets$MT > thresholds$subsets$MT[b]))

# (1.3) The filter does not penalise the AD samples (6.8% vs 8.7% discarded).
round(prop.table(table(diagnosis, discarded = !qc.keep), 1), 3)

# (1.4) Visual diagnostics per block (1.4).
sce$sum <- qc$sum; sce$detected <- qc$detected; sce$MT <- qc$subsets$MT
sce$keep <- qc.keep; sce$block <- donor
png(fig("02b_qc_metricas_por_bloque.png"), 2000, 1300, res = 150)
grid.arrange(
  plotColData(sce, y = "sum",      x = "block", colour_by = "keep") + scale_y_log10() + ggtitle("UMIs"),
  plotColData(sce, y = "detected", x = "block", colour_by = "keep") + scale_y_log10() + ggtitle("Genes"),
  plotColData(sce, y = "MT",       x = "block", colour_by = "keep") + ggtitle("MT"),
  ncol = 3)
dev.off()

filtered    <- counts[, qc.keep, drop = FALSE]
donor.f     <- droplevels(donor[qc.keep])
diagnosis.f <- droplevels(diagnosis[qc.keep])


# (2) Normalisation
# Source: scrapbook ch. 2. Size factors per library + log (2.3.1, 2.4).

size.factors <- centerSizeFactors(qc$sum[qc.keep], block = donor.f)
normalized   <- scrapper::normalizeCounts(filtered, size.factors)   # scrapper:: avoids the warning

png(fig("03_sizefactors_por_diagnostico.png"), 1400, 1200, res = 150)
boxplot(size.factors ~ diagnosis.f, main = "Size factors by diagnosis")
dev.off()


# (3) Feature selection (HVGs)
# Source: scrapbook ch. 3. Variance~mean residual; block only for nuisance factors (3.3).

gene.var <- modelGeneVariances(normalized, block = donor.f, num.threads = nthreads)
st.orig  <- gene.var$statistics
st       <- modelGeneVariances(normalized, block = donor.f,
                               use.min.width = TRUE, num.threads = nthreads)$statistics  # 3.4

top.hvgs <- chooseHighlyVariableGenes(st$residuals, top = 4000)
is.hvg   <- logical(nrow(normalized)); is.hvg[top.hvgs] <- TRUE

# Fitted curve drawn with lines() over the ordered means (curve() threw an error).
png(fig("04_varianza_media.png"), 1400, 1200, res = 150)
plot(st$means, st$variances, col = ifelse(is.hvg, "black", "grey"), pch = 16, cex = 0.5,
     xlab = "Mean (log)", ylab = "Variance", main = "Mean-variance relationship")
ord.m <- order(st$means); lines(st$means[ord.m], st$fitted[ord.m], col = "dodgerblue", lwd = 2)
legend("topright", col = c("black", "grey"), pch = 16, legend = c("HVG", "not HVG"))
dev.off()

# Biological sanity check: astrocyte markers should sit near the top of the ranking.
ord <- order(st$residuals, decreasing = TRUE)
astro.markers <- c("GFAP","AQP4","SLC1A2","SLC1A3","CLU","MT2A","CST3","APOE")
data.frame(gene = astro.markers, rank = match(astro.markers, rownames(normalized)[ord]))


# (4) PCA
# Source: scrapbook ch. 4. 'number' is how many PCs to keep (4.3).

pca <- runPca(normalized[top.hvgs, ], number = 25, num.threads = nthreads)
percent.var <- pca$variance.explained / pca$total.variance * 100

png(fig("06_scree_plot.png"), 1400, 1200, res = 150)
plot(percent.var, log = "y", type = "b", xlab = "PC", ylab = "Variance explained (%)")
dev.off()
png(fig("07_pc1_por_donante.png"), 1400, 1200, res = 150)
boxplot(pca$components[1, ] ~ donor.f, las = 2, ylab = "PC1", main = "PC1 by donor")
dev.off()
png(fig("08_pc1_por_diagnostico.png"), 1400, 1200, res = 150)
boxplot(pca$components[1, ] ~ diagnosis.f, ylab = "PC1", main = "PC1 by diagnosis")
dev.off()

# The size of the problem: per PC, how much variance library explains vs diagnosis (4.4).
r2.donor <- apply(pca$components, 1, function(pc) summary(lm(pc ~ donor.f))$r.squared)
r2.diag  <- apply(pca$components, 1, function(pc) summary(lm(pc ~ diagnosis.f))$r.squared)
round(head(cbind(library = r2.donor, diagnosis = r2.diag), 6), 3)


# (5) Visualisation (uncorrected)
# Source: scrapbook ch. 5. For looking only; distances are not quantitative (5.3).

set.seed(1234)
tsne.uncorr <- runTsne(pca$components, num.threads = nthreads)
umap.uncorr <- runUmap(pca$components, num.threads = nthreads)

sce.f <- sce[, qc.keep]
sizeFactors(sce.f) <- size.factors
assay(sce.f, "logcounts") <- normalized
reducedDim(sce.f, "PCA")  <- t(pca$components)
reducedDim(sce.f, "TSNE") <- tsne.uncorr
reducedDim(sce.f, "UMAP") <- umap.uncorr

ggsave(fig("09_umap_por_diagnostico.png"), plotReducedDim(sce.f, "UMAP", colour_by = COL_DIAG),  width = 7, height = 6, dpi = 150)
ggsave(fig("10_umap_por_donante.png"),     plotReducedDim(sce.f, "UMAP", colour_by = COL_DONOR), width = 7, height = 6, dpi = 150)
ggsave(fig("10b_tsne_por_donante.png"),    plotReducedDim(sce.f, "TSNE", colour_by = COL_DONOR), width = 7, height = 6, dpi = 150)


# (6) Clustering (uncorrected)
# Source: scrapbook ch. 6. SNN graph + community detection (6.2). The target is
# substates within astrocytes, not cell types.

snn.uncorr   <- buildSnnGraph(pca$components, num.threads = nthreads)
clust.uncorr <- clusterGraph(snn.uncorr)
sce.f$cluster.uncorr <- factor(clust.uncorr$membership)
ggsave(fig("11_umap_por_cluster.png"), plotReducedDim(sce.f, "UMAP", colour_by = "cluster.uncorr"), width = 7, height = 6, dpi = 150)

# Mathys' control (6.5): a cluster dominated by one library is batch, not biology.
tab.uncorr <- table(cluster = clust.uncorr$membership, donor = donor.f)
sum(apply(prop.table(tab.uncorr, 1), 1, max) > 0.70)


# (7) Markers (uncorrected)
# Source: scrapbook ch. 7. Per-gene effect sizes; block = donor.f (7.5).

markers.uncorr <- scoreMarkers(normalized, groups = clust.uncorr$membership,
                               block = donor.f, num.threads = nthreads)
df.uncorr  <- reportGroupMarkerStatistics(markers.uncorr, "1")
top.uncorr <- head(rownames(df.uncorr)[order(df.uncorr$cohens.d.mean, decreasing = TRUE)], 15)

ggsave(fig("15_marcadores_violines.png"),
       plotExpression(sce.f, features = top.uncorr[1:6], x = "cluster.uncorr", colour_by = "cluster.uncorr"),
       width = 8, height = 6, dpi = 150)
png(fig("16_marcadores_heatmap.png"), 1600, 1400, res = 150)
plotGroupedHeatmap(sce.f, features = top.uncorr, group = "cluster.uncorr", center = TRUE)
dev.off()
# 7.7: no p-values between clusters (they were defined on these very same data).


# (8) Batch correction
# Source: scrapbook ch. 8. MNN from Haghverdi et al. 2018.

# (8.1) MNN: mutual nearest neighbours across libraries; the batch difference is
#       subtracted (8.2).
mnn.out   <- correctMnn(pca$components, block = donor.f, num.threads = nthreads)
corrected <- mnn.out$corrected

# Check: the per-library PC1 medians should move closer together.
round(cbind(before = tapply(pca$components[1, ], donor.f, median),
            after  = tapply(corrected[1, ],      donor.f, median)), 2)

# (8.2) Warning (8.3-8.4): here block = diagnosis, so MNN erases part of the
#       disease signal. 'corrected' is used only to SEE and to CLUSTER;
#       expression is always measured on 'normalized'.

# (8a) Redo visualisation on the corrected space.
set.seed(1234)
tsne.corr <- runTsne(corrected, num.threads = nthreads)
umap.corr <- runUmap(corrected, num.threads = nthreads)
reducedDim(sce.f, "PCA.corr")  <- t(corrected)
reducedDim(sce.f, "TSNE.corr") <- tsne.corr
reducedDim(sce.f, "UMAP.corr") <- umap.corr

ggsave(fig("12_umap_corregido_por_donante.png"),     plotReducedDim(sce.f, "UMAP.corr", colour_by = COL_DONOR), width = 7, height = 6, dpi = 150)
ggsave(fig("13_umap_corregido_por_diagnostico.png"), plotReducedDim(sce.f, "UMAP.corr", colour_by = COL_DIAG),  width = 7, height = 6, dpi = 150)

png(fig("13b_umap_antes_vs_despues.png"), 2200, 1000, res = 150)
grid.arrange(
  plotReducedDim(sce.f, "UMAP",      colour_by = COL_DONOR) + ggtitle("Uncorrected"),
  plotReducedDim(sce.f, "UMAP.corr", colour_by = COL_DONOR) + ggtitle("Corrected"),
  ncol = 2)
dev.off()

# (8b) Redo clustering on the corrected space.
snn      <- buildSnnGraph(corrected, num.threads = nthreads)
clusters <- clusterGraph(snn)
colLabels(sce.f) <- factor(clusters$membership)
ggsave(fig("14_umap_corregido_por_cluster.png"), plotReducedDim(sce.f, "UMAP.corr", colour_by = "label"), width = 7, height = 6, dpi = 150)

# The number that matters: clusters dominated by one library, before vs after.
tab.corr <- table(cluster = clusters$membership, donor = donor.f)
c(before = sum(apply(prop.table(tab.uncorr, 1), 1, max) > 0.70),
  after  = sum(apply(prop.table(tab.corr,   1), 1, max) > 0.70))

# (8c) Redo markers on 'normalized' (rule from 8.4); only the labels change.
markers <- scoreMarkers(normalized, groups = clusters$membership,
                        block = donor.f, num.threads = nthreads)
for (k in levels(factor(clusters$membership))) {
  d <- reportGroupMarkerStatistics(markers, k)
  write.csv(d, fig(paste0("22_marcadores_corregido_cluster", k, ".csv")))
}
df <- reportGroupMarkerStatistics(markers, "1")
top.corr <- head(rownames(df)[order(df$cohens.d.mean, decreasing = TRUE)], 15)
png(fig("17_marcadores_corregido_heatmap.png"), 1600, 1400, res = 150)
plotGroupedHeatmap(sce.f, features = top.corr, group = "label", center = TRUE)
dev.off()

# (8d) Substate annotation with a panel of astrocyte markers.
astro.panel <- c("GFAP","VIM","CD44","CHI3L1","S100B",
                 "SLC1A2","SLC1A3","AQP4","GPC5","NRXN1",
                 "HSPA1A","HSPB1","CRYAB","DDIT4","FOS",
                 "MT2A","APOE","CLU","CST3")
astro.panel <- astro.panel[astro.panel %in% rownames(sce.f)]
ggsave(fig("18_panel_astrocito_dotplot.png"),
       plotDots(sce.f, features = astro.panel, group = "label", center = TRUE), width = 8, height = 7, dpi = 150)

# (8.5) Multi-condition: differential expression and abundance (edgeR) on the
#       ORIGINAL counts aggregated per sample (never on 'corrected').
#       Source: 8.5.1-8.5.2.
library(edgeR)

# Pseudobulk via matrix product (aggregateAcrossCells does not take a dgCMatrix).
ind <- Matrix::sparse.model.matrix(~ 0 + donor.f)
colnames(ind) <- levels(donor.f)
sums <- as.matrix(filtered %*% ind)                 # 10850 genes x 6 samples
info <- data.frame(library = colnames(sums),
                   diagnosis = diagnosis.f[match(colnames(sums), as.character(donor.f))])

y <- DGEList(sums, samples = info)
y <- y[filterByExpr(y, group = y$samples$diagnosis), , keep.lib.sizes = FALSE]
y <- normLibSizes(y)
design <- model.matrix(~ factor(diagnosis), y$samples)
y <- estimateDisp(y, design)
png(fig("19_bcv_pseudobulk.png"), 1400, 1200, res = 150); plotBCV(y); dev.off()

fit <- glmQLFit(y, design, robust = TRUE)
res <- glmQLFTest(fit, coef = 2)
summary(decideTests(res))                           # 47 down / 79 up
write.csv(topTags(res, n = Inf)$table, fig("20_DE_pseudobulk.csv"))

# Batch = diagnosis: library cannot enter as a covariate, the model loses rank.
m <- model.matrix(~ factor(diagnosis) + library, y$samples)
c(rank = qr(m)$rank, columns = ncol(m))

# Differential abundance of the substates.
abund  <- table(cluster = clusters$membership, library = donor.f)
y.ab   <- DGEList(unclass(abund), samples = info)
y.ab   <- estimateDisp(y.ab, design, trend = "none")
fit.ab <- glmQLFit(y.ab, design, robust = TRUE, abundance.trend = FALSE)
res.ab <- glmQLFTest(fit.ab, coef = 2)
write.csv(topTags(res.ab, n = Inf)$table, fig("21_abundancia_diferencial.csv"))

# (8.6) Replication: true n is 3 vs 3 and batch = condition. Grubman validates
#       the pipeline; concluding about AD needs cohorts where the two are not
#       confounded (Mathys, Gazestani).


# (9) Protein multi-omics
# Source: scrapbook ch. 9. Not applicable: Grubman is RNA only (no ADTs).
altExpNames(sce)     # character(0)


# Save
saveRDS(sce.f, file.path(project, "results/objects/grubman_astro_corregido.rds"))
sessionInfo()
