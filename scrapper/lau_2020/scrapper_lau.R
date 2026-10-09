# Lau et al. 2020: snRNA-seq preprocessing in a single script
#
# MSc thesis "Transcriptional Noise in Alzheimer's Disease Progression" (CBMSO)
# Alba Pintado. MSc in Bioinformatics and Data Science (ISCIII)
#
# Before running: set the "project" variable to your local path.
#
# Steps:
#   (1) Quality control
#   (2) Ambient RNA with decontX and contamination filter
#   (3) Doublets with scDblFinder on the original counts
#   (4) Normalisation on the decontX counts
#   (5) Gene filtering
#   (6) Summary and saving
# Everything is computed here; the final object is saved for the noise script.
#
# General sources:
#   scrapper book: https://bioconductor.org/books/release/scrapbook/
#   OSCA (Amezquita et al. 2020, Nat Methods, doi:10.1038/s41592-019-0654-x):
#     https://bioconductor.org/books/release/OSCA/
#   Data: Lau et al. 2020, PNAS, doi:10.1073/pnas.2008762117



# (0) SETUP

library(SingleCellExperiment)  # https://bioconductor.org/packages/SingleCellExperiment
library(scrapper)              # https://bioconductor.org/packages/scrapper
library(celda)                 # decontX: https://bioconductor.org/packages/celda
library(scDblFinder)           # https://bioconductor.org/packages/scDblFinder

set.seed(1234)

project <- "/home/vant/Documentos/proyectos/alzheimer_tn"
path    <- file.path(project, "data/raw/Lau_32989152/Lau2020_sce.rds")
outdir  <- file.path(project, "results/01_preprocessing/Lau_32989152")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

sce <- readRDS(path)
dim(sce)                         # 33538 genes x 179392 nuclei
table(sce$donor, sce$condition)  # each donor is either AD or control
n_original <- ncol(sce)



# (1) QUALITY CONTROL (scrapbook ch. 1)
#     Median +/- 3 MAD thresholds within each donor. The metrics are already
#     stored in the object (genes are Ensembl IDs, so mitochondrial genes
#     cannot be found with "MT-").
#     Source: ?suggestRnaQcThresholds and ?filterRnaQcMetrics (scrapper).

qc <- list(sum      = sce$sum,
           detected = sce$detected,
           subsets  = list(MT = sce$subsets_mito_percent / 100))
thresholds <- suggestRnaQcThresholds(qc, block = sce$donor)
sce$pass_qc <- filterRnaQcMetrics(thresholds, qc, block = sce$donor)

table(sce$pass_qc)
table(sce$condition, sce$pass_qc)

sce <- sce[, sce$pass_qc]
n_qc <- ncol(sce)

# (1.1) To save memory: remove genes with no counts in any nucleus.
sce <- sce[rowSums(counts(sce)) > 0, ]
dim(sce)
anyDuplicated(colnames(sce))     # must be 0 (unique nucleus names)
gc()



# (2) AMBIENT RNA WITH decontX
#     decontX estimates which fraction of each nucleus' counts comes from the
#     free RNA in the sample and creates a new matrix, "decontXcounts", with
#     those counts removed. The original "counts" matrix is not modified.
#     It is run donor by donor: this is the same as decontX(sce, batch = donor),
#     which also treats each donor separately, but keeps only one donor in
#     memory at a time. With all 163,000 nuclei at once the computer ran out
#     of RAM.
#     Source: Yang et al. 2020, Genome Biol 21:57, doi:10.1186/s13059-020-1950-6
#     Vignette: https://bioconductor.org/packages/release/bioc/vignettes/celda/inst/doc/decontX.html

contamination <- numeric(ncol(sce))
corrected     <- list()

for (d in unique(sce$donor)) {
  print(d)                                       # shows which donor is running
  idx <- which(sce$donor == d)
  res <- decontX(counts(sce)[, idx], verbose = FALSE)
  contamination[idx] <- res$contamination
  corrected[[d]] <- round(res$decontXcounts)     # integer counts
  rm(res)
  gc()
}

# Join the donors in the same nucleus order as the object.
corrected <- do.call(cbind, corrected)
sce$decontX_contamination <- contamination
assay(sce, "decontXcounts") <- corrected[, colnames(sce)]
rm(corrected)
gc()
assayNames(sce)                  # "counts" "decontXcounts"

# (2.1) Contamination distribution, to choose the threshold.
summary(sce$decontX_contamination)

png(file.path(outdir, "01_contamination_histogram.png"), 1400, 1000, res = 150)
hist(sce$decontX_contamination, breaks = 50,
     xlab = "Estimated contamination", main = "decontX: all nuclei")
abline(v = 0.5, col = "red")
dev.off()

# (2.2) Nuclei that would be lost with different thresholds.
sum(sce$decontX_contamination > 0.3)
sum(sce$decontX_contamination > 0.4)
sum(sce$decontX_contamination > 0.5)

# (2.3) Filter: contamination > 0.5. Also nuclei left with fewer than 100
#       corrected counts, because normalisation fails with 0 counts.
sce$discard_contam <- sce$decontX_contamination > 0.5 |
                      colSums(assay(sce, "decontXcounts")) < 100
table(sce$discard_contam)
table(sce$condition, sce$discard_contam)

sce <- sce[, !sce$discard_contam]
n_contam <- ncol(sce)
gc()



# (3) DOUBLETS WITH scDblFinder
#     Doublets are searched on the original counts ("counts"), which is the
#     scDblFinder default: the corrected matrix may remove the mixture of
#     cell types the method needs to see. samples = "donor" analyses each
#     donor separately.
#     Source: Germain et al. 2021, F1000Research 10:979, doi:10.12688/f1000research.73600.2
#     Vignette: https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html

sce <- scDblFinder(sce, samples = "donor", BPPARAM = BiocParallel::SerialParam())

table(sce$scDblFinder.class)
table(sce$condition, sce$scDblFinder.class)

sce <- sce[, sce$scDblFinder.class == "singlet"]
n_doublets <- ncol(sce)
gc()



# (4) NORMALISATION (scrapbook ch. 2)
#     On the decontX-corrected counts: size factors proportional to the total
#     counts, centred within each donor, and log2.
#     scrapper:: is written in front because celda has another function with
#     the same name.
#     Source: ?centerSizeFactors and ?normalizeCounts (scrapper).

sf <- centerSizeFactors(colSums(assay(sce, "decontXcounts")), block = sce$donor)
summary(sf)                      # the minimum must be above 0
logcounts(sce) <- scrapper::normalizeCounts(assay(sce, "decontXcounts"), sf)



# (5) GENE FILTERING
#     Remove genes detected in fewer than 3 nuclei.

keep_genes <- rowSums(counts(sce) > 0) >= 3
table(keep_genes)
sce <- sce[keep_genes, ]



# (6) SUMMARY AND SAVING

summary_table <- data.frame(
  step   = c("original", "after QC", "after contamination filter", "after doublet removal"),
  nuclei = c(n_original, n_qc, n_contam, n_doublets))
summary_table
write.csv(summary_table, file.path(outdir, "02_filtering_summary.csv"), row.names = FALSE)

dim(sce)
saveRDS(sce, file.path(project, "results/objects/lau_sce_preprocessed.rds"))
sessionInfo()