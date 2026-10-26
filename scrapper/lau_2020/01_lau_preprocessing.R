# Lau et al. 2020: preprocessing in a single script

# MSc thesis "Transcriptional Noise in Alzheimer's Disease Progression" (CBMSO)
# MSc in Bioinformatics and Data Science (ISCIII)

# Order of the steps:
#   (1) Quality control
#   (2) Ambient RNA with decontX and contamination filter
#   (3) Doublets with scDblFinder on the original matrix
#   (4) Normalisation on the decontX counts
#   (5) Gene filtering and check of the nucleus order
#   (6) Summary and saving
# Everything is computed here; the final object is saved for the noise script.

# Run from a terminal (not from the editor), because of memory use:
#   Rscript scrapper_lau_preprocessing.R > log_preprocessing.txt 2>&1

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

proyecto <- "/path/to/alzheimer_tn"   # change to the local project path
ruta     <- file.path(proyecto, "data/raw/Lau_32989152/Lau2020_sce.rds")
outdir   <- file.path(proyecto, "results/01_preprocesado/Lau_32989152")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

sce <- readRDS(ruta)
dim(sce)                         # 33538 genes x 179392 nuclei
table(sce$donor, sce$condition)  # each donor is either AD or control only
n_original <- ncol(sce)

# Copy of the original metadata of each nucleus, used at the end (step 5.1)
# to check that no nucleus has changed donor, condition or cell type.
meta_original <- as.data.frame(colData(sce))

# Original counts of 500 random nuclei (very little memory). Some will be
# lost in the filters; step 5.1 checks the ones that remain.
muestra <- sample(colnames(sce), 500)
cuentas_muestra <- counts(sce)[, muestra]
nucleos_por_donante <- data.frame(donante = names(table(sce$donor)),
                                  original = as.vector(table(sce$donor)))



# (1) QUALITY CONTROL (scrapbook ch. 1)
#     Thresholds at median +/- NUM_MADS MAD within each donor (default 3). The metrics are
#     already in the object (genes are Ensembl IDs, so mitochondrial genes
#     cannot be found by the "MT-" prefix).
#     Source: ?suggestRnaQcThresholds and ?filterRnaQcMetrics (scrapper).

qc <- list(sum      = sce$sum,
           detected = sce$detected,
           subsets  = list(MT = sce$subsets_mito_percent / 100))
NUM_MADS <- 2.5   # 3 = default; stricter (2.5 or 2) to remove more low-quality nuclei
umbrales <- suggestRnaQcThresholds(qc, block = sce$donor, num.mads = NUM_MADS)
sce$pasa_qc <- filterRnaQcMetrics(umbrales, qc, block = sce$donor)

table(sce$pasa_qc)
table(sce$condition, sce$pasa_qc)

sce <- sce[, sce$pasa_qc]
n_qc <- ncol(sce)
nucleos_por_donante$tras_qc <- as.vector(table(factor(sce$donor, levels = nucleos_por_donante$donante)))

# (1.1) To save memory: genes with no counts in any nucleus are removed
#       (they carry no information).
sce <- sce[rowSums(counts(sce)) > 0, ]
dim(sce)
anyDuplicated(colnames(sce))     # must be 0 (unique nucleus names)
gc()



# (2) AMBIENT RNA WITH decontX
#     decontX estimates which part of the counts of each nucleus comes from
#     free RNA in the sample and creates a new matrix, "decontXcounts", with
#     those counts removed. The original "counts" matrix is not modified.
#     It is run donor by donor: this is the same as decontX(sce, batch = donor),
#     which also treats each donor separately, but only one donor is held in
#     memory at a time. With all 163,000 nuclei at once the computer runs out
#     of RAM and the editor crashes.
#     Source: Yang et al. 2020, Genome Biol 21:57, doi:10.1186/s13059-020-1950-6
#     Vignette: https://bioconductor.org/packages/release/bioc/vignettes/celda/inst/doc/decontX.html

contaminacion <- numeric(ncol(sce))
corregidas    <- list()

for (d in unique(sce$donor)) {
  print(d)                                       # shows which donor is running
  cuales <- which(sce$donor == d)
  res <- decontX(counts(sce)[, cuales], verbose = FALSE)
  contaminacion[cuales] <- res$contamination
  corregidas[[d]] <- round(res$decontXcounts)    # integer counts
  rm(res)
  gc()
}

# The donors are joined back in the same nucleus order as the object.
corregidas <- do.call(cbind, corregidas)
sce$decontX_contamination <- contaminacion
assay(sce, "decontXcounts") <- corregidas[, colnames(sce)]
rm(corregidas)
gc()
assayNames(sce)                  # "counts" "decontXcounts"

# (2.1) Distribution of the contamination, to choose the threshold.
summary(sce$decontX_contamination)

png(file.path(outdir, "01_histograma_contaminacion.png"), 1400, 1000, res = 150)
hist(sce$decontX_contamination, breaks = 50,
     xlab = "Estimated contamination", main = "decontX: all nuclei")
abline(v = 0.5, col = "red")
dev.off()

# (2.2) Nuclei that would be lost with different thresholds.
sum(sce$decontX_contamination > 0.3)
sum(sce$decontX_contamination > 0.4)
sum(sce$decontX_contamination > 0.5)

# (2.3) Filter: contamination > 0.5. Nuclei left with fewer than 100
#       corrected counts are also removed, because normalisation fails
#       with 0 counts.
sce$discard_contam <- sce$decontX_contamination > 0.5 |
                      colSums(assay(sce, "decontXcounts")) < 100
table(sce$discard_contam)
table(sce$condition, sce$discard_contam)

sce <- sce[, !sce$discard_contam]
n_contam <- ncol(sce)
nucleos_por_donante$tras_contam <- as.vector(table(factor(sce$donor, levels = nucleos_por_donante$donante)))
gc()


# (2.4) LOSS PER DONOR IN THE CONTAMINATION FILTER
#       In an earlier version AD6 lost 24 % of its nuclei in this filter and
#       was excluded. With the current preprocessing it loses 10 %, like AD4
#       (9.8 %): it no longer stands out, so all 21 donors are kept (12 AD,
#       9 control). The loss per donor is still printed to check that no
#       donor separates from the rest.

perdida <- 1 - nucleos_por_donante$tras_contam / nucleos_por_donante$tras_qc
names(perdida) <- nucleos_por_donante$donante
round(sort(perdida, decreasing = TRUE), 3)
table(sce$condition)



# (3) DOUBLETS WITH scDblFinder
#     Detected on the original matrix (counts), which is scDblFinder's
#     default: the corrected matrix can remove the mixture between cell
#     types that the method needs to see. samples = "donor" analyses each
#     donor separately.
#     Source: Germain et al. 2021, F1000Research 10:979, doi:10.12688/f1000research.73600.2
#     Vignette: https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html

sce <- scDblFinder(sce, samples = "donor", BPPARAM = BiocParallel::SerialParam())

table(sce$scDblFinder.class)
table(sce$condition, sce$scDblFinder.class)

sce <- sce[, sce$scDblFinder.class == "singlet"]
n_dobletes <- ncol(sce)
nucleos_por_donante$tras_dobletes <- as.vector(table(factor(sce$donor, levels = nucleos_por_donante$donante)))
gc()



# (4) NORMALISATION (scrapbook ch. 2)
#     On the decontX-corrected counts: size factors proportional to the
#     total counts, centred within each donor, and log2.
#     scrapper:: is written in front because celda has another function
#     with the same name.
#     Source: ?centerSizeFactors and ?normalizeCounts (scrapper).

sf <- centerSizeFactors(colSums(assay(sce, "decontXcounts")), block = sce$donor)
summary(sf)                      # the minimum must be greater than 0
logcounts(sce) <- scrapper::normalizeCounts(assay(sce, "decontXcounts"), sf)



# (5) GENE FILTERING
#     Genes detected in fewer than 3 nuclei are removed.

keep_genes <- rowSums(counts(sce) > 0) >= 3
table(keep_genes)
sce <- sce[keep_genes, ]



# (5.1) CHECK OF THE NUCLEUS ORDER
#       After splitting by donor (decontX, scDblFinder) and joining back,
#       each column must still be the same nucleus, with its metadata and
#       its counts; matching dimensions is not enough.
#       If anything fails, stopifnot() stops the script with an error.

# a) Nucleus names are unique and all present in the original object.
stopifnot(anyDuplicated(colnames(sce)) == 0)
stopifnot(all(colnames(sce) %in% rownames(meta_original)))

# b) The metadata of each nucleus matches the original.
idx <- match(colnames(sce), rownames(meta_original))
for (v in intersect(colnames(meta_original), colnames(colData(sce)))) {
  stopifnot(identical(as.vector(sce[[v]]), as.vector(meta_original[[v]][idx])))
}
print("Metadata: OK")

# c) The count matrices have the nuclei in the same order.
stopifnot(identical(colnames(assay(sce, "counts")), colnames(sce)))
stopifnot(identical(colnames(assay(sce, "decontXcounts")), colnames(sce)))
stopifnot(identical(colnames(logcounts(sce)), colnames(sce)))

# d) decontX only removes counts: in each nucleus the corrected counts
#    cannot exceed the original ones (this would fail if the nuclei had
#    been shuffled when joining the donors).
stopifnot(all(colSums(assay(sce, "decontXcounts")) <= colSums(counts(sce)) + 1))
print("decontX: OK")

# e) The counts of the sampled nuclei still in the object are identical to
#    the original ones (in the remaining genes). The original object is not
#    read again, to keep the RAM peak low.
quedan <- intersect(muestra, colnames(sce))
length(quedan)                   # how many of the 500 passed the filters
stopifnot(all(counts(sce)[, quedan] ==
              cuentas_muestra[rownames(sce), quedan]))
rm(cuentas_muestra)
print("Cuentas: OK")

# f) Nuclei per donor at each step, to see whether any donor loses many
#    nuclei or disappears.
nucleos_por_donante$condicion <- meta_original$condition[
  match(nucleos_por_donante$donante, meta_original$donor)]
nucleos_por_donante
write.csv(nucleos_por_donante, file.path(outdir, "03_nucleos_por_donante.csv"),
          row.names = FALSE)
table(sce$donor, sce$condition)



# (6) SUMMARY AND SAVING

resumen <- data.frame(
  paso    = c("originales", "tras QC", "tras filtrar contaminacion",
              "tras quitar dobletes"),
  nucleos = c(n_original, n_qc, n_contam, n_dobletes))
resumen
write.csv(resumen, file.path(outdir, "02_resumen_filtrado.csv"), row.names = FALSE)

dim(sce)
saveRDS(sce, file.path(proyecto, "results/objects/lau_2020/lau_sce_preprocesado.rds"))
sessionInfo()