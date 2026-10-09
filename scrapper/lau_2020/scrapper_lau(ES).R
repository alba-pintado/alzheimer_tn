# Lau et al. 2020: preprocesado de snRNA-seq en un unico script
#
# TFM "Transcriptional Noise in Alzheimer's Disease Progression" (CBMSO)
# Alba Pintado. Master en Bioinformatica y Ciencia de Datos (ISCIII)
#
# Antes de ejecutar: cambiar la variable "proyecto" por la ruta local.
#
# Pasos:
#   (1) Control de calidad
#   (2) ARN ambiental con decontX y filtrado por contaminacion
#   (3) Dobletes con scDblFinder sobre la matriz original
#   (4) Normalizacion sobre las cuentas de decontX
#   (5) Filtrado de genes
#   (6) Resumen y guardado
# Todo se calcula aqui; al final se guarda el objeto para el script de ruido.
#
# Fuentes generales:
#   Libro de scrapper: https://bioconductor.org/books/release/scrapbook/
#   OSCA (Amezquita et al. 2020, Nat Methods, doi:10.1038/s41592-019-0654-x):
#     https://bioconductor.org/books/release/OSCA/
#   Datos: Lau et al. 2020, PNAS, doi:10.1073/pnas.2008762117



# (0) PREPARACION

library(SingleCellExperiment)  # https://bioconductor.org/packages/SingleCellExperiment
library(scrapper)              # https://bioconductor.org/packages/scrapper
library(celda)                 # decontX: https://bioconductor.org/packages/celda
library(scDblFinder)           # https://bioconductor.org/packages/scDblFinder

set.seed(1234)

proyecto <- "/home/vant/Documentos/proyectos/alzheimer_tn"
ruta     <- file.path(proyecto, "data/raw/Lau_32989152/Lau2020_sce.rds")
outdir   <- file.path(proyecto, "results/01_preprocesado/Lau_32989152")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

sce <- readRDS(ruta)
dim(sce)                         # 33538 genes x 179392 nucleos
table(sce$donor, sce$condition)  # cada donante es solo AD o solo control
n_original <- ncol(sce)



# (1) CONTROL DE CALIDAD (scrapbook cap. 1)
#     Umbrales mediana +/- 3 MAD dentro de cada donante. Las metricas ya
#     vienen calculadas en el objeto (los genes estan en Ensembl, asi que
#     no se pueden buscar los mitocondriales por "MT-").
#     Fuente: ?suggestRnaQcThresholds y ?filterRnaQcMetrics (scrapper).

qc <- list(sum      = sce$sum,
           detected = sce$detected,
           subsets  = list(MT = sce$subsets_mito_percent / 100))
umbrales <- suggestRnaQcThresholds(qc, block = sce$donor)
sce$pasa_qc <- filterRnaQcMetrics(umbrales, qc, block = sce$donor)

table(sce$pasa_qc)
table(sce$condition, sce$pasa_qc)

sce <- sce[, sce$pasa_qc]
n_qc <- ncol(sce)

# (1.1) Para ahorrar memoria: se quitan los genes que no tienen ninguna
#       cuenta en ningun nucleo (no aportan nada).
sce <- sce[rowSums(counts(sce)) > 0, ]
dim(sce)
anyDuplicated(colnames(sce))     # tiene que dar 0 (nombres de nucleo unicos)
gc()



# (2) ARN AMBIENTAL CON decontX
#     decontX estima que parte de las cuentas de cada nucleo viene del ARN
#     libre de la muestra y crea una matriz nueva, "decontXcounts", con esas
#     cuentas restadas. La matriz original "counts" no se toca.
#     Se calcula donante a donante: es lo mismo que decontX(sce, batch = donor),
#     que tambien trata cada donante por separado, pero solo tiene un donante
#     en memoria cada vez. Con los 163.000 nucleos a la vez el ordenador se
#     queda sin RAM y el editor se cierra.
#     Fuente: Yang et al. 2020, Genome Biol 21:57, doi:10.1186/s13059-020-1950-6
#     Vinieta: https://bioconductor.org/packages/release/bioc/vignettes/celda/inst/doc/decontX.html

contaminacion <- numeric(ncol(sce))
corregidas    <- list()

for (d in unique(sce$donor)) {
  print(d)                                       # para ver por que donante va
  cuales <- which(sce$donor == d)
  res <- decontX(counts(sce)[, cuales], verbose = FALSE)
  contaminacion[cuales] <- res$contamination
  corregidas[[d]] <- round(res$decontXcounts)    # cuentas enteras
  rm(res)
  gc()
}

# Se juntan los donantes en el mismo orden de nucleos que el objeto.
corregidas <- do.call(cbind, corregidas)
sce$decontX_contamination <- contaminacion
assay(sce, "decontXcounts") <- corregidas[, colnames(sce)]
rm(corregidas)
gc()
assayNames(sce)                  # "counts" "decontXcounts"

# (2.1) Distribucion de la contaminacion para decidir el umbral.
summary(sce$decontX_contamination)

png(file.path(outdir, "01_histograma_contaminacion.png"), 1400, 1000, res = 150)
hist(sce$decontX_contamination, breaks = 50,
     xlab = "Contaminacion estimada", main = "decontX: todos los nucleos")
abline(v = 0.5, col = "red")
dev.off()

# (2.2) Nucleos que se perderian con distintos umbrales.
sum(sce$decontX_contamination > 0.3)
sum(sce$decontX_contamination > 0.4)
sum(sce$decontX_contamination > 0.5)

# (2.3) Filtrado: contaminacion > 0,5. Tambien los nucleos que se quedan
#       con menos de 100 cuentas corregidas, porque con 0 cuentas la
#       normalizacion falla.
sce$discard_contam <- sce$decontX_contamination > 0.5 |
                      colSums(assay(sce, "decontXcounts")) < 100
table(sce$discard_contam)
table(sce$condition, sce$discard_contam)

sce <- sce[, !sce$discard_contam]
n_contam <- ncol(sce)
gc()



# (3) DOBLETES CON scDblFinder
#     Se buscan sobre la matriz original (counts), que es la que usa
#     scDblFinder por defecto: la corregida puede borrar la mezcla entre
#     tipos celulares que el metodo necesita ver. samples = "donor" hace que
#     cada donante se analice por separado.
#     Fuente: Germain et al. 2021, F1000Research 10:979, doi:10.12688/f1000research.73600.2
#     Vinieta: https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html

sce <- scDblFinder(sce, samples = "donor", BPPARAM = BiocParallel::SerialParam())

table(sce$scDblFinder.class)
table(sce$condition, sce$scDblFinder.class)

sce <- sce[, sce$scDblFinder.class == "singlet"]
n_dobletes <- ncol(sce)
gc()



# (4) NORMALIZACION (scrapbook cap. 2)
#     Sobre las cuentas corregidas por decontX: size factors proporcionales
#     al total de cuentas, centrados dentro de cada donante, y log2.
#     Se escribe scrapper:: delante porque celda tiene otra funcion con el
#     mismo nombre.
#     Fuente: ?centerSizeFactors y ?normalizeCounts (scrapper).

sf <- centerSizeFactors(colSums(assay(sce, "decontXcounts")), block = sce$donor)
summary(sf)                      # el minimo tiene que ser mayor que 0
logcounts(sce) <- scrapper::normalizeCounts(assay(sce, "decontXcounts"), sf)



# (5) FILTRADO DE GENES
#     Se quitan los genes detectados en menos de 3 nucleos.

keep_genes <- rowSums(counts(sce) > 0) >= 3
table(keep_genes)
sce <- sce[keep_genes, ]



# (6) RESUMEN Y GUARDADO

resumen <- data.frame(
  paso    = c("originales", "tras QC", "tras filtrar contaminacion", "tras quitar dobletes"),
  nucleos = c(n_original, n_qc, n_contam, n_dobletes))
resumen
write.csv(resumen, file.path(outdir, "02_resumen_filtrado.csv"), row.names = FALSE)

dim(sce)
saveRDS(sce, file.path(proyecto, "results/objects/lau_sce_preprocesado.rds"))
sessionInfo()