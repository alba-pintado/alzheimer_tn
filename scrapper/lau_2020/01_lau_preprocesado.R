# Lau et al. 2020: preprocesado en un unico script

# TFM "Transcriptional Noise in Alzheimer's Disease Progression" (CBMSO)
# Master en Bioinformatica y Ciencia de Datos (ISCIII)

# Orden de los pasos:
#   (1) Control de calidad
#   (2) ARN ambiental con decontX y filtrado por contaminacion
#   (3) Dobletes con scDblFinder sobre la matriz original
#   (4) Normalizacion sobre las cuentas de decontX
#   (5) Filtrado de genes y comprobacion del orden de los nucleos
#   (6) Resumen y guardado
# Todo se calcula aqui; al final se guarda el objeto para el script de ruido.

# Usa mucha memoria (unos 20 GB en el preprocesado): si el editor se cierra,
# lanzarlo desde la terminal con Rscript.

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

proyecto <- "~/Documentos/proyectos/alzheimer_tn"   # carpeta del proyecto
ruta     <- file.path(proyecto, "data/raw/Lau_32989152/Lau2020_sce.rds")
outdir   <- file.path(proyecto, "results/01_preprocesado/Lau_32989152")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

sce <- readRDS(ruta)
dim(sce)                         # 33538 genes x 179392 nucleos
table(sce$donor, sce$condition)  # cada donante es solo AD o solo control
n_original <- ncol(sce)

# Copia de la metadata original de cada nucleo, para comprobar al final
# (paso 5.1) que ningun nucleo ha cambiado de donante, condicion o tipo.
meta_original <- as.data.frame(colData(sce))

# Cuentas originales de 500 nucleos al azar (ocupa muy poco). Algunos se
# perderan en los filtros; en el paso 5.1 se comprueban los que queden.
muestra <- sample(colnames(sce), 500)
cuentas_muestra <- counts(sce)[, muestra]
nucleos_por_donante <- data.frame(donante = names(table(sce$donor)),
                                  original = as.vector(table(sce$donor)))



# (1) CONTROL DE CALIDAD (scrapbook cap. 1)
#     Umbrales mediana +/- NUM_MADS MAD dentro de cada donante (por defecto 3). Las metricas ya
#     vienen calculadas en el objeto (los genes estan en Ensembl, asi que
#     no se pueden buscar los mitocondriales por "MT-").
#     Fuente: ?suggestRnaQcThresholds y ?filterRnaQcMetrics (scrapper).

qc <- list(sum      = sce$sum,
           detected = sce$detected,
           subsets  = list(MT = sce$subsets_mito_percent / 100))
NUM_MADS <- 2.5   # 3 = por defecto; mas estricto (2.5 o 2) si el UMAP se resuelve
umbrales <- suggestRnaQcThresholds(qc, block = sce$donor, num.mads = NUM_MADS)
sce$pasa_qc <- filterRnaQcMetrics(umbrales, qc, block = sce$donor)

table(sce$pasa_qc)
table(sce$condition, sce$pasa_qc)

sce <- sce[, sce$pasa_qc]
n_qc <- ncol(sce)
nucleos_por_donante$tras_qc <- as.vector(table(factor(sce$donor, levels = nucleos_por_donante$donante)))

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
nucleos_por_donante$tras_contam <- as.vector(table(factor(sce$donor, levels = nucleos_por_donante$donante)))
gc()


# (2.4) PERDIDA POR DONANTE EN EL FILTRO DE CONTAMINACION
#       En una version anterior AD6 perdia un 24 % de sus nucleos en este
#       filtro y se excluyo. Al repetir el preprocesado pierde un 10 %, igual
#       que AD4 (9,8 %): ya no destaca, asi que se mantienen los 21 donantes
#       (12 AD, 9 control). Se sigue calculando la perdida por donante para
#       vigilar que ninguno se separe del resto.

perdida <- 1 - nucleos_por_donante$tras_contam / nucleos_por_donante$tras_qc
names(perdida) <- nucleos_por_donante$donante
round(sort(perdida, decreasing = TRUE), 3)
table(sce$condition)



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
nucleos_por_donante$tras_dobletes <- as.vector(table(factor(sce$donor, levels = nucleos_por_donante$donante)))
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



# (5.1) COMPROBACION DEL ORDEN DE LOS NUCLEOS
#       Tras separar por donante (decontX, scDblFinder) y volver a unir, hay
#       que confirmar que cada columna sigue siendo el mismo nucleo, con su
#       metadata y sus cuentas, y no solo que las dimensiones cuadran.
#       Si algo falla, stopifnot() detiene el script con un error.

# a) Nombres de nucleo unicos y todos presentes en el objeto original.
stopifnot(anyDuplicated(colnames(sce)) == 0)
stopifnot(all(colnames(sce) %in% rownames(meta_original)))

# b) La metadata de cada nucleo coincide con la original.
idx <- match(colnames(sce), rownames(meta_original))
for (v in intersect(colnames(meta_original), colnames(colData(sce)))) {
  stopifnot(identical(as.vector(sce[[v]]), as.vector(meta_original[[v]][idx])))
}
print("Metadata: OK")

# c) Las matrices de cuentas tienen los nucleos en el mismo orden.
stopifnot(identical(colnames(assay(sce, "counts")), colnames(sce)))
stopifnot(identical(colnames(assay(sce, "decontXcounts")), colnames(sce)))
stopifnot(identical(colnames(logcounts(sce)), colnames(sce)))

# d) decontX solo resta cuentas: en cada nucleo, las cuentas corregidas
#    no pueden superar a las originales (si se hubieran desordenado los
#    nucleos al unir los donantes, esto fallaria).
stopifnot(all(colSums(assay(sce, "decontXcounts")) <= colSums(counts(sce)) + 1))
print("decontX: OK")

# e) Las cuentas de los nucleos de la muestra que siguen en el objeto son
#    identicas a las originales (en los genes que quedan). No se vuelve a
#    leer el objeto original para no subir el pico de RAM.
quedan <- intersect(muestra, colnames(sce))
length(quedan)                   # cuantos de los 500 han pasado los filtros
stopifnot(all(counts(sce)[, quedan] ==
              cuentas_muestra[rownames(sce), quedan]))
rm(cuentas_muestra)
print("Cuentas: OK")

# f) Nucleos por donante en cada paso. Sirve para ver si algun donante
#    pierde muchos nucleos o desaparece.
nucleos_por_donante$condicion <- meta_original$condition[
  match(nucleos_por_donante$donante, meta_original$donor)]
nucleos_por_donante
write.csv(nucleos_por_donante, file.path(outdir, "03_nucleos_por_donante.csv"),
          row.names = FALSE)
table(sce$donor, sce$condition)



# (6) RESUMEN Y GUARDADO

resumen <- data.frame(
  paso    = c("originales", "tras QC", "tras filtrar contaminacion",
              "tras quitar dobletes"),
  nucleos = c(n_original, n_qc, n_contam, n_dobletes))
resumen
write.csv(resumen, file.path(outdir, "02_resumen_filtrado.csv"), row.names = FALSE)

dim(sce)
saveRDS(sce, file.path(proyecto, "results/objects/lau_2020/lau_sce_preprocesado.rds"))
sessionInfo()