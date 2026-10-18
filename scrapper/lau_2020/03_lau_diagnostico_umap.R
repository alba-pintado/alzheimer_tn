# Lau et al. 2020: diagnóstico de las poblaciones pequeñas del UMAP
#
# TFM "Transcriptional Noise in Alzheimer's Disease Progression" (CBMSO)
#
# El UMAP muestra muchas islas pequeñas repartidas por el espacio. Este script
# no modifica nada: describe cada isla para saber si son subtipos reales,
# núcleos de baja calidad, dobletes residuales o un efecto de la corrección MNN.
#
#   (1) Islas del UMAP: grupos de núcleos separados en el dibujo (DBSCAN)
#   (2) Tabla por isla: tamaño, tipo celular, donantes y métricas de calidad
#   (3) UMAP coloreado por las métricas de calidad
#   (4) UMAP sin corregir (sobre el PCA, sin MNN) para comparar

library(SingleCellExperiment)
library(scrapper)
library(scater)
library(ggplot2)
library(dbscan)                    # islas del UMAP: https://cran.r-project.org/package=dbscan

set.seed(1234)
proyecto <- "~/Documentos/proyectos/alzheimer_tn"   # carpeta del proyecto
outdir <- file.path(proyecto, "results/02_analisis/Lau_32989152/diagnostico_umap")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

sce <- readRDS(file.path(proyecto, "results/objects/lau_2020/lau_sce_analisis.rds"))
dim(sce)
names(colData(sce))                # columnas disponibles



# (1) ISLAS DEL UMAP
#     DBSCAN agrupa los puntos que están juntos en el dibujo: cada grupo es una
#     isla. eps es la distancia máxima entre vecinos (en unidades del UMAP).

umap <- reducedDim(sce, "UMAP")
islas <- dbscan(umap, eps = 0.3, minPts = 20)
sce$isla <- factor(islas$cluster)  # 0 = puntos sueltos que no forman isla
table(sce$isla == 0)
length(unique(islas$cluster[islas$cluster > 0]))   # número de islas



# (2) TABLA POR ISLA
#     Para cada isla: núcleos, tipo celular mayoritario y su fracción, número de
#     donantes, fracción del donante principal, fracción AD, y medianas de genes
#     detectados, % mitocondrial, contaminación (decontX) y puntuación de doblete.

# Algunas columnas pueden no existir; se usan solo las que haya.
metrica <- function(v) if (v %in% names(colData(sce))) sce[[v]] else rep(NA, ncol(sce))
mito  <- metrica("subsets_mito_percent")
contam <- metrica("decontX_contamination")
dbl   <- metrica("scDblFinder.score")
genes <- metrica("detected")

filas <- list()
for (k in levels(sce$isla)) {
  sel <- which(sce$isla == k)
  tipos <- table(sce$celltype[sel], useNA = "ifany")
  donantes <- table(sce$donor[sel])
  filas[[k]] <- data.frame(
    isla           = k,
    nucleos        = length(sel),
    tipo           = names(tipos)[which.max(tipos)],
    frac_tipo      = round(max(tipos) / length(sel), 2),
    clusters       = paste(names(sort(table(sce$cluster[sel]), decreasing = TRUE))[1:min(3, length(unique(sce$cluster[sel])))], collapse = ","),
    n_donantes     = sum(donantes > 0),
    frac_donante   = round(max(donantes) / length(sel), 2),
    frac_AD        = round(mean(sce$condition[sel] == "AD"), 2),
    genes_mediana  = round(median(genes[sel])),
    mito_mediana   = round(median(mito[sel]), 2),
    contam_mediana = round(median(contam[sel]), 3),
    doblete_mediana = round(median(dbl[sel]), 3))
}
tabla <- do.call(rbind, filas)
tabla <- tabla[order(tabla$nucleos), ]

# Referencia: valores medianos de todo el objeto, para comparar cada isla.
print(c(genes = median(genes), mito = median(mito),
        contam = median(contam), doblete = median(dbl)))
print(tabla, row.names = FALSE)
write.csv(tabla, file.path(outdir, "01_tabla_islas.csv"), row.names = FALSE)

# Señales de alarma: islas pequeñas (< 500 núcleos) con algo raro.
pequenas <- tabla[tabla$isla != "0" & tabla$nucleos < 500, ]
print("Islas pequeñas (< 500 núcleos):")
print(pequenas, row.names = FALSE)
print(paste("Con un donante > 50 %:", sum(pequenas$frac_donante > 0.5)))
print(paste("Con % mito > 2 veces la mediana:", sum(pequenas$mito_mediana > 2 * median(mito), na.rm = TRUE)))
print(paste("Con genes < mitad de la mediana:", sum(pequenas$genes_mediana < median(genes) / 2, na.rm = TRUE)))



# (3) UMAP COLOREADO POR CALIDAD Y POR ISLA
orden <- sample(ncol(sce))
sce_plot <- sce[, orden]

ggsave(file.path(outdir, "02_umap_islas.png"),
       plotReducedDim(sce_plot, "UMAP", colour_by = "isla", text_by = "isla",
                      point_size = 0.3) + theme(legend.position = "none"),
       width = 9, height = 8)
for (v in c("subsets_mito_percent", "decontX_contamination", "scDblFinder.score", "detected")) {
  if (v %in% names(colData(sce))) {
    ggsave(file.path(outdir, paste0("03_umap_", v, ".png")),
           plotReducedDim(sce_plot, "UMAP", colour_by = v, point_size = 0.3),
           width = 8, height = 7)
  }
}



# (4) UMAP SIN CORREGIR
#     Mismo UMAP pero sobre el PCA sin MNN. Si las islas pequeñas desaparecen o
#     se unen a las grandes, las crea la corrección; si siguen, están en los datos.

reducedDim(sce, "UMAP_sin_corregir") <- runUmap(t(reducedDim(sce, "PCA")), num.threads = 2)
sce_plot <- sce[, orden]
ggsave(file.path(outdir, "04_umap_sin_corregir_celltype.png"),
       plotReducedDim(sce_plot, "UMAP_sin_corregir", colour_by = "celltype", point_size = 0.3),
       width = 8, height = 7)
ggsave(file.path(outdir, "05_umap_sin_corregir_donante.png"),
       plotReducedDim(sce_plot, "UMAP_sin_corregir", colour_by = "donor", point_size = 0.3),
       width = 8, height = 7)

# Cuántas islas tiene el UMAP sin corregir, con el mismo criterio.
islas_sc <- dbscan(reducedDim(sce, "UMAP_sin_corregir"), eps = 0.3, minPts = 20)
print(paste("Islas con MNN:", length(unique(islas$cluster[islas$cluster > 0])),
            "| sin corregir:", length(unique(islas_sc$cluster[islas_sc$cluster > 0]))))

sessionInfo()