
# 01_exploration.R

# Análisis exploratorio de objetos SingleCellExperiment
# Dataset: Grubman_31768052 — 7 tipos celulares


# Este script carga 7 objetos SingleCellExperiment (SCE) en formato
# .rds, uno por tipo celular del dataset Grubman et al. 2019
# (snRNA-seq de corteza entorrinal humana, AD vs Control).
# Realiza un análisis exploratorio completo: estructura de los
# objetos, metadatos de células y genes, métricas de calidad (QC),
# y genera visualizaciones en PDF + tablas resumen en CSV.


# 0. CARGA DE LIBRERÍAS 
# suppressPackageStartupMessages() evita que se impriman mensajes
# de carga de paquetes en la consola, manteniendo el output limpio.
suppressPackageStartupMessages({

  # SingleCellExperiment: clase S4 de Bioconductor que representa
  # datos de single-cell RNA-seq. Los objetos .rds de este proyecto
  # son instancias de esta clase. Proporciona acceso a matrices de
  # expresión (assays), metadatos de células (colData), metadatos
  # de genes (rowData) y reducciones dimensionales (reducedDims).
  
  library(SingleCellExperiment)

  # SummarizedExperiment: clase base de la que hereda
  # SingleCellExperiment. Necesaria para funciones como assay(),
  # colData() y rowData() que acceden a los componentes del objeto.
  
  library(SummarizedExperiment)

  # scuttle: paquete de Bioconductor para control de calidad en
  # datos single-cell. Usamos perCellQCMetrics() para calcular
  # por célula: nº total de UMIs (sum), nº de genes detectados
  # (detected) y porcentaje de reads mitocondriales (%MT).
  
  library(scuttle)

  # ggplot2: sistema de visualización basado en la gramática de
  # gráficos. Se usa para generar todos los plots del script:
  # violines, barplots, scatter plots.
  
  library(ggplot2)

  # dplyr: manipulación de data.frames de forma legible.
  # Se usa para agrupar células por tipo (group_by) y tomar
  # muestras aleatorias (slice_sample) para el scatter plot.
  
  library(dplyr)

  # patchwork: permite combinar múltiples objetos ggplot en una
  # sola figura usando operadores | (lado a lado) y / (uno encima
  # de otro). Se usa para componer los paneles del PDF final.
  
  library(patchwork)

  # Matrix: soporte para matrices dispersas (sparse matrices) en
  # formato dgCMatrix, que es el formato de almacenamiento de las
  # matrices de conteos en los objetos SCE. Necesaria para nnzero()
  # (cuenta elementos no-cero) y operaciones eficientes sobre la
  # matriz de expresión sin convertirla a formato denso.
  
  library(Matrix)
})


# 1. DEFINICIÓN DE RUTAS 

# Ruta relativa al directorio de datos del dataset Grubman.
# Contiene los 7 archivos .rds (uno por tipo celular).
# Las rutas son relativas a la raíz del proyecto; el script debe
# ejecutarse desde alzheimer_tn/ para que funcionen correctamente.

DATA_DIR <- "data/organized_data/human/brain/snRNA/Grubman_31768052"

# Directorio de salida donde se guardarán los resultados:
# summary_cell_types.csv, sparsity_summary.csv y 01_qc_overview.pdf

OUT_DIR  <- "results/01_exploration/Grubman_31768052"

# Crea el directorio de salida si no existe. recursive=TRUE crea
# todos los directorios intermedios necesarios. showWarnings=FALSE
# suprime el aviso si el directorio ya existía.

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


# SECCIÓN 2: CARGA DE LOS ARCHIVOS .RDS 

# list.files() busca en DATA_DIR todos los archivos que terminan
# en .rds. full.names=TRUE devuelve rutas completas (necesario para
# que readRDS() funcione independientemente del directorio de trabajo).

rds_files <- list.files(DATA_DIR, pattern = "\\.rds$", full.names = TRUE)

# Extrae el nombre del tipo celular de cada archivo eliminando el
# sufijo "_data_arranged...rds" del nombre base del archivo.
# Ejemplo: "Astro_data_arranged_idUpdated_final.rds" → "Astro"

ct_names  <- sub("_data_arranged.*\\.rds$", "", basename(rds_files))

# Carga todos los archivos .rds en una lista nombrada.
# lapply() itera sobre los índices, imprime el progreso con cat()
# y carga cada objeto con readRDS().
# setNames() asigna los nombres de tipos celulares a la lista.
# El resultado es: sce_list$Astro, sce_list$MG, etc.

sce_list  <- setNames(
  lapply(seq_along(rds_files), function(i) {
    cat(sprintf("  [%d/%d] %s\n", i, length(rds_files), ct_names[i]))
    readRDS(rds_files[i])
  }),
  ct_names
)

# SECCIÓN 3: FUNCIÓN DE EXPLORACIÓN POR TIPO CELULAR 
# explore_sce() recibe el nombre del tipo celular (ct) y su objeto
# SCE (sce), y realiza un análisis completo imprimiendo resultados
# en consola. Devuelve una lista con los datos calculados para uso
# posterior en la tabla resumen y los plots.

explore_sce <- function(ct, sce) {

  # 3.1 CLASE Y DIMENSIONES 
  # class() devuelve la clase S4 del objeto (SingleCellExperiment).
  # nrow() = número de genes (filas de la matriz de expresión).
  # ncol() = número de células (columnas de la matriz).
  cat(sprintf("\nClase   : %s\n", paste(class(sce), collapse = " / ")))
  cat(sprintf("Dims    : %d genes  x  %d células\n", nrow(sce), ncol(sce)))

  # 3.2 ASSAYS DISPONIBLES 
  # assayNames() lista las capas de expresión almacenadas en el
  # objeto. En este dataset hay una sola: "counts" (conteos crudos
  # sin normalizar, enteros no-negativos). Otros datasets pueden
  # tener además "logcounts" (log-normalizado) o "normcounts".
  
  an <- assayNames(sce)
  cat(sprintf("Assays  : %s\n", paste(an, collapse = ", ")))

  # 3.3 REDUCCIONES DIMENSIONALES 
  # reducedDimNames() lista las reducciones dimensionales
  # pre-calculadas almacenadas en el objeto (PCA, UMAP, tSNE...).
  # En este dataset están precalculadas las coordenadas UMAP de
  # alineamiento (UMAP_alignment_1 y UMAP_alignment_2).
  
  rd <- reducedDimNames(sce)
  cat(sprintf("ReducDim: %s\n", if (length(rd)) paste(rd, collapse=", ") else "ninguna"))

  # 3.4 colData: METADATOS DE CÉLULAS 
  # colData() extrae el DataFrame con una fila por célula y una
  # columna por variable de metadatos. as.data.frame() lo convierte
  # a formato estándar de R para facilitar su manipulación.
  # En este dataset hay 57 columnas que incluyen: condición AD/Control
  # (disease.state.ch1), ID de donante (individuals.ch1), batch,
  # subclusters, métricas QC pre-calculadas y coordenadas UMAP.
  
  cd <- as.data.frame(colData(sce))
  cat(sprintf("\n--- colData (%d columnas, %d células) ---\n", ncol(cd), nrow(cd)))

  if (ncol(cd) > 0) {
    for (col in names(cd)) {
      vals <- cd[[col]]

      # Para columnas numéricas: imprime mediana, rango y nº de NAs.
      # La mediana es más robusta que la media ante valores extremos.
    
        if (is.numeric(vals)) {
        cat(sprintf("  %-25s [numeric]  median=%.1f  range=[%.0f, %.0f]  NAs=%d\n",
                    col, median(vals, na.rm=TRUE),
                    min(vals, na.rm=TRUE), max(vals, na.rm=TRUE),
                    sum(is.na(vals))))
      } else {
        
        # Para columnas categóricas: cuenta cuántas células hay por nivel.
        # useNA="ifany" incluye NAs en el recuento si los hay.
        # Si hay ≤10 niveles únicos, los muestra todos (ej: AD=142 Control=307).
        # Si hay >10 niveles (ej: IDs de donante), solo muestra el total.
        
        tbl <- table(vals, useNA="ifany")
        n_unique <- length(tbl)
        if (n_unique <= 10) {
          cat(sprintf("  %-25s [factor/chr, %d niveles]: %s\n",
                      col, n_unique,
                      paste(names(tbl), tbl, sep="=", collapse="  ")))
        } else {
          cat(sprintf("  %-25s [factor/chr, %d niveles únicos]\n", col, n_unique))
        }
      }
    }
  } else {
    cat("  (colData vacío)\n")
  }

  # 3.5 rowData: METADATOS DE GENES 
  # rowData() extrae el DataFrame con una fila por gen y una columna
  # por variable de anotación. En este dataset los genes se
  # identifican solo por su ID ENSEMBL (ENSG...) en los rownames,
  # sin columnas adicionales de anotación (rowData vacío).
 
  rwd <- as.data.frame(rowData(sce))
  cat(sprintf("\n--- rowData (%d columnas) ---\n", ncol(rwd)))
  if (ncol(rwd) > 0) {
    cat(sprintf("  Columnas: %s\n", paste(names(rwd), collapse=", ")))
    cat("  Primeras 5 filas:\n")
    print(head(rwd, 5))
  } else {
    # Si rowData está vacío, mostramos directamente los rownames
    # (IDs ENSEMBL) para confirmar el formato de los identificadores.
    cat("  Solo IDs de fila (ENSG...)\n")
    cat(sprintf("  Primeros genes: %s\n",
                paste(head(rownames(sce), 6), collapse=", ")))
  }

  # 3.6 metadata(): METADATOS GLOBALES DEL OBJETO 
  # metadata() devuelve una lista con información sobre el objeto
  # completo (no por célula ni por gen): versiones de software
  # usadas, parámetros de procesamiento, etc.
  
  meta <- metadata(sce)
  if (length(meta) > 0) {
    cat(sprintf("\n--- metadata() [%d entradas] ---\n", length(meta)))
    cat(sprintf("  Claves: %s\n", paste(names(meta), collapse=", ")))
  }

  # 3.7 MÉTRICAS DE CONTROL DE CALIDAD (QC)
  
  cat("\n--- QC metrics ---\n")

  # Extrae la primera (y única) capa de expresión: la matriz de
  # conteos crudos. Es una matriz sparse (dgCMatrix) de dimensiones
  # n_genes × n_células con valores enteros no-negativos.
  
  count_mat <- assay(sce, an[1])

  # Detecta genes mitocondriales buscando el prefijo "MT-" en los
  # nombres de fila. Los genes mitocondriales (MT-CO1, MT-ND1, etc.)
  # se usan como indicador de calidad: un alto % de reads
  # mitocondriales indica daño celular o muerte. En snRNA-seq los
  # valores son generalmente bajos (<5%) porque los núcleos
  # contienen poco RNA mitocondrial (las mitocondrias están en el
  # citoplasma). Si los rownames son IDs ENSEMBL, busca también
  # en la columna "Symbol" de rowData como fallback.
  
  rn       <- rownames(sce)
  is_mito  <- grepl("^MT-", rn, ignore.case = FALSE)
  if (!any(is_mito) && "Symbol" %in% names(rowData(sce)))
    is_mito <- grepl("^MT-", rowData(sce)$Symbol)
  cat(sprintf("  Genes mitocondriales detectados: %d\n", sum(is_mito)))

  # perCellQCMetrics() de scuttle calcula para cada célula:
  #   - sum: total de UMIs (profundidad de secuenciación)
  #   - detected: número de genes con al menos 1 count (complejidad)
  #   - subsets_MT_percent: porcentaje de UMIs en genes mitocondriales
  # Si falla (ej: versión incompatible), el tryCatch calcula
  # manualmente sum y detected con colSums sobre la matriz sparse.
  
  qc <- tryCatch(
    as.data.frame(perCellQCMetrics(
      sce, subsets = list(MT = which(is_mito))
    )),
    error = function(e) {
      data.frame(
        sum      = colSums(count_mat),
        detected = colSums(count_mat > 0),
        subsets_MT_percent = rep(0, ncol(count_mat))
      )
    }
  )

  # Imprime estadísticas de distribución de UMIs:
  # mediana (valor central), percentil 5 (cola baja), percentil 95
  # (cola alta) y máximo. El formato %7.0f reserva 7 caracteres sin
  # decimales para alinear columnas correctamente en la consola.
 
  cat(sprintf("  UMI      — mediana: %7.0f  p5: %6.0f  p95: %7.0f  max: %7.0f\n",
              median(qc$sum), quantile(qc$sum, .05),
              quantile(qc$sum, .95), max(qc$sum)))

  # Misma distribución para el número de genes detectados por célula.
  
  cat(sprintf("  Genes    — mediana: %7.0f  p5: %6.0f  p95: %7.0f  max: %7.0f\n",
              median(qc$detected), quantile(qc$detected, .05),
              quantile(qc$detected, .95), max(qc$detected)))

  # Mediana y percentil 95 del porcentaje mitocondrial.
  # %% en sprintf produce un literal %. na.rm=TRUE ignora NAs.
  
  if ("subsets_MT_percent" %in% names(qc)) {
    cat(sprintf("  %%MT      — mediana: %7.2f  p95: %.2f\n",
                median(qc$subsets_MT_percent, na.rm=TRUE),
                quantile(qc$subsets_MT_percent, .95, na.rm=TRUE)))
  }

  # Cuenta células con UMI < 500: posibles empty droplets (gotas
  # sin célula que capturan RNA ambiental) o células muertas.
  
  low_umi  <- sum(qc$sum < 500)

  # Cuenta outliers extremos de UMI alto: células cuyo recuento
  # supera 3x el percentil 99, que podrían ser doublets (dos células
  # capturadas en la misma gota).
  
  high_umi <- sum(qc$sum > quantile(qc$sum, .99) * 3)
  cat(sprintf("  Células con UMI < 500  : %d (%.1f%%)\n",
              low_umi,  100*low_umi/nrow(qc)))
  cat(sprintf("  Células outlier UMI alto: %d (%.1f%%)\n",
              high_umi, 100*high_umi/nrow(qc)))

  # Sparsity: proporción de ceros en la matriz de expresión.
  # nnzero() cuenta eficientemente los elementos no-cero en la
  # matriz sparse sin convertirla a formato denso (que consumiría
  # mucha memoria). 1 - (no_ceros / total_entradas) = fracción ceros.
  # En snRNA-seq valores del 90-95% son normales y esperados.
  
  sparsity <- 1 - nnzero(count_mat) / (nrow(count_mat) * ncol(count_mat))
  cat(sprintf("  Sparsity de la matriz  : %.2f%%\n", 100 * sparsity))

  # Genes con 0 conteos en TODAS las células del tipo celular.
  # rowSums() suma por filas (genes); == 0 identifica genes
  # completamente ausentes en este tipo celular. Son candidatos
  # a filtrar antes del análisis de ruido con RaceID/VarID2.
  
  zero_genes <- sum(rowSums(count_mat) == 0)
  cat(sprintf("  Genes con 0 counts     : %d / %d (%.1f%%)\n",
              zero_genes, nrow(count_mat),
              100 * zero_genes / nrow(count_mat)))

  # Devuelve una lista con todos los datos calculados para este
  # tipo celular. Esta lista se usará en la sección 5 para construir
  # la tabla resumen y en la sección 6 para los plots.
  
  list(ct = ct, n_cells = ncol(sce), n_genes = nrow(sce),
       assays = an, rd_cols = names(rwd), cd_cols = names(cd),
       qc = qc, is_mito = is_mito)
}


# 4. APLICAR LA FUNCIÓN A TODOS LOS TIPOS CELULARES ────
# lapply() aplica explore_sce() a cada elemento de sce_list.
# names(sce_list) proporciona tanto el nombre (ct) como el objeto
# (sce_list[[ct]]). El resultado es una lista de listas con los
# resultados de cada tipo celular.

results <- lapply(names(sce_list), function(ct) explore_sce(ct, sce_list[[ct]]))

# Asigna nombres a la lista de resultados para acceder por tipo
# celular: results$Astro, results$MG, etc.

names(results) <- names(sce_list)


# 5. TABLA RESUMEN GLOBAL
# Construye un data.frame con una fila por tipo celular y columnas
# con las métricas más importantes extraídas de results.
# sapply() aplica una función a cada elemento y simplifica el
# resultado a un vector (a diferencia de lapply que devuelve lista).

summary_df <- data.frame(
  # Nombre del tipo celular
  CellType      = names(results),
  # Número de células
  N_cells       = sapply(results, `[[`, "n_cells"),
  # Número de genes en la matriz
  N_genes       = sapply(results, `[[`, "n_genes"),
  # Nombre(s) de los assays disponibles, concatenados con "/"
  Assays        = sapply(results, function(x) paste(x$assays, collapse="/")),
  # Nombres de todas las columnas de colData, separados por coma
  ColData_cols  = sapply(results, function(x)
                    if (length(x$cd_cols)) paste(x$cd_cols, collapse=", ") else "(vacío)"),
  # Mediana de UMIs por célula, redondeada al entero más cercano
  Median_UMI    = sapply(results, function(x) round(median(x$qc$sum))),
  # Mediana de genes detectados por célula
  Median_genes  = sapply(results, function(x) round(median(x$qc$detected)))
)

# Imprime en consola solo las columnas más relevantes (sin ColData_cols
# que sería muy larga) para una revisión rápida de los resultados.
print(summary_df[, c("CellType","N_cells","N_genes","Median_UMI","Median_genes","Assays")])

# Guarda la tabla completa (incluyendo ColData_cols) en CSV.
# row.names=FALSE evita que R añada una columna numérica de índices.
write.csv(summary_df, file.path(OUT_DIR, "summary_cell_types.csv"), row.names = FALSE)


# 6. VISUALIZACIONES 

# Construye un data.frame "largo" con las métricas QC de todas las
# células de todos los tipos celulares en un único data.frame.
# do.call(rbind, ...) une por filas las 7 tablas QC individuales.
# Se añade la columna "ct" con el nombre del tipo celular para poder
# separar por color/faceta en los plots.

qc_long <- do.call(rbind, lapply(names(results), function(ct) {
  qc       <- results[[ct]]$qc
  qc$ct    <- ct
  qc
}))

# Paleta de colores manual con un color distinto por tipo celular.
# Los nombres deben coincidir exactamente con los valores de la
# columna "ct" para que scale_fill/color_manual funcione.
pal <- c(Astro="#4E9AF1", Endo="#F4A460", ExN="#6BCB77",
         InN="#FF6B6B", MG="#C77DFF", Oligo="#FFD93D", OPC="#00C9A7")


# PLOT 1: NÚMERO DE CÉLULAS POR TIPO 
# Barplot mostrando cuántas células hay en cada tipo celular.
# reorder(CellType, -N_cells) ordena las barras de mayor a menor.
# geom_text() añade el valor numérico encima de cada barra.
# expansion(mult=c(0, 0.12)) añade espacio arriba para el texto.

p1 <- ggplot(summary_df, aes(x = reorder(CellType, -N_cells),
                              y = N_cells, fill = CellType)) +
  geom_bar(stat="identity", width=0.7) +
  geom_text(aes(label = format(N_cells, big.mark=",")),
            vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = pal) +
  scale_y_continuous(expand = expansion(mult=c(0, 0.12))) +
  theme_minimal(base_size = 12) +
  labs(title = "Número de células por tipo celular",
       subtitle = "Grubman_31768052", x = NULL, y = "N células") +
  theme(legend.position="none",
        axis.text.x = element_text(angle=35, hjust=1))


# PLOT 2: DISTRIBUCIÓN DE UMIs (VIOLIN)
# Violin plot de log10(UMI + 1) por tipo celular.
# Se usa log10 porque los UMIs tienen distribución log-normal.
# El +1 evita log10(0) = -Inf en células con 0 UMIs.
# geom_violin() muestra la distribución de densidad completa.
# geom_boxplot(width=0.08) superpone un boxplot estrecho para ver
# mediana y cuartiles con más precisión. outlier.shape=NA oculta
# los outliers del boxplot (ya visibles en el violin).

p2 <- ggplot(qc_long, aes(x = ct, y = log10(sum + 1), fill = ct)) +
  geom_violin(scale="width", alpha=0.85, trim=FALSE, linewidth=0.3) +
  geom_boxplot(width=0.08, fill="white", outlier.shape=NA, linewidth=0.4) +
  scale_fill_manual(values = pal) +
  theme_minimal(base_size = 12) +
  labs(title = "Distribución de UMIs",
       x = NULL, y = "log10(UMI + 1)") +
  theme(legend.position="none",
        axis.text.x = element_text(angle=35, hjust=1))


# PLOT 3: DISTRIBUCIÓN DE GENES DETECTADOS (VIOLIN)
# Mismo formato que p2 pero para el número de genes detectados
# por célula. Complementa el análisis de UMIs: baja detección de
# genes con muchos UMIs puede indicar amplificación sesgada.

p3 <- ggplot(qc_long, aes(x = ct, y = log10(detected + 1), fill = ct)) +
  geom_violin(scale="width", alpha=0.85, trim=FALSE, linewidth=0.3) +
  geom_boxplot(width=0.08, fill="white", outlier.shape=NA, linewidth=0.4) +
  scale_fill_manual(values = pal) +
  theme_minimal(base_size = 12) +
  labs(title = "Genes detectados por célula",
       x = NULL, y = "log10(genes + 1)") +
  theme(legend.position="none",
        axis.text.x = element_text(angle=35, hjust=1))


# PLOT 4: UMIs vs GENES DETECTADOS (SCATTER) 
# set.seed() fija la semilla aleatoria para reproducibilidad.
# slice_sample() toma una muestra aleatoria de 300 células por tipo
# (o todas si hay menos de 300) para no sobrecargar el plot.
# El scatter plot muestra la correlación entre profundidad de
# secuenciación (UMI) y complejidad transcriptómica (genes).
# Una correlación lineal en escala log-log es la esperada y normal.
# Puntos fuera de la tendencia pueden indicar células anómalas.
# facet_wrap() crea un panel separado por tipo celular para ver
# cada distribución sin solapamiento.

set.seed(42)
qc_sub <- qc_long %>%
  group_by(ct) %>%
  slice_sample(n = min(300, n())) %>%
  ungroup()

p4 <- ggplot(qc_sub, aes(x = log10(sum+1), y = log10(detected+1),
                          color = ct)) +
  geom_point(alpha=0.35, size=0.7) +
  scale_color_manual(values = pal) +
  facet_wrap(~ct, ncol=4) +
  theme_minimal(base_size = 10) +
  labs(title = "UMIs vs genes detectados",
       x = "log10(UMI + 1)", y = "log10(genes + 1)") +
  theme(legend.position="none",
        strip.text = element_text(face="bold"))


# PLOT 5: PORCENTAJE DE READS MITOCONDRIALES
# Si la columna subsets_MT_percent fue calculada por perCellQCMetrics,
# genera un violin plot del %MT por tipo celular.
# geom_hline() dibuja una línea discontinua en y=20 como referencia
# del umbral típico de QC (células con >20% MT se filtran en la
# mayoría de pipelines como potencialmente dañadas).
# En snRNA-seq los valores son generalmente bajos (<5%) por la
# ausencia de RNA mitocondrial en los núcleos aislados.

if ("subsets_MT_percent" %in% names(qc_long)) {
  p5 <- ggplot(qc_long, aes(x = ct, y = subsets_MT_percent, fill = ct)) +
    geom_violin(scale="width", alpha=0.85, trim=FALSE, linewidth=0.3) +
    geom_boxplot(width=0.08, fill="white", outlier.shape=NA, linewidth=0.4) +
    geom_hline(yintercept=20, linetype="dashed", color="red", linewidth=0.5) +
    scale_fill_manual(values = pal) +
    theme_minimal(base_size = 12) +
    labs(title = "% reads mitocondriales",
         subtitle = "línea roja = umbral 20%",
         x = NULL, y = "% MT") +
    theme(legend.position="none",
          axis.text.x = element_text(angle=35, hjust=1))
} else {
  # Si no se detectaron genes MT- (rownames son IDs ENSEMBL sin
  # prefijo MT-), genera un plot vacío con mensaje informativo.
  p5 <- ggplot() +
    annotate("text", x=0.5, y=0.5, label="% MT no disponible\n(no se detectaron genes MT-)") +
    theme_void()
}


# PLOT 6: SPARSITY DE LA MATRIZ DE EXPRESIÓN 
# Calcula la sparsity (% de ceros) y los genes sin ningún count
# para cada tipo celular directamente desde los objetos SCE.
# Se accede a la matriz con assay(sce, 1) (primer assay = counts).
# nnzero() cuenta eficientemente los no-ceros sin descomprimir.
# rowSums(m)==0 identifica genes con cero conteos en TODAS las
# células de ese tipo: son genes sin expresión detectable en ese
# tipo celular y candidatos a ser filtrados en análisis posteriores.

sparsity_df <- data.frame(
  CellType = names(sce_list),
  Sparsity = sapply(names(sce_list), function(ct) {
    m <- assay(sce_list[[ct]], 1)
    round(100 * (1 - nnzero(m) / (nrow(m) * ncol(m))), 2)
  }),
  ZeroGenes_pct = sapply(names(sce_list), function(ct) {
    m <- assay(sce_list[[ct]], 1)
    round(100 * sum(rowSums(m) == 0) / nrow(m), 2)
  })
)

# Barplot de sparsity con etiquetas de porcentaje encima de cada
# barra. scale_y_continuous(limits=c(0,105)) asegura que las
# etiquetas no se corten en el borde superior del gráfico.

p6 <- ggplot(sparsity_df, aes(x = CellType, y = Sparsity, fill = CellType)) +
  geom_bar(stat="identity", width=0.7) +
  geom_text(aes(label = paste0(Sparsity, "%")),
            vjust=-0.4, size=3.2) +
  scale_fill_manual(values = pal) +
  scale_y_continuous(limits=c(0,105), expand=c(0,0)) +
  theme_minimal(base_size = 12) +
  labs(title = "Sparsity de la matriz de conteos",
       x = NULL, y = "% ceros") +
  theme(legend.position="none",
        axis.text.x = element_text(angle=35, hjust=1))


# PLOT 7: PROPORCIÓN AD vs CONTROL POR TIPO CELULAR 
# Busca en el colData de cada tipo celular una columna que contenga
# información de condición (AD/Control) buscando palabras clave
# con grep() en los nombres de columna (ignore.case para robustez).
# Construye una tabla de frecuencias y la almacena en condition_plots.

condition_plots <- list()
for (ct in names(sce_list)) {
  cd   <- as.data.frame(colData(sce_list[[ct]]))
  cols <- grep("condition|diagnosis|disease|group|status",
               names(cd), ignore.case=TRUE, value=TRUE)
  if (length(cols) > 0) {
    # Usa la primera columna que coincide con las palabras clave.
    # table() cuenta cuántas células hay por nivel (AD, Control).
    col <- cols[1]
    tbl <- as.data.frame(table(Condition = cd[[col]]))
    condition_plots[[ct]] <- tbl
  }
}

if (length(condition_plots) > 0) {
  # Combina las tablas de todos los tipos celulares en un data.frame
  # largo y genera un barplot apilado normalizado (position="fill")
  # que muestra proporciones en lugar de conteos absolutos.
  # scales::percent formatea el eje Y como porcentajes (0%, 25%...).
  cond_df <- do.call(rbind, lapply(names(condition_plots), function(ct) {
    df     <- condition_plots[[ct]]
    df$ct  <- ct
    df
  }))
  p7 <- ggplot(cond_df, aes(x = ct, y = Freq, fill = Condition)) +
    geom_bar(stat="identity", position="fill", width=0.7) +
    scale_y_continuous(labels = scales::percent) +
    theme_minimal(base_size=12) +
    labs(title = "Proporción AD vs Control por tipo celular",
         x = NULL, y = "Proporción") +
    theme(axis.text.x = element_text(angle=35, hjust=1))
} else {
  p7 <- ggplot() +
    annotate("text", x=0.5, y=0.5,
             label="Columna de condición (AD/Control)\nno encontrada en colData") +
    theme_void()
}


# 7. EXPORTAR PLOTS A PDF 
# pdf() abre un dispositivo gráfico que redirige todos los plots
# al archivo especificado hasta que se llame dev.off().
# width=14, height=9 define el tamaño en pulgadas de cada página.
# Se generan 4 páginas:
#   Pág 1: barplot células (p1) + barplot sparsity (p6) lado a lado
#   Pág 2: violin UMIs (p2) + violin genes (p3) + violin %MT (p5)
#   Pág 3: scatter UMI vs genes por tipo celular (p4) en facetas
#   Pág 4: barplot apilado AD/Control (p7)
# El operador | de patchwork coloca plots lado a lado.
# plot_annotation() añade un título global a la composición.

pdf(file.path(OUT_DIR, "01_qc_overview.pdf"), width=14, height=9)
print((p1 | p6) + plot_annotation(title="Grubman_31768052 — Overview"))
print((p2 | p3 | p5) + plot_annotation(title="Distribuciones QC por tipo celular"))
print(p4 + plot_annotation(title="UMIs vs Genes detectados"))
print(p7 + plot_annotation(title="Condición AD/Control"))

# dev.off() cierra el dispositivo gráfico y finaliza la escritura
# del PDF. Es obligatorio llamarlo; sin él el PDF queda incompleto.

dev.off()


# 8. GUARDAR TABLAS RESUMEN EN CSV 
# Guarda sparsity_summary.csv con los valores de sparsity y genes
# con cero conteos por tipo celular. row.names=FALSE evita la
# columna de índices numéricos que R añade por defecto.

write.csv(sparsity_df, file.path(OUT_DIR, "sparsity_summary.csv"), row.names=FALSE)

# Imprime en consola la ruta de todos los archivos generados.

cat(sprintf(" Resultados guardados en: %s\n", OUT_DIR))
cat(sprintf("   · summary_cell_types.csv\n"))
cat(sprintf("   · sparsity_summary.csv\n"))
cat(sprintf("   · 01_qc_overview.pdf\n"))
