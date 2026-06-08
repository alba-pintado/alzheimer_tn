# ================================================================
# convert_qs_to_rds.R — Plan B: convertir .qs a .rds sin conda
# ================================================================
# Si el entorno conda sigue fallando, ejecuta esto con el R del 
# SISTEMA (no con el de conda):
#
#   /usr/bin/Rscript scripts/convert_qs_to_rds.R
#
# El R del sistema (/usr/bin/R) instala qs de forma independiente
# y no tiene el conflicto de TBB del entorno conda.
# ================================================================

cat("R del sistema:", R.home(), "\n")
cat("Versión:", R.version$major, ".", R.version$minor, "\n\n")

# Instalar qs si no está (solo necesario la primera vez)
if (!requireNamespace("qs", quietly = TRUE)) {
  install.packages("qs", repos = "https://cran.r-project.org")
}

library(qs)
cat("qs versión:", as.character(packageVersion("qs")), "\n\n")

# ── Configuración ──────────────────────────────────────────────
DATASET_DIR <- "data/organized_data/human/brain/snRNA/Grubman_31768052"

# Encontrar todos los .qs del dataset
qs_files <- list.files(DATASET_DIR, pattern = "\\.qs$", full.names = TRUE)
cat(sprintf("Archivos .qs encontrados: %d\n\n", length(qs_files)))

# ── Conversión ─────────────────────────────────────────────────
for (f in qs_files) {
  rds_path <- sub("\\.qs$", ".rds", f)
  
  if (file.exists(rds_path)) {
    cat(sprintf("  YA EXISTE: %s, saltando\n", basename(rds_path)))
    next
  }
  
  cat(sprintf("  Convirtiendo: %s ...", basename(f)))
  tryCatch({
    obj <- qread(f)
    saveRDS(obj, rds_path, compress = FALSE)  # compress=FALSE es más rápido
    cat(sprintf(" OK (%s → %s)\n", 
                format(file.size(f), big.mark="."), 
                format(file.size(rds_path), big.mark=".")))
  }, error = function(e) {
    cat(sprintf(" ERROR: %s\n", e$message))
  })
}

cat("\n¡Conversión completada!\n")
cat("Ahora puedes usar readRDS() en tu entorno conda sin necesitar qs.\n")

# ── Verificación ───────────────────────────────────────────────
rds_files <- list.files(DATASET_DIR, pattern = "\\.rds$", full.names = TRUE)
cat(sprintf("\nArchivos .rds generados: %d\n", length(rds_files)))
for (f in rds_files) {
  cat(sprintf("  %s (%s bytes)\n", basename(f), format(file.size(f), big.mark=".")))
}
