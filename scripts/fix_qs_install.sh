#!/usr/bin/env bash
# ================================================================
# fix_qs_install.sh
# Soluciona el error stringfish TBB en un entorno conda de R
# Uso: bash fix_qs_install.sh [nombre_del_entorno]
# Ejemplo: bash fix_qs_install.sh alzheimer_R
# ================================================================

ENV_NAME="${1:-alzheimer_R}"

echo "======================================================"
echo "  Reparando qs/stringfish en el entorno: $ENV_NAME"
echo "======================================================"

# Activar entorno
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

echo ""
echo "[1/5] Versión actual de TBB en el entorno:"
conda list tbb

echo ""
echo "[2/5] Eliminando paquetes rotos..."
R -e "
pkgs_to_remove <- c('qs', 'qs2', 'stringfish')
installed <- rownames(installed.packages())
to_remove <- intersect(pkgs_to_remove, installed)
if (length(to_remove) > 0) {
  remove.packages(to_remove)
  cat('Eliminados:', paste(to_remove, collapse=', '), '\n')
} else {
  cat('No había paquetes que eliminar\n')
}
"

echo ""
echo "[3/5] Reinstalando TBB compatible desde conda-forge..."
conda install -c conda-forge r-stringfish -y

echo ""
echo "[4/5] Verificando que stringfish carga correctamente..."
R -e "library(stringfish); cat('stringfish OK, version:', as.character(packageVersion('stringfish')), '\n')"

STRINGFISH_OK=$?

if [ $STRINGFISH_OK -eq 0 ]; then
    echo ""
    echo "[5/5] Instalando qs desde conda-forge (binario compatible)..."
    conda install -c conda-forge r-qs -y

    echo ""
    echo "Verificando qs..."
    R -e "
    library(qs)
    cat('qs OK, version:', as.character(packageVersion('qs')), '\n')
    # Test de escritura/lectura
    tmp <- tempfile(fileext='.qs')
    qsave(list(a=1:10, b='test'), tmp)
    obj <- qread(tmp)
    cat('Test lectura/escritura: OK\n')
    cat('Objeto leido:', class(obj), 'con', length(obj), 'elementos\n')
    file.remove(tmp)
    "
else
    echo ""
    echo "stringfish sigue fallando. Intentando compilación desde fuente..."
    echo "Esto puede tardar 5-10 minutos."
    
    # Opción 2: compilar desde fuente desactivando TBB
    R -e "
    Sys.setenv(PKG_LIBS = '-L$(conda info --base)/envs/$ENV_NAME/lib')
    Sys.setenv(PKG_CPPFLAGS = '-I$(conda info --base)/envs/$ENV_NAME/include')
    install.packages('stringfish', type='source', repos='https://cran.r-project.org')
    "
    
    R -e "
    install.packages('qs', type='source', repos='https://cran.r-project.org')
    library(qs)
    cat('qs (desde fuente) OK, version:', as.character(packageVersion('qs')), '\n')
    "
fi

echo ""
echo "======================================================"
echo "  Diagnóstico final"
echo "======================================================"
R -e "
cat('R version:', R.version\$major, '.', R.version\$minor, '\n')
for (pkg in c('stringfish', 'qs')) {
  status <- tryCatch({
    library(pkg, character.only=TRUE)
    paste0('OK (', packageVersion(pkg), ')')
  }, error = function(e) paste0('FALLO: ', e\$message))
  cat(sprintf('  %-15s %s\n', pkg, status))
}
"
