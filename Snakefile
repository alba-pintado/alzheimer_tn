import os

# -----------------------------
# (1) RUTAS BASE DEL PROYECTO
# -----------------------------
DATA = "Data"
OUT = "Output"

VARID2_DIR = "Scripts/VarID2/VarID2-main"
DECIBELL_DIR = "Scripts/decibell/decibel-main/module"


# ---------------------------------
# (2) OBJETIVO FINAL DEL PIPELINE
# --------------------------------
rule all:
    input:
        f"{OUT}/varid2/done.txt",
        f"{OUT}/scallop/done.txt",
        f"{OUT}/decibell/done.txt"


# -------------------------
# VARID2 (R - FIGURA 1)
# -------------------------
rule varid2_fig1:
    input:
        script=f"{VARID2_DIR}/VarID2_Figure1.R",
        data=f"{DATA}/NPH_DE_res/NPH_cellType_DEres.qs"
    output:
        f"{OUT}/varid2/done.txt"
    shell:
        """
        mkdir -p {OUT}/varid2
        Rscript {input.script}
        touch {output}
        """


# -------------------------
# SCALLOP (PLACEHOLDER)
# -------------------------
rule scallop:
    output:
        f"{OUT}/scallop/done.txt"
    shell:
        """
        mkdir -p {OUT}/scallop
        echo "Scallop pipeline not implemented yet" > {output}
        """


# -------------------------
# DECIBELL (PYTHON MODULE)
# -------------------------
rule decibell:
    output:
        f"{OUT}/decibell/done.txt"
    shell:
        """
        mkdir -p {OUT}/decibell
        python -c "import sys; sys.path.append('{DECIBELL_DIR}'); import decibel"
        echo "Decibell loaded successfully" > {output}
        """
```
