# Lau et al. 2020: analisis en Python (scanpy)

# TFM "Transcriptional Noise in Alzheimer's Disease Progression" (CBMSO)

# Version en Python de 02_lau_analisis.R. Parte del mismo objeto que R,
# exportado a .h5ad con 00_exportar_h5ad.R + 00b_construir_h5ad.py, de modo que
# los nucleos son los mismos y los resultados se pueden comparar con los de R.

# Version del 27-10-2026: objeto definitivo (2,5 MADs, sin nucleos NA, 21
# donantes: 12 AD y 9 control; 137.303 nucleos x 29.447 genes). Como en R se
# mantiene MNN, si el .h5ad trae el espacio MNN de R (obsm["X_mnn_R"], lo
# exporta 00_exportar_h5ad.R desde el 27-10-2026) el UMAP y el clustering se
# hacen sobre el; si no lo trae, se usa Harmony calculado aqui, como en la
# version del 20-10-2026. Se puede forzar con CORRECCION.

# Equivalencias R -> Python:
#   genes variables  modelGeneVariances      -> highly_variable_genes
#   PCA              runPca(block = donante) -> pca (sin bloqueo; ver (2))
#   correccion       MNN                     -> MNN de R o Harmony
#   clustering       SNN + clusterGraph      -> kNN + Leiden
#   marcadores       scoreMarkers (AUC)      -> rank_genes_groups (Wilcoxon)
#   pseudobulk       edgeR                   -> PyDESeq2
#   abundancia       edgeR                   -> GLM quasi-binomial

# Ejecutar desde la carpeta raiz del proyecto (-u para que el log se escriba
# en tiempo real):
#   conda activate alzheimer_snRNA
#   nohup python -u scripts/SCANPY/lau_2020/01_analisis_lau.py > scripts/SCANPY/lau_2020/log_python.txt 2>&1 &

# Fuentes: scanpy (Wolf et al. 2018, https://scanpy.readthedocs.io); Harmony
# (Korsunsky et al. 2019); Leiden (Traag et al. 2019); PyDESeq2 (Muzellec et al.
# 2023); pseudobulk con el donante como replica (Squair et al. 2021);
# single-cell best practices (Heumos et al. 2023, https://www.sc-best-practices.org).


# (0) PREPARACION
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")                      # guardar figuras sin pantalla
import matplotlib.pyplot as plt
import scanpy as sc
from sklearn.metrics import adjusted_rand_score
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

SEMILLA = 1234
CORRECCION = None                          # None = automatico; o "mnn_R" / "harmony"
np.random.seed(SEMILLA)
sc.settings.verbosity = 1

proyecto = os.getcwd()                     # ejecutar desde la carpeta raiz del proyecto
entrada  = os.path.join(proyecto, "results/objects/lau_2020/lau_sce_analisis.h5ad")
outdir   = os.path.join(proyecto, "results/05_scanpy/Lau_32989152")
os.makedirs(outdir, exist_ok=True)
sc.settings.figdir = outdir

adata = sc.read_h5ad(entrada)
print(adata)                                # 137303 nucleos x 29447 genes

# (0.1) COMPROBACIONES DE LA IMPORTACION
# Mismos nucleos y donantes que en R, sin nucleos NA y con las cuentas crudas.
assert adata.obs_names.is_unique and adata.var_names.is_unique
assert "counts" in adata.layers
assert adata.obs["donor"].nunique() == 21
assert adata.obs["celltype"].notna().all()
print(pd.crosstab(adata.obs["donor"], adata.obs["condition"]))
if CORRECCION is None:
    CORRECCION = "mnn_R" if "X_mnn_R" in adata.obsm else "harmony"
if CORRECCION == "mnn_R":
    assert "X_mnn_R" in adata.obsm, "Falta X_mnn_R: vuelve a exportar con 00_exportar_h5ad.R (27-10-2026)"
print("Correccion usada para UMAP y clustering:", CORRECCION)

for c in ["donor", "condition", "celltype", "cluster"]:
    adata.obs[c] = adata.obs[c].astype(str).astype("category")
adata.obs["cluster_R"] = adata.obs["cluster"]          # clusters de R, para comparar
simbolo = adata.var["simbolo"].astype(str).where(adata.var["simbolo"].notna(),
                                                adata.var_names.to_series())

# (0.2) Nucleos por donante y tipo celular. Es la tabla que decide cuantos
#       nucleos se pueden submuestrear en el analisis de ruido (DECIBEL).
n_donante_tipo = pd.crosstab(adata.obs["donor"], adata.obs["celltype"])
print(n_donante_tipo)
print("Minimo de nucleos por donante en cada tipo:")
print(n_donante_tipo.min())
n_donante_tipo.to_csv(os.path.join(outdir, "00_nucleos_por_donante_y_tipo.csv"))


# (1) GENES VARIABLES
# Varianza frente a media sobre la matriz logaritmica, calculada dentro de cada
# donante y combinada (batch_key), y los 3000 genes mas variables, como en R.
# Fuente: ?sc.pp.highly_variable_genes; Heumos et al. 2023, "Feature selection".
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=3000,
                            batch_key="donor")
print(adata.var["highly_variable"].sum())               # 3000


# (2) PCA Y R² FRENTE A TIPO CELULAR Y DONANTE
# scanpy no tiene un PCA bloqueado por donante como runPca(block=) de scrapper,
# asi que el PCA de Python es el normal. Se calcula tambien el R² sobre el PCA
# bloqueado de R (X_pca_R) y sobre el espacio corregido que se use.
sc.tl.pca(adata, n_comps=25, mask_var="highly_variable", random_state=SEMILLA)
var_pct = adata.uns["pca"]["variance_ratio"] * 100
print(np.round(var_pct[:6], 2))

def r2_por_grupo(componentes, grupos):
    """R² de cada componente frente a un factor (igual que lm(pc ~ factor) en R)."""
    df = pd.DataFrame(componentes)
    medias = df.groupby(np.asarray(grupos)).transform("mean")
    ss_entre = ((medias - df.mean()) ** 2).sum()
    ss_total = ((df - df.mean()) ** 2).sum()
    return (ss_entre / ss_total).values

def tabla_r2(clave_obsm):
    pcs = adata.obsm[clave_obsm][:, :6]
    return pd.DataFrame({
        "pc": range(1, 7),
        "tipo": r2_por_grupo(pcs, adata.obs["celltype"].values),
        "donante": r2_por_grupo(pcs, adata.obs["donor"].values)}).round(3)

r2 = {"pca_py": tabla_r2("X_pca")}
if "X_pca_R" in adata.obsm:
    r2["pca_R_bloqueado"] = tabla_r2("X_pca_R")
for nombre, tabla in r2.items():
    print(nombre); print(tabla)

plt.figure(figsize=(5, 4))
plt.semilogy(range(1, 26), var_pct, "o-")
plt.xlabel("PC"); plt.ylabel("Varianza explicada (%)")
plt.tight_layout(); plt.savefig(os.path.join(outdir, "01_scree_plot.png"), dpi=150); plt.close()


# (3) ESPACIO CORREGIDO POR DONANTE
# Cada donante es solo AD o solo control: corregir por donante elimina tambien
# parte del efecto de la enfermedad. El espacio corregido se usa SOLO para UMAP
# y clustering; marcadores, pseudobulk, abundancia y ruido usan datos sin corregir.
# - "mnn_R": el espacio MNN de R (correctMnn sobre el PCA bloqueado).
# - "harmony": harmonypy sobre el PCA de Python. Se llama a harmonypy
#   directamente; segun la version devuelve la matriz en una u otra orientacion.
if CORRECCION == "mnn_R":
    clave_corr = "X_mnn_R"
else:
    import harmonypy
    ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, "donor",
                               random_state=SEMILLA, verbose=False)
    z = np.asarray(ho.Z_corr)
    adata.obsm["X_pca_harmony"] = z if z.shape[0] == adata.n_obs else z.T
    clave_corr = "X_pca_harmony"
r2["corregido_" + CORRECCION] = tabla_r2(clave_corr)
print(r2["corregido_" + CORRECCION])
pd.concat([t.add_suffix("_" + n) for n, t in r2.items()], axis=1) \
  .to_csv(os.path.join(outdir, "02_r2_pca.csv"), index=False)


# (4) UMAP
sc.pp.neighbors(adata, use_rep=clave_corr, n_neighbors=15, random_state=SEMILLA)
sc.tl.umap(adata, random_state=SEMILLA)

# Nucleos en orden aleatorio para que ningun donante ni condicion tape al resto.
orden = np.random.permutation(adata.n_obs)
ad_plot = adata[orden]
for col, nombre in [("celltype", "03_umap_por_celltype.png"),
                    ("donor", "04_umap_por_donante.png"),
                    ("condition", "05_umap_por_diagnostico.png")]:
    fig = sc.pl.umap(ad_plot, color=col, size=2, show=False, return_fig=True)
    fig.savefig(os.path.join(outdir, nombre), dpi=150, bbox_inches="tight"); plt.close(fig)

fig, ejes = plt.subplots(1, 2, figsize=(14, 6), sharex=True, sharey=True)
for eje, cond in zip(ejes, ["AD", "Ctrl"]):
    sc.pl.umap(ad_plot[ad_plot.obs["condition"] == cond], color="celltype",
               size=2, ax=eje, show=False, title=cond,
               legend_loc="right margin" if cond == "Ctrl" else None)
fig.savefig(os.path.join(outdir, "05b_umap_por_diagnostico_paneles.png"),
            dpi=150, bbox_inches="tight"); plt.close(fig)
del ad_plot


# (5) CLUSTERING (LEIDEN) Y VALIDACION
# Con "mnn_R" el grafo se construye sobre el mismo espacio que en R, asi que el
# ARI con los clusters de R compara solo el algoritmo (SNN + clusterGraph
# frente a kNN + Leiden).
sc.tl.leiden(adata, resolution=1.0, key_added="cluster_py",
             flavor="igraph", n_iterations=2, random_state=SEMILLA)
print(adata.obs["cluster_py"].value_counts().sort_index())

tab = pd.crosstab(adata.obs["cluster_py"], adata.obs["celltype"])
pureza = tab.max(axis=1) / tab.sum(axis=1)
tipo_cluster = tab.idxmax(axis=1)
print("Clusters con pureza >= 0,99:", (pureza >= 0.99).sum(), "de", len(pureza))
print(pureza.sort_values().head().round(3))

ari_tipo = adjusted_rand_score(adata.obs["celltype"], adata.obs["cluster_py"])
ari_R    = adjusted_rand_score(adata.obs["cluster_R"], adata.obs["cluster_py"])
print(f"ARI frente a la anotacion: {ari_tipo:.3f}")
print(f"ARI frente a los clusters de R: {ari_R:.3f}")

max_donante = pd.crosstab(adata.obs["cluster_py"], adata.obs["donor"], normalize="index").max(axis=1)
print("Clusters dominados (> 70 %) por un donante:", list(max_donante[max_donante > 0.7].index))
mito = adata.obs.groupby("cluster_py", observed=True)["subsets_mito_percent"].median()
print(mito.sort_values(ascending=False).head().round(2))

fig = sc.pl.umap(adata, color="cluster_py", legend_loc="on data", size=2,
                 show=False, return_fig=True)
fig.savefig(os.path.join(outdir, "06_umap_por_cluster.png"), dpi=150, bbox_inches="tight"); plt.close(fig)


# (6) MARCADORES
# Sobre la matriz logaritmica sin corregir (X). Wilcoxon, cada cluster frente al resto.
sc.tl.rank_genes_groups(adata, groupby="cluster_py", method="wilcoxon",
                        use_raw=False, pts=False)
resumen_marcadores = []
for k in adata.obs["cluster_py"].cat.categories:
    top = sc.get.rank_genes_groups_df(adata, group=k).head(10)["names"]
    resumen_marcadores.append({"cluster": k, "tipo": tipo_cluster.get(k),
                               "pureza": round(pureza.get(k, np.nan), 3),
                               "top10": ", ".join(simbolo[top].values)})
resumen_marcadores = pd.DataFrame(resumen_marcadores)
print(resumen_marcadores.to_string())
resumen_marcadores.to_csv(os.path.join(outdir, "07_marcadores_top10.csv"), index=False)


# (7) AD FRENTE A CONTROL: PSEUDOBULK CON PYDESEQ2
# Cuentas crudas sumadas por donante dentro de cada tipo celular: la replica es
# el donante, no el nucleo (Squair et al. 2021). Control como referencia.
def pseudobulk(ad, capa, grupo):
    """Suma las cuentas de cada grupo (donante). Devuelve donantes x genes."""
    g = pd.Categorical(ad.obs[grupo])
    indicadora = sp.csr_matrix((np.ones(ad.n_obs), (g.codes, np.arange(ad.n_obs))),
                               shape=(len(g.categories), ad.n_obs))
    suma = indicadora @ ad.layers[capa]
    suma = suma.toarray() if sp.issparse(suma) else np.asarray(suma)
    return pd.DataFrame(suma.astype(int), index=g.categories, columns=ad.var_names)

resumen_DE = []
for tipo in sorted(adata.obs["celltype"].unique()):
    sub = adata[adata.obs["celltype"] == tipo]
    cuentas = pseudobulk(sub, "counts", "donor")
    cuentas = cuentas.loc[cuentas.sum(axis=1) > 0]          # donantes sin nucleos de este tipo fuera
    meta = (sub.obs.groupby("donor", observed=True)["condition"].first()
            .reindex(cuentas.index).to_frame())
    meta["condition"] = pd.Categorical(meta["condition"].astype(str),
                                       categories=["Ctrl", "AD"])

    # Filtro equivalente a filterByExpr: al menos 10 cuentas en tantos donantes
    # como tiene el grupo mas pequeno.
    n_min = meta["condition"].value_counts().min()
    cuentas = cuentas.loc[:, (cuentas >= 10).sum(axis=0) >= n_min]

    dds = DeseqDataSet(counts=cuentas, metadata=meta, design="~condition",
                       quiet=True)
    dds.deseq2()
    res = DeseqStats(dds, contrast=["condition", "AD", "Ctrl"], quiet=True)
    res.summary()
    tabla = res.results_df.sort_values("padj")
    tabla["simbolo"] = simbolo.reindex(tabla.index).values
    tabla.to_csv(os.path.join(outdir, f"08_DE_{tipo}.csv"))
    resumen_DE.append({"tipo": tipo, "donantes": len(cuentas),
                       "genes": cuentas.shape[1],
                       "significativos": int((tabla["padj"] < 0.05).sum())})
resumen_DE = pd.DataFrame(resumen_DE)
print(resumen_DE)
resumen_DE.to_csv(os.path.join(outdir, "09_resumen_DE.csv"), index=False)


# (8) ABUNDANCIA DIFERENCIAL
# Para cada tipo celular, proporcion de nucleos por donante ~ condicion con un
# GLM quasi-binomial (pesos = nucleos del donante). Equivale a la pregunta del
# analisis de abundancia de R: ¿cambia la fraccion de cada tipo?
cuenta_tipo = n_donante_tipo
total = cuenta_tipo.sum(axis=1)
cond = adata.obs.groupby("donor", observed=True)["condition"].first().reindex(cuenta_tipo.index)
X = sm.add_constant((cond.astype(str) == "AD").astype(float).rename("AD"))
filas = []
for tipo in cuenta_tipo.columns:
    prop = cuenta_tipo[tipo] / total
    ajuste = sm.GLM(prop, X, family=sm.families.Binomial(),
                    var_weights=total).fit(scale="X2")
    filas.append({"tipo": tipo, "logFC": ajuste.params["AD"] / np.log(2),
                  "p": ajuste.pvalues["AD"]})
abundancia = pd.DataFrame(filas)
abundancia["FDR"] = multipletests(abundancia["p"], method="fdr_bh")[1]
abundancia = abundancia.sort_values("p")
print(abundancia.round(4))
abundancia.to_csv(os.path.join(outdir, "10_abundancia_diferencial.csv"), index=False)
# Nota: aqui logFC es el log2 del odds ratio de la proporcion, no el logFC de edgeR;
# se compara con R la direccion y la significacion, no el valor exacto.


# (9) GUARDADO
adata.write_h5ad(os.path.join(proyecto, "results/objects/lau_2020/lau_adata_python.h5ad"),
                 compression="gzip")
print("Hecho")