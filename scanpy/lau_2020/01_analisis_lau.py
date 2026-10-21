# Lau et al. 2020: análisis en Python (scanpy)

# TFM "Transcriptional Noise in Alzheimer's Disease Progression" (CBMSO)

# Versión en Python de scrapper_lau_analisis.R. Parte del mismo objeto
# preprocesado en R, exportado a .h5ad con 00_exportar_h5ad.R, de modo que los
# núcleos son los mismos y los resultados se pueden comparar con los de R.

# Equivalencias R -> Python:
#   genes variables  modelGeneVariances      -> highly_variable_genes
#   corrección       MNN                     -> Harmony
#   clustering       SNN + clusterGraph      -> kNN + Leiden
#   marcadores       scoreMarkers (AUC)      -> rank_genes_groups (Wilcoxon)
#   pseudobulk       edgeR                   -> PyDESeq2
#   abundancia       edgeR                   -> GLM quasi-binomial

# Ejecutar desde la carpeta raíz del proyecto:
#   python scripts/SCANPY/lau_2020/01_analisis_lau.py > scripts/SCANPY/lau_2020/log_python.txt 2>&1

# Fuentes: scanpy (Wolf et al. 2018); Harmony (Korsunsky et al. 2019);
# Leiden (Traag et al. 2019); PyDESeq2 (Muzellec et al. 2023);
# single-cell best practices (Heumos et al. 2023, https://www.sc-best-practices.org)


# (0) PREPARACIÓN
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
np.random.seed(SEMILLA)
sc.settings.verbosity = 1

proyecto = os.getcwd()                     # ejecutar desde la carpeta raíz del proyecto
entrada  = os.path.join(proyecto, "results/objects/lau_2020/lau_sce_analisis.h5ad")
outdir   = os.path.join(proyecto, "results/05_scanpy/Lau_32989152")
os.makedirs(outdir, exist_ok=True)
sc.settings.figdir = outdir

adata = sc.read_h5ad(entrada)
print(adata)                                # 140727 núcleos x 29456 genes

# (0.1) COMPROBACIONES DE LA IMPORTACIÓN
# Mismos núcleos y donantes que en R, y X debe ser la matriz logarítmica.
assert adata.obs_names.is_unique and adata.var_names.is_unique
assert "counts" in adata.layers
print(pd.crosstab(adata.obs["donor"], adata.obs["condition"]))
assert adata.obs["donor"].nunique() == 20 and "AD6" not in set(adata.obs["donor"])

for c in ["donor", "condition", "celltype", "cluster"]:
    adata.obs[c] = adata.obs[c].astype("category")
adata.obs["cluster_R"] = adata.obs["cluster"]          # clusters de R, para comparar
anotada = adata.obs["celltype"].notna().values
simbolo = adata.var["simbolo"].astype(str).where(adata.var["simbolo"].notna(),
                                                adata.var_names.to_series())

# (1) GENES VARIABLES
# Varianza frente a media sobre la matriz logarítmica, calculada dentro de cada
# donante y combinada (batch_key), y los 3000 genes más variables, como en R.
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=3000,
                            batch_key="donor")
print(adata.var["highly_variable"].sum())               # 3000

# (2) PCA Y R² FRENTE A TIPO CELULAR Y DONANTE
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
    pcs = adata.obsm[clave_obsm][anotada, :6]
    return pd.DataFrame({
        "pc": range(1, 7),
        "tipo": r2_por_grupo(pcs, adata.obs["celltype"].values[anotada]),
        "donante": r2_por_grupo(pcs, adata.obs["donor"].values[anotada])}).round(3)

r2_pca = tabla_r2("X_pca")
print(r2_pca)

plt.figure(figsize=(5, 4))
plt.semilogy(range(1, 26), var_pct, "o-")
plt.xlabel("PC"); plt.ylabel("Varianza explicada (%)")
plt.tight_layout(); plt.savefig(os.path.join(outdir, "01_scree_plot.png"), dpi=150); plt.close()

# (3) CORRECCIÓN POR DONANTE (HARMONY)
# Como en R, cada donante es solo AD o solo control: la corrección elimina
# también parte del efecto de la enfermedad. El espacio corregido se usa SOLO
# para UMAP y clustering; marcadores y AD frente a control usan datos sin corregir.
# Se llama a harmonypy directamente (más estable entre versiones que el envoltorio
# de scanpy); según la versión devuelve la matriz en una u otra orientación.
import harmonypy
ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, "donor",
                           random_state=SEMILLA, verbose=False)
z = np.asarray(ho.Z_corr)
adata.obsm["X_pca_harmony"] = z if z.shape[0] == adata.n_obs else z.T
r2_harmony = tabla_r2("X_pca_harmony")
print(r2_harmony)
pd.concat([r2_pca.add_suffix("_pca"), r2_harmony.add_suffix("_harmony")], axis=1) \
  .to_csv(os.path.join(outdir, "02_r2_pca_harmony.csv"), index=False)

# (4) UMAP
sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=15,
                random_state=SEMILLA)
sc.tl.umap(adata, random_state=SEMILLA)

# Núcleos en orden aleatorio para que ningún donante ni condición tape al resto.
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

# (5) CLUSTERING (LEIDEN) Y VALIDACIÓN
sc.tl.leiden(adata, resolution=1.0, key_added="cluster_py",
             flavor="igraph", n_iterations=2, random_state=SEMILLA)
print(adata.obs["cluster_py"].value_counts().sort_index())

tab = pd.crosstab(adata.obs["cluster_py"][anotada], adata.obs["celltype"][anotada])
pureza = tab.max(axis=1) / tab.sum(axis=1)
tipo_cluster = tab.idxmax(axis=1)
print("Clusters con pureza >= 0,99:", (pureza >= 0.99).sum(), "de", len(pureza))
print(pureza.sort_values().head().round(3))

ari_tipo = adjusted_rand_score(adata.obs["celltype"][anotada], adata.obs["cluster_py"][anotada])
ari_R    = adjusted_rand_score(adata.obs["cluster_R"], adata.obs["cluster_py"])
print(f"ARI frente a la anotación: {ari_tipo:.3f}")
print(f"ARI frente a los clusters de R: {ari_R:.3f}")

max_donante = pd.crosstab(adata.obs["cluster_py"], adata.obs["donor"], normalize="index").max(axis=1)
print("Clusters dominados (> 70 %) por un donante:", list(max_donante[max_donante > 0.7].index))
mito = adata.obs.groupby("cluster_py", observed=True)["subsets_mito_percent"].median()
print(mito.sort_values(ascending=False).head().round(2))

fig = sc.pl.umap(adata, color="cluster_py", legend_loc="on data", size=2,
                 show=False, return_fig=True)
fig.savefig(os.path.join(outdir, "06_umap_por_cluster.png"), dpi=150, bbox_inches="tight"); plt.close(fig)

# (6) MARCADORES
# Sobre la matriz logarítmica sin corregir (X). Wilcoxon, cada cluster frente al resto.
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
# Cuentas crudas sumadas por donante dentro de cada tipo celular: la réplica es
# el donante, no el núcleo (Squair et al. 2021). Control como referencia.
def pseudobulk(ad, capa, grupo):
    """Suma las cuentas de cada grupo (donante). Devuelve donantes x genes."""
    g = pd.Categorical(ad.obs[grupo])
    indicadora = sp.csr_matrix((np.ones(ad.n_obs), (g.codes, np.arange(ad.n_obs))),
                               shape=(len(g.categories), ad.n_obs))
    suma = indicadora @ ad.layers[capa]
    suma = suma.toarray() if sp.issparse(suma) else np.asarray(suma)
    return pd.DataFrame(suma.astype(int), index=g.categories, columns=ad.var_names)

resumen_DE = []
for tipo in sorted(adata.obs["celltype"].dropna().unique()):
    sub = adata[adata.obs["celltype"] == tipo]
    cuentas = pseudobulk(sub, "counts", "donor")
    meta = (sub.obs.groupby("donor", observed=True)["condition"].first()
            .reindex(cuentas.index).to_frame())
    meta["condition"] = pd.Categorical(meta["condition"].astype(str),
                                       categories=["Ctrl", "AD"])

    # Filtro equivalente a filterByExpr: al menos 10 cuentas en tantos donantes
    # como tiene el grupo más pequeño.
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
# Para cada tipo celular, proporción de núcleos por donante ~ condición con un
# GLM quasi-binomial (pesos = núcleos anotados del donante). Equivale a la
# pregunta del análisis de abundancia de R: ¿cambia la fracción de cada tipo?
cuenta_tipo = pd.crosstab(adata.obs["donor"][anotada], adata.obs["celltype"][anotada])
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
# Nota: aquí logFC es el log2 del odds ratio de la proporción, no el logFC de edgeR;
# se compara con R la dirección y la significación, no el valor exacto.

# (9) GUARDADO
adata.write_h5ad(os.path.join(proyecto, "results/objects/lau_2020/lau_adata_python.h5ad"),
                 compression="gzip")
print("Hecho")