Tutorial
==================

In order to use ``decibel`` you need a basic knowledge of how ``scanpy`` and 
``annData`` objects work.

Computing transcriptional noise
-------------------------------------------------

In order to work with decibel you need to first run ``scanpy`` and load 
a dataset in an ``annData`` object. ::

    adata = sc.read('path_to_file/filename.h5ad')
    

	
Then process data (normalization, log-transformation, QC filtering on cells and genes). ::

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_cells(adata, min_genes=100)


  
Run PCA, feature selection (triku), batch effect correction (harmony) and dimensionality reduction (UMAP). ::

    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    tk.tl.triku(adata)
    sc.pp.pca(adata)
    sce.pp.harmony_integrate(adata, 'batch')
    sc.pp.neighbors(adata, use_rep='X_pca_harmony')
    sc.tl.umap(adata)
	


Compute transcriptional noise as in Enge et al, (2017). ::

    dcb.enge_transcriptional_noise(adata, 'batch')



    
