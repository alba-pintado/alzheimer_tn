# Enge methods: 
    # 1) Biological variation/technical variation
    # 2) Euclidean distance to cell type mean (additionally: we compute the correlation distance)
    # 3) Euclidean distance to average expression across cell types, using set of invariant genes
    # 4) Pairwise Euclidean distance of N randomly selected cells from the same cell type
    
# Enge 1: Biological variation/technical variation

def get_bins(adata):
        """
        Sort genes according to their mean expression, create equally sized bins, 
        remove first and last bins. 
        """
        import pandas as pd
        ## Get bin size (n//10) and bin indices
        n_bins, n_genes = 10, len(adata.var_names)
        bins = [(i, i + n_genes//10) for i in range(0, n_genes-n_genes//10, n_genes//10)]
    #     bins[-1] = (bins[-1][0], n_genes) # not necessary as we remove first and last bins
        
        ## Remove first and last bin (most and least expressed genes)
        bins = bins[1:-1]
        return bins
    
def get_least_variable(adata, bins):
        """
        adata: annData object
        bins: list of (idx0, idx1) first and last indices that correspond to that bin
        Returns list of least variable genes within each bin.
        """
        from scipy.stats import variation
        import pandas as pd
        ## Sort genes according to their mean expression
        gene_order = adata.var.sort_values(by='means', ascending=False).index
        least_variable = []
        for b in bins:
            adata_x = adata[:, gene_order[b[0]:b[1]]]
            n = adata_x.shape[1]
            least_var = adata_x.var.sort_index(key=lambda gene: variation(adata_x[:,adata_x.var_names == gene].X.toarray()), ascending=True)[0:(n//10)].index.tolist()
            least_variable.extend(least_var)
        return least_variable

def enge_transcriptional_noise(adata, batch):

    """
    Compute the transcriptional noise as the biological variation over the technical 
    variation. It can only be computed in datasets with ERCC spike-ins. The biological
    variation is computed as the correlation distance between each cell and the average 
    gene expression of all the cells of the same cell type and the same batch. The
    technical variation is computed as the correlation distance between each cell and the 
    mean ERCC spike-in expression of all the cell of the same cell type and batch. 
    The transcriptional noise is computed as the biological variation over the technical
    variation. 
    
    Parameters
    ----------
    adata: :class:`annData`
        annData object with gene expression data. It must contain a slot with the 
        batch identity of each cell in adata.obs
    batch: ``str``
        batch label (Examples: 'donor', 'patient', 'mouse')
    Returns 
    -------
    adata: :class:`annData`
        annData object with gene expression data. Euclidean distances are
        stored in adata.obs['cordist_bio'], adata.obs['cordist_tech']
        and adata.obs['noise']
        
    """
    import pandas as pd
    import numpy as np
    from tqdm import tqdm
    from scipy.spatial.distance import correlation


    ercc = [var for var in adata.var_names if var.startswith('ERCC-')]
    if len(ercc) == 0:
        print('Technical noise cannot be quantified as there are no ERCC spike-in controls.')
        return adata
#     bio = [var for var in adata.var_names if not var.startswith('ERCC-')]

    # Create categories (combinations of batch x cell type)
    batches = list(set(adata.obs[batch]))
    cell_types = list(set(adata.obs['cell_type']))
    combs = [f'{batch}_{ct}' for batch in batches for ct in cell_types]
    
    # Compute mean expression vector per category
    mean_expr_biological = pd.DataFrame(index=adata.var_names, columns=combs)
    # mean_expr_biological = pd.DataFrame(index=bio, columns=combs)
    mean_expr_technical = pd.DataFrame(index=ercc, columns=combs)
    for b in tqdm(batches):
        for ct in cell_types: 
            mean_expr_biological.loc[:, f'{b}_{ct}'] = adata[(adata.obs[batch] == b) & (adata.obs['cell_type'] == ct), :].X.toarray().mean(axis=0)
#             mean_expr_biological.loc[:, f'{b}_{ct}'] = adata[(adata.obs[batch] == b) & (adata.obs['cell_type'] == ct), bio].X.mean(axis=0)
            mean_expr_technical.loc[:, f'{b}_{ct}'] = adata[(adata.obs[batch] == b) & (adata.obs['cell_type'] == ct), ercc].X.toarray().mean(axis=0)
    
    # Compute correlation distance between each cell and its mean expression vector (per category)
    df_cordist = pd.DataFrame(index=adata.obs_names, columns=['cordist_bio', 'cordist_tech'])
    for b in tqdm(batches):
        for ct in cell_types:
            cells = adata[(adata.obs[batch] == b) & (adata.obs['cell_type'] == ct)].obs_names
#             print(donor, ct, len(cells))
            for cell in cells: 
                df_cordist.loc[cell, 'cordist_bio'] = correlation(adata[cell, :].X.toarray().ravel(), mean_expr_biological.loc[:, f'{b}_{ct}'])
    #             df_cordist.loc[cell, 'cordist_bio'] = correlation(adata[cell, bio].X.toarray(), mean_expr_biological.loc[bio, f'{b}_{ct}'])
                df_cordist.loc[cell, 'cordist_tech'] = correlation(adata[cell, ercc].X.toarray().ravel(), mean_expr_technical.loc[ercc, f'{b}_{ct}'])
    df_cordist.loc[:, 'transcriptional_noise'] = df_cordist.cordist_bio.astype('float64')/df_cordist.cordist_tech.astype('float64')

    # Update annData object
    adata.obs[f'bio_var'] = df_cordist.loc[:, 'cordist_bio'].astype('float64')
    adata.obs[f'tech_var'] = df_cordist.loc[:, 'cordist_tech'].astype('float64')
    adata.obs[f'noise'] = df_cordist.loc[:, 'transcriptional_noise'].astype('float64')
    
    return adata

# Enge 2: Euclidean distance to cell type mean (additionally: we compute the correlation and manhattan distances)
def distance_to_celltype_mean(adata, batch):

    """
    Compute the distance between each cell and the mean expression of its cell type in the 
    same batch (donor/mouse). It computes three distances: euclidean, correlation and manhattan. 
    
    Parameters
    ----------
    adata: :class:`annData`
        annData object with gene expression data. It must contain a slot with the 
        batch identity of each cell in adata.obs
    batch: ``str``
        batch label (Examples: 'donor', 'patient', 'mouse')
    Returns 
    -------
    adata: :class:`annData`
        annData object with gene expression data. Euclidean distances are
        stored in adata.obs['cor_dist'], adata.obs['euc_dist'] and adata.obs['man_dist']
        
    """
    import pandas as pd
    import numpy as np
    from tqdm import tqdm
    from scipy.spatial.distance import euclidean, correlation, cityblock

    
    # Create categories (combinations of batch x cell type)
    batches = list(set(adata.obs[batch]))
    cell_types = list(set(adata.obs['cell_type']))
    combs = [f'{batch}_{ct}' for batch in batches for ct in cell_types]
    
    # Compute mean expression vector per category
    mean_expr = pd.DataFrame(index=adata.var_names, columns=combs)
    for b in tqdm(batches):
        for ct in cell_types: 
            mean_expr.loc[:, f'{b}_{ct}'] = adata[(adata.obs[batch] == b) & (adata.obs['cell_type'] == ct), :].X.toarray().mean(axis=0)
    
    # Compute correlation distance between each cell and its mean expression vector (per category)
    df_dist = pd.DataFrame(index=adata.obs_names, columns=['cor_dist', 'euc_dist'])
    for b in batches:
        for ct in cell_types:
            cells = adata[(adata.obs[batch] == b) & (adata.obs['cell_type'] == ct)].obs_names
            for cell in tqdm(cells): 
                df_dist.loc[cell, 'cor_dist'] = correlation(adata[cell, :].X.toarray().ravel(), mean_expr.loc[:, f'{b}_{ct}'])
                df_dist.loc[cell, 'euc_dist'] = euclidean(adata[cell, :].X.toarray().ravel(), mean_expr.loc[:, f'{b}_{ct}'])
                df_dist.loc[cell, 'man_dist'] = cityblock(adata[cell, :].X.toarray().ravel(), mean_expr.loc[:, f'{b}_{ct}'])

    # Update annData object
    adata.obs[f'cor_dist'] = df_dist.loc[:, 'cor_dist'].astype('float64')
    adata.obs[f'euc_dist'] = df_dist.loc[:, 'euc_dist'].astype('float64')
    adata.obs[f'man_dist'] = df_dist.loc[:, 'man_dist'].astype('float64')
    
    return adata


# Enge 3: Euclidean distance to average expression across cell types, using set of invariant genes
def enge_euclidean_dist(adata):
    """
    Compute the Euclidean distance to the average expression across cell types
    using a set of invariant genes. The invariant genes are selected as follows: 
    1) Create equally sized bins of genes according to their mean expression
    2) Discard the two most extreme bins (lowest and highest mean expression)
    3) Select the 10% with the lowest coefficient of variation within each of the 
    remaining bins
        
    Parameters
    ----------
    adata: :class:`annData`
        annData object with gene expression data. 
        
    Returns 
    -------
    adata: :class:`annData`
        annData object with gene expression data. Euclidean distances are
        stored in adata.obs['euc_dist_tissue_invar']
        
    """
    from tqdm import tqdm
    import pandas as pd
    import numpy as np
    from scipy.spatial.distance import euclidean
    
    ## Sort genes according to their mean expression. 
    gene_order = adata.var.sort_values(by='means', ascending=False).index
    
    # Euclidean distance method as described in their methods section
    bins = get_bins(adata)
    least_variable = get_least_variable(adata, bins)
    # Mean expression across all cells
    mean_expr = adata[:, least_variable].X.mean(axis=0)
        
    adata.obs['euc_dist_tissue_invar'] = [euclidean(mean_expr, adata[i, least_variable].X.toarray().ravel()) for i in tqdm(range(adata.shape[0]))]
    
    return adata


def pairwise_euclidean_sample(adata, cell_type, n):
    from tqdm import tqdm

    """
    Computes the euclidean distance between n randomly selected cells of the same cell type. 
    
    Parameters
    ----------
    adata: :class:`annData`
        annData object with gene expression data. Cell type annotations must be 
        stored in adata.obs['cell_type']
    cell_type: ``str``
        cell type label
    n: ``int``
        number of cell pairs to 
    
    
    Returns 
    -------
    list_distances: ``list``
        list of the n euclidean pairwise distances
    
    """
    
    def euclidean_random_pair(adata, cell_type):
        
        """
        Takes an annData object and a cell type label
        and returns the euclidean distance between a pair of randomly
        selected cells of that cell type. 
        """
        from scipy.spatial.distance import euclidean
        import numpy as np

        # Pairwise distances of randomly selected 10,000 pairs of cells
        if cell_type == 'all':
            adata_x = adata
        else:
            adata_x = adata[adata.obs['cell_type'] == cell_type]
        n_cells = adata_x.shape[0]        
        vec = np.random.choice(np.arange(n_cells), 2, replace=False)
       
        return euclidean(adata_x[vec[0], :].X[0].toarray().ravel(), adata_x[vec[1], :].X[0].toarray().ravel())
    
    list_distances = []
    for i in tqdm(range(n)): 
        list_distances.append(euclidean_random_pair(adata, cell_type))
    
    return list_distances


def hernando_herraez(adata, batch):
    """
    Computes the correlation distance of each cell to the cell type median using 
    the 500 most variably expressed genes. 
    
    Parameters
    ----------
    adata: :class:`annData`
        annData object with gene expression data. Cell type annotations must be 
        stored in adata.obs['cell_type']
    cell_type: ``str``
        cell type label
    batch: ``str``
        batch label (Examples: 'donor', 'patient', 'mouse')
  
    Returns 
    -------
    adata: :class:`annData`
        annData object with gene expression data. Euclidean distances are
        stored in adata.obs['cor_dist_median']
    """
    import pandas as pd
    import numpy as np
    from scipy.spatial.distance import correlation
    from tqdm import tqdm

    # Create categories (combinations of batch x cell type)
    batches = list(set(adata.obs[batch]))
    cell_types = list(set(adata.obs['cell_type']))
    combs = [f'{batch}_{ct}' for batch in batches for ct in cell_types]
    
    # Take 500 genes with the greatest dispersion
    hvgs = adata.var.sort_values(by='dispersions_norm', ascending=False).index[0:500]
    
    # Compute median expression vector per category
    median_expr = pd.DataFrame(index=hvgs, columns=combs)
    for b in tqdm(batches):
        for ct in cell_types: 
            median_expr.loc[:, f'{b}_{ct}'] = np.median(adata[(adata.obs[batch] == b) & (adata.obs['cell_type'] == ct), hvgs].X.toarray(), axis=0)
    

    # Compute correlation distance between each cell and its median expression vector (per category)
    df_dist = pd.DataFrame(index=adata.obs_names, columns=['cor_dist_median'])
    for b in tqdm(batches):
        for ct in cell_types:
            cells = adata[(adata.obs[batch] == b) & (adata.obs['cell_type'] == ct)].obs_names
            for cell in cells: 
                df_dist.loc[cell, 'cor_dist_median'] = correlation(adata[cell, hvgs].X.toarray(), median_expr.loc[:, f'{b}_{ct}'])
               
    # Update annData object
    adata.obs[f'cor_dist_median'] = df_dist.loc[:, 'cor_dist_median'].astype('float64')
    return adata

    
def distance_to_celltype_mean_invariant(adata, batch):

    """
    Compute the distance between each cell and the mean expression of its cell type in the 
    same batch (donor/mouse), using a set of invariant genes as in Enge (2017). 
    It computes three distances: euclidean, correlation and manhattan. 
    
    Parameters
    ----------
    adata: :class:`annData`
        annData object with gene expression data. It must contain a slot with the 
        batch identity of each cell in adata.obs
    batch: ``str``
        batch label (Examples: 'donor', 'patient', 'mouse')
    Returns 
    -------
    adata: :class:`annData`
        annData object with gene expression data. Euclidean distances are
        stored in adata.obs['cor_dist_invar'], adata.obs['euc_dist_invar'] and adata.obs['man_dist_invar']
        
    """
    import pandas as pd
    import numpy as np
    from tqdm import tqdm
    from scipy.spatial.distance import euclidean, correlation, cityblock

    
    # Create categories (combinations of batch x cell type)
    batches = list(set(adata.obs[batch]))
    cell_types = list(set(adata.obs['cell_type']))
    combs = [f'{batch}_{ct}' for batch in batches for ct in cell_types]
    
   
    # Euclidean distance method as described in their methods section
    bins = get_bins(adata)
    lvgs = get_least_variable(adata, bins)
    
    # Compute mean expression vector per category
    mean_expr = pd.DataFrame(index=lvgs, columns=combs)
    for b in tqdm(batches):
        for ct in cell_types: 
            mean_expr.loc[:, f'{b}_{ct}'] = adata[(adata.obs[batch] == b) & (adata.obs['cell_type'] == ct), lvgs].X.toarray().mean(axis=0)
    
    # Compute correlation distance between each cell and its mean expression vector (per category)
    df_dist = pd.DataFrame(index=adata.obs_names, columns=['cor_dist', 'euc_dist'])
    for b in tqdm(batches):
        for ct in cell_types:
            cells = adata[(adata.obs[batch] == b) & (adata.obs['cell_type'] == ct)].obs_names
            for cell in cells: 
                df_dist.loc[cell, 'cor_dist_invar'] = correlation(adata[cell, lvgs].X.toarray().ravel(), mean_expr.loc[:, f'{b}_{ct}'])
                df_dist.loc[cell, 'euc_dist_invar'] = euclidean(adata[cell, lvgs].X.toarray().ravel(), mean_expr.loc[:, f'{b}_{ct}'])
                df_dist.loc[cell, 'man_dist_invar'] = cityblock(adata[cell, lvgs].X.toarray().ravel(), mean_expr.loc[:, f'{b}_{ct}'])

    # Update annData object
    adata.obs[f'cor_dist_invar'] = df_dist.loc[:, 'cor_dist_invar'].astype('float64')
    adata.obs[f'euc_dist_invar'] = df_dist.loc[:, 'euc_dist_invar'].astype('float64')
    adata.obs[f'man_dist_invar'] = df_dist.loc[:, 'man_dist_invar'].astype('float64')
    
    return adata
    
    
    
def gcl(adata, num_divisions):
    """
    Following the original GCL.m script provided by the authors
    (https://github.com/guy531/gcl). 
    """
    from sklearn.metrics.pairwise import euclidean_distances
    import pandas as pd
    import numpy as np

    n, num_genes = adata.shape #changed
    data = adata.X.todense().transpose() #changed
    gcl_output = np.zeros(num_divisions)

    def Vn(Aij, Bij):                    
        return 1/(n*(n-3)) * (np.sum(Aij*Bij) - n/(n-2)*np.sum(np.diag(Aij)*np.diag(Bij)))

    def Rn(Aij, Bij):
        return Vn(Aij, Bij)/np.sqrt(Vn(Aij, Aij)*Vn(Bij, Bij))
    
    def compute_Aij(X, n):
        d1 = euclidean_distances(X, X)
        m1 = np.mean(d1, 0)
        M1 = np.mean(d1)
        Aij = d1 - m1.T * np.ones((1, n))
        Aij = Aij.T - np.ones((n, 1)) * m1 
        Aij += M1
        Aij = Aij - d1 / n
        Aij[np.diag_indices(len(Aij))] = m1 - M1
        Aij = (n / (n-1)) * Aij
        return Aij


    for i in range(num_divisions):
        random_genes = np.arange(num_genes)
        np.random.shuffle(random_genes)

        X1 = data[random_genes[:num_genes//2],:].T
        X2 = data[random_genes[num_genes//2:],:].T
    
        Aij1 = compute_Aij(X1, n)
        Aij2 = compute_Aij(X2, n)
        
        gcl_output[i] = Rn(Aij2, Aij1)
        
    return gcl_output

def gcl_per_cell_type_and_batch(adata, num_divisions, batch): 
    """
    Compute GCL for each cell type and batch in adata.obs['batch'].
    
    Parameters
    ----------
    adata: :class:`annData`
        annData object with gene expression data. It must contain a slot with the 
        batch identity of each cell in adata.obs
    num_divisions: `int`
        number of iterations to use in  gcl()
    
    Returns 
    -------
    output: :class:`pd.DataFrame`
        Pandas dataframe with the GCL per cell type x batch x iteration
        
    """
    import pandas as pd
    from tqdm import tqdm
    cell_types = list(set(adata.obs['cell_type']))
    batches = list(set(adata.obs[batch]))
    # ages = list(set(adata_small.obs['age']))
    combinations = [f'{ct}_{b}_{i}' for ct in cell_types for b in batches for i in range(num_divisions)]
    output = pd.DataFrame(index=combinations, columns=['GCL', 'condition', 'cell_type', batch])
    for ct in tqdm(cell_types):
        for b in batches:
            adata_sel = adata[(adata.obs['cell_type']==ct) & (adata.obs[batch]==b)]
#             sc.pp.filter_genes(adata_sel, min_cells=3)
#             sc.pp.filter_cells(adata_sel, min_genes=5)
            if adata_sel.shape[0] > 10:
                output.loc[[f'{ct}_{b}_{i}' for i in range(num_divisions)], 'GCL'] = gcl(adata_sel, num_divisions=num_divisions)
                output.loc[[f'{ct}_{b}_{i}' for i in range(num_divisions)], 'condition'] = adata_sel.obs['condition'][0]
                output.loc[[f'{ct}_{b}_{i}' for i in range(num_divisions)], 'cell_type'] = ct
                output.loc[[f'{ct}_{b}_{i}' for i in range(num_divisions)], batch] = b
    return output



def rerun_preprocessing(adata, batch_key):
    """
    Re-runs preprocessing steps: filter lowly expressed genes, compute HVGs, 
    run batch-effect corrected PCA (harmony), neighbors. 
    
    Parameters
    ----------
    adata: :class:`annData`
        annData object with gene expression data. It must contain a slot with the 
        batch identity of each cell in adata.obs
    batch_key: ``str``
        batch label to run batch-effect correction (Examples: 'donor', 'patient', 'mouse')
    Returns 
    -------
    adata: :class:`annData`
        updated annData object with gene expression data. 
        
    """
    
    import scanpy as sc
    import scanpy.external as sce
    
    sc.pp.filter_genes(adata, min_cells=5)
    sc.pp.highly_variable_genes(adata)
    sce.pp.harmony_integrate(adata, key=batch_key)
    sc.tl.umap(adata)
    
    return adata
    
    
def scallop_pipeline(adata, res_vals=None): 
    
    """
    Compute transcriptional noise as 1 - membership score (averaged over aa range of resolution values). 
    It runs the whole Scallop pipeline: 
    1) Create separate annData object per condition (young/old, smoker/non-smoker) and cell type. 
    2) Re-run preprocessing on annData
    3) Run Scallop over range of resolution values
    4) Compute average membership score across resolutions
    5) Transcriptional noise = 1 - mean membership score
    
    Parameters
    ----------
    adata: :class:`annData`
        annData object with gene expression data. It must contain a slot with the 
        batch identity of each cell in adata.obs
    batch: ``str``
        batch label (Examples: 'donor', 'patient', 'mouse')
    Returns 
    -------
    adata: :class:`annData`
        annData object with gene expression data. Euclidean distances are
        stored in adata.obs['cor_dist_invar'], adata.obs['euc_dist_invar'] and adata.obs['man_dist_invar']
        
    """
    import pandas as pd
    import numpy as np
    from tqdm import tqdm
    import scallop as sl
    import scanpy as sc
    
    sc.tl.leiden(adata, resolution=0.1)
    if res_vals == None: 
        res_vals = [round(i*0.1, 2) for i in range(1, 11)]

    big_clusters = list(set(adata.obs['leiden']))
#     individuals = list(set(adata.obs['batch'])) # mouse/donor
    conditions = list(set(adata.obs['condition'])) # smoker/non-smoker, young/old
    adata.obs['scallop_noise'] = np.nan

    for cond in conditions:
        for clus in big_clusters: 
            adata_sample = adata[(adata.obs['condition'] == cond) & (adata.obs['leiden'] == clus)]
            adata_sample.obs.loc[:, [f'freq_score_{res}' for res in res_vals]] = np.nan
            sc.pp.filter_cells(adata_sample, min_genes=100)
            if adata_sample.shape[0] > 50:
            	sc.pp.filter_genes(adata_sample, min_cells=3)
            	sc.pp.pca(adata_sample)
            	sc.pp.neighbors(adata_sample)
            	# Run Scallop
            	scal = sl.Scallop(adata_sample)
            	for i, res in enumerate(res_vals):
                	adata_sample.obs[f'freq_score_{res}'] = sl.tl.getScore(scal, score_type='freq', res=res, n_trials=30, frac_cells=0.8, do_return=True)


            mean_membership = adata_sample.obs.loc[:, [f'freq_score_{res}' for res in res_vals]].mean(axis=1)
            adata.obs.loc[mean_membership.index, 'scallop_noise'] = 1 - mean_membership

    return adata
