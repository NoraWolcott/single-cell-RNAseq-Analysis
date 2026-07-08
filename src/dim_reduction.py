import scanpy as sc


def run_pca(adata, n_comps=None):
    """
    Perform PCA.

    Parameters
    ----------
    adata : AnnData
        Annotated single-cell dataset.
    n_comps : int, optional
        Number of principal components to compute.

    Returns
    -------
    AnnData
        AnnData object with PCA results stored.
    """

    sc.tl.pca(
        adata,
        use_highly_variable=True,
        n_comps=n_comps
    )

    return adata


def get_pca_variance(adata):
    """
    Calculate variance explained by each PC.
    """

    return adata.uns["pca"]["variance_ratio"]


def run_umap(
    adata,
    n_neighbors=15,
    n_pcs=35
):
    """
    Compute neighborhood graph and UMAP embedding.

    Parameters
    ----------
    adata : AnnData
        Annotated single-cell dataset.
    n_neighbors : int
        Number of neighbors for graph construction.
    n_pcs : int
        Number of PCs used for neighborhood graph.

    Returns
    -------
    AnnData
        AnnData object with neighbors and UMAP results.
    """

    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs
    )

    sc.tl.umap(adata)

    return adata


def leiden_clustering(
    adata,
    resolution=0.5
):
    """
    Perform Leiden clustering.

    Parameters
    ----------
    adata : AnnData
        Annotated single-cell dataset.
    resolution : float
        Leiden clustering resolution parameter.

    Returns
    -------
    AnnData
        AnnData object with Leiden clusters.
    """

    sc.tl.leiden(
        adata,
        resolution=resolution
    )

    return adata