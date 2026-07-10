import scanpy as sc


def run_pca(
    adata,
    n_comps=None
):
    """
    Perform PCA using highly variable genes.
    """

    sc.tl.pca(
        adata,
        use_highly_variable=True,
        n_comps=n_comps
    )

    return adata


def run_umap(
    adata,
    n_neighbors=15,
    n_pcs=35
):
    """
    Compute neighborhood graph and UMAP embedding.
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
    """

    sc.tl.leiden(
        adata,
        resolution=resolution
    )

    return adata

def run_dimensionality_reduction(
    adata,
    n_top_genes=2000,
    n_pcs=20,
    leiden_resolution=None
):
    """
    Run HVG selection, PCA, neighborhood graph, UMAP,
    and optionally Leiden clustering.
    """

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes
    )

    sc.tl.pca(
        adata,
        use_highly_variable=True
    )

    sc.pp.neighbors(
        adata,
        n_pcs=n_pcs
    )

    sc.tl.umap(adata)

    if leiden_resolution is not None:
        sc.tl.leiden(
            adata,
            resolution=leiden_resolution
        )

    return adata