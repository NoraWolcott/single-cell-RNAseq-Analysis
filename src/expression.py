import numpy as np
import scanpy as sc


def rank_cluster_markers(
    adata,
    groupby="leiden",
    method="wilcoxon"
):
    """
    Identify marker genes for each cluster.

    Parameters
    ----------
    adata : AnnData
        Annotated single-cell dataset.
    groupby : str
        Metadata column defining groups.
    method : str
        Statistical test used for ranking genes.

    Returns
    -------
    AnnData
        AnnData object with rank_genes_groups results stored.
    """

    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method
    )

    return adata


def get_top_marker_genes(
    adata,
    n_genes=1
):
    """
    Extract top-ranked marker genes from each group.

    Returns
    -------
    dict
        Dictionary mapping groups to top marker genes.
    """

    ranked_genes = adata.uns["rank_genes_groups"]["names"]

    top_genes = {}

    for group in ranked_genes.dtype.names:
        top_genes[group] = list(
            ranked_genes[group][:n_genes]
        )

    return top_genes


def calculate_gene_counts(
    adata,
    genes,
    layer="counts"
):
    """
    Calculate total expression counts for selected genes per cell.

    Parameters
    ----------
    adata : AnnData
        Annotated single-cell dataset.
    genes : list
        Genes to summarize.
    layer : str
        AnnData layer containing counts.

    Returns
    -------
    array
        Total counts per cell.
    """

    counts = np.asarray(
        adata[:, genes].layers[layer].sum(axis=1)
    ).flatten()

    return counts