import numpy as np
import scanpy as sc
import pandas as pd
import os


def rank_cluster_markers(
    adata,
    groupby="leiden",
    method="wilcoxon"
):
    """
    Identify marker genes for each cluster.
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
    """

    counts = np.asarray(
        adata[:, genes].layers[layer].sum(axis=1)
    ).flatten()

    return counts


def save_top_genes_csv(
    top_genes,
    filename="top_DE_genes.csv",
    results_dir="../results"
):
    """
    Save top ranked genes for each group to CSV.
    """

    os.makedirs(results_dir, exist_ok=True)

    gene_table = pd.DataFrame(top_genes)

    filepath = os.path.join(
        results_dir,
        filename
    )

    gene_table.to_csv(
        filepath,
        index_label="rank"
    )

    return filepath