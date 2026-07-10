import scanpy as sc

def score_gene_sets(
    adata,
    marker_sets
):
    """
    Calculate Scanpy module scores for multiple gene sets.
    """

    for score_name, genes in marker_sets.items():

        sc.tl.score_genes(
            adata,
            genes,
            score_name=score_name
        )

    return adata

def identify_mature_osn_clusters(
    adata,
    cluster_key="leiden"
):
    """
    Identify clusters enriched for mature OSN signatures.
    """

    cluster_scores = (
        adata.obs.groupby(cluster_key)[
            [
                "mature_OSN_score",
                "immature_OSN_score",
                "non_neuronal_score",
                "total_counts",
                "pct_counts_mito"
            ]
        ]
        .mean()
    )

    print("Cluster score summary:")
    print(cluster_scores)

    likely_mature_osns = cluster_scores[
        (cluster_scores["mature_OSN_score"] >
         cluster_scores["mature_OSN_score"].median()) &
        (cluster_scores["immature_OSN_score"] <
         cluster_scores["immature_OSN_score"].median()) &
        (cluster_scores["non_neuronal_score"] <
         cluster_scores["non_neuronal_score"].median())
    ].index.tolist()

    print("\nLikely mature OSN clusters:")
    print(likely_mature_osns)

    return likely_mature_osns, cluster_scores

def screen_osn_doublets(
    adata,
    cluster_key="leiden",
    count_multiplier=1.75,
    multiple_or_fraction=0.25
):
    """
    Identify clusters with features consistent with OSN-OSN doublets.
    """

    cluster_doublet_screen = (
        adata.obs.groupby(cluster_key)[
            [
                "total_counts",
                "n_genes_by_counts",
                "n_ORs"
            ]
        ]
        .mean()
    )

    # Fraction of cells with multiple expressed OR genes
    or_frequency = (
        adata.obs.groupby(cluster_key)["n_ORs"]
        .apply(lambda x: (x > 1).mean())
    )

    cluster_doublet_screen["pct_cells_multiple_ORs"] = or_frequency

    # Identify high-count clusters enriched for multiple OR expression
    mean_counts_all_clusters = (
        cluster_doublet_screen["total_counts"].mean()
    )

    candidate_clusters = cluster_doublet_screen[
        (cluster_doublet_screen["total_counts"] >
         count_multiplier * mean_counts_all_clusters) &
        (cluster_doublet_screen["pct_cells_multiple_ORs"] >
         multiple_or_fraction)
    ].index.tolist()

    return candidate_clusters, cluster_doublet_screen

def screen_osn_nonneuronal_doublets(
    adata,
    cluster_key="leiden"
):
    """
    Identify clusters with both OSN and non-neuronal signatures.
    """

    cluster_screen = (
        adata.obs.groupby(cluster_key)[
            [
                "mature_OSN_score",
                "non_neuronal_score",
                "immature_OSN_score",
                "total_counts",
                "pct_counts_mito"
            ]
        ]
        .mean()
    )

    # Identify clusters with both OSN and non-neuronal signatures
    candidate_clusters = cluster_screen[
        (cluster_screen["mature_OSN_score"] >
         cluster_screen["mature_OSN_score"].median()) &
        (cluster_screen["non_neuronal_score"] >
         cluster_screen["non_neuronal_score"].median())
    ].index.tolist()

    return candidate_clusters, cluster_screen

def filter_doublet_clusters(
    mature_osn_clusters,
    doublet_clusters
):
    """
    Remove candidate doublet clusters from mature OSN candidates.
    """

    filtered_clusters = [
        cluster
        for cluster in mature_osn_clusters
        if cluster not in doublet_clusters
    ]

    return filtered_clusters

def remove_clusters(
    adata,
    clusters_to_remove,
    cluster_key="leiden"
):
    """
    Remove selected clusters from an AnnData object.
    """

    adata_filtered = adata[
        ~adata.obs[cluster_key].isin(clusters_to_remove)
    ].copy()

    return adata_filtered

def subset_clusters(
    adata,
    clusters,
    cluster_key="leiden"
):
    """
    Subset AnnData object to only mature OSNs, no doublets.
    """

    adata_subset = adata[
        adata.obs[cluster_key].isin(clusters)
    ].copy()

    return adata_subset

def identify_clusters_by_score(
    adata,
    score_key,
    threshold,
    cluster_key="leiden",
    comparison=">"
):
    """
    Identify clusters based on the mean value of a module score.
    """

    cluster_scores = (
        adata.obs
        .groupby(cluster_key)[
            [score_key, "total_counts", "pct_counts_mito"]
        ]
        .mean()
    )

    if comparison == ">":
        selected_clusters = cluster_scores[
            cluster_scores[score_key] > threshold
        ].index.tolist()
    elif comparison == "<":
        selected_clusters = cluster_scores[
            cluster_scores[score_key] < threshold
        ].index.tolist()
    else:
        raise ValueError("comparison must be '>' or '<'")

    return selected_clusters, cluster_scores