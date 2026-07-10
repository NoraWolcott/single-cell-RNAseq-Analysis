import os
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from adjustText import adjust_text


# Repository root (one level above src/)
REPO_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..")
)

RESULTS_DIR = os.path.join(REPO_ROOT, "results")

def plot_umap_metadata(
    adata,
    colors=["sample", "condition"],
    results_dir=RESULTS_DIR,
    filename="umap_metadata.png"
):
    """
    Plot UMAP colored by sample/condition metadata.
    """

    sc.pl.umap(
        adata,
        color=colors,
        wspace=0.4,
        show=False
    )

    plt.savefig(
        os.path.join(results_dir, filename),
        bbox_inches="tight",
        dpi=300
    )

    plt.close()


def plot_umap_clusters(
    adata,
    cluster_key="leiden",
    results_dir=RESULTS_DIR,
    filename="umap_leiden_clusters.png"
):
    """
    Plot UMAP colored by clustering assignments.
    """

    sc.pl.umap(
        adata,
        color=cluster_key,
        legend_loc="on data",
        show=False
    )

    plt.savefig(
        os.path.join(results_dir, filename),
        bbox_inches="tight",
        dpi=300
    )

    plt.close()

def plot_marker_umap(
    adata,
    genes,
    cmap="cubehelix_r",
    size=2,
    vmax="p99.75",
    ncols=4,
    layer="counts",
    legend_loc="on data",
    save_path=None
):
    """
    Plot selected marker genes on UMAP.
    """

    sc.pl.umap(
        adata,
        color=genes,
        cmap=cmap,
        s=size,
        legend_loc=legend_loc,
        vmax=vmax,
        ncols=ncols,
        frameon=False,
        layer=layer,
        save=save_path
    )


def plot_score_umap(
    adata,
    scores=None,
    cmap="cubehelix_r",
    results_dir=RESULTS_DIR,
    filename="OSN_scores_umap.png"
):
    """
    Plot module scores and QC metrics on UMAP and save figure.
    """

    if scores is None:
        scores = [
            "mature_OSN_score",
            "immature_OSN_score",
            "non_neuronal_score",
            "total_counts",
            "pct_counts_mito"
        ]


    sc.pl.umap(
        adata,
        color=scores,
        cmap=cmap,
        show=False
    )


    plt.savefig(
        os.path.join(results_dir, filename),
        bbox_inches="tight",
        dpi=300
    )

    plt.close()


def plot_marker_umap(
    adata,
    genes,
    layer="counts",
    cmap="cubehelix_r",
    size=2,
    vmax="p99.75",
    ncols=4,
    results_dir=RESULTS_DIR,
    filename="marker_umap.png"
):
    """
    Plot UMAPs colored by expression of marker genes.
    """

    sc.pl.umap(
        adata,
        color=genes,
        layer=layer,
        cmap=cmap,
        s=size,
        vmax=vmax,
        ncols=ncols,
        frameon=False,
        show=False
    )

    plt.savefig(
        os.path.join(results_dir, filename),
        bbox_inches="tight",
        dpi=300
    )

    plt.close()


def plot_or_gene_counts(
    adata,
    prefix="Olfr",
    cmap="cubehelix_r",
    vmax=10,
    results_dir=RESULTS_DIR,
    filename="OR_gene_counts_umap.png"
):

    # Identify OR genes
    or_genes = [
        g for g in adata.var_names
        if g.startswith(prefix)
    ]


    # Count expressed OR genes per cell
    n_ORs = (adata[:, or_genes].X > 0).sum(axis=1)

    # Convert sparse matrix / matrix output to 1D array
    if hasattr(n_ORs, "A1"):
        n_ORs = n_ORs.A1
    else:
        n_ORs = np.asarray(n_ORs).flatten()

    adata.obs["n_ORs"] = n_ORs

    sc.pl.umap(
        adata,
        color="n_ORs",
        cmap=cmap,
        vmax=vmax,
        show=False
    )


    plt.savefig(
        os.path.join(results_dir, filename),
        bbox_inches="tight",
        dpi=300
    )

    plt.close()

    return adata


def plot_volcano(
    de_results,
    filename="volcano_plot.png",
    results_dir=RESULTS_DIR,
    fc_threshold=0.25,
    padj_threshold=0.05,
    label_genes=True,
    n_labels=20
):
    """
    Create volcano plot from differential expression results.
    """

    volcano = de_results.copy()

    # Avoid log(0)
    volcano["neg_log10_padj"] = -np.log10(
        volcano["pvals_adj"].replace(0, 1e-300)
    )

    volcano["significant"] = (
        (volcano["pvals_adj"] < padj_threshold) &
        (abs(volcano["logfoldchanges"]) > fc_threshold)
    )

    plt.figure(figsize=(8, 6))

    # All genes
    plt.scatter(
        volcano["logfoldchanges"],
        volcano["neg_log10_padj"],
        s=5,
        alpha=0.5
    )

    # Significant genes
    sig = volcano[volcano["significant"]]

    plt.scatter(
        sig["logfoldchanges"],
        sig["neg_log10_padj"],
        s=10
    )

    # Label top significant genes
    if label_genes:
        top_sig = (
            sig.sort_values(
                "pvals_adj"
            )
            .head(n_labels)
        )

        texts = []

        for _, row in top_sig.iterrows():
            texts.append(
                plt.text(
                    row["logfoldchanges"],
                    row["neg_log10_padj"],
                    row["names"]
                )
            )

        adjust_text(texts)

    # Threshold lines
    plt.axvline(
        fc_threshold,
        linestyle="--"
    )

    plt.axvline(
        -fc_threshold,
        linestyle="--"
    )

    plt.axhline(
        -np.log10(padj_threshold),
        linestyle="--"
    )

    plt.xlabel(
        "Log fold change"
    )

    plt.ylabel(
        "-log10 adjusted p-value"
    )

    plt.title(
        "Differential expression volcano plot"
    )

    plt.savefig(
        os.path.join(results_dir, filename),
        bbox_inches="tight",
        dpi=300
    )

    plt.close()