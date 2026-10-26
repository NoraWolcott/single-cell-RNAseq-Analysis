import os
import scanpy as sc
import numpy as np
import seaborn as sns
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
from adjustText import adjust_text


# Repository root (one level above src/)
REPO_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..")
)

RESULTS_DIR = os.path.join(REPO_ROOT, "results")

# Define GEP naming conventions
GEPS_NAMED = ["GEP_Dorsal", "GEP_Ventral", "GEP_Anterior", "GEP_Posterior", "GEP_Cd36", "GEP_High", "GEP_Low"]
GEPS_NUMBERED = ["GEP_8", "GEP_9", "GEP_10", "GEP_11_UPR", "GEP_12", "GEP_13", "GEP_14", "GEP_15", "GEP_16"]

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

def gep_label(gep):
    """
    Plot label for a GEP obs column (GEP_Dorsal -> Dorsal, GEP_11_UPR -> GEP 11 (UPR)).
    """
 
    label = gep.replace("GEP_", "").replace("_UPR", " (UPR)")
 
    return f"GEP {label}" if label[0].isdigit() else label
 
 
def or_means(
    adata,
    columns,
    groups,
    group_key="condition",
    or_key="top_Olfr",
    min_cells=10
):
    """
    Mean of obs columns per OR within each group, keeping ORs with enough cells in every group.
    """
 
    df = adata.obs.loc[adata.obs[group_key].isin(groups), [group_key, or_key] + list(columns)]
    grouped = df.groupby([or_key, group_key], observed=True)
    keep = grouped.size().unstack().ge(min_cells).all(axis=1)
 
    return grouped.mean().loc[keep[keep].index]
 
 
def plot_score_scatter(
    adata,
    groups,
    labels=None,
    scores=["DV_score", "AP_score", "ES_score"],
    title=None,
    group_key="condition",
    or_key="top_Olfr",
    min_cells=10,
    results_dir=RESULTS_DIR,
    filename="score_scatter.png"
):
    """
    Scatter per-OR mean scores in groups[0] (x) vs groups[1] (y), with Spearman rho.
    """
 
    labels = labels or groups
    means = or_means(adata, scores, groups, group_key, or_key, min_cells)
 
    fig, axes = plt.subplots(1, len(scores), figsize=(4.5 * len(scores), 4.5), squeeze=False)
 
    for ax, score in zip(axes[0], scores):
        df = means[score].unstack(group_key)
        lo, hi = df.values.min(), df.values.max()
        pad = 0.05 * (hi - lo)
 
        ax.scatter(df[groups[0]], df[groups[1]], s=12, alpha=0.5, lw=0)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "--", color="gray", zorder=0)
        ax.set_xlim(lo - pad, hi + pad)
        ax.set_ylim(lo - pad, hi + pad)
        ax.text(
            0.98, 0.02,
            rf"$\rho$:{spearmanr(df[groups[0]], df[groups[1]])[0]:.3f}",
            transform=ax.transAxes, ha="right", va="bottom"
        )
        ax.set_title(score.replace("_", " "))
        ax.set_xlabel(labels[0])
 
    axes[0, 0].set_ylabel(labels[1])
    fig.suptitle(title or "")
    sns.despine()
 
    plt.savefig(
        os.path.join(results_dir, filename),
        bbox_inches="tight",
        dpi=300
    )
 
    plt.close()
 
 
def plot_score_pairwise(
    adata,
    groups,
    score="ES_score",
    group_key="condition",
    or_key="top_Olfr",
    min_cells=10,
    results_dir=RESULTS_DIR,
    filename="score_pairwise.png"
):
    """
    Lower-triangle grid of per-OR mean score scatters for every pair of the given groups (2 or more).
    """
 
    df = or_means(adata, [score], groups, group_key, or_key, min_cells)[score].unstack(group_key)
    lo, hi = df.values.min(), df.values.max()
    pad = 0.05 * (hi - lo)
    n = len(groups) - 1
 
    fig, axes = plt.subplots(n, n, figsize=(3.6 * n, 3.6 * n), sharex=True, sharey=True, squeeze=False)
 
    for i in range(n):
        for j in range(n):
            ax = axes[i, j]
 
            if j > i:
                ax.axis("off")
                continue
 
            x, y = groups[j], groups[i + 1]
 
            ax.scatter(df[x], df[y], s=8, alpha=0.5, lw=0)
            ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "--", color="gray", zorder=0)
            ax.set_xlim(lo - pad, hi + pad)
            ax.set_ylim(lo - pad, hi + pad)
            ax.text(0.98, 0.02, rf"$\rho$:{spearmanr(df[x], df[y])[0]:.3f}", transform=ax.transAxes, ha="right", va="bottom")
 
            if j == 0:
                ax.set_ylabel(y)
 
            if i == n - 1:
                ax.set_xlabel(x)
 
    fig.suptitle(score.replace("_", " "), y=0.92)
    sns.despine()
 
    plt.savefig(
        os.path.join(results_dir, filename),
        bbox_inches="tight",
        dpi=300
    )
 
    plt.close()
 
 
def plot_gep_deltas(
    adata,
    groups,
    geps=GEPS_NAMED,
    labels=None,
    title=None,
    ylim=70,
    group_key="condition",
    or_key="top_Olfr",
    min_cells=10,
    results_dir=RESULTS_DIR,
    filename="gep_deltas.png"
):
    """
    Boxen plot of the per-OR change in mean GEP usage (groups[1] minus groups[0]).
    """
 
    labels = labels or groups
    means = or_means(adata, geps, groups, group_key, or_key, min_cells)
    delta = means.xs(groups[1], level=group_key) - means.xs(groups[0], level=group_key)
    delta.columns = [gep_label(g) for g in geps]
 
    fig, ax = plt.subplots(figsize=(0.6 * len(geps) + 1.5, 5))
 
    sns.boxenplot(
        data=delta,
        color="tab:blue",
        k_depth="full",
        linewidth=0,
        showfliers=False,
        line_kws=dict(linestyle="--", linewidth=0.6, color="k"),
        ax=ax
    )
    ax.axhline(0, ls="--", color="gray")
    ax.set_ylim(-ylim, ylim)
    ax.set_ylabel(f"Δ GEP usage ({labels[1]} vs {labels[0]})")
    ax.tick_params(axis="x", rotation=90)
    ax.set_title(title or "")
    sns.despine()
 
    plt.savefig(
        os.path.join(results_dir, filename),
        bbox_inches="tight",
        dpi=300
    )
 
    plt.close()
 