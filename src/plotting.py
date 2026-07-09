import scanpy as sc


def plot_marker_umap(
    adata,
    genes,
    cmap="viridis_r",
    size=2,
    vmax="p99.9",
    ncols=4,
    layer="counts",
    legend_loc="on data"
):
    """
    Plot selected marker genes on UMAP.

    Parameters
    ----------
    adata : AnnData
        Annotated single-cell dataset.
    genes : list
        Genes to visualize.
    cmap : str
        Colormap.
    size : float
        Marker size.
    vmax : str or float
        Upper expression limit for scaling.
    ncols : int
        Number of columns in multi-panel plot.
    layer : str
        AnnData layer to plot from.
    legend_loc : str
        UMAP legend location.
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
        layer=layer
    )