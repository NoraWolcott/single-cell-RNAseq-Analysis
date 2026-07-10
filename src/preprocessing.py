import scanpy as sc


def store_raw_counts(adata):
    """
    Store raw counts before normalization.
    """
    adata.layers["counts"] = adata.X.copy()
    return adata


def normalize_data(adata, target_sum=1e4):
    """
    Normalize counts and log transform.
    """
    sc.pp.normalize_total(
        adata,
        target_sum=target_sum
    )

    sc.pp.log1p(adata)

    return adata


def identify_hvgs(
    adata,
    n_top_genes=3000,
    exclude_prefixes=("mt-", "Rpl", "Rps", "Olfr", "Gm"),
    exclude_suffixes=("Rik",)
):
    """
    Identify highly variable genes, excluding riobosomal, mitochondrial, and OR genes.
    """

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        subset=False
    )

    genes_to_exclude = (
        adata.var_names.str.startswith(exclude_prefixes)
    )

    for suffix in exclude_suffixes:
        genes_to_exclude |= adata.var_names.str.endswith(suffix)

    adata.var.loc[
        genes_to_exclude,
        "highly_variable"
    ] = False

    return adata


def filter_cells_by_qc(
    adata,
    min_counts=1000,
    max_pct_mt=10
):
    """
    Filter cells based on total counts and mitochondrial gene expression.
    """

    # Identify mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("mt-")

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True
    )

    # Filter
    adata = adata[
        (adata.obs["total_counts"] > min_counts) &
        (adata.obs["pct_counts_mt"] < max_pct_mt)
    ].copy()

    return adata


def normalize_and_select_hvgs(
    adata,
    n_top_genes=3000,
    target_sum=1e4,
    exclude_prefixes=None,
    exclude_suffixes=None
):
    """
    Normalize counts, log transform, and identify highly variable genes.
    """

    if exclude_prefixes is None:
        exclude_prefixes = [
            "mt-",
            "Rpl",
            "Rps",
            "Olfr",
            "Gm"
        ]

    if exclude_suffixes is None:
        exclude_suffixes = [
            "Rik"
        ]

    # Store raw counts
    adata.layers["counts"] = adata.X.copy()

    # Normalize and log transform
    sc.pp.normalize_total(
        adata,
        target_sum=target_sum
    )

    sc.pp.log1p(adata)

    # Identify genes to exclude
    genes_to_exclude = False

    for prefix in exclude_prefixes:
        genes_to_exclude |= adata.var_names.str.startswith(prefix)

    for suffix in exclude_suffixes:
        genes_to_exclude |= adata.var_names.str.endswith(suffix)

    # Find HVGs
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        subset=False
    )

    # Remove excluded genes
    adata.var.loc[
        genes_to_exclude,
        "highly_variable"
    ] = False

    return adata