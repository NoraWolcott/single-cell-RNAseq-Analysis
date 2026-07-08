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