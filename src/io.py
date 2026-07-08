import os
import glob
import scanpy as sc


def load_samples(parent_dir, samples):
    """
    Load h5ad files and concatenate.

    Parameters
    ----------
    parent_dir : str
        Directory containing one subdirectory per sample.
    samples : list[str]
        Sample names (e.g. ["D1", "D2", "P1", "P2"]).

    Returns
    -------
    adata_all : AnnData
        Concatenated AnnData object.
    """

    adata_list = []

    for sample in samples:

        sample_dir = os.path.join(parent_dir, sample)

        h5ad_files = glob.glob(
            os.path.join(sample_dir, "*.h5ad")
        ) # load in first h5ad file in directory

        if len(h5ad_files) == 0:
            raise FileNotFoundError(
                f"No .h5ad file found in {sample_dir}"
            )

        if len(h5ad_files) > 1:
            print(
                f"Warning: multiple .h5ad files found in {sample_dir}. "
                f"Using {h5ad_files[0]}"
            )

        h5ad_path = h5ad_files[0]

        print(f"Loading {h5ad_path}")

        adata = sc.read_h5ad(h5ad_path)

        adata.obs["sample"] = sample

        adata_list.append(adata)

    adata_all = sc.concat(
        adata_list,
        join="inner",
        label="sample",
        keys=samples
    )

    return adata_all