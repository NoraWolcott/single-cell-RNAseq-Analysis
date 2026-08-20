import os
import glob
import scanpy as sc

def load_samples(parent_dir, samples):

    adata_list = []

    for sample in samples:

        sample_dir = os.path.join(parent_dir, sample)

        # Find h5ad and h5 files
        h5ad_files = glob.glob(os.path.join(sample_dir, "*.h5ad"))
        h5_files = glob.glob(os.path.join(sample_dir, "*.h5"))

        all_files = h5ad_files + h5_files

        if len(all_files) == 0:
            raise FileNotFoundError(
                f"No .h5ad or .h5 file found in {sample_dir}"
            )

        if len(all_files) > 1:
            print(
                f"Warning: multiple .h5ad/.h5 files found in {sample_dir}. "
                f"Using {all_files[0]}"
            )

        file_path = all_files[0]

        print(f"Loading {file_path}")

        # Load the appropriate file type
        if file_path.endswith(".h5ad"):
            adata = sc.read_h5ad(file_path)

        elif file_path.endswith(".h5"):
            adata = sc.read_10x_h5(file_path)

        else:
            raise ValueError(
                f"Unsupported file type: {file_path}"
            )

        # Make gene names unique
        adata.var_names_make_unique()

        # Make cell/barcode names unique within each sample
        adata.obs_names_make_unique()

        # Add sample identifier
        adata.obs["sample"] = sample

        print(
            f"  Loaded {adata.n_obs:,} cells × "
            f"{adata.n_vars:,} genes"
        )

        adata_list.append(adata)

    # Concatenate samples
    adata_all = sc.concat(
        adata_list,
        join="inner",
        merge="same"
    )

    return adata_all