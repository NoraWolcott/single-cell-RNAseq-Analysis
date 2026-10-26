import re
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import non_negative_factorization

URL = "https://raw.githubusercontent.com/dattalab/Tsukahara_Brann_OSN/main/data/tables/"
FILES = ["GSE173947_Table_S2_nmf_loadings.csv", "hvg_std_for_cnmf.csv"]
OR_URL = "https://raw.githubusercontent.com/dattalab/Brann_olfactory_dorsoventral/main/data/tables/Olfr_biomart_mm39_Ensembl_105.csv"
HIGH_EXPR = ["Malat1", "mt-Rnr1", "mt-Rnr2", "mt-Cytb", "Gm42418", "Tmsb4x", "Xist"]

GEP_NAMES = {"GEP_" + n: n for n in ["Dorsal", "Ventral", "Anterior", "Posterior", "Cd36", "High", "Low"]}

for i in range(8, 17):
    tag = "_UPR" if i == 11 else ""
    GEP_NAMES[f"GEP_{i}{tag}"] = f"GEP{i}{tag}" + ("_" if tag else "")


def load_tables(tables_dir):
    """
    Load GEP loadings and per-gene SDs (Tsukahara et al. 2021), downloading if needed.
    """

    tables_dir = Path(tables_dir).expanduser()
    tables_dir.mkdir(parents=True, exist_ok=True)

    for f in FILES:
        if not (tables_dir / f).exists():
            urllib.request.urlretrieve(URL + f, tables_dir / f)

    loadings = pd.read_csv(tables_dir / FILES[0], index_col=0)
    sd = pd.read_csv(tables_dir / FILES[1], index_col=0)

    return loadings, sd


def calculate_geps(
    adata,
    counts_layer="counts",
    tables_dir="~/.cache/gep_tables",
    min_gene_fraction=0.95,
    chunk=20000
):
    """
    Calculate the 16 GEP usages (obs columns GEP_*) and DV, AP, ES scores
    (DV = Dorsal - Ventral, AP = Anterior - Posterior, ES = High - Low)
    from raw counts by refitting NMF with the published loadings held fixed.
    """
 
    loadings, sd = load_tables(tables_dir)
    names = adata.var_names
    assert names.is_unique, "var_names must be unique gene symbols"
 
    counts = sparse.csr_matrix(adata.layers[counts_layer])
    assert np.allclose(counts.data, np.round(counts.data)), "counts layer must be raw integers"
 
    # TPT normalization (no log); per-cell total excludes ribosomal, mito, high-expression genes
    excluded = np.asarray(names.str.contains("^Rpl|^Rps|^mt-")) | np.asarray(names.isin(HIGH_EXPR))
    tot = np.asarray(counts[:, ~excluded].sum(axis=1)).ravel()
    scale = 1e4 / np.maximum(tot, 1)
 
    # loading genes present in the data
    is_hvg = np.asarray(names.isin(loadings.index))
    hvg = names[is_hvg]
    print(f"{len(hvg)}/{len(loadings)} GEP loading genes found")
    assert len(hvg) / len(loadings) >= min_gene_fraction, "too many loading genes missing; check gene names"
 
    L = loadings.loc[hvg]
    H = (L / L.sum(axis=0)).T.values
    inv_sd = sparse.diags(1 / sd.loc[hvg].values.mean(axis=1))
    X = counts[:, np.flatnonzero(is_hvg)]
 
    # usages: NMF refit per cell with the loadings fixed
    usage = np.vstack([
        non_negative_factorization(
            (sparse.diags(scale[i:i + chunk]) @ X[i:i + chunk] @ inv_sd).toarray(),
            H=H, n_components=H.shape[0], update_H=False, init="random", solver="cd",
            beta_loss="frobenius", l1_ratio=0.0, max_iter=400, tol=1e-4, random_state=0
        )[0]
        for i in range(0, X.shape[0], chunk)
    ])
 
    for name, values in zip(loadings.columns, usage.T):
        adata.obs["GEP_" + re.sub(r"\W+", "_", name.replace("GEP ", "")).strip("_")] = values
 
    adata.obs["DV_score"] = adata.obs["GEP_Dorsal"] - adata.obs["GEP_Ventral"]
    adata.obs["AP_score"] = adata.obs["GEP_Anterior"] - adata.obs["GEP_Posterior"]
    adata.obs["ES_score"] = adata.obs["GEP_High"] - adata.obs["GEP_Low"]
 
    return adata

def assign_dominant_or(
    adata,
    counts_layer="counts",
    tables_dir="~/.cache/gep_tables",
    expr_thresh=3,
    diff_thresh=2,
    chunk=20000
):
    """
    Assign each cell its dominant OR in obs["top_Olfr"] (NaN if no single OR). A cell has one OR if
    exactly one OR has >= expr_thresh UMIs and it leads the next OR by >= diff_thresh (Tsukahara et al. 2021).
    """
 
    tables_dir = Path(tables_dir).expanduser()
    tables_dir.mkdir(parents=True, exist_ok=True)
    path = tables_dir / "Olfr_biomart_mm39_Ensembl_105.csv"
 
    if not path.exists():
        urllib.request.urlretrieve(OR_URL, path)
 
    is_or = np.asarray(adata.var_names.isin(pd.read_csv(path, index_col=0).index))
    names = adata.var_names[is_or]
    counts = sparse.csr_matrix(adata.layers[counts_layer])[:, np.flatnonzero(is_or)]
    top = np.full(adata.n_obs, None, dtype=object)
 
    for i in range(0, adata.n_obs, chunk):
        D = counts[i:i + chunk].toarray()
        two = np.partition(D, -2, axis=1)[:, -2:]
        has = ((D >= expr_thresh).sum(axis=1) == 1) & (two[:, 1] - two[:, 0] >= diff_thresh)
        top[i:i + chunk][has] = names[D.argmax(axis=1)[has]]
 
    adata.obs["top_Olfr"] = pd.Categorical(top)
    print(f"{adata.obs['top_Olfr'].notna().mean():.1%} of cells have a single dominant OR")
 
    return adata

def save_cell_table(
    adata,
    results_dir,
    filename="cell_expression_geps.parquet",
    layer=None,
    min_cells=3,
    sample_key="sample",
    group_key="condition",
    cell_col="cell",
    chunk=2000
):
    """
    Save one row per cell to parquet: normalized expression of genes detected in >= min_cells cells,
    then the 16 GEPs, DV, AP, ES, top_Olfr, sample, group (and the cell barcode as cell_col; None to omit).
    Written in chunks of cells so the dense table never sits in memory.
    """
 
    import pyarrow as pa
    import pyarrow.parquet as pq
 
    assert "top_Olfr" in adata.obs, "run assign_dominant_or first"
 
    X = sparse.csr_matrix(adata.X if layer is None else adata.layers[layer])
    keep = np.flatnonzero(np.asarray((X > 0).sum(axis=0)).ravel() >= min_cells)
    X = X[:, keep]
    genes = list(adata.var_names[keep])
 
    obs = adata.obs
    numeric = {GEP_NAMES[k]: obs[k] for k in GEP_NAMES}
    numeric.update({"DV": obs["DV_score"], "AP": obs["AP_score"], "ES": obs["ES_score"]})
    text = {"top_Olfr": obs["top_Olfr"], "sample": obs[sample_key], "group": obs[group_key]}
 
    if cell_col:
        text[cell_col] = obs.index.to_series()
 
    # a GEP can share a name with a gene (Cd36), so suffix any clash
    gene_set = set(genes)
    names = [c + "_GEP" if c in gene_set else c for c in list(numeric) + list(text)]
    renamed = [n for n in names if n.endswith("_GEP")]
 
    path = Path(results_dir) / filename
    writer = None
 
    for i in range(0, adata.n_obs, chunk):
        dense = np.asfortranarray(X[i:i + chunk].toarray(), dtype=np.float32)
        arrays = [pa.array(dense[:, j]) for j in range(dense.shape[1])]
        arrays += [pa.array(v.iloc[i:i + chunk].to_numpy(np.float32)) for v in numeric.values()]
        arrays += [pa.array(v.iloc[i:i + chunk].astype(object), type=pa.string(), from_pandas=True) for v in text.values()]
        table = pa.Table.from_arrays(arrays, names=genes + names)
 
        if writer is None:
            writer = pq.ParquetWriter(path, table.schema, compression="zstd")
 
        writer.write_table(table)
 
    writer.close()
    print(f"saved {adata.n_obs:,} cells x {len(genes):,} genes (+{len(names)} annotation columns) to {path}")

 