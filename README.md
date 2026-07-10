# MOE Single-Cell RNA-seq Analysis Pipeline

A modular Python pipeline for analyzing single-cell RNA-seq data from the mammalian olfactory epithelium (MOE), focusing on mature olfactory sensory neurons (mOSNs). 

## Overview

This repository contains analysis tools compiled for MOE scRNA-seq.

## Repository Structure

```
single-cell-RNAseq-Analysis/
│
├── notebooks/
│   └── singleCellAnalysis_masterNotebook.ipynb
│
├── src/
│   ├── __init__.py
│   ├── io.py
│   ├── preprocessing.py
│   ├── dim_reduction.py
│   ├── expression.py
│   ├── scoring.py
│   ├── markers.py
│   ├── gep.py
│   └── plotting.py
│
├── data/
├── results/
├── requirements.txt
└── README.md
```

## Installation

Clone the repository:

```bash
git clone https://github.com/USERNAME/single-cell-RNAseq-Analysis.git
cd single-cell-RNAseq-Analysis
```

Create a Python environment:

```bash
conda create -n scrna_analysis python=3.11
conda activate scrna_analysis
```

Install dependencies:

```bash
pip install -r requirements.txt
```

## Data

Sequencing data are not included in this repository. Data directory is provided. 

## Usage

The main analysis workflow is contained in:

```
notebooks/singleCellAnalysis_masterNotebook.ipynb
```

The notebook calls reusable functions from the `src/` directory.

Example:

```python
from src.io import load_data
from src.expression import calculate_gene_expression
from src.plotting import plot_expression

data = load_data("data/raw/example.parquet")

expression = calculate_gene_expression(data)

plot_expression(
    expression,
    genes=["Gnal", "Omp"]
)
```
## Citation

Based on analyses developed in:

Tsukahara T, Brann DH, Pashkovski SL, Guitchounts G, Bozza T, Datta SR. (2021). A transcriptional rheostat couples past activity to future sensory responses. *Cell*, 184(26), 6326–6343.e32. https://doi.org/10.1016/j.cell.2021.11.022

