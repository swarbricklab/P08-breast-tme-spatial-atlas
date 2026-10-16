#!/usr/bin/env python3

# --- Aim
# Import xenium data into anndata object & and do some basic qc
# Save results as .h5ad file


# --- Define environment
# env: xen-hev

# --- Load packages
import argparse
import os

import scanpy as sc
import pandas as pd
import numpy as np


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Xenium data & save as .h5ad file",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--project_dir", help="Path to the root of the project directory")
parser.add_argument("--sample_sheet", help="Path to sample sheet csv file")
parser.add_argument("--res_dir", help="Path to results directory")

args = vars(parser.parse_args())


# --- Arguments
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
res_dir=args["res_dir"]
project_dir=args["project_dir"]

# Print the arguments
print(f'.h5ad dir: {h5ad_dir}')
print(f'Sample id: {sample_id}')

# Make directories if they don't exist yet
os.makedirs(h5ad_dir, exist_ok=True)


# --- Extract some info from the sample sheet

# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_xenium.csv')

sample_sheet = pd.read_csv(sample_sheet_path)


# --- Load data

raw_data_dir_sample = os.path.join(project_dir, 'data', 'processed', 'xenium_outs', sample_id)

# Load the cell feature matrix as an anndata object

cell_feature_matrix = 'cell_feature_matrix.h5'

adata = sc.read_10x_h5(
    filename=os.path.join(raw_data_dir_sample, cell_feature_matrix)
)


# -- Load cell info

# Note: make sure pyarrow is installed in the environment to read parquet files (e.g. pip install pyarrow)
cells_file = 'cells.parquet'
cells_df = pd.read_parquet(os.path.join(raw_data_dir_sample, cells_file))
cells_df.set_index(adata.obs_names, inplace=True)

# Add to anndata obs
adata.obs = cells_df.copy()


# --- Copy spatial coordinates to obsm['spatial']

adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].copy().to_numpy()


# --- Basic qc

# Percentage of control probes and control codewords
cprobes = (
    adata.obs["control_probe_counts"].sum() / adata.obs["total_counts"].sum() * 100
)
cwords = (
    adata.obs["control_codeword_counts"].sum() / adata.obs["total_counts"].sum() * 100
)
print(f"Negative DNA probe count % : {cprobes}")
print(f"Negative decoding count % : {cwords}")


# --- Filter low quality cells (low counts) before log normalization

min_count=10
adata.obs['low_qc'] = np.where((adata.obs['total_counts'] <= min_count), 'low counts', 'pass')

n_low_qc = sum(adata.obs['low_qc'] == 'low counts')

print(f'Low quality cells: {n_low_qc}/{len(adata.obs_names)} ({round(n_low_qc/len(adata.obs_names)*100, 2)}%)')

# Filter before log normalization
adata.layers['filtered_counts'] = adata.X.copy()
sc.pp.filter_cells(adata, min_counts=min_count)

# Filter genes expressed in fewer than 5 cells
sc.pp.filter_genes(adata, min_cells=5)


# --- Normalize
# layers['filtered_counts']: counts scaled to 10,000 per cell, used for bin-level
# sums (TLS zonation, chemokine densities)
sc.pp.normalize_total(adata, layer='filtered_counts', target_sum=1e4, inplace=True)

# X: log1p of the counts.
sc.pp.log1p(adata)


# --- Save the anndata object

# save AnnData object to *.h5ad file
outs_dir = os.path.join(h5ad_dir, 'xenium')

os.makedirs(outs_dir, exist_ok=True)

print(f"Saving .h5ad files to {outs_dir}")

# save processed data object
adata.write(filename = os.path.join(outs_dir, f"{sample_id}.h5ad"), compression='gzip')

print("Script succesfully completed")