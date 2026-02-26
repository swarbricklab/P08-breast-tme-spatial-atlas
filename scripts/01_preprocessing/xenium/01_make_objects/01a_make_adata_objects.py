#!/usr/bin/env python3

# --- Aim
# Import xenium data into anndata object & and do some basic qc
# Save results as .h5ad file

# --- Define environment
# env: xen-hev

# --- Load packages
import argparse
import os
import glob
import re

import scanpy as sc
import pandas as pd
import numpy as np

# --- Argparse arguments
parser = argparse.ArgumentParser(description="Import Visium data & save as .h5ad file",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# Required
parser.add_argument("--raw_data_dir", help="Path to raw data")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--project_dir", help="Path to the root of the project directory")
parser.add_argument("--outs_path", help="Path to spaceranger 'outs'")
parser.add_argument("--use_data", help="Which data to use (log_norm or SCT)")

args = vars(parser.parse_args())

# --- Arguments
raw_data_dir = args["raw_data_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
outs_path=args["outs_path"]

# Print the arguments
print(f'Raw data dir: {raw_data_dir}')
print(f'.h5ad dir: {h5ad_dir}')
print(f'Sample id: {sample_id}')
print(f'Spatial : {outs_path}')

# Make directories if they don't exist yet
os.makedirs(h5ad_dir, exist_ok=True)

# --- Load data
raw_data_dir_sample = os.path.join(raw_data_dir, sample_id)
# Load the cell feature matrix as an anndata object
cell_feature_matrix = 'cell_feature_matrix.h5'
adata = sc.read_10x_h5(
    filename=os.path.join(raw_data_dir_sample, cell_feature_matrix)
)

# -- Load cell info
cells_file = 'cells.csv.gz'
cells_df = pd.read_csv(os.path.join(raw_data_dir_sample, cells_file))
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

# --- Mark cells as low quality (low counts) but do not filter from object
min_count=10
adata.obs['low_qc'] = np.where((adata.obs['total_counts'] <= min_count), 'low counts', 'pass')
n_low_qc = sum(adata.obs['low_qc'] == 'low counts')
print(f'Low quality cells: {n_low_qc}/{len(adata.obs_names)} ({round(n_low_qc/len(adata.obs_names)*100, 2)}%)')
# Filter before log normalization
adata.layers['filtered_counts'] = adata.X
sc.pp.filter_cells(adata, min_counts=min_count)
# Filter genes expressed in fewer than 5 cells
sc.pp.filter_genes(adata, min_cells=5)

# --- Normalize
sc.pp.normalize_total(adata, layer='filtered_counts', target_sum=1e4, inplace=True)
# log transform
sc.pp.log1p(adata)

# --- Import & add the on-board clustering results
cluster_dir = os.path.join(raw_data_dir_sample, 'analysis', 'clustering')
# List the directories here
sub_dirs = glob.glob(os.path.join(cluster_dir, '*'))
sub_dirs_base_names = []
for dir in sub_dirs:
    sub_dirs_base_names.append(os.path.basename(dir))
# Sort the list (so that it's graphclust, kmeans_2, kmeans_3 etc)
def extract_numbers(filename):
    # Use regular expression to extract numbers from the filename
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', filename)]
sub_dirs_base_names = sorted(sub_dirs_base_names, key=extract_numbers)
# -- Extract clustering resolutions (i.e. 'graphclust', 'kmeans_2')
def trim_gene_and_cluster(s):
    s = s.replace("gene_expression_", "")  # Remove "gene_" prefix
    s = s.replace("_clusters", "")  # Remove "_cluster" suffix
    return s
cluster_res = [trim_gene_and_cluster(s) for s in sub_dirs_base_names]
# -- For each resolution: read in and add to a dataframe
for n, resolution in enumerate(sub_dirs_base_names):
    if n == 0:
        cluster_res_df = pd.read_csv(os.path.join(cluster_dir, resolution, 'clusters.csv'), index_col=0)
        new_col_name = {cluster_res_df.columns[0]: cluster_res[n]}
        cluster_res_df.rename(columns=new_col_name, inplace=True)
    else:
        tmp_df = pd.read_csv(os.path.join(cluster_dir, resolution, 'clusters.csv'), index_col=0)
        new_col_name = {tmp_df.columns[0]: cluster_res[n]}
        tmp_df.rename(columns=new_col_name, inplace=True)
        cluster_res_df = pd.concat([cluster_res_df, tmp_df], axis=1)
# -- Save the clustering results in obsm
# Filter the cells that were filtered from adata
cluster_res_df_filtered = cluster_res_df.loc[adata.obs_names, :]
adata.obsm['on_board_clustering'] = cluster_res_df_filtered

# --- Save the anndata object
# save AnnData object to *.h5ad file
outs_dir = os.path.join(h5ad_dir, 'log_norm')
os.makedirs(outs_dir, exist_ok=True)
print(f"Saving .h5ad files to {outs_dir}")
# save processed data object
adata.write(filename = os.path.join(outs_dir, f"{sample_id}_log_norm.h5ad"))
print("Script succesfully completed")
