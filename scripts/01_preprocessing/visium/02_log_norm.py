#!/usr/bin/env python3

# --- Aim
# Use scanpy's sc.pp.normalize_total & sc.pp.log1p to normalize the data


# --- Define environment
# env: preprocessing

# --- Load packages
import argparse
import os

import scanpy as sc
import pandas as pd
import anndata as ad

import numpy as np

sc.settings.verbosity = 3

# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium data & save as .h5ad file",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required
parser.add_argument("--raw_data_dir", help="Path to raw data")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--project_dir", help="Path to the root of the project directory")
parser.add_argument("--sample_sheet", help="Path to sample sheet csv file")
parser.add_argument("--res_dir", help="Path to results directory")

args = vars(parser.parse_args())


# --- Arguments
raw_data_dir = args["raw_data_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
project_dir=args["project_dir"]
res_dir=args["res_dir"]


# --- Extract some info from the sample sheet

# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_dir = os.path.join(project_dir, 'config')
    sample_sheet_file = 'sample_sheet_visium.csv'
    sample_sheet_path = os.path.join(sample_dir, sample_sheet_file)

sample_sheet = pd.read_csv(sample_sheet_path)


# --- Load the raw anndata object

input_adata_dir = os.path.join(h5ad_dir, 'visium', 'raw')
adata = ad.read(os.path.join(input_adata_dir, f'{sample_id}_raw.h5ad'))


# --- Do some very conservative filtering: only spots with 0 counts

sc.pp.filter_cells(adata, min_counts=1)

adata.layers['filtered_counts'] = adata.X


# --- Normalization

adata.layers['norm'] = sc.pp.normalize_total(adata, layer='filtered_counts', target_sum=1e4, inplace=False)['X']

# log transform
adata.layers['log_norm'] = sc.pp.log1p(adata.layers['norm'], copy=True)

# Remove intermediate layers to save space (optional)
# Uncomment these lines if you want to reduce file size:
del adata.layers['norm']
del adata.layers['filtered_counts']

# Scale gene expression to unit variance (clip values that exceed 10 standard deviations)
adata.layers['scaled'] = sc.pp.scale(adata.layers['log_norm'], max_value=10, copy=True)

# Detect highly variable genes as well
sc.pp.highly_variable_genes(adata, layer='log_norm', flavor='seurat', n_top_genes=2000)


# --- PCA, UMAP + clustering
# Optional

# seed = 42

# # create a new object
# pca_obj = ad.AnnData(adata.layers['scaled'])
# pca_obj.var = adata.var
# sc.tl.pca(pca_obj, n_comps = 50, zero_center = True, svd_solver='arpack', random_state = seed, return_info = True, use_highly_variable = True)

# adata.uns['pca'] = pca_obj.uns['pca']
# adata.obsm['X_pca'] = pca_obj.obsm['X_pca']
# adata.varm['PCs'] = pca_obj.varm['PCs']

# # neighbourhood graph
# sc.pp.neighbors(adata, n_pcs=50, random_state=seed)

# # calculate umap
# sc.tl.umap(adata, random_state=seed)

# # Leiden clustering
# sc.tl.leiden(adata, key_added='leiden', resolution=1)


# --- Save output

outs_dir = os.path.join(h5ad_dir, 'visium')

os.makedirs(outs_dir, exist_ok=True)

print(f"Saving .h5ad files to {outs_dir}")

ad.AnnData.write(adata, filename=os.path.join(outs_dir, f"{sample_id}.h5ad"), compression='gzip')

print("Script succesfully completed")
