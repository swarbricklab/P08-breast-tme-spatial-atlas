#!/usr/bin/env python3

# --- Aim
# Import visium data into anndata object & and do some (very basic) qc
# Save results as .h5ad file

# --- Define environment
# env: preprocessing

# --- Load packages
import argparse
import scanpy as sc
import pandas as pd
import anndata as ad

import numpy as np
import os

import matplotlib.pyplot as plt

sc.settings.verbosity = 3

# --- Argparse arguments
parser = argparse.ArgumentParser(description="Import Visium data & save as .h5ad file",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# Required
parser.add_argument("--raw_data_dir", help="Path to raw data")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--project_dir", help="Path to the root of the project directory")
parser.add_argument("--min_count", help="Minimal number of counts per spot", type=int)
parser.add_argument("--min_gene", help="Minimal number of genes per spot", type=int)
args = vars(parser.parse_args())

# --- Arguments
raw_data_dir = args["raw_data_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
project_dir=args["project_dir"]
min_count=args["min_count"]
min_gene=args["min_gene"]

# Make directories if they don't exist yet
os.makedirs(h5ad_dir, exist_ok=True)

# --- Extract some info from the sample sheet
sample_dir = os.path.join(project_dir, 'config')
sample_sheet_file = 'sample_sheet.csv'
sample_sheet = pd.read_csv(os.path.join(sample_dir, sample_sheet_file))
external_id = sample_sheet[sample_sheet['section_id'] == sample_id]['external_id'].values[0]
metadata_dict = {}
metadata_dict['section_id'] = sample_id
metadata_dict['external_id'] = external_id
metadata_dict['donor_id'] = sample_sheet[sample_sheet['section_id'] == sample_id]['donor_id'].values[0].astype(str)
metadata_dict['vis_exp_id'] = sample_sheet[sample_sheet['section_id'] == sample_id]['vis_exp_id'].values[0]

# --- Specify the path to spaceranger outs
joined_id = f'{external_id}-{sample_id}'
outs_path = os.path.join(raw_data_dir, joined_id, 'outs')
v2_spatial_file = outs_path + '/spatial/tissue_positions.csv'
v1_spatial_file = outs_path + '/spatial/tissue_positions_list.csv'
if os.path.exists(v1_spatial_file) == False:
    v2_spatial_df = pd.read_csv(v2_spatial_file)
    v2_spatial_df.to_csv(path_or_buf = v1_spatial_file, header=None,index=False)
adata = sc.read_visium(
    path=outs_path,
    source_image_path=os.path.join(outs_path, 'spatial'),
    library_id=sample_id
)
adata.var_names_make_unique()
adata.uns['metadata'] = metadata_dict
# --- Initial QC
adata.var['mt'] = adata.var_names.str.startswith(('MT-'))
print(f"Sample {sample_id}: Mitochondrial genes:", adata.var['mt'].sum())
adata.var['hb'] = adata.var_names.str.contains(("^HB[AB]"))
print(f"Sample {sample_id}: Hemoglobin genes:", adata.var['hb'].sum())
sc.pp.calculate_qc_metrics(
    adata, qc_vars=['mt','hb'], 
    percent_top=None, log1p=False, inplace=True)
# --- Flag spots with low QC values (but don't filter)
adata.obs['low_qc'] = np.where((adata.obs['total_counts'] <= min_count) | (adata.obs['n_genes_by_counts'] <= min_gene), 'low counts/genes', 'pass')
adata.obs['low_qc'] = pd.Categorical(adata.obs['low_qc'], categories=['pass', 'low counts/genes'], ordered=True)
adata.obs['low_counts'] = np.where(adata.obs['total_counts'] <= min_count, 'low_counts', 'pass')
adata.obs['low_counts'] = pd.Categorical(adata.obs['low_counts'], categories=['pass', 'low_counts'], ordered=True)
adata.obs['low_gene'] = np.where(adata.obs['n_genes_by_counts'] <= min_gene, 'low_gene', 'pass')
adata.obs['low_gene'] = pd.Categorical(adata.obs['low_gene'], categories=['pass', 'low_gene'], ordered=True)
adata.obs['no_counts'] = np.where(adata.obs.total_counts == 0.0, 'no counts', 'pass')
adata.obs['no_counts'] = pd.Categorical(adata.obs['total_counts'], categories=['pass', 'no counts'], ordered=True)
# --- For convenience: add the cytassist image to the object (if there is one)
cytassist_path = os.path.join(outs_path, 'spatial/cytassist_image.tiff')
if os.path.exists(cytassist_path):
    img = plt.imread(os.path.join(outs_path, 'spatial/cytassist_image.tiff'))
    adata.uns['spatial'][sample_id]['images']['cytassist'] = img
    adata.uns['spatial'][sample_id]['scalefactors']['tissue_cytassist_scalef'] = 0.042
# --- Save output
outs_dir = os.path.join(h5ad_dir, 'raw')
os.makedirs(outs_dir, exist_ok=True)
print(f"Saving .h5ad files to {outs_dir}")
ad.AnnData.write(adata, filename=os.path.join(outs_dir, f"{joined_id}_raw.h5ad"))
# --- Also make an object from the raw feature bc matrix
adata = sc.read_visium(
    path=outs_path,
    count_file='raw_feature_bc_matrix.h5',
    source_image_path=os.path.join(outs_path, 'spatial'),
    library_id=sample_id
)
adata.var_names_make_unique()
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
os.makedirs(os.path.join(h5ad_dir, "raw_feature_bc_matrix"), exist_ok=True)
print(f"Saving .h5ad files to {h5ad_dir}")
ad.AnnData.write(adata, filename=os.path.join(h5ad_dir, "raw_feature_bc_matrix", f"{joined_id}_raw_feature_bc_matrix.h5ad"))
print("Script succesfully completed")
