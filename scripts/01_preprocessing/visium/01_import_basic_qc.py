#!/usr/bin/env python3

# --- Aim
# Import visium data into anndata object & and do some (very basic) qc
# Save results as .h5ad file

# --- Define environment
# env: spatial-pilot

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
parser.add_argument("--sample_sheet", help="Path to sample sheet csv file")
parser.add_argument("--res_dir", help="Path to results directory")

args = vars(parser.parse_args())


# --- Arguments
raw_data_dir = args["raw_data_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
project_dir=args["project_dir"]
res_dir=args["res_dir"]


# Make directories if they don't exist yet
os.makedirs(h5ad_dir, exist_ok=True)


# --- Extract some info from the sample sheet

# Extract the section_id
# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_dir = os.path.join(project_dir, 'config')
    sample_sheet_file = 'sample_sheet_visium.csv'
    sample_sheet_path = os.path.join(sample_dir, sample_sheet_file)


# --- Specify the path to spaceranger outs

outs_path = os.path.join(project_dir, 'data', 'processed', 'spaceranger_outs', 'visium', sample_id)

# -- Load data

v2_spatial_file = outs_path + '/spatial/tissue_positions.csv'
v1_spatial_file = outs_path + '/spatial/tissue_positions_list.csv'

# if V1 spatial position file does not exist, create a new one
created_v1_file = False
if os.path.exists(v1_spatial_file) == False:
    v2_spatial_df = pd.read_csv(v2_spatial_file)
    v2_spatial_df.to_csv(path_or_buf = v1_spatial_file, header=None,index=False)
    created_v1_file = True

adata = sc.read_visium(
    path=outs_path,
    source_image_path=os.path.join(outs_path, 'spatial'),
    library_id=sample_id
)

adata.var_names_make_unique()

# Clean up the temporary V1 spatial file if we created it
if created_v1_file:
    os.remove(v1_spatial_file)


# --- Initial QC

# Annotate the group of mitochondrial genes as 'mt'
adata.var['mt'] = adata.var_names.str.startswith(('MT-'))
print(f"Sample {sample_id}: Mitochondrial genes:", adata.var['mt'].sum())

# Annotate hemoglobin genes as hb
adata.var['hb'] = adata.var_names.str.contains(("^HB[AB]"))
print(f"Sample {sample_id}: Hemoglobin genes:", adata.var['hb'].sum())

# There are no ribosomal genes in the probe set

sc.pp.calculate_qc_metrics(
    adata, qc_vars=['mt','hb'], 
    percent_top=None, log1p=False, inplace=True)


# --- Flag spots with low QC values (but don't filter)

min_count = 1000
min_gene = 200

adata.obs['low_qc'] = np.where((adata.obs['total_counts'] <= min_count) | (adata.obs['n_genes_by_counts'] <= min_gene), 'low counts/genes', 'pass')

# Need to do this to make sure the pass isn't colored red if there is none below pass
adata.obs['low_qc'] = pd.Categorical(adata.obs['low_qc'], categories=['pass', 'low counts/genes'], ordered=True)

# Flag spots with low counts
adata.obs['low_counts'] = np.where(adata.obs['total_counts'] <= min_count, 'low_counts', 'pass')
adata.obs['low_counts'] = pd.Categorical(adata.obs['low_counts'], categories=['pass', 'low_counts'], ordered=True)

# Flag spots with low genes
adata.obs['low_gene'] = np.where(adata.obs['n_genes_by_counts'] <= min_gene, 'low_gene', 'pass')
adata.obs['low_gene'] = pd.Categorical(adata.obs['low_gene'], categories=['pass', 'low_gene'], ordered=True)

# Flag spots with 0 counts
adata.obs['no_counts'] = np.where(adata.obs.total_counts == 0.0, 'no counts', 'pass')
adata.obs['no_counts'] = pd.Categorical(adata.obs['total_counts'], categories=['pass', 'no counts'], ordered=True)


# --- Save output

outs_dir = os.path.join(h5ad_dir, 'visium', 'raw')

os.makedirs(outs_dir, exist_ok=True)

print(f"Saving .h5ad files to {outs_dir}")

ad.AnnData.write(adata, filename=os.path.join(outs_dir, f"{sample_id}_raw.h5ad"))


print("Script succesfully completed")