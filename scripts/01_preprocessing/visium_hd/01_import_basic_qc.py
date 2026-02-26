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

args = vars(parser.parse_args())


# --- Arguments
raw_data_dir = args["raw_data_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
project_dir=args["project_dir"]


# Make directories if they don't exist yet
os.makedirs(h5ad_dir, exist_ok=True)


# --- Extract some info from the sample sheet

# Extract the section_id
sample_dir = os.path.join(project_dir, 'config')
sample_sheet_file = 'sample_sheet.csv'

sample_sheet = pd.read_csv(os.path.join(sample_dir, sample_sheet_file))
external_id = sample_sheet[sample_sheet['section_id'] == sample_id]['external_id'].values[0]

# Create dictionary to store metadata
metadata_dict = {}

metadata_dict['section_id'] = sample_id
metadata_dict['external_id'] = external_id

# Patient id
metadata_dict['donor_id'] = sample_sheet[sample_sheet['section_id'] == sample_id]['donor_id'].values[0].astype(str)

# Visium experiment id
metadata_dict['vis_exp_id'] = sample_sheet[sample_sheet['section_id'] == sample_id]['vis_exp_id'].values[0]

# Visium experiment id
metadata_dict['slide_type'] = sample_sheet[sample_sheet['section_id'] == sample_id]['slide_type'].values[0]


# --- Specify the path to spaceranger outs

outs_path = os.path.join(raw_data_dir, sample_id, 'outs')

# Deal with visium HD resolutions
binned_outputs_path = os.path.join(outs_path, 'binned_outputs')

if os.path.exists(binned_outputs_path):
    contents = os.listdir(binned_outputs_path)
    resolutions = [item for item in contents if os.path.isdir(os.path.join(binned_outputs_path, item))]

# -- Load data per resolution

for resolution in resolutions:
    outs_path_res = os.path.join(binned_outputs_path, resolution)

    # Make a fake v1 spatial file so scanpy doesn't complain
    hd_spatial_file = outs_path_res + '/spatial/tissue_positions.parquet'
    v1_spatial_file = outs_path_res + '/spatial/tissue_positions_list.csv'

    # if V1 spatial position file does not exist, create a new one
    if os.path.exists(v1_spatial_file) == False:
        v2_spatial_df = pd.read_parquet(hd_spatial_file)
        v2_spatial_df.to_csv(path_or_buf = v1_spatial_file, header=None,index=False)

    adata = sc.read_visium(
        path=outs_path_res,
        source_image_path=os.path.join(outs_path_res, 'spatial'),
        library_id=sample_id
    )

    adata.var_names_make_unique()

    # Add the metadata to the anndata object for convenience
    adata.uns['metadata'] = metadata_dict


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

    outs_dir = os.path.join(h5ad_dir, resolution, 'raw')

    os.makedirs(outs_dir, exist_ok=True)

    print(f"Saving .h5ad files to {outs_dir}")

    ad.AnnData.write(adata, filename=os.path.join(outs_dir, f"{sample_id}_raw.h5ad"))


    # --- Also make an object from the raw feature bc matrix
    # Has spots that were detected to not be 'on tissue'
    # Good to check for diffusion (counts in these spots should be low)

    adata = sc.read_visium(
        path=outs_path_res,
        count_file='raw_feature_bc_matrix.h5',
        source_image_path=os.path.join(outs_path_res, 'spatial'),
        library_id=sample_id
    )

    adata.var_names_make_unique()

    # Calculate qc metrics
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)


    # --- For convenience: add the cytassist image to the object (if there is one)

    cytassist_path = os.path.join(outs_path_res, 'spatial/cytassist_image.tiff')

    if os.path.exists(cytassist_path):
        # Read in the cytassist image
        img = plt.imread(os.path.join(outs_path_res, 'spatial/cytassist_image.tiff'))
        adata.uns['spatial'][sample_id]['images']['cytassist'] = img # Save the cytassist image
        adata.uns['spatial'][sample_id]['scalefactors']['tissue_cytassist_scalef'] = 0.042 # Add scale factors for plotting with scanpy (note: this is an approximation, not exact)

    # --- Save output

    os.makedirs(os.path.join(h5ad_dir, resolution, 'raw_feature_bc_matrix'), exist_ok=True)
    print(f"Saving .h5ad files to {h5ad_dir}")

    ad.AnnData.write(adata, filename=os.path.join(h5ad_dir, resolution, "raw_feature_bc_matrix", f"{sample_id}_raw_feature_bc_matrix.h5ad"))

    # Remove the fake v1 spatial file after making it (so that dvc doesn't think the data has changed)
    if os.path.exists(v1_spatial_file):
        os.remove(v1_spatial_file)
        print(f"File '{v1_spatial_file}' deleted successfully.")


print("Script succesfully completed")