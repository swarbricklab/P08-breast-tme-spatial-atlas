#!/usr/bin/env python3

# --- Aim
# Import transcript-level xenium data into anndata object
# To be able to plot the location of transcripts in space

# --- Define environment
# env: xen-hev

# --- Load packages
import argparse
import os
import anndata as ad

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
parser.add_argument("--min_count", help="Minimal number of counts per spot", type=int)
parser.add_argument("--min_gene", help="Minimal number of genes per spot", type=int)
parser.add_argument("--outs_path", help="Path to spaceranger 'outs'")
parser.add_argument("--use_data", help="Which data to use (log_norm or SCT)")
args = vars(parser.parse_args())
# --- Arguments
raw_data_dir = args["raw_data_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
min_count=args["min_count"]
min_gene=args["min_gene"]
outs_path=args["outs_path"]
# Make directories if they don't exist yet
os.makedirs(h5ad_dir, exist_ok=True)
# --- Load data
raw_data_dir_sample = os.path.join(raw_data_dir, sample_id)
transcripts_file = 'transcripts.parquet'
transcripts_df = pd.read_parquet(os.path.join(raw_data_dir_sample, transcripts_file))
transcripts_of_interest = ['CCL19', 'CCL21', 'TSPAN7', 'MEOX2', 'RAMP3', 'SPNS2', 'IL33']
transcripts_sub = transcripts_df[transcripts_df['feature_name'].isin(transcripts_of_interest)].copy()
transcripts_sub.set_index('transcript_id', inplace=True)
adata_tx = ad.AnnData(
    X=None,
    obs=transcripts_sub,
    var=pd.DataFrame(index=["gene_info"])
)
adata_tx.obsm["spatial"] = adata_tx.obs[["x_location", "y_location"]].copy().to_numpy()
# --- Save the anndata object
outs_dir = os.path.join(h5ad_dir, 'transcript_objects')
os.makedirs(outs_dir, exist_ok=True)
print(f"Saving .h5ad files to {outs_dir}")
adata_tx.write(filename = os.path.join(outs_dir, f"{sample_id}_HEV_transcripts.h5ad"), compression = "gzip")
print("Script succesfully completed")
