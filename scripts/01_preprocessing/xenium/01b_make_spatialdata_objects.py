#!/usr/bin/env python3

# --- Aim
# Use spatialdata package to import xenium data into anndata objects with segmentation mask & transcript locations

# --- Define environment
# env: spatdat

# --- Load packages
import argparse
import os
import shutil

import scanpy as sc
import pandas as pd
import numpy as np
from spatialdata_io import xenium

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
project_dir = args["project_dir"]
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
sdata = xenium(path=raw_data_dir_sample, n_jobs=8, cells_as_circles=False, morphology_mip=False, morphology_focus=False, transcripts=True)
# --- Update the ids of the anndata object to match the cell_id
sdata.tables['table'].obs.index = sdata.shapes['cell_boundaries'].index
sdata.tables['table'].obs.index.name = None
# --- Basic qc
cprobes = (
    sdata.tables['table'].obs["control_probe_counts"].sum() / sdata.tables['table'].obs["total_counts"].sum() * 100
)
cwords = (
    sdata.tables['table'].obs["control_codeword_counts"].sum() / sdata.tables['table'].obs["total_counts"].sum() * 100
)
print(f"Negative DNA probe count % : {cprobes}")
print(f"Negative decoding count % : {cwords}")
min_count=10
sdata.tables['table'].obs['low_qc'] = np.where((sdata.tables['table'].obs['total_counts'] <= min_count), 'low counts', 'pass')
n_low_qc = sum(sdata.tables['table'].obs['low_qc'] == 'low counts')
print(f'Low quality cells: {n_low_qc}/{len(sdata.tables['table'].obs_names)} ({round(n_low_qc/len(sdata.tables['table'].obs_names)*100, 2)}%)')
sdata.tables['table'].layers['filtered_counts'] = sdata.tables['table'].X
sc.pp.filter_cells(sdata.tables['table'], min_counts=min_count)
sc.pp.filter_genes(sdata.tables['table'], min_cells=5)
sc.pp.normalize_total(sdata.tables['table'], layer='filtered_counts', target_sum=1e4, inplace=True)
sc.pp.log1p(sdata.tables['table'])
# --- Save the spatialdata object
outs_dir = os.path.join(project_dir, 'data', 'processed', 'spatialdata_objects', 'log_norm')
os.makedirs(outs_dir, exist_ok=True)
print(f"Saving .h5ad files to {outs_dir}")
sdata.write(os.path.join(outs_dir, f"{sample_id}.zarr"))
shutil.make_archive(os.path.join(outs_dir, f"{sample_id}"), 'zip', os.path.join(outs_dir, f"{sample_id}.zarr"))
shutil.rmtree(os.path.join(outs_dir, f"{sample_id}.zarr"))
print("Script succesfully completed")
