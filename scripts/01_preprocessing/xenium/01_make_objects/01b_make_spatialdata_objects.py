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

parser = argparse.ArgumentParser(description="Create Xenium SpatialData objects with segmentation masks and transcript locations",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required
parser.add_argument("--raw_data_dir", help="Path to raw data")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--project_dir", help="Path to the root of the project directory")
parser.add_argument("--outs_path", help="Path to spaceranger 'outs'")
parser.add_argument("--sample_sheet", help="Path to sample sheet csv file")
parser.add_argument("--res_dir", help="Path to results directory")


args = vars(parser.parse_args())


# --- Arguments
project_dir = args["project_dir"]
raw_data_dir = args["raw_data_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
res_dir=args["res_dir"]

outs_path=args["outs_path"]

# Make directories if they don't exist yet
os.makedirs(h5ad_dir, exist_ok=True)


# --- Load data

raw_data_dir_sample = os.path.join(project_dir, 'data', 'processed', 'xenium_outs', sample_id)

if os.path.exists(os.path.join(raw_data_dir_sample, 'cells.zarr.zip')):
    sdata = xenium(path=raw_data_dir_sample, n_jobs=8, cells_as_circles=False, morphology_mip=False, morphology_focus=False, transcripts=True)
else:
    # The GEO deposit has no cells.zarr.zip, which holds the cell and nucleus label images and the z_level
    # column of the table (none of these are used downstream). Read without the label images, take
    # nucleus_count from cells.parquet instead, leave z_level empty, and link the table to the cell
    # boundary polygons instead of the cell label image.
    import spatialdata_io.readers.xenium as xenium_reader

    def cells_metadata_from_parquet(path, file, specs):
        cells = pd.read_parquet(path / 'cells.parquet', columns=['cell_id', 'nucleus_count'])
        cells['cell_id'] = xenium_reader._decode_cell_id_column(cells['cell_id'])
        cells['z_level'] = np.nan
        cells.index = cells.index.astype(str)
        return cells

    xenium_reader._get_cells_metadata_table_from_zarr = cells_metadata_from_parquet
    sdata = xenium(path=raw_data_dir_sample, n_jobs=8, cells_as_circles=False, morphology_mip=False, morphology_focus=False, transcripts=True,
                   cells_labels=False, nucleus_labels=False)
    sdata.tables['table'].obs['region'] = 'cell_boundaries'
    sdata.set_table_annotates_spatialelement('table', region='cell_boundaries', region_key='region', instance_key='cell_id')


# --- Update the ids of the anndata object to match the cell_id
# anndata.obs.index comes through as a sequence of numbers rather than the cell_ids
# Copy them over from the segmentation mask
# If you don't do this you won't be able to plot anything on the segmentation mask because it won't know how to match
# data from the anndata object to the segmentation mask

sdata.tables['table'].obs.index = sdata.shapes['cell_boundaries'].index
sdata.tables['table'].obs.index.name = None


# --- Basic qc

# Percentage of control probes and control codewords
cprobes = (
    sdata.tables['table'].obs["control_probe_counts"].sum() / sdata.tables['table'].obs["total_counts"].sum() * 100
)
cwords = (
    sdata.tables['table'].obs["control_codeword_counts"].sum() / sdata.tables['table'].obs["total_counts"].sum() * 100
)
print(f"Negative DNA probe count % : {cprobes}")
print(f"Negative decoding count % : {cwords}")


# --- Mark cells as low quality (low counts) but do not filter from object

min_count=10
sdata.tables['table'].obs['low_qc'] = np.where((sdata.tables['table'].obs['total_counts'] <= min_count), 'low counts', 'pass')

n_low_qc = sum(sdata.tables['table'].obs['low_qc'] == 'low counts')

print(f'Low quality cells: {n_low_qc}/{len(sdata.tables['table'].obs_names)} ({round(n_low_qc/len(sdata.tables['table'].obs_names)*100, 2)}%)')

# Filter before log normalization
sdata.tables['table'].layers['filtered_counts'] = sdata.tables['table'].X.copy()
sc.pp.filter_cells(sdata.tables['table'], min_counts=min_count)

# Filter genes expressed in fewer than 5 cells
sc.pp.filter_genes(sdata.tables['table'], min_cells=5)


# --- Normalize
# layers['filtered_counts']: counts scaled to 10,000 per cell, used for bin-level
# sums (TLS zonation, chemokine densities)
sc.pp.normalize_total(sdata.tables['table'], layer='filtered_counts', target_sum=1e4, inplace=True)

# X: log1p of the counts.
sc.pp.log1p(sdata.tables['table'])


# --- Save the spatialdata object

# save AnnData object to *.h5ad file
outs_dir = os.path.join(project_dir, 'data', 'spatialdata_objects', 'xenium', 'without_image')
os.makedirs(outs_dir, exist_ok=True)

print(f"Saving .h5ad files to {outs_dir}")

# save processed data object
sdata.write(os.path.join(outs_dir, f"{sample_id}.zarr"))

# Zip the folder to avoid having thousands of file in the dvc object
shutil.make_archive(os.path.join(outs_dir, f"{sample_id}"), 'zip', os.path.join(outs_dir, f"{sample_id}.zarr"))

# Remove the unzipped folder to save space
shutil.rmtree(os.path.join(outs_dir, f"{sample_id}.zarr"))


print("Script successfully completed")