#!/usr/bin/env python3

# --- Aim
# Extract raw per-cell transcript counts for a list of genes (CXCL13, CCL19, CCL21, TSPAN7)
# from the Xenium SpatialData object. The AnnData `filtered_counts` layer is library-normalised
# to 10000 / cell (CP10K), NOT raw integer counts, so real integer counts come from the
# SpatialData transcripts table. Outputs one CSV per (sample, gene) with the per-cell integer
# transcript count. These CSVs feed the downstream source-attribution scripts (04-07).

# --- Define environment
# env: spatdat

# --- Load packages

import argparse
import os
import shutil
import tempfile

import spatialdata as sd


# --- Argparse arguments
# Single-sample: dispatches through the array runner
# (submit_array_jobs_nci.sh -> run_single_sample_job.sh), which passes
# --h5ad_dir --sample_id --project_dir --res_dir --sample_sheet.

parser = argparse.ArgumentParser(
    description="Extract raw per-cell transcript counts from a Xenium SpatialData object",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--sample_id", required=True, help="Sample ID (section_id)")
parser.add_argument("--project_dir", required=True)
parser.add_argument("--res_dir", required=True)
parser.add_argument("--h5ad_dir", default=None, help="(inert) accepted for the runner contract")
parser.add_argument("--sample_sheet", default=None, help="(inert) accepted for the runner contract")
parser.add_argument("--spatdata_dir", default=None,
                    help="Path to spatialdata .zarr/.zip files")
parser.add_argument("--genes", nargs="*", default=['CXCL13', 'CCL19', 'CCL21', 'TSPAN7'],
                    help="Genes to extract")

args = vars(parser.parse_args())

sample_id    = args["sample_id"]
project_dir  = args["project_dir"]
res_dir      = args["res_dir"]
spatdata_dir = args["spatdata_dir"] or os.path.join(project_dir, 'data',
                                                    'spatialdata_objects', 'xenium')
genes = args["genes"]


# --- Output dir (plots / csv_files split; per-sample intermediates -> csv_files only)

tx_out_dir = os.path.join(res_dir, '05_xenium_mechanisms', '09_chemokine_sources',
                          '02_extract_transcripts', 'csv_files')
os.makedirs(tx_out_dir, exist_ok=True)


# --- Read the SpatialData object (use the persistent .zarr if present, else unzip a .zip to tmp)

print(f'\n=== {sample_id} ===')
zip_path = os.path.join(spatdata_dir, f'{sample_id}.zip')
zarr_persist = os.path.join(spatdata_dir, f'{sample_id}.zarr')
use_persistent = (os.path.isdir(zarr_persist) and
                  sum(os.path.getsize(os.path.join(dp, f))
                      for dp, _, fns in os.walk(zarr_persist) for f in fns) > 1e9)

tmp_root = None
try:
    if use_persistent:
        zarr_path = zarr_persist
        print(f'  Reading existing zarr ({zarr_persist})')
    else:
        tmp_root = tempfile.mkdtemp(prefix=f'tx_extract_{sample_id}_')
        zarr_path = os.path.join(tmp_root, f'{sample_id}.zarr')
        print(f'  Unzipping {zip_path} -> {zarr_path}')
        shutil.unpack_archive(zip_path, zarr_path)

    sdata = sd.read_zarr(zarr_path)
    tx = sdata.points['transcripts']

    for gene in genes:
        cx = tx[tx['feature_name'] == gene]
        cx_df = cx.compute().reset_index(drop=True)
        assigned = cx_df[cx_df['cell_id'].astype(str) != 'UNASSIGNED']
        per_cell = (assigned.groupby('cell_id', observed=True).size()
                    .rename(f'n_{gene}_transcripts').reset_index())
        out_path = os.path.join(tx_out_dir, f'{sample_id}_{gene}_per_cell.csv')
        per_cell.to_csv(out_path, index=False)
        print(f'  {gene}: {len(cx_df):,} transcripts ({len(cx_df) - len(assigned):,} UNASSIGNED); '
              f'{len(per_cell):,} cells -> {out_path}')

finally:
    if tmp_root and os.path.isdir(tmp_root):
        shutil.rmtree(tmp_root, ignore_errors=True)
        print(f'  Cleaned up tmp')


print('\nScript successfully completed')
