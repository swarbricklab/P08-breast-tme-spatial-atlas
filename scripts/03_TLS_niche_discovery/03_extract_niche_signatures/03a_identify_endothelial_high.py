#!/usr/bin/env python3

# --- Aim
# Refine HEV (high endothelial venule) gene signature by filtering out genes broadly expressed in endothelial-dominant spots.
# Integrate clustering and deconvolution results, identify endothelial-high regions, and visualize spatial distributions.
# Support downstream HEV-specific marker discovery.

# --- Define environment
# env: spatial-pilot

# --- Load packages
import argparse
import scanpy as sc
import anndata as ad
import os
import sys
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from datetime import date
from datetime import datetime

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium data & do basic QC",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required

parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)

parser.add_argument("--min_count", help="Min count for QC", type=int)
parser.add_argument("--min_gene", help="Min gene for QC", type=int)
parser.add_argument("--min_spots", help="Minimal spots per gene", type=int)

args = vars(parser.parse_args())


# --- Arguments

res_dir = args["res_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
project_dir=args["project_dir"]

sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import auto_crop_scanpy, save_multi_image


# --- Get some info from the sample sheet
# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet.csv')

sample_sheet = pd.read_csv(sample_sheet_path)
tissue_id = sample_sheet[sample_sheet['section_id'] == sample_id]['tissue_id'].values[0]
slide_size = sample_sheet[sample_sheet['section_id'] == sample_id]['slide_size'].values[0]


# --- Also get some tissue metadata

resource_dir = os.path.join(project_dir, 'resources')
metadata_file = 'tissue_metadata.csv'

tissue_metadata = pd.read_csv(os.path.join(resource_dir, metadata_file))

subtype = tissue_metadata[tissue_metadata['tissue_id'] == str(tissue_id)]['clinical_subtype'].values[0]


# --- Load the adata object

adata = ad.read_h5ad(os.path.join(h5ad_dir, 'visium', f'{sample_id}.h5ad'))
use_data = 'log_norm'

# --- Calculate crop coords for optimal presentation

coords = auto_crop_scanpy(adata, border_size=0.1)


# --- Set spot size for large & small capture area
# Plotting the spots larger to make it easier to see

if slide_size == '11 mm':
    spot_size = 1.8
else:
    spot_size = 1.3


# --- Load clustering results

# Function to find the most recent results
def find_most_recent(dir):
    # Find most recent results
    subdirs = [d for d in os.listdir(dir) if os.path.isdir(os.path.join(dir, d))]
    # Filter subdirectories to keep only those with date format YYYYmmdd
    date_format = '%Y%m%d'
    valid_subdirs = []
    for subdir in subdirs:
        try:
            datetime.strptime(subdir, date_format)
            valid_subdirs.append(subdir)
        except ValueError:
            pass
    # Sort the valid subdirectories by date and get the most recent one
    if valid_subdirs:
        most_recent_dir = max(valid_subdirs, key=lambda x: datetime.strptime(x, date_format))
        most_recent_dir_path = os.path.join(dir, most_recent_dir)
        print("Most recent results:", most_recent_dir_path)
    else:
        print("No valid directories found.")
    return most_recent_dir_path


# Load harmony clustering results
harmony_dir = os.path.join(project_dir, 'results', '03_TLS_niche_discovery', '02_identify_TLS_niches', '01_cluster_TLS_high')
most_recent = find_most_recent(harmony_dir)

filename = 'All_samples_leiden_clusters.csv'
cluster_res = pd.read_csv(os.path.join(most_recent, 'csv_files', filename), index_col=0)

# Add cluster results to anndata 
results_sub = cluster_res[cluster_res['sample_id'] == sample_id]
adata.obs['leiden'] = results_sub['leiden'].astype('category')

# Also load the extracted signature list for HEV
signature_dir = os.path.join(project_dir, 'results', '03_TLS_niche_discovery', '03_extract_niche_signatures', '01_extract_niche_signatures')
filename = 'hevs_top_genes_filtered.csv'

hev_signature = pd.read_csv(os.path.join(signature_dir, filename))
hev_cluster = hev_signature['group'].values[0]


# --- Load the decon results

# Decide which level to plot
level = 'major'

decon_dir = os.path.join(project_dir, 'results', '02_cell_type_mapping', 'deconvolution', 'visium', f'{level}_level')

# Read in the results & add to anndata
decon_filename = f'{sample_id}_{level}_q05_cell_abundance_w_sf.csv'

# Load decon results
decon_results = pd.read_csv(os.path.join(decon_dir, decon_filename), index_col=0)

cell_types = list(decon_results.columns)

adata.obs[cell_types] = decon_results


# --- Prepare for saving

qc_dir = os.path.join(res_dir, '03_TLS_niche_discovery', '03_extract_niche_signatures', '03_filter_HEV_signature', 'identify_endothelial_high')

os.makedirs(qc_dir, exist_ok=True)

today = date.today()
today = today.strftime("%Y%m%d")


out_dir = os.path.join(qc_dir, today)
os.makedirs(out_dir, exist_ok=True)


# --- Define some colors

# Define some colors per score
color_dict = {
    'germinal_center': mpl.colors.cnames["cyan"],
    'B_cell_zone': mpl.colors.cnames["cornflowerblue"],
    'T_cell_zone': mpl.colors.cnames["crimson"],
    'mixed': mpl.colors.cnames["indigo"]
}


# --- Plot 1: Plot decon results & clustering results

endothelial_key = 'Endothelial'

sc.set_figure_params(scanpy=True, fontsize=18)

nrows = 1
ncols = 3

# Make figure
fig, axs = plt.subplots(nrows, ncols, figsize=(7.5*ncols, 6*nrows))

# First plot H&E
n = 0
ax = axs[n]

sc.pl.spatial(adata, color=None, size=spot_size, ax=ax, bw=False, alpha_img=1,
    crop_coord=coords, title=f'{sample_id} ({subtype})\nH&E',)

ax.set_title(ax.get_title(), fontweight='bold')


# Then plot clustering results
n += 1
ax = axs[n]

sc.pl.spatial(adata, color='leiden', groups=[hev_cluster], size=spot_size, ax=ax, bw=False, alpha_img=0.7,
    crop_coord=coords, title=f'{sample_id} ({subtype})\nHEV cluster', palette=['#1f77b4'])

ax.set_title(ax.get_title(), fontweight='bold')

# Then plot score classifications
n += 1
ax = axs[n]

sc.pl.spatial(adata, color=endothelial_key, size=spot_size, ax=ax, bw=False, alpha_img=0.7, cmap='magma_r',
    crop_coord=coords, title=f'{sample_id} ({subtype})\nEndothelial cells', vmax='p99.2')

ax.set_title(ax.get_title(), fontweight='bold')

plt.suptitle(f'{sample_id} ({subtype}) - clustering results & decon results', y=1.1, fontweight='bold', x=0.5)


# --- Find endothelial-high spots

total_per_spot = decon_results.sum(axis=1)
decon_prop = decon_results.div(total_per_spot, axis=0)

# Make a new df with just the endothelial cell info
endo_df = pd.DataFrame()
endo_df['endo_cells'] = decon_results[endothelial_key]
endo_df['endo_prop'] = decon_prop[endothelial_key]

# Set some thresholds
quantile = 0.8
cells_thresh = decon_results[endothelial_key].quantile(quantile)
prop_thresh = 0.1

endo_df['pass'] = np.where((endo_df['endo_cells'] > cells_thresh) & (endo_df['endo_prop'] > prop_thresh),
                           'pass', None)

adata.obs['endo_high'] = endo_df['pass']

ncols = 5

# Make figure
fig, axs = plt.subplots(nrows, ncols, figsize=(7.5*ncols, 6*nrows))

# First plot H&E
n = 0
ax = axs[n]

sc.pl.spatial(adata, color=None, size=spot_size, ax=ax, bw=False, alpha_img=1,
    crop_coord=coords, title=f'{sample_id} ({subtype})\nH&E',)

ax.set_title(ax.get_title(), fontweight='bold')

# Then plot clustering results
n += 1
ax = axs[n]

sc.pl.spatial(adata, color='leiden', groups=[hev_cluster], size=spot_size, ax=ax, bw=False, alpha_img=0.7,
    crop_coord=coords, title=f'{sample_id} ({subtype})\nHEV cluster', palette=['#1f77b4'])

ax.set_title(ax.get_title(), fontweight='bold')

# Then plot endothelial high
n += 1
ax = axs[n]

if len(adata.obs['endo_high'].dropna()) > 0:
    plot_key = 'endo_high'
else:
    plot_key = None

sc.pl.spatial(adata, color=plot_key, size=spot_size, ax=ax, bw=False, alpha_img=0.7,
    crop_coord=coords, title=f'{sample_id} ({subtype})\nEndothelial-high', palette=['#1f77b4'])

ax.set_title(ax.get_title(), fontweight='bold')

# Then plot decon results
n += 1
ax = axs[n]

sc.pl.spatial(adata, color=endothelial_key, size=spot_size, ax=ax, bw=False, alpha_img=0.7, cmap='magma_r',
    crop_coord=coords, title=f'{sample_id} ({subtype})\nEndothelial cells', vmax='p99.2')

ax.set_title(ax.get_title(), fontweight='bold')


# Then plot decon proportion
adata.obs['endothelial_prop'] = endo_df['endo_prop']

n += 1
ax = axs[n]

sc.pl.spatial(adata, color='endothelial_prop', size=spot_size, ax=ax, bw=False, alpha_img=0.7, cmap='magma_r',
    crop_coord=coords, title=f'{sample_id} ({subtype})\nEndothelial cells (proportion)', vmax=None)

ax.set_title(ax.get_title(), fontweight='bold')

plt.suptitle(f'{sample_id} ({subtype}) - clustering results & decon results (endothelial-high: min_prop={prop_thresh}, min_cells={cells_thresh.round(2)} (quantile {quantile}))', y=1.1, fontweight='bold', x=0.5)


# --- Save results

# Save plots
out_dir_plots = os.path.join(out_dir, 'plots')
os.makedirs(out_dir_plots, exist_ok=True)

filename = f'{sample_id}_identify_endothelial_high.pdf'
save_multi_image(os.path.join(out_dir_plots, filename))

# Save csv files
out_dir_csv = os.path.join(out_dir, 'csv_files')
os.makedirs(out_dir_csv, exist_ok=True)

# Also save endothelial spots
adata.obs['endo_sub'] = np.where(adata.obs['endo_high'].notna(), 'Endothelial-high', None)
adata.obs['endo_sub'] = np.where(adata.obs['leiden'] == hev_cluster, 'HEV', adata.obs['endo_sub'])

filename = f'{sample_id}_endothelial_high_spots.csv'
to_save = adata.obs[['leiden', 'endo_high', 'endo_sub']]

to_save.to_csv(os.path.join(out_dir_csv, filename))


# --- Finish

print("Script successfully completed")
