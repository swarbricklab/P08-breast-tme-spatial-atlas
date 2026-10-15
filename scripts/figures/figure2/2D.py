#!/usr/bin/env python3

# --- Aim
# Plot Figure 2D: heatmap of top marker genes per TLS cluster

# --- Define environment
# env: spatial-pilot

# --- Load packages
import argparse
import os

import scanpy as sc
import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import seaborn as sns

mpl.rcParams['pdf.fonttype'] = 42

sc.set_figure_params(facecolor="white", figsize=(8, 8), fontsize=10)
sc.settings.verbosity = 3


# --- Argparse arguments
parser = argparse.ArgumentParser(description="Import Visium unfiltered data & plot",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-s", "--samples",
    nargs="*",  # 0 or more values expected => creates a list
    type=str,
    help="List of samples")
parser.add_argument("-r", "--res_dir", help="Path to results directory")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)

args = vars(parser.parse_args())

# Arguments
project_dir = args["project_dir"]
h5ad_dir = args["h5ad_dir"]
res_dir = args["res_dir"] or os.path.join(project_dir, 'results')  # default: <project_dir>/results
sample_list = args["samples"]

print(f'Samples in python: {sample_list}')


# --- Extract some info from the sample sheet

# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_dir = os.path.join(project_dir, 'config')
    sample_sheet_file = 'sample_sheet_visium.csv'
    sample_sheet_path = os.path.join(sample_dir, sample_sheet_file)

sample_sheet = pd.read_csv(sample_sheet_path)
sample_lut = sample_sheet.set_index('section_id')


# --- Also get some tissue metadata

resource_dir = os.path.join(project_dir, 'resources')
metadata_file = 'tissue_metadata.csv'

tissue_metadata = pd.read_csv(os.path.join(resource_dir, metadata_file))

# Merge the two DataFrames on the common column 'tissue_id'
merged_df = pd.merge(sample_sheet, tissue_metadata, on='tissue_id')

# Create a sample - subtype dictionary
subtype_dict = dict(zip(merged_df['section_id'], merged_df['clinical_subtype']))

# And make a sample patient LUT
sample_donor_lut = {}

for sample in sample_list:
    sample_donor_lut[sample] = sample_lut.loc[sample, 'donor_id']


# --- Load the clustering results

# Load harmony clustering results
harmony_dir = os.path.join(res_dir, '03_TLS_niche_discovery', '02_identify_TLS_niches', '01_cluster_TLS_high')
most_recent = harmony_dir

# Load the adata containing results
all_adata = ad.read_h5ad(os.path.join(most_recent, 'csv_files', 'All_samples_TLS_umap.h5ad'))

use_data = 'log_norm'
use_layer = 'log_norm'


# --- Also load top genes

filename = 'All_samples_leiden_clusters_DEGS.csv'
top_genes_table = pd.read_csv(os.path.join(most_recent, 'csv_files', filename))

# Cluster numbers may shift a bit but should be pretty clear from gene expression which is which
top_marker_dict = {
    'T cell zone': ['CCL19', 'LAMP3'],
    'B cell zone': ['CXCL13', 'MS4A1'],
    'HEVs': ['CCL21', 'ACKR1'],
    'Germinal center': ['RGS13', 'MEF2B'],
}

# Get all unique cluster labels as strings
all_clusters = sorted(top_genes_table['group'].astype(str).unique(), key=lambda x: int(x) if x.isdigit() else x)

clusters_of_interest = {}

for region, markers in top_marker_dict.items():
    for cluster in all_clusters:
        # Get the top 20 genes for this cluster
        cluster_genes = top_genes_table[top_genes_table['group'].astype(str) == cluster]['genes'].values[:20]
        if all(marker in cluster_genes for marker in markers):
            clusters_of_interest[cluster] = region
            break

print("Clusters of interest based on marker gene expression:")
for cluster, region in clusters_of_interest.items():
    print(f"{region}: Cluster {cluster}")


# --- Plot 1: Dot plots for top 20 genes per cluster - separated by cluster

leiden_clusters = all_adata.obs['leiden'].cat.categories

# Prepare a dict with top 20 genes per cluster
top = 20

top_genes_dict = {}

for cluster in sorted(clusters_of_interest.keys(), key=lambda x: int(x) if x.isdigit() else x):
    top_genes = top_genes_table[top_genes_table['group'].astype(str) == cluster].head(top)
    top_genes_dict[clusters_of_interest[cluster]] = list(top_genes['genes'])


# --- Add new column to adata with leiden clusters of interest

all_adata.obs['leiden_adj'] = all_adata.obs['leiden'].astype(str)
# Replace other clusters with NaN
all_adata.obs['leiden_adj'] = all_adata.obs['leiden_adj'].where(all_adata.obs['leiden_adj'].isin(clusters_of_interest.keys()), other=np.nan)

# Rename the clusters according to the dict
for cluster, name in clusters_of_interest.items():
    all_adata.obs['leiden_adj'] = all_adata.obs['leiden_adj'].replace(cluster, f'c{cluster}: {name}')


# --- Plot 1: Heatmap genes of interest

# Disable the grid
plt.rcParams['axes.grid'] = False

dp = sc.pl.matrixplot(all_adata, var_names=top_genes_dict, groupby=['leiden_adj'], layer=use_data, return_fig=True,
                   standard_scale='var', cmap='viridis')
dp.add_totals().show()

# Adjust the size of the axes
axs = plt.gcf().get_axes()

heatmap_ax = axs[0]

# Set the desired ratio
heatmap_ax.set_aspect(1.75)  # Adjust the aspect ratio as needed

# Adjust the other axs
new_heatmap_pos = heatmap_ax.get_position()

bar_ax = axs[1]
bar_ax_pos = bar_ax.get_position()
bar_ax.set_position([new_heatmap_pos.x1, new_heatmap_pos.y0, bar_ax_pos.width, new_heatmap_pos.height])

# Move the cbar
cbar_ax = axs[-1]
bar_ax_pos = bar_ax.get_position()
cbar_ax_pos = cbar_ax.get_position()
cbar_ax.set_position([bar_ax_pos.x1, bar_ax_pos.y0, cbar_ax_pos.width*0.75, cbar_ax_pos.height])

# Make the top ax invisible
top_ax = axs[2]
top_ax.set_visible(False)

# Add cluster labels above the heatmap
x_offset = 0  # adjust if your heatmap has an offset
y = 1.05  # relative to axis, adjust as needed

current_x = 0
for cluster_name, gene_list in top_genes_dict.items():
    n_genes = len(gene_list)
    x_center = current_x + n_genes / 2 - 0.5  # center of the gene block
    heatmap_ax.text(
        x_center, y, cluster_name, 
        ha='center', va='bottom', 
        fontsize=12, fontweight='bold', 
        transform=heatmap_ax.get_xaxis_transform()
    )
    current_x += n_genes

# Calculate position for plot title
title_x = (new_heatmap_pos.x0 + new_heatmap_pos.x1) / 2
title_y = new_heatmap_pos.y1 + 0.2

plt.suptitle(f'All samples - Top {top} (use_data={use_data})', x=title_x, y=title_y, fontweight='bold')


# --- Save the output

# Prepare for saving
out_dir_plots = os.path.join(res_dir, 'figures')
os.makedirs(out_dir_plots, exist_ok=True)

#  Save plots
filename = 'Fig_2D.pdf'
plt.savefig(os.path.join(out_dir_plots, filename), bbox_inches='tight')
plt.close('all')


# --- Finish

print("Script successfully completed")