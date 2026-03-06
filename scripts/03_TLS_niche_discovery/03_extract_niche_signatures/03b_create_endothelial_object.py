#!/usr/bin/env python3

# --- Aim
# Aggregate and integrate endothelial-high spots across all samples into a unified object.

# --- Define environment
# env: spatial-pilot

# --- Load packages
import argparse
import scanpy as sc
import scanpy.external as sce
import anndata as ad
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
from datetime import date
from datetime import datetime
import numpy as np
import matplotlib.lines as mlines

import math
import scanpy.plotting.palettes as sc_palettes

sc.set_figure_params(facecolor="white", figsize=(8, 8), fontsize=20)
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
res_dir = args["res_dir"]
sample_list = args["samples"]

print(f'Samples in python: {sample_list}')

sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import auto_crop_scanpy, save_multi_image


# --- Extract some info from the sample sheet
# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet.csv')

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


# -- Load data
# Multiple samples as a dictionary seems like an okay idea
adata_dict = {}

for sample in sample_list:
    adata_dict[sample] = ad.read_h5ad(os.path.join(h5ad_dir, 'visium', f'{sample}.h5ad'))

use_layer = 'log_norm'

# Add some info for plotting
for sample in sample_list:
    adata_dict[sample].obs['sample_id'] = sample
    adata_dict[sample].obs['subtype'] = subtype_dict[sample]
    adata_dict[sample].obs['donor_id'] = sample_lut.loc[sample, 'donor_id']


# --- Calculate crop coords for optimal presentation

coords = {}
for sample in adata_dict.keys():
    coords[sample] = auto_crop_scanpy(adata_dict[sample], border_size=0.1)


# --- Order based on sample ID

ordered_keys = sorted(adata_dict.keys())


# --- Set spot size for large & small capture area
# Plotting the spots larger to make it easier to see

spot_size={}

for sample in adata_dict.keys():
    if sample_lut.loc[sample, 'slide_size'] == '11 mm':
        spot_size[sample]=1.8
    else:
        spot_size[sample]=1.6


# --- Load the endothelial high spots

class_res_dir = os.path.join(res_dir, '03_TLS_niche_discovery', '03_extract_niche_signatures', '03_filter_HEV_signature', 'identify_endothelial_high',)

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

most_recent_dir_path = find_most_recent(class_res_dir)

# Read in the results & add to anndata
for sample in adata_dict.keys():
    TLS_file_name = f'{sample}_endothelial_high_spots.csv'
    TLS_results = pd.read_csv(os.path.join(most_recent_dir_path, 'csv_files', TLS_file_name), index_col=0)
    adata_dict[sample].obs[TLS_results.columns] = TLS_results
    columns = TLS_results.columns

classification_key = 'endo_sub'


# --- Define a colordict

# Define some colors per score
color_dict = {
    'Endothelial-high': '#1f77b4',
    'HEV': '#ff7f0e',
}


# --- Prepare for saving

# Put in a timestamped folder to avoid overwriting older plots
qc_dir = os.path.join(res_dir, '03_TLS_niche_discovery', '03_extract_niche_signatures', '03_filter_HEV_signature', 'create_endothelial_object',)
os.makedirs(qc_dir, exist_ok=True)

today = date.today()
today = today.strftime("%Y%m%d")

out_dir = os.path.join(qc_dir, today)
os.makedirs(out_dir, exist_ok=True)

# Save plots
out_dir_plots = os.path.join(out_dir, 'plots')
os.makedirs(out_dir_plots, exist_ok=True)


# --- Prepare for plotting

# Prepare for plotting the result per subtype: define a plotting order based on subtype

# Extract all subtypes
all_subtypes = set([subtype_dict[sample] for sample in sample_list if sample in subtype_dict])

# Sort according to preferred plotting order
subtype_order = ['Luminal A', 'Luminal B', 'Luminal B (HER2+)', 'HER2+', 'TNBC']
sorted_subtypes = sorted(all_subtypes, key=lambda x: subtype_order.index(x) if x in subtype_order else float('inf'))

# Sort samples first by type (using sorted_types order) and then by sample ID
list_all_samples = sorted(sample_list, key=lambda s: (sorted_subtypes.index(subtype_dict[s]), s))


# --- Plot 1: Spatial overview of endothelial high

if len(adata_dict.keys()) < 5:
        ncols = len(adata_dict.keys())
else:
    ncols = 5

nrows = math.ceil((len(adata_dict.keys()))/ncols)

plt.figure(figsize=(8*ncols, 7*nrows))

for n, sample in enumerate(list_all_samples):
    # add a new subplot iteratively
    ax = plt.subplot(nrows, ncols, n + 1)
    # filter df and plot ticker on the new subplot axis
    subtype = subtype_dict[sample]
    if len(adata_dict[sample].obs[classification_key].dropna()) > 0:
        plot_key = classification_key
    else:
        plot_key = None
    # Add legend to last plot on first row
    sc.pl.spatial(adata_dict[sample], color=plot_key, size=spot_size[sample], bw=False, alpha_img=0.5, palette=color_dict,
        crop_coord=coords[sample], title=f'{sample} ({subtype})\nEndothelial classification', ax=ax, legend_loc=None)
    ax.set_title(ax.get_title(), fontweight='bold')

# Add a legend
values = ['Endothelial-high', 'HEV']
legend_handles = [
    mlines.Line2D([], [], color=color_dict[val], marker='o', linestyle='None', markersize=8, label=val)
    for val in values
]

plt.legend(handles=legend_handles, title='Endothelial classification', title_fontproperties={'weight':'bold'}, bbox_to_anchor=(1, 0.5), loc='center left',
           frameon=False, alignment='left')

plt.suptitle(f'All samples - Endothelial classification', y=1, fontweight='bold')
plt.subplots_adjust(hspace=0.3)


# --- Save the spatial plots

filename = 'All_samples_endothelial_classification_spatial.pdf'
save_multi_image(os.path.join(out_dir_plots, filename))

plt.close('all')

# Set fontsize
sc.set_figure_params(fontsize=10)

# --- Extract a single adata and subset to just the endothelial-high spots

# Convert to adata list & concatenate
adata_list = list(adata_dict.values())
all_adata = ad.concat(adata_list, merge='same')

all_adata = all_adata[all_adata.obs[classification_key].notna()]


# --- Get colors per subtype

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
color_file = 'SCSubtype_colors.csv'

color_lut = pd.read_csv(os.path.join(color_dir, color_file))
color_lut = color_lut.set_index('subtype')

subtype_colors = list(color_lut.reindex(sorted_subtypes)['color'])


# --- Now calculate LFC for genes in HEV score

# Go back to the unfiltered gene list
# Bc we want to filter against vascular background rather than TLS background
harmony_dir = os.path.join(project_dir, 'results', '03_TLS_niche_discovery', '02_identify_TLS_niches', '01_cluster_TLS_high')

gene_dir = find_most_recent(harmony_dir)

filename = 'csv_files/All_samples_leiden_clusters_DEGS.csv'

gene_table = pd.read_csv(os.path.join(gene_dir, filename))

# Filter for just the HEV genes
# Find the HEV cluster based on markers
hev_markers = ['CCL21', 'ACKR1']
# Find the cluster with these markers in top 20 genes
hev_cluster = None
for cluster in gene_table['group'].unique():
    cluster_genes = gene_table[gene_table['group'] == cluster]['genes'].tolist()[:20]
    if all(marker in cluster_genes for marker in hev_markers):
        hev_cluster = cluster
        break

print(f'HEV cluster: {hev_cluster}')
print(f'HEV cluster genes: {gene_table[gene_table["group"] == hev_cluster]["genes"].tolist()[:20]}')

gene_table = gene_table[gene_table['group'] == hev_cluster]

# Filter for p value & LFC
gene_table = gene_table[gene_table['pvals_adj'] < 0.05]
gene_table = gene_table.sort_values(by='logfoldchanges', ascending=False)
gene_list = gene_table['genes'].to_list()

# Remove any potential nan values
gene_list = [x for x in gene_list if isinstance(x, str) or not np.isnan(x)]

# Select only the genes that are present in the anndata objects
to_plot = np.intersect1d(gene_list , all_adata.var.index).tolist()

print(f'Total genes in HEV cluster after filtering: {len(gene_list)}')

missing = list(set(gene_list) - set(to_plot))

print(f'Genes not found: {missing}') 

# Cut off at 100 if gene list exceeds 100 so the plot doesn't become too big
cut = 100
if cut > len(gene_list):
    cut = len(gene_list)

dp = sc.pl.dotplot(all_adata, var_names=gene_list[:cut], groupby=[classification_key], layer=use_layer, swap_axes=False, 
                return_fig=True, )
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5, cmap='Reds',).show()

plt.suptitle(f'All samples - marker genes endothelial-high & HEV cluster spots (top {cut})', y=1.1, fontweight='bold', x=0.5)


# --- Calculate LFCs to filter

# Subset the data for each group
adata_A = all_adata[all_adata.obs['endo_sub'] == 'Endothelial-high']
adata_B = all_adata[all_adata.obs['endo_sub'] == 'HEV']

# Compute mean expression per group
mean_A = np.mean(adata_A[:, gene_list].to_df(layer=use_layer), axis=0)
mean_B = np.mean(adata_B[:, gene_list].to_df(layer=use_layer), axis=0)

# Compute log2 fold change
lfc = np.log2((mean_B) / (mean_A))

pos_lfc_genes = list(lfc[lfc > 0].index)
neg_lfc_genes = list(lfc[lfc < 0].index)

gene_dict = {'Positive': pos_lfc_genes,
            'Negative': neg_lfc_genes}

# Limit the number of genes per category
cut = 50
for key, gene_list in gene_dict.items():
    if cut > len(gene_list):
        cut = len(gene_list)
    gene_dict[key] = gene_list[:cut]

dp = sc.pl.dotplot(all_adata, var_names=gene_dict, groupby=['endo_sub'], layer=use_layer, swap_axes=False, 
                return_fig=True, )
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5, cmap='Reds',).show()

plt.suptitle(f'All samples - marker genes endothelial-high & HEV cluster spots (top {cut} per category)', y=1.1, fontweight='bold', x=0.5)


# --- Save output

filename = f'All_samples_endothelial_vs_HEV_genes.pdf'
save_multi_image(os.path.join(out_dir_plots, filename))

out_dir_csv = os.path.join(out_dir, 'csv_files')
os.makedirs(out_dir_csv, exist_ok=True)

# Save genes & whether they had pos or neg LFC
lfc = pd.DataFrame(lfc)
lfc.columns = ['LFC_HEV_vs_endothelial']
lfc['direction'] = np.where(lfc['LFC_HEV_vs_endothelial'] > 0, 'positive', 'negative')

# Add the original score info
merged_df = gene_table.merge(lfc, right_index=True, left_on='genes')

filename = f'All_samples_HEV_vs_endothelial_genes.csv'
merged_df.to_csv(os.path.join(out_dir_csv, filename), index=False)

# Also save a df that has just the genes with positive LFC
merged_df = merged_df[merged_df['direction'] == 'positive']

filename = f'All_samples_HEV_vs_endothelial_positive_LFC.csv'
merged_df.to_csv(os.path.join(out_dir_csv, filename), index=False)

# Save the whole object as well to have all the umap stuff
ad.AnnData.write(all_adata, filename=os.path.join(out_dir_csv, f"All_samples_endothelial_object.h5ad"))


# --- Finish

print("Script successfully completed")