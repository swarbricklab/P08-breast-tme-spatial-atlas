#!/usr/bin/env python3

# --- Aim
# Plot Figure 4A: gene expression comparison between ACKR1+ (venous) ECs in HEV vs non-HEV niches

# --- Define environment
# env: xen-hev

# --- Load packages
import argparse
import os

import anndata as ad

import matplotlib as mpl
import matplotlib.pyplot as plt

import scanpy as sc

import pandas as pd
import numpy as np

mpl.rcParams['pdf.fonttype'] = 42


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


# --- Extract some info from the sample sheet

# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_dir = os.path.join(project_dir, 'config')
    sample_sheet_file = 'sample_sheet_xenium.csv'
    sample_sheet_path = os.path.join(sample_dir, sample_sheet_file)

sample_sheet = pd.read_csv(sample_sheet_path)

# Also load tissue metadata
metadata_dir = os.path.join(project_dir, 'resources',)
tissue_metadata_file = 'tissue_metadata.csv'
tissue_metadata = pd.read_csv(os.path.join(metadata_dir, tissue_metadata_file))

# Make a section_id to subtype mapping
section_to_subtype = tissue_metadata.set_index('tissue_id')['clinical_subtype'].to_dict()
# Map section_id to subtype
sample_sheet['subtype'] = sample_sheet['tissue_id'].map(section_to_subtype)

subtype_dict = sample_sheet.set_index('section_id')['subtype'].to_dict()


# --- Load the anndata objects

adata_dict = {}

for sample in sample_list:
    adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample}.h5ad'))
    adata_dict[sample] = adata


# --- Load the HEV annotations

distance_dir = os.path.join(res_dir, '05_xenium_mechanisms', '02_identify_HEV_regions')
most_recent_distance = distance_dir

for sample in sample_list:
    # Also load HEV niche barcode and add to adata
    filename = f'{sample}_HEV_niches_cell_to_niche.csv'
    niche_barcode_df = pd.read_csv(os.path.join(most_recent_distance, 'csv_files', filename), index_col=0)
    niche_barcode_df.dropna(subset=['niche_id'], inplace=True)
    # Prepend sample name to niche_id
    niche_barcode_df['niche_id'] = sample + '_' + niche_barcode_df['niche_id'].astype(str)
    adata_dict[sample].obs['niche_id'] = niche_barcode_df['niche_id']
    # Also load the distance calculations
    distance_file = f'{sample}_HEV_niches_info.csv'
    distance_df = pd.read_csv(os.path.join(most_recent_distance, 'csv_files', distance_file))
    # Prepend sample name to niche_id
    distance_df['niche_id'] = sample + '_' + distance_df['niche_id'].astype(str)
    # Make a mapping dict of distance to TLS
    distance_dict = distance_df.set_index('niche_id')['close_to_TLS'].to_dict()
    # Map to adata
    adata_dict[sample].obs['close_to_TLS'] = adata_dict[sample].obs['niche_id'].map(distance_dict)


# --- Also load cell annotations

# Decide which annotation levels to use
level = 'subset'

for sample in sample_list:
    annotation_file = f'{sample}_{level}_pred_cell_types.csv'
    annotation_dir = os.path.join(res_dir, '02_cell_type_mapping', 'label_transfer', 'xenium')
    annotations = pd.read_csv(os.path.join(annotation_dir, f'{level}_level', annotation_file), index_col=0)
    # Add to adata
    adata_dict[sample].obs[level] = annotations['pred_cell_type']


# --- Look at gene expression in HEV niches vs non-HEV niches

# Subset adatas to endothelial ACKR1+ cells only
for sample in sample_list:
    adata = adata_dict[sample]
    adata_ACKR1 = adata[adata.obs[level] == 'Endothelial ACKR1'].copy()
    adata_dict[sample] = adata_ACKR1

# Extract a single adata object with all samples
adata_list = list(adata_dict.values())
all_adata = ad.concat(adata_list, merge='same')


# --- Load the endothelial colors

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
filename = 'endothelial_colors.csv'

endothelial_colors_df = pd.read_csv(os.path.join(color_dir, filename), index_col=0)
color_dict = endothelial_colors_df['color'].to_dict()


# --- Also add the HEV colors

filename = 'TLS_HEV_colors.csv'

HEV_colors_df = pd.read_csv(os.path.join(color_dir, filename), index_col=0)

# Add to color_dict
for cell_type, row in HEV_colors_df.iterrows():
    color_dict[cell_type] = row['color']


# --- Define a dict to rename the labels for consistency

# Rename for consistency
rename_dict = {
    'HEV far from TLS': 'HEV-like (stroma)',
    'HEV close to TLS': 'HEV (TLS)',
}


# --- Define gene groups for dotplot

genes_dict = {
    'HEV markers': ['TSPAN7', 'RAMP3', 'MEOX2', 'SPNS2', 'IL33'],
    'LRI receptors': ['ACKR1', 'SELE', 'SELP',],
    'Other': ['TMEM173'],
}

# Filter any genes not in the dataset
for gene_group in genes_dict.keys():
    genes_dict[gene_group] = [gene for gene in genes_dict[gene_group] if gene in all_adata.var_names]

genes_list = []
for gene_group in genes_dict.values():
    genes_list.extend(gene_group)


# --- Plot 1: Make a dotplot of the HEV markers and LRI receptors in HEV vs non-HEV niches

# Add niche type back to all_adata
all_adata.obs['niche_type'] = np.where(all_adata.obs['niche_id'].notna(), 'HEV', 'Venous')

# Add niche type tls proximity back to all_adata
all_adata.obs['niche_TLS_proximity'] = 'Venous'
hev_close_mask = (all_adata.obs['niche_id'].notna()) & (all_adata.obs['close_to_TLS'] == True)
hev_far_mask = (all_adata.obs['niche_id'].notna()) & (all_adata.obs['close_to_TLS'] == False)
all_adata.obs.loc[hev_close_mask, 'niche_TLS_proximity'] = 'HEV close to TLS'
all_adata.obs.loc[hev_far_mask, 'niche_TLS_proximity'] = 'HEV far from TLS' 

plot_order_TLS = ['Venous', 'HEV far from TLS', 'HEV close to TLS']

# Set categorical
all_adata.obs['niche_TLS_proximity'] = pd.Categorical(all_adata.obs['niche_TLS_proximity'], categories=plot_order_TLS, ordered=True)

dp = sc.pl.dotplot(
    all_adata,
    var_names=genes_dict,
    groupby='niche_TLS_proximity',
    show=False,
    swap_axes=True,
    return_fig=True,
    standard_scale='var',
    figsize=(3, 3)
)

dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5, cmap='Reds').show()

# Extract the dotplot ax
dp_ax = dp.get_axes()['mainplot_ax']

# Use rename dict to rename the xticklabels (swapped axes)
xticklabels = [tick.get_text() for tick in dp_ax.get_xticklabels()]
xticklabels = [rename_dict.get(label, label) for label in xticklabels]
dp_ax.set_xticklabels(xticklabels, rotation=90)

# Adjust the plot dimensions
n_xticks = len(dp_ax.get_xticks())
n_yticks = len(dp_ax.get_yticks())

ax_pos = dp_ax.get_position()
dp_ax.set_position([ax_pos.x0, ax_pos.y0, 0.08 * n_xticks, 0.06 * n_yticks])

# Adjust xlim
dp_ax.set_xlim(-0.1, n_xticks + 0.1)
dp_ax.set_ylim(n_yticks + 0.3, -0.3)

# Hide the gene annotation ax
gene_ax = dp.get_axes()['gene_group_ax']
gene_ax.set_visible(False)

# Move the count bar ax
count_ax = dp.get_axes()['group_extra_ax']
ax_pos = dp_ax.get_position()
count_ax.set_position([ax_pos.x0, ax_pos.y1, ax_pos.width, ax_pos.height / n_yticks * 2])

# Move the legends closer
cbar_ax = dp.get_axes()['color_legend_ax']
cbar_pos = cbar_ax.get_position()
cbar_ax.set_position([ax_pos.x1 + 0.05, cbar_pos.y0, cbar_pos.width*0.75, cbar_pos.height])

size_ax = dp.get_axes()['size_legend_ax']
size_pos = size_ax.get_position()
size_ax.set_position([ax_pos.x1 + 0.05, cbar_ax.get_position().y1 + 0.2, size_pos.width*0.75, size_pos.height])

# Position the title
title_x = (dp_ax.get_position().x0 + dp_ax.get_position().x1) / 2
title_y = count_ax.get_position().y1 + 0.1
plt.suptitle(f'All samples - Endothelial ACKR1+ expression split by niche (HEV or non-HEV, close vs far from TLS)', fontweight='bold', y=title_y, x=title_x)


# --- Save the output

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_4A.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')


# --- Finish

print("Script successfully completed")