#!/usr/bin/env python3

# --- Aim
# After classifying vessels based on decon results: investigate which vessels are likely to be arterial or venous
# Using marker genes identified in Kumar 2023 normal breast study
# This script: all samples together

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import os

from datetime import datetime
import anndata as ad
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


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

print(f'Samples in python: {sample_list}')


# --- Extract some info from the sample sheet

# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_dir = os.path.join(project_dir, 'config')
    sample_sheet_file = 'sample_sheet_visium_hd.csv'
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


# --- Load miniatlas classifications

decon_level = 'subset'

# To make it easier to visualize trends by lineage
class_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
class_file = 'miniatlas_subset_classifications.csv'

class_df = pd.read_csv(os.path.join(class_dir, class_file))

# Get endothelial classes
endothelial = list(class_df[class_df['miniatlas_major'] == 'Endothelial'][f'miniatlas_{decon_level}'].unique())
endothelial = sorted(endothelial)


# --- Load the adata objects

resolution = 'square_016um'

# Multiple samples as a dictionary seems like an okay idea
adata_dict = {}

use_data = 'log_norm'

for sample in sample_list:
    adata_dict[sample] = ad.read_h5ad(os.path.join(h5ad_dir, 'visium_hd', resolution, f'{sample}.h5ad'))

# Add the sample_id to the obs
for sample in adata_dict.keys():
    adata_dict[sample].obs['sample_id'] = sample


# --- Load the vessel classifications

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

# Specify path
TLS_group_dir = os.path.join(project_dir, 'results', '04_vascular_microenvironment', '02_identify_all_vasculature', '01_classify_vessels', 'csv_files')
most_recent_dir_path = find_most_recent(TLS_group_dir)

# Load per sample
for sample in sample_list:
    filename = f'{sample}_vessel_classification.csv'
    vessel_classes = pd.read_csv(os.path.join(most_recent_dir_path, filename), index_col=0)
    adata_dict[sample].obs[vessel_classes.columns] = vessel_classes
    # Collapse to single vessel class per spot
    collapsed_col = vessel_classes[endothelial].apply(lambda row: row.first_valid_index() if row.first_valid_index() is not None else None, axis=1)
    adata_dict[sample].obs['vessel_class'] = collapsed_col


# --- Load endothelial marker genes

gene_dir = os.path.join(project_dir, 'resources', 'gene_lists', 'marker_genes')
gene_file = 'kumar_2023_endothelial_genes.csv'

gene_table = pd.read_csv(os.path.join(gene_dir, gene_file))

# Create a dict of genes per class
gene_dict = gene_table.groupby('class')['gene'].apply(list).to_dict()

# Remove any genes not in the adata object
filtered_gene_dict = {}
missing_genes = {}

# Check which genes are missing
sample_id = sample_list[0]  # Use the first sample as reference
adata = adata_dict[sample_id]

for cls, genes in gene_dict.items():
    valid_genes = [g for g in genes if g in adata.var_names]
    missing = [g for g in genes if g not in adata.var_names]
    filtered_gene_dict[cls] = valid_genes
    if missing:
        missing_genes[cls] = missing
        print(f"Class '{cls}' is missing {len(missing)} genes: {missing}")

# Optional: see overall missing stats
total_missing = sum(len(v) for v in missing_genes.values())
print(f"Total missing genes: {total_missing}")


# --- Extract a single adata object for plotting

# Subset  to only vessels with a classification
adata_dict = {sample: adata_dict[sample][adata_dict[sample].obs['vessel_class'].notna()] for sample in adata_dict.keys()}

adata_list = list(adata_dict.values())
all_adata = ad.concat(adata_list, merge='same')

# Make obs unique
all_adata.obs_names_make_unique()

sc.set_figure_params(fontsize=12)


# --- Look at gene expression per vessel group

# Create a copy to manipulate
tmp = all_adata.obs[endothelial + ['sample_id']].copy()

# Prepend the sample_id to the values in each column
# Only if the value is not NaN
for col in endothelial:  # Iterate only over endothelial columns
    tmp[col] = tmp.apply(lambda row: f"{row['sample_id']}-{row[col]}" if pd.notna(row[col]) else row[col], axis=1)

# Drop the sample_id column after prepending
tmp = tmp.drop(columns=['sample_id'])

# Replace each value with "colname-value" if not NaN
for col in endothelial:
    tmp[col] = tmp[col].apply(lambda x: f"{col}-{x}" if pd.notna(x) else x)

# Collapse across columns by taking the first non-null entry per row
all_adata.obs['final_group'] = tmp.apply(lambda row: next((val for val in row if pd.notna(val)), None), axis=1)

# Get an unnested list of genes
genes_of_interest = [gene for sublist in gene_dict.values() for gene in sublist if gene in all_adata.var_names]
expr_df = all_adata.to_df(layer=use_data)[genes_of_interest]

# Add composite group column
expr_df['group'] = all_adata.obs['final_group'].values

# Compute mean expression per composite group
grouped_expr = expr_df.groupby('group').mean()

adata_grouped = ad.AnnData(X=grouped_expr.values)
adata_grouped.obs_names = grouped_expr.index
adata_grouped.var_names = grouped_expr.columns
vessel_classes_df = grouped_expr.index.to_series().str.split('-', expand=True)[0]

adata_grouped.obs['vessel_class'] = vessel_classes_df


# --- Plot 1: Heatmap of marker gene expression per vessel group

# Define a plotting order
vessel_class_order = ['arterial', 'capillary', 'venous', 'lymphatic']
# Sort the gene_dict keys to desired order
gene_dict = {k: gene_dict[k] for k in vessel_class_order if k in gene_dict}

# Also sort the adata_grouped obs
plot_order = ['CXCL12', 'RGS5', 'ACKR1', 'LYVE1']

# Sort endothelial based on substrings in plot_order
sorted_endothelial = sorted(
    endothelial,
    key=lambda x: plot_order.index(next((g for g in plot_order if g in x), None))
)

# Reorder the adata_grouped obs based on sorted_endothelial
adata_grouped.obs['vessel_class'] = pd.Categorical(adata_grouped.obs['vessel_class'], categories=sorted_endothelial, ordered=True)

# Plot the heatmap
dp = sc.pl.matrixplot(adata_grouped, var_names=gene_dict, groupby=['vessel_class'], return_fig=True,
                     standard_scale='var', )
dp.add_totals().show()

# Adjust the size of the axes
axs = plt.gcf().get_axes()

heatmap_ax = axs[0]

# Set the desired ratio
heatmap_ax.set_aspect(1.75)  # Adjust the aspect ratio as needed

# Update ytick labels
ytick_labels = heatmap_ax.get_yticklabels()
new_labels = []

for label in ytick_labels:
    txt = label.get_text()
    new_text = txt.split(' ')[-1]
    label.set_text(new_text)
    new_labels.append(label)

heatmap_ax.set_yticklabels(new_labels)

# Adjust the other axs
new_heatmap_pos = heatmap_ax.get_position()

bar_ax = axs[1]
bar_ax_pos = bar_ax.get_position()
bar_ax.set_position([new_heatmap_pos.x1, new_heatmap_pos.y0, bar_ax_pos.width, new_heatmap_pos.height])

# Move the cbar
cbar_ax = axs[-1]
bar_ax_pos = bar_ax.get_position()
cbar_ax_pos = cbar_ax.get_position()
cbar_ax.set_position([bar_ax_pos.x1+0.02, bar_ax_pos.y0, cbar_ax_pos.width*0.8, cbar_ax_pos.height])

# Make the top ax invisible
top_ax = axs[2]
top_ax_pos = top_ax.get_position()
top_ax.set_position([new_heatmap_pos.x0, top_ax_pos.y0, new_heatmap_pos.width, top_ax_pos.height])

# Rotate the top_ax labels
for label in top_ax.texts:
    label.set_text(label.get_text().capitalize())
    label.set_rotation(0)
    label.set_ha('center')  # Optional: align center

# Calculate position for plot title
title_x = (new_heatmap_pos.x0 + new_heatmap_pos.x1) / 2
title_y = top_ax.get_position().y1 + 0.2

plt.suptitle(f'All samples - Endothelial marker genes (Kumar 2023) - mean expression per vessel group (use_data={use_data})', y=title_y, x=title_x, fontweight='bold')


# --- Save the output

out_dir = os.path.join(res_dir, 'figures',)
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_3D.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')


# --- Finish

print("Script successfully completed")
