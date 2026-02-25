#!/usr/bin/env python3

# --- Aim
# Plot Figure 3I: differential gene expression in HEV-overlapping vs non-overlapping vessels
# Shows volcano plot and dotplot of genes enriched in vessels with HEV signatures
# Stratified by vessel type based on deconvolution results

# --- Define environment
# env: shapely (has adjustText)

# --- Load packages

import argparse
import os
from datetime import  datetime

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from adjustText import adjust_text


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

decon_level = 'minor'

# To make it easier to visualize trends by lineage
class_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
class_file = 'miniatlas_subset_classifications.csv'

class_df = pd.read_csv(os.path.join(class_dir, class_file))

# Get endothelial classes
endothelial = list(class_df[class_df['miniatlas_major'] == 'Endothelial'][f'miniatlas_{decon_level}'].unique())
endothelial = sorted(endothelial)


# -- Load data

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
        return most_recent_dir_path
    else:
        print("No valid directories found.")

# Define which endothelial type want to look at
endo_of_interest = 'Endothelial ACKR1'
out_name = endo_of_interest.split(' ')[-1]

deg_dir = os.path.join(project_dir, 'results', '04_vascular_microenvironment', '03_compare_HEV_to_other_venous', '02_compare_gene_expression', 'csv_files')
most_recent_dir_path = find_most_recent(deg_dir)

use_data = 'log_norm'

# Load the adata object with the DEG analysis
filename = f'All_samples_HEV_{out_name}.h5ad'
adata = ad.read_h5ad(os.path.join(most_recent_dir_path, filename))

# Also load the DEG list
deg_file = 'All_samples_rank_genes_HEV_overlap.csv'
deg_df = pd.read_csv(os.path.join(most_recent_dir_path, deg_file))

# Take the the overlap HEV group only
deg_df = deg_df[deg_df['group'] == 'overlap']


# --- Identify genes that are not expressed

# Remove the genes that are not expressed in many cell types
# Filter genes that are expressed in at least 50 cells
min_cells = 50
test = sc.pp.filter_genes(adata, min_cells=min_cells, copy=True)
genes_to_keep = list(test.var_names)


# --- Functions for plotting

def adjust_text_func(ax, sorted_df, x_val, y_val, min_distance=20, labels_to_annotate=None, lw=0.2, color='gray', gene_list=None,
                     label_size=5):
    texts = []
    annotated_coords = []
    for i, row in sorted_df.iterrows():
        label = row['genes']
        # If a list of labels is provided, only annotate those
        if labels_to_annotate is not None and label not in labels_to_annotate:
            continue
        x = row[x_val]
        y = row[y_val]
        # Check distance to all previous annotated points
        too_close = False
        for x0, y0 in annotated_coords:
            dist = np.sqrt((x - x0)**2 + (y - y0)**2)
            if dist < min_distance:
                too_close = True
                break
        if not too_close:
            if gene_list is not None and label in gene_list:
                # If the label is in the gene list, use a different color
                text_color = 'black'
                fontweight = 'bold'
            else:
                text_color = 'black'
                fontweight = 'normal'
            text = ax.text(x, y, label, fontsize=label_size, ha='right', va='bottom', color=text_color, fontweight=fontweight)
            texts.append(text)
            annotated_coords.append((x, y))
    # Now bring text in front of the arrows
    for text in texts:
        text.set_zorder(200)  # Bring text to front
    # Adjust remaining text to avoid overlap
    if texts:
        adjust_text(
            texts,
            ax=ax,
            arrowprops=dict(arrowstyle="-", color=color, lw=lw),
            expand_points=(1.2, 1.4),
            expand_text=(1.2, 1.4),
            force_text=(0.5, 0.5)
        )


# --- Load the known HEV markers from sawada

gene_dir = os.path.join(project_dir, 'resources', 'gene_lists', 'HEV_markers')
filename =  'HEV_sawada_2022.csv'

gene_table = pd.read_csv(os.path.join(gene_dir, filename))

# Extract the genes
HEV_genes = []
if 'gene' in gene_table.columns:
    HEV_genes = list(gene_table.gene)


# --- Plot 1:  Volcano plot of HEV vs non-HEV

# Prepare for plotting
# Log transform the p-values
deg_df['-np.log10(pvals_adj)'] = -np.log10(deg_df['pvals_adj'])

# Filter any genes that aren't expressed in the HEVs / vessels
deg_df_filtered = deg_df[deg_df['genes'].isin(genes_to_keep)]

# Set grid behind plot
plt.rcParams['axes.axisbelow'] = True

# Set grid color
plt.rcParams['grid.color'] = 'lightgray'

# -- Plotting

# Thresholds for plotting
LFC_threshold = 0.5

# Annotation threshold (LFC)
annotate_threshold = 1

# Additional genes to always annotate
genes_of_interest = 'TMEM173'

# Annotate genes with very low pvals even if LFC is low
p_add_threshold = 10

fig, ax = plt.subplots(figsize=(3, 3))

# Force ax to be square
ax_pos = ax.get_position()
ax.set_position([ax_pos.x0, ax_pos.y0, ax_pos.width, ax_pos.width])

# First plot non significant in gray
deg_df_non_sig = deg_df_filtered[(deg_df_filtered['-np.log10(pvals_adj)'] < -np.log10(0.05)) | (deg_df_filtered['logfoldchanges'].abs() < LFC_threshold)]
sns.scatterplot(data=deg_df_non_sig, x='logfoldchanges', y='-np.log10(pvals_adj)', ax=ax, s=5, alpha=1, linewidth=0, color='gray')

# Add a horizontal line for significance threshold
ax.axhline(y=-np.log10(0.05), color='red', linestyle='--', linewidth=0.8, zorder=1)

# Add a vertical line for fold change threshold
ax.axvline(x=LFC_threshold, color='red', linestyle='--', linewidth=0.8, zorder=1)
ax.axvline(x=-LFC_threshold, color='red', linestyle='--', linewidth=0.8, zorder=1)

# Then plot the significant ones in color
deg_sig = deg_df_filtered[(deg_df_filtered['-np.log10(pvals_adj)'] > -np.log10(0.05)) & (abs(deg_df_filtered['logfoldchanges']) > LFC_threshold)]
sns.scatterplot(data=deg_sig, x='logfoldchanges', y='-np.log10(pvals_adj)', ax=ax, s=10, alpha=1, linewidth=0)

# Set the x axis equal
xlim = ax.get_xlim()
xmax = max(abs(xlim[0]), abs(xlim[1]))
ax.set_xlim(-xmax, xmax)

# Create a df with genes that should be annotated
annotate_df = deg_df_filtered[(deg_df_filtered['-np.log10(pvals_adj)'] > -np.log10(0.05)) & (abs(deg_df_filtered['logfoldchanges']) > annotate_threshold)]

# Also add genes of interest
annotate_df = pd.concat([annotate_df, deg_df_filtered[deg_df_filtered['genes'].isin([genes_of_interest])]])

# Also annotate genes with very low pvals even if LFC is low
annotate_df = pd.concat([annotate_df, deg_df_filtered[(deg_df_filtered['-np.log10(pvals_adj)'] > p_add_threshold) & (abs(deg_df_filtered['logfoldchanges']) > LFC_threshold)]])

min_distance = 0.1
label_size = 5
adjust_text_func(ax, annotate_df, x_val='logfoldchanges', y_val='-np.log10(pvals_adj)', min_distance=min_distance, lw=0.5, color='gray',
                  gene_list=HEV_genes, label_size=label_size)

# Give genes of interest a black outline around dots to make them more visible
for i, row in annotate_df.iterrows():
    label = row['genes']
    if label in [genes_of_interest]:
        x = row['logfoldchanges']
        y = row['-np.log10(pvals_adj)']
        ax.scatter(x, y, s=10, facecolors='none', edgecolors='#0f3c5a', linewidths=0.5, zorder=3)

# Adjust the labels
ax.set_xlabel('LFC(HEV vs non-HEV)')
ax.set_ylabel('-log10(p-val)')

# Add arrows just below the text
arrowprops = dict(arrowstyle='-|>', color='black', lw=1, shrinkA=0, shrinkB=0, joinstyle='miter')
arrow_y = 1.03
ax.annotate('', xy=(1, arrow_y), xytext=(0.7, arrow_y), xycoords='axes fraction', textcoords='axes fraction', arrowprops=arrowprops, clip_on=False)
ax.annotate('', xy=(0, arrow_y), xytext=(0.3, arrow_y), xycoords='axes fraction', textcoords='axes fraction', arrowprops=arrowprops, clip_on=False)

# Add text just above and centered the arrows
midpoint = 0.3/2
label_size = 10
buffer = 0.03
ax.text(1-midpoint, 1.04, 'HEV', fontsize=label_size, ha='center', va='bottom', color='black', clip_on=False, transform=ax.transAxes)    
ax.text(0+midpoint, 1.04, 'Venous', fontsize=label_size, ha='center', va='bottom', color='black', clip_on=False, transform=ax.transAxes)

# Increase the xlim 10% to accommodate labels
xlim = ax.get_xlim()
adj = 1.05
xlim = (xlim[0] * adj, xlim[1] * adj)
ax.set_xlim(xlim)
ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1] * 1.02)

# Enable grid
ax.grid(True, linewidth=0.5, alpha=0.7)

# Add a title
plt.suptitle(f'HEV vs non-HEV gene expression - {out_name} (min_cells={min_cells}, bold = Sawada genes, LFC_threshold={LFC_threshold})', y=1.1, fontweight='bold')


# --- Save the output

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_3I.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')


# --- Finish

print("Script successfully completed")