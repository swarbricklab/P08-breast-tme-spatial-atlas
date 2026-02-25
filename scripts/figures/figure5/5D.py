#!/usr/bin/env python3

# --- Aim
# Plot Figure 5D: ligand-receptor interactions for HEV-upregulated genes
# Shows dotplot of LR pairs enriched in HEVs compared to ACKR1+ ECs across HEV compartments
# Uses bin_size=0 (HEV bins only)

# --- Define environment
# env: vis-hd

# --- Load packages

import argparse
import os
import sys
from datetime import datetime

import anndata as ad
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Plot ligand-receptor interactions for HEV DEGs",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--sample_id", help="Sample ID (not used but required for compatibility)")
parser.add_argument("--project_dir", help="Path to the root of the project directory")
parser.add_argument("--min_count", help="Minimal number of counts per spot", type=int)
parser.add_argument("--min_gene", help="Minimal number of genes per spot", type=int)
parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)
parser.add_argument("--samples", help="List of sample IDs to analyze", nargs='+', default=None)

args = vars(parser.parse_args())

# Arguments
project_dir = args['project_dir']
res_dir = args['res_dir']
sample_list = args['samples']

print(f'Samples: {sample_list}')


# --- Import custom plotting functions

# Import function script as a module
sys.path.append(os.path.join(project_dir, 'scripts', '04_vascular_microenvironment'))
import plotting_functions


# --- Parameters

threshold = 0.5
threshold_adj = str(threshold).split('.')[-1]
bin_size = 0  # um - perivascular bins
resolution = 'square_016um'


# --- Function to find most recent results

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
        return None


# --- Load the adata objects with LR results

# Find the directory with individual sample results
lr_sample_dir = os.path.join(res_dir, '04_vascular_microenvironment', '04_ligand_receptor_analysis',
                              '01_run_LR_analysis_stlearn', f'{bin_size}um_bin', 'csv_files')

most_recent_sample_path = find_most_recent(lr_sample_dir)

# Load adata objects for each sample
adata_dict = {}

for sample in sample_list:
    adata_name = f'{sample}_HEV_LR_analysis_stlearn.h5ad'
    adata_dict[sample] = ad.read_h5ad(os.path.join(most_recent_sample_path, adata_name))


# --- Load LR summary data

lr_dir = os.path.join(res_dir, '04_vascular_microenvironment', '04_ligand_receptor_analysis', 
                        '02_analyze_LR_all', f'{bin_size}um_bin', 'csv_files')

most_recent_path = find_most_recent(lr_dir)

# Load LR summary
lr_summary_file = 'All_samples_LR_summary.csv'
lr_summary = pd.read_csv(os.path.join(most_recent_path, lr_summary_file), index_col=0)

# Load LFC data
lr_lfc_file = 'All_samples_group_counts.csv'
lr_lfc = pd.read_csv(os.path.join(most_recent_path, lr_lfc_file), index_col=0)

# Sort the lr_list by LFC
lr_list = lr_lfc.sort_values(by='LFC', ascending=False).index.tolist()

# Extract list of lrs found in all samples
lr_shared = lr_summary[lr_summary['n_samples_sig'] == len(sample_list)].index.tolist()

# Subset to shared LRs
lr_list = [lr for lr in lr_list if lr in lr_shared]

print(f"Found {len(lr_shared)} LR pairs significant in all {len(sample_list)} samples")


# --- Build pivot_df from adata objects

# Extract cell type counts per sample
all_rows = []

for sample_id, adata in adata_dict.items():
    for lr_pair, df in adata.uns['per_lr_cci_cell_type'].items():
        for source in df.index:
            for target in df.columns:
                count = df.loc[source, target]
                if count > 0:  # skip zeroes if desired
                    all_rows.append({
                        'LR_pair': lr_pair,
                        'source_cell_type': source,
                        'target_cell_type': target,
                        'sample_id': sample_id,
                        'n_interactions': count
                    })

flat_df = pd.DataFrame(all_rows)

# Pivot to wide format
pivot_df = flat_df.pivot_table(
    index=['LR_pair', 'source_cell_type', 'target_cell_type'],
    columns='sample_id',
    values='n_interactions',
    fill_value=0
)

# Optional: flatten the multi-index columns
pivot_df.columns.name = None
pivot_df = pivot_df.reset_index()

# Sum across all sample columns
sample_cols = [col for col in pivot_df.columns if col not in ['LR_pair', 'source_cell_type', 'target_cell_type']]
pivot_df['total_n_interactions'] = pivot_df[sample_cols].sum(axis=1)


# --- Load HEV colors

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
color_file = 'TLS_HEV_colors.csv'

color_df = pd.read_csv(os.path.join(color_dir, color_file), index_col=0)

# Make a color_dict
color_dict = {}
for label, color in zip(color_df.index, color_df['color']):
    color_dict[label] = color


# --- Prepare rename dict for labels

rename_dict = {
    'near_TLS_HEV': 'HEV (TLS)',
    'non_TLS_HEV': 'HEV-like (stroma)',
    'nearby_spots': f'Perivascular bins' if bin_size > 0 else 'Non-HEV',
}


# --- Load DE results (HEV vs ACKR1 ECs)

de_dir = os.path.join(res_dir, '04_vascular_microenvironment', '03_compare_HEV_to_other_venous', 
                      '02_compare_gene_expression', 'csv_files')

most_recent_de_path = find_most_recent(de_dir)

de_filename = 'All_samples_rank_genes_HEV_overlap.csv'
de_df = pd.read_csv(os.path.join(most_recent_de_path, de_filename))

# Filter genes enriched in overlap HEV group
min_lfc = 1.25
p_adj_min = 0.05

# Subset de_df to just genes up in HEV
de_df = de_df[de_df['group'] == 'overlap']
de_df_filt = de_df[(de_df['logfoldchanges'] > min_lfc) & (de_df['pvals_adj'] < p_adj_min)]
de_genes = de_df_filt['genes'].tolist()

print(f"Found {len(de_genes)} DE genes enriched in HEV")


# --- Filter LR pairs involving DE genes

lr_interest = []

for lr in lr_list:
    if any(gene in lr for gene in de_genes):
        lr_interest.append(lr)

print(f"Found {len(lr_interest)} LR pairs involving HEV DE genes")


# --- Prepare data for plotting

df_sub = pivot_df[pivot_df['LR_pair'].isin(lr_interest)]

# Create the combined source→target column
df_sub['source_target'] = df_sub['source_cell_type'] + '→' + df_sub['target_cell_type']

# Get all unique sources and targets
sources = df_sub['source_cell_type'].unique()
targets = sorted(df_sub['target_cell_type'].unique())

# Build the desired order
ordered_source_targets = []
for source in sources:
    for target in targets:
        ordered_source_targets.append(f"{source}→{target}")

# Set as categorical with the desired order
# Subset to only categories that are in the df_sub
ordered_source_targets = [st for st in ordered_source_targets if st in df_sub['source_target'].unique()]
# Set the categorical order
df_sub['source_target'] = pd.Categorical(df_sub['source_target'], categories=ordered_source_targets, ordered=True)

merged = df_sub

# Reorder based on LR pair list
merged['LR_pair'] = pd.Categorical(merged['LR_pair'], categories=lr_interest, ordered=True)


# --- Setup plotting parameters

size_per_spot = 0.1  # Size of each dot in the dotplot
cmap = 'magma_r'
xval = 'source_target'
yval = 'LR_pair'

size_max = merged['total_n_interactions'].max()

# Disable grid
plt.rcParams['axes.grid'] = False


# --- Create the plot

nrows = 1
ncols = 1

fig, ax = plt.subplots(ncols=ncols, nrows=nrows, figsize=(3, 3))

# Extract number of cell types and LRs
n_cell_types = len(merged[xval].unique())
n_lrs = len(merged[yval].unique())

ax_pos = ax.get_position()
ax_pos = ax.set_position([ax_pos.x0, ax_pos.y0, n_cell_types * size_per_spot, n_lrs * size_per_spot*0.8])

dotplot_ax, _ = plotting_functions.seaborn_dotplot(
    merged,
    x=xval,
    y=yval,
    hue='total_n_interactions',
    size='total_n_interactions',
    cmap=cmap,
    dotplot_ax=ax,
    max_size_px=300,
    size_max=size_max,
)

# Adjust limits of plot
ax.set_xlim(-0.75, len(merged[xval].unique()) - 0.25)
ax.set_ylim(len(merged[yval].unique()) - 0.25, -0.75)

# Remove the legend if it exists
if dotplot_ax.legend_ is not None:
    dotplot_ax.legend_.remove()

# Add labels
ax.set_xlabel('Direction', fontweight='bold')
ax.set_ylabel('LR pair', fontweight='bold')

# Set xticklabels
dotplot_ax.set_xticklabels(dotplot_ax.get_xticklabels(), rotation=90, ha='center')


# --- Add colorbar

norm = mpl.colors.Normalize(vmin=merged['total_n_interactions'].min(), vmax=merged['total_n_interactions'].max())
sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

cbar_ax = fig.add_axes([0.91, 0.15, 0.02, 0.7])
cbar_ax.grid(False)
cbar = plt.colorbar(sm, cax=cbar_ax)
cbar.set_label('# of interactions', rotation=270, labelpad=15)

# Set colorbar position
main_ax_pos = ax.get_position()
cbar_ax.set_position([main_ax_pos.x1 + size_per_spot*0.2, main_ax_pos.y0, size_per_spot*0.5, main_ax_pos.height])


# --- Add rectangles showing source→target relationships

# Create a new axis above the main axis for rectangles
rect_ax = fig.add_axes([ax.get_position().x0, ax.get_position().y1 + 0.02, 
                        ax.get_position().width, ax.get_position().height/len(ax.get_yticklabels())*1.25])
rect_ax.set_xlim(ax.get_xlim())

# Add rectangles according to xticklabels
for i, label in enumerate(ax.get_xticklabels()):
    label_text = label.get_text()
    # Split on the arrow
    if '→' in label_text:
        source, target = label_text.split('→')
    else:
        continue
    
    # Top rectangle (source)
    if source in color_dict:
        color = color_dict[source]
        rect_ax.add_patch(
            patches.Rectangle((i - 0.35, 0.7), 0.75, 0.3, facecolor=color, 
                            transform=rect_ax.transData, edgecolor='black', linewidth=0.8)
        )
    
    # Bottom rectangle (target)
    if target in color_dict:
        color = color_dict[target]
        rect_ax.add_patch(
            patches.Rectangle((i - 0.35, 0), 0.75, 0.3, facecolor=color, 
                            transform=rect_ax.transData, edgecolor='black', linewidth=0.8)
        )
    
    # Add an arrow from top rectangle to bottom rectangle
    rect_ax.annotate('', 
        xy=(i, 0.35), xytext=(i, 0.7),
        arrowprops=dict(facecolor='black', width=0.5, headwidth=5, headlength=4),
        ha='center'
    )

# Set title for rectangle axis
rect_ax.set_title('LRs involving HEV DEGs', fontweight='bold')

# Make the spines of rect_ax invisible
rect_ax.spines['top'].set_visible(False)
rect_ax.spines['right'].set_visible(False)
rect_ax.spines['left'].set_visible(False)
rect_ax.spines['bottom'].set_visible(False)

# Remove ticks and labels from the rectangle axis
rect_ax.set_xticks([])
rect_ax.set_yticks([])


# --- Create legend

# Figure out unique labels from xticklabels
legend_labels = []
for label in ax.get_xticklabels():
    label_text = label.get_text()
    if '→' in label_text:
        source, target = label_text.split('→')
        if source not in legend_labels:
            legend_labels.append(source)
        if target not in legend_labels:
            legend_labels.append(target)

# Make handles based on color_dict
handles = []
for label in legend_labels:
    if label in color_dict:
        color = color_dict[label]
        handles.append(patches.Patch(facecolor=color, label=label, edgecolor='black', linewidth=0.8))

# Add legend
x_offset = 1/len(ax.get_xticklabels())
rect_ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1+x_offset, 0.5),
               frameon=False, handlelength=0.7, handletextpad=0.5, borderpad=0)

# Update the labels of the legend
for text in rect_ax.get_legend().get_texts():
    if text.get_text() in rename_dict:
        text.set_text(rename_dict[text.get_text()])


# --- Update xticklabels with renamed labels

new_xticklabels = []
for label in dotplot_ax.get_xticklabels():
    label_text = label.get_text()
    # Labels to update will be substrings
    for key, value in rename_dict.items():
        if key in label_text:
            label_text = label_text.replace(key, value)
    # Also make sure to put space around the arrow
    if '→' in label_text:
        source, target = label_text.split('→')
        label_text = f'{source} → {target}'
    new_xticklabels.append(label_text)

dotplot_ax.set_xticklabels(new_xticklabels)


# --- Add main title

# Calculate the midpoint to position the title
title_x = (ax.get_position().x0 + ax.get_position().x1)/2
title_y = ax.get_position().y1 + 0.4

bin_text = f'{bin_size}µm' if bin_size > 0 else 'no bins'
plt.suptitle(f'LR interactions involving HEV vs ACKR1 DEGs ({bin_text})', 
             y=title_y, x=title_x, fontsize=20, fontweight='bold')


# --- Save plot

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_5D.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')


# --- Finish

print("Script successfully completed")
print(f"Saved plot to: {os.path.join(out_dir, filename)}")
