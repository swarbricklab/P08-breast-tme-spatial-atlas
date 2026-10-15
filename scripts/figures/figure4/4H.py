#!/usr/bin/env python3

# --- Aim
# Plot dotplot showing top LR pairs enriched in TLS-HEVs (Figure 4H)
# Shows the source→target cell type combinations for each LR pair

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import os
import sys

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments
parser = argparse.ArgumentParser(description="Plot top LR pairs dotplot for Figure 4H",
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


# --- Import custom plotting functions

# Import function script as a module
sys.path.append(os.path.join(project_dir, 'scripts', '04_vascular_microenvironment'))
import plotting_functions


# --- Extract some info from the sample sheet
# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_visium_hd.csv')

sample_sheet = pd.read_csv(sample_sheet_path)
sample_lut = sample_sheet.set_index('section_id')


# --- Parameters

threshold = 0.5
threshold_adj = str(threshold).split('.')[-1]
bin_size = 25  # um
resolution = 'square_016um'


# --- Load the adata objects with LR results

lr_sample_dir = os.path.join(res_dir, '04_vascular_microenvironment', '04_ligand_receptor_analysis',
                              '01_run_LR_analysis_stlearn', f'{bin_size}um_bin', 'csv_files')

most_recent_sample_path = lr_sample_dir

adata_dict = {}

for sample in sample_list:
    filename = f'{sample}_HEV_LR_analysis_stlearn.h5ad'
    adata_dict[sample] = ad.read_h5ad(os.path.join(most_recent_sample_path, filename))


# --- Load LR summary data

lr_dir = os.path.join(res_dir, '04_vascular_microenvironment', '04_ligand_receptor_analysis', 
                        '02_analyze_LR_all', f'{bin_size}um_bin', 'csv_files')

most_recent_path = lr_dir

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
    'nearby_spots': f'Perivascular bins',
}


# --- Select top LR pairs for plotting

top = 15
lr_sub = lr_list[:top]

df_sub = pivot_df[pivot_df['LR_pair'].isin(lr_sub)]

# Create the combined source→target column
df_sub['source_target'] = df_sub['source_cell_type'] + '→' + df_sub['target_cell_type']

# Manually define the desired order
ordered_source_targets = [
    'near_TLS_HEV→near_TLS_HEV',
    'near_TLS_HEV→nearby_spots',
    'nearby_spots→near_TLS_HEV',
    'nearby_spots→nearby_spots',
    'nearby_spots→non_TLS_HEV',
    'non_TLS_HEV→non_TLS_HEV',
]

# Subset to only categories that are actually in the df_sub
ordered_source_targets = [st for st in ordered_source_targets if st in df_sub['source_target'].unique()]

# Set the categorical order
df_sub['source_target'] = pd.Categorical(df_sub['source_target'], categories=ordered_source_targets, ordered=True)

merged = df_sub

# Reorder based on LR pair list
merged['LR_pair'] = pd.Categorical(merged['LR_pair'], categories=lr_sub, ordered=True)

size_max = merged['total_n_interactions'].max()


# --- Setup plotting parameters

size_per_spot = 0.1  # Size of each dot in the dotplot
cmap = 'magma_r'
xval = 'source_target'
yval = 'LR_pair'

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
ax_pos = ax.set_position([ax_pos.x0, ax_pos.y0, top * size_per_spot, n_cell_types * size_per_spot*0.8])

dotplot_ax, _ = plotting_functions.seaborn_dotplot(
    merged,
    x=yval,
    y=xval,
    hue='total_n_interactions',
    size='total_n_interactions',
    cmap=cmap,
    dotplot_ax=ax,
    max_size_px=300,
    size_max=size_max,
)

# Adjust limits of plot
ax.set_xlim(-0.75, len(merged[yval].unique()) - 0.25)
ax.set_ylim(len(merged[xval].unique()) - 0.25, -0.75)

# Remove the legend
dotplot_ax.legend_.remove()

# Add labels
ax.set_xlabel('LR pair', fontweight='bold')
ax.set_ylabel('', fontweight='bold')

# Set xticklabels
ax.set_xticklabels(dotplot_ax.get_xticklabels(), rotation=90, ha='center')


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

# Create a new axis on the left side of the main axis
rect_ax = fig.add_axes([ax.get_position().x0 - 0.02 - size_per_spot*1, ax.get_position().y0, 
                        size_per_spot*1, ax.get_position().height])
rect_ax.set_ylim(ax.get_ylim())

# Add rectangles according to yticklabels
# Add small rectangles on the left to indicate the group labels
for i, label in enumerate(ax.get_yticklabels()):
    label_text = label.get_text()
    # Split on the arrow
    if '→' in label_text:
        source, target = label_text.split('→')
    else:
        continue
    
    # Left rectangle (source)
    if source in color_dict:
        color = color_dict[source]
        rect_ax.add_patch(
            patches.Rectangle((0, i - 0.35), 0.3, 0.75, facecolor=color, 
                            transform=rect_ax.transData, edgecolor='black', linewidth=0.8)
        )
    
    # Right rectangle (target)
    if target in color_dict:
        color = color_dict[target]
        rect_ax.add_patch(
            patches.Rectangle((0.7, i - 0.35), 0.3, 0.75, facecolor=color, 
                            transform=rect_ax.transData, edgecolor='black', linewidth=0.8)
        )
    
    # Add an arrow from left rectangle to right rectangle
    rect_ax.annotate('', 
        xy=(0.7, i), xytext=(0.35, i),
        arrowprops=dict(facecolor='black', width=0.5, headwidth=5, headlength=4),
        ha='center'
    )

# Make the spines of rect_ax invisible
rect_ax.spines['top'].set_visible(False)
rect_ax.spines['right'].set_visible(False)
rect_ax.spines['left'].set_visible(False)
rect_ax.spines['bottom'].set_visible(False)

# Remove ticks and labels from the rectangle axis
rect_ax.set_xticks([])

# Take the y-ticks from the main axis
rect_ax.set_yticks(ax.get_yticks())
rect_ax.set_yticklabels(ax.get_yticklabels())

# Remove y-tick labels from the main axis
ax.set_yticks([])


# --- Create legend

# Figure out the labels from the yticklabels
legend_labels = []

for label in rect_ax.get_yticklabels():
    label_text = label.get_text()
    # Split on the arrow
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

# Add legend to the axis
x_offset = 1/len(ax.get_xticklabels())*3  # Adjust x_offset as needed
ax.legend(handles=handles, loc='lower left', bbox_to_anchor=(0, 1),
          frameon=False, handlelength=0.7, handletextpad=0.5, borderpad=0)

# Update the labels of the legend
for text in ax.get_legend().get_texts():
    if text.get_text() in rename_dict:
        text.set_text(rename_dict[text.get_text()])

# Update the yticklabels of the rect_ax
new_yticklabels = []
for label in rect_ax.get_yticklabels():
    label_text = label.get_text()
    # Labels to update will be substrings
    for key, value in rename_dict.items():
        if key in label_text:
            label_text = label_text.replace(key, value)
    # Also make sure to put space around the arrow
    if '→' in label_text:
        source, target = label_text.split('→')
        label_text = f'{source} → {target}'
    new_yticklabels.append(label_text)

rect_ax.set_yticklabels(new_yticklabels)

# Calculate the midpoint to position the title
title_x = (ax.get_position().x0 + ax.get_position().x1)/2
title_y = ax.get_position().y1 + 0.3

plt.suptitle(f'Top {top} LR pairs enriched in TLS-HEVs', 
             y=title_y, x=title_x, fontweight='bold')


# --- Save plot

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_4H.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')


# --- Finish

print("Script successfully completed")
print(f"Saved plot to: {os.path.join(out_dir, filename)}")
