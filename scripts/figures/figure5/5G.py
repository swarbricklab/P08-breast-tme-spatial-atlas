#!/usr/bin/env python3

# --- Aim
# Plot Figure 5G: gene expression in cells inside vs outside TLS regions
# Shows cell type-specific expression patterns comparing TLS interior with surrounding tissue

# --- Define environment
# env: xen-hev

# --- Load packages

import argparse
import os
import anndata as ad
import scanpy as sc

import pandas as pd

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches

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


# --- Load the anndata objects

adata_dict = {}

for sample in sample_list:
    adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample}.h5ad'))
    adata_dict[sample] = adata


# --- Load the tangram annotations

# Decide which annotation level to use
level = 'subset'

# Initialize a df to store all annotations
anno_df = pd.DataFrame()

for sample in sample_list:
    tmp_df = pd.DataFrame()
    annotation_file = f'{sample}_{level}_pred_cell_types.csv'
    annotation_dir = os.path.join(project_dir, 'results', '02_cell_type_mapping', 'label_transfer', 'xenium')
    annotations = pd.read_csv(os.path.join(annotation_dir, f'{level}_level', annotation_file), index_col=0)
    # Add to adata
    adata_dict[sample].obs[level] = annotations['pred_cell_type']
    # Prepend sample id to index
    annotations.index = [f'{sample}_{idx}' for idx in annotations.index]
    tmp_df[level] = annotations['pred_cell_type']
    anno_df = pd.concat([anno_df, tmp_df])

# Add sample_id column
# Split on the last _ and take the first part to get the sample id
anno_df['sample_id'] = [idx.rsplit('_', 1)[0] for idx in anno_df.index]


# --- Load plotting colors

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')

# Use the major level for colors
color_level = 'major'
color_file = f'miniatlas_{color_level}_colors.csv'

color_dict = pd.read_csv(os.path.join(color_dir, color_file), index_col=0).to_dict()['color']


# --- Load the miniatlas classifications

# These are needed to check the alignment between major and minor levels
miniatlas_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
miniatlas_file = 'miniatlas_subset_classifications.csv'

miniatlas_df = pd.read_csv(os.path.join(miniatlas_dir, miniatlas_file),)

# Add the other levels to the color_dict
for level in ['minor', 'subset']:
    # For each cell type in the miniatlas at this level, check if it is in the color_dict
    # If not, get the corresponding major cell type and use its color
    for cell_type in miniatlas_df[f'miniatlas_{level}'].unique():
        if cell_type not in color_dict.keys():
            # Check that the value is not nan
            if pd.isna(cell_type):
                continue
            print(f'Adding color for {cell_type}')
        # Get the corresponding major cell type
        major_type = miniatlas_df[miniatlas_df[f'miniatlas_{level}'] == cell_type][f'miniatlas_major'].values
        # Find which value is in the color_dict
        major_type_in_dict = [m for m in major_type if m in color_dict.keys()]
        if len(major_type_in_dict) > 0:
            color_dict[cell_type] = color_dict[major_type_in_dict[0]]
        else:
            color_dict[cell_type] = 'lightgrey'

# Also define a dict with shorter names for the subset level cell types
miniatlas_sub_dict = miniatlas_df.set_index('miniatlas_subset')[['miniatlas_subset_short']].to_dict()


# --- Load the TLS annotations

tls_dir = os.path.join(res_dir, '05_xenium_mechanisms', '01_identify_TLS_regions', 'csv_files')

tls_df = pd.DataFrame()

for sample in sample_list:
    tls_file = f'{sample}_cell_TLS_annotations.csv'
    tmp_df = pd.read_csv(os.path.join(tls_dir, tls_file), index_col=0)
    # Add to adata
    adata_dict[sample].obs['in_TLS'] = tmp_df['in_TLS']
    # Prepend sample id to index
    tmp_df.index = [f'{sample}_{idx}' for idx in tmp_df.index]
    tls_df = pd.concat([tls_df, tmp_df])

# Add sample_id column
tls_df['sample_id'] = [idx.rsplit('_', 1)[0] for idx in tls_df.index]

# Add the annotations to the tls_df
tls_df = tls_df.merge(anno_df[[level]], left_index=True, right_index=True, how='left')


# --- Sort by lineage, then by major

lineage_order = ['Epithelial', 'Immune', 'Stromal']

major_cell_types = miniatlas_df[f'miniatlas_major'].dropna().unique().tolist()
minor_cell_types = miniatlas_df[f'miniatlas_minor'].dropna().unique().tolist()
subset_cell_types = miniatlas_df[f'miniatlas_subset'].dropna().unique().tolist()


# Sort based on lineage first (recorded in the miniatlas df)
def get_lineage(cell_type, level='major'):
    if cell_type in miniatlas_df[f'miniatlas_{level}'].values:
        lineage = miniatlas_df[miniatlas_df[f'miniatlas_{level}'] == cell_type]['lineage'].values[0]
    else:
        lineage = 'other'
    return lineage

major_cell_types_sorted = sorted(major_cell_types, key=lambda x: (lineage_order.index(get_lineage(x)) if get_lineage(x) in lineage_order else len(lineage_order), x))

# Now sort minor level by corresponding major level
def get_major(cell_type, level, upper_level='major'):
    level_clean = level.replace('_stepwise', '')
    if cell_type in miniatlas_df[f'miniatlas_{level_clean}'].values:
        major = miniatlas_df[miniatlas_df[f'miniatlas_{level_clean}'] == cell_type][f'miniatlas_{upper_level}'].values[0]
    else:
        major = 'other'
    return major

minor_cell_types_sorted = sorted(minor_cell_types, key=lambda x: (major_cell_types_sorted.index(get_major(x, level='minor')) if get_major(x, level='minor') in major_cell_types_sorted else len(major_cell_types_sorted), x))
subset_cell_types_sorted = sorted(subset_cell_types, key=lambda x: (major_cell_types_sorted.index(get_major(x, level='subset')) if get_major(x, level='subset') in major_cell_types_sorted else len(major_cell_types_sorted), x))

cell_type_orders = {
    'major': major_cell_types_sorted,
    'minor': minor_cell_types_sorted,
    'subset': subset_cell_types_sorted
}


# --- Extract a single adata object with all samples

# Extract a single adata object with all samples
adata_list = list(adata_dict.values())
all_adata = ad.concat(adata_list, merge='same')

all_adata.obs_names_make_unique()


# --- Some functions for plotting

# Automatic resizing function
def auto_resize_dotplot(ax, space_per_dot=0.1, width_ratio=1):
    """
    Automatically resize the dotplot axes based on the number of genes.
    """
    # Get the current position of the axes
    pos = ax.get_position()
    # Resize based on xticks & yticks
    xticks = len(ax.get_xticks())
    yticks = len(ax.get_yticks())
    # Calculate new width based on number of genes
    new_width = xticks * space_per_dot * width_ratio
    new_height = yticks * space_per_dot
    # Set the new position
    ax.set_position([pos.x0, pos.y0, new_width, new_height])

def reposition_dotplot_axs(axs, y_offset=0.1, y_offset_within_group=0):
    """
    Reposition dotplot axes so that groups of axes (lists) are separated by y_offset,
    and axes within the same group are separated by y_offset_within_group.
    Accepts axs as a list of axes or a list of lists of axes.
    """
    # Ensure axs is a list of lists
    if not isinstance(axs[0], (list, tuple)):
        axs = [[ax] for ax in axs]
    y_pos = None
    for group in axs:
        for i, ax in enumerate(group):
            pos = ax.get_position()
            if y_pos is None:
                y_new = pos.y0
            else:
                y_new = y_pos  # Start at previous position
            ax.set_position([pos.x0, y_new, pos.width, pos.height])
            # Move down for next ax in group
            y_pos = ax.get_position().y0 - pos.height - (y_offset_within_group if i < len(group)-1 else 0)
        # After last in group, apply y_offset for next group
        y_pos -= y_offset

# Add legends
def position_legend(ax):
    """
    Given a single axis, return the bottom, top, left, and right edges
    in figure coordinates, for robust legend placement.
    Returns (bottom_edge, top_edge, left_edge, right_edge)
    """
    fig = ax.figure
    fig.canvas.draw()
    bbox = ax.get_tightbbox(fig.canvas.get_renderer())
    inv_trans = fig.transFigure.inverted()
    bbox_fig = inv_trans.transform(bbox)
    # bbox_fig: [[x0, y0], [x1, y1]]
    left_edge = bbox_fig[0][0]
    right_edge = bbox_fig[1][0]
    bottom_edge = bbox_fig[0][1]
    top_edge = bbox_fig[1][1]
    return bottom_edge, top_edge, left_edge, right_edge

# Function to add legend
def add_color_legend(fig, labels, color_dict, title='', x=0, y=0, offset=0.01, spacing_per_label=0.05, loc='center left', fontsize=10):
    # Create legend handles (colored boxes with labels)
    handles = [
        mpatches.Patch(facecolor=color_dict[label], label=label, edgecolor='black', linewidth=0.8)
        for label in labels if label in color_dict
    ]
    # Estimate height needed
    est_height = spacing_per_label * len(handles)
    bbox_anchor = (x + offset, y)
    # Add the legend
    legend = fig.legend(
        handles=handles,
        title=title,
        loc=loc,
        bbox_to_anchor=bbox_anchor,
        borderaxespad=0,
        frameon=False,
        handlelength=0.8,
        bbox_transform=fig.transFigure,  # <-- This is the key fix!
        fontsize=fontsize,
    )
    legend.get_title().set_fontweight('bold')
    legend.get_title().set_fontsize(fontsize)  # Set the title fontsize
    legend._legend_box.align = "left"
    # Return updated vertical position for next stacked legend
    return legend

# Prepare for plotting
def prepare_plot(tls_df, annotation_key):
    # Calculate proportions in vs outside of TLS per sample
    sample_prop_all = pd.DataFrame()
    for sample_id in sample_list:
        sample_tls_df = tls_df[tls_df['sample_id'] == sample_id]
        sample_cell_counts = sample_tls_df.groupby([annotation_key, 'in_TLS']).size().unstack(fill_value=0)
        # Calculate proportions
        sample_cell_props = sample_cell_counts.div(sample_cell_counts.sum(axis=1), axis=0)
        # Reverse the column order
        sample_cell_props.columns = sample_cell_props.columns.astype(str)
        sample_cell_props = sample_cell_props[['True', 'False']]
        # Add sample column
        sample_cell_props['sample_id'] = sample_id
        # Reset index
        sample_cell_props = sample_cell_props.reset_index()
        # Concat
        sample_prop_all = pd.concat([sample_prop_all, sample_cell_props])
    # Compute mean and std for each major cell type
    mean_true = sample_prop_all.groupby(annotation_key)['True'].mean()
    std_true = sample_prop_all.groupby(annotation_key)['True'].std()
    # Calculate false fractions based on mean_true
    mean_false = 1 - mean_true
    return mean_true, std_true, mean_false

def plot_heatmap(adata, ax, gene, annotation_key, vmin=0, vmax=1, last_ax=False):
    dp = sc.pl.matrixplot(
        adata,
        var_names=gene,
        groupby=annotation_key,
        ax=ax,
        cmap='viridis',
        vmin=vmin,
        vmax=vmax,
        show=False,
        swap_axes=True,
        standard_scale='var',
    )
    # Adjust the size of the axes
    auto_resize_dotplot(dp['mainplot_ax'], space_per_dot=0.04, width_ratio=0.7)
    heatmap_ax = dp['mainplot_ax']
    # Adjust the ylims
    buffer = 0
    dp['mainplot_ax'].set_ylim(len(dp['mainplot_ax'].get_yticks()) + buffer, -buffer)
    # If not the last plot, make all other axes invisible
    if not last_ax:
        for name, dp_ax in dp.items():
            if not name == 'mainplot_ax':
                dp_ax.set_visible(False)
    cbar_ax = dp['color_legend_ax']
    return heatmap_ax, cbar_ax


# --- Plot 1: Gene expression inside vs outside TLS

# Set plotting fontsize
sc.settings.set_figure_params(fontsize=10)

# Disable the grid
mpl.rcParams['axes.grid'] = False

color_dict_TLS = {
    True: 'cornflowerblue',
    False: 'lightgrey'
}

# Set vmin and vmax for the heatmaps
vmin=0
vmax=1

# Define genes of interest
top_genes_dict = {
    'SPP1-CD44': ['SPP1', 'CD44'],
}

# Plot each gene on a separate ax
total_genes = sum(len(genes) for genes in top_genes_dict.values())

# Also extract a flat list of genes
all_genes = [gene for genes in top_genes_dict.values() for gene in genes]

# Sort cells according to enrichment within TLS
mean_true, std_true, mean_false = prepare_plot(tls_df, level)

# Sort by highest mean_true
sorted_cells = mean_true.sort_values(ascending=False).index
mean_true = mean_true.loc[sorted_cells]
std_true = std_true.loc[sorted_cells]
mean_false = mean_false.loc[sorted_cells]

# Set adata as sorted categorical
all_adata.obs[level] = pd.Categorical(all_adata.obs[level], categories=sorted_cells, ordered=True)

fig, axs = plt.subplots(ncols=1, nrows=total_genes*2, figsize=(5, 5))

# Group the axs in groups of 2 (inside vs outside TLS)
axs = [axs[i:i + 2] for i in range(0, len(axs), 2)]

heatmap_axs = []

for i, gene in enumerate(all_genes):
    last_ax = False
    gene_axs = axs[i]
    # First plot outside TLS
    adata_sub = all_adata[all_adata.obs['in_TLS'] == False].copy()
    heatmap_ax_1, cbar_ax = plot_heatmap(adata_sub, gene_axs[0], gene, annotation_key=level, vmin=vmin, vmax=vmax, last_ax=last_ax)
    # Remove the xticks
    heatmap_ax_1.set_xticks([])
    # Then plot inside TLS
    # Check if it's the last gene
    if i == total_genes - 1:
        last_ax = True
    else:
        last_ax = False
    adata_sub = all_adata[all_adata.obs['in_TLS'] == True].copy()
    heatmap_ax_2, cbar_ax = plot_heatmap(adata_sub, gene_axs[1], gene, annotation_key=level, vmin=vmin, vmax=vmax, last_ax=last_ax)
    # Remove the xticks if not last_ax
    if not last_ax:
        heatmap_ax_2.set_xticks([])
    # Append heatmap axs as a pair
    heatmap_axs.append([heatmap_ax_1, heatmap_ax_2])

# Reposition the dotplot axes
reposition_dotplot_axs(heatmap_axs, y_offset=0.02, y_offset_within_group=0.01)

# Position the cbar
ax_pos = heatmap_axs[-1][-1].get_position()
space_per_dot = 0.04
cbar_height = space_per_dot * 0.8
cbar_width = space_per_dot * 5
cbar_ax.set_position([ax_pos.x1 + 0.03, ax_pos.y0, cbar_width*0.75, cbar_height*0.9])

# Replace mean by scaled in the cbar ax title
cbar_title = cbar_ax.get_title()
cbar_title = cbar_title.replace('Mean', 'Scaled')
cbar_ax.set_title(cbar_title, fontsize=8)

# Add an ax for the enrichment barplot on top
top_ax = heatmap_axs[0][0]
bar_ax = plt.axes([0.1, 0.1, 0.2, 0.2])
top_ax_pos = top_ax.get_position()
bar_ax.set_position([top_ax_pos.x0, top_ax_pos.y1 + 0.02, top_ax_pos.width, top_ax_pos.height/len(top_ax.get_yticks())*3])

# Extract colors from color_dict for the cell types in the order of sorted_cells
colors = [color_dict[cell_type] for cell_type in sorted_cells]
palette = color_dict_TLS[True]
bar_ax.bar(sorted_cells, mean_true, color=palette, yerr=std_true, ecolor='black',
        error_kw=dict(lw=0.8, capsize=3, capthick=1), width=0.9)

# Plot 'Not corresponding' (False_fraction) on top
bar_ax.bar(sorted_cells, mean_false, bottom=mean_true, label='Outside TLS', color='lightgray', width=0.9)

# Adjust the xlims
bar_ax.set_xlim(-0.5, len(sorted_cells)-0.5)
bar_ax.set_ylim(0, 1)

# Remove xticks
bar_ax.set_xticks([])
bar_ax.set_ylabel('Proportion', fontsize=8)

# Set only 2 yticks
bar_ax.set_yticks([0, 1])
bar_ax.set_yticklabels(['0', '1'], fontsize=8)

# Add a colorlegend for in TLS true/false
bar_ax.legend(handles=[
    mpatches.Patch(facecolor=color_dict_TLS[False], label='Outside TLS', edgecolor='black', linewidth=0.8),
    mpatches.Patch(facecolor=color_dict_TLS[True], label='Inside TLS', edgecolor='black', linewidth=0.8),
], loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8, frameon=False, handlelength=0.8)

# Add another ax on top with colors by major type
top_ax = bar_ax
top_ax_pos = top_ax.get_position()

color_list = [color_dict[get_major(cell, level=level)] for cell in sorted_cells]
color_ax = plt.axes([0.1, 0.1, 0.2, 0.02])
color_ax.set_position([top_ax_pos.x0, top_ax_pos.y1 + 0.01, top_ax_pos.width, 0.02])
color_ax.bar(sorted_cells, [1]*len(sorted_cells), color=color_list, width=0.9)

# Remove xticks & yticks
color_ax.set_xticks([])
color_ax.set_yticks([])

# Adjust the lims
color_ax.set_xlim(top_ax.get_xlim()[0], top_ax.get_xlim()[1])

# Add a legend
bottom_edge, top_edge, left_edge, right_edge = position_legend(cbar_ax)
x_legend = right_edge + 0.01
bottom_edge, top_edge, left_edge, right_edge = position_legend(color_ax)
y_legend = top_edge

add_color_legend(
    fig,
    labels=major_cell_types_sorted,
    color_dict=color_dict,
    title='Major class',
    x=x_legend,
    y=y_legend,
    offset=0,
    spacing_per_label=0.025,
    loc='upper left',
    fontsize=8,
)

# Set a title
color_ax.set_title(f'Gene expression inside vs outside TLS regions', fontweight='bold')

ax_width = 0.04*0.7
rect_axs = []

# Add small axes just to the right of each heatmap ax colored by whether in TLS or not
for i, gene_axs in enumerate(heatmap_axs):
    tmp_axs = []
    for j, ax in enumerate(gene_axs):
        in_tls = True if j == 1 else False
        side_ax = plt.axes([0.1, 0.1, 0.02, 0.02])
        ax_pos = ax.get_position()
        side_ax.set_position([ax_pos.x0 - ax_width - 0.01, ax_pos.y0, ax_width, ax_pos.height])
        side_ax.bar([0], [1], color=color_dict_TLS[in_tls], width=1)
        side_ax.set_xticks([])
        side_ax.set_yticks([])
        side_ax.set_xlim(-0.5, 0.5)
        side_ax.set_ylim(0, 1)
        tmp_axs.append(side_ax)
    rect_axs.append(tmp_axs)

# Draw a bracket connecting each pair of rect axs
for i, tmp_axs in enumerate(rect_axs):
    if len(tmp_axs) == 2:
        ax1, ax2 = tmp_axs
        fig = ax1.figure
        # Get left edge and vertical centers in figure coordinates
        pos1 = ax1.get_position()
        pos2 = ax2.get_position()
        x = pos1.x0 - 0.01  # bracket just to the left
        y1 = pos1.y0 + pos1.height / 2
        y2 = pos2.y0 + pos2.height / 2
        # Bracket shape: horizontal, vertical, horizontal
        x_bracket = [x + 0.01, x, x, x + 0.01]
        y_bracket = [y1, y1, y2, y2]
        plt.plot(x_bracket, y_bracket, color='black', linewidth=0.8, transform=fig.transFigure, clip_on=False)
        # Add a small spine to the middle
        plt.plot([x, x - 0.01], [(y1 + y2) / 2, (y1 + y2) / 2], color='black', linewidth=0.8, transform=fig.transFigure, clip_on=False)
        # Add gene name to the middle
        y_text = (y1 + y2) / 2
        gene = all_genes[i]
        fig.text(
            x - 0.02, y_text,
            gene,
            ha='right', va='center',
            fontsize=9,
            color='black',
        )

# Remove the yticklabels from the heatmap axs
for gene_axs in heatmap_axs:
    for ax in gene_axs:
        ax.set_yticklabels([])

# If subset level, replace xticklabels of the last heatmap ax by the short names
if 'subset' in level:
    last_heatmap_ax = heatmap_axs[-1][-1]
    xticklabels = last_heatmap_ax.get_xticklabels()
    new_labels = []
    for label in xticklabels:
        cell_type = label.get_text()
        if cell_type in miniatlas_sub_dict['miniatlas_subset_short']:
            new_label = miniatlas_sub_dict['miniatlas_subset_short'][cell_type]
        else:
            new_label = cell_type
        new_labels.append(new_label)
    last_heatmap_ax.set_xticklabels(new_labels)

# Position the title
top_ax = color_ax
title_x = (top_ax.get_position().x0 + top_ax.get_position().x1) / 2
title_y = top_ax.get_position().y1 + 0.1

plt.suptitle(f'All samples - Gene expression inside vs outside TLS regions (annotation_key={level}, scale=var)', x=title_x, y=title_y, fontsize=12, fontweight='bold')


# --- Save the output

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_5G.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')


# --- Finish

print("Script successfully completed")