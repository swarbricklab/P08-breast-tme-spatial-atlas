#!/usr/bin/env python3

# --- Aim
# Plot members of the lymphotoxin system in TLS HEVs vs stromal HEVs (Figure 4I)
# Shows expression of lymphotoxin signaling components in perivascular bins

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import os
import re

import seaborn as sns

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments
parser = argparse.ArgumentParser(description="Plot lymphotoxin system expression for Figure 4I",
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
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_visium_hd.csv')

sample_sheet = pd.read_csv(sample_sheet_path)
sample_lut = sample_sheet.set_index('section_id')


# --- Parameters

bin_size = 25  # um
resolution = 'square_016um'
use_data = 'log_norm'


# --- Load the adata objects

adata_dict = {}

for sample in sample_list:
    adata_dict[sample] = ad.read_h5ad(os.path.join(h5ad_dir, 'visium_hd', resolution, f'{sample}.h5ad'))
    # Add the sample_id to the obs
    adata_dict[sample].obs['sample_id'] = sample


# --- Load the spot labels from LR analysis and add to adata

lr_label_dir = os.path.join(res_dir, '04_vascular_microenvironment', '04_ligand_receptor_analysis',
                             '03_plot_LR_per_sample', f'{bin_size}um_bin', 'csv_files')

most_recent_label_path = lr_label_dir

for sample in sample_list:
    spot_label_file = f'{sample}_HEV_spot_grouping.csv'
    label_df = pd.read_csv(os.path.join(most_recent_label_path, spot_label_file), index_col=0)
    adata_dict[sample].obs[label_df.columns] = label_df


# --- Extract a single adata object for plotting

adata_list = list(adata_dict.values())
all_adata = ad.concat(adata_list, merge='same')
all_adata.obs_names_make_unique()


# --- Define a color_dict for the spot labels

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
color_file = 'TLS_HEV_colors.csv'

# Load the color file
color_df = pd.read_csv(os.path.join(color_dir, color_file), index_col=0)

# Make a color_dict
color_dict = {}
for label, color in zip(color_df.index, color_df['color']):
    color_dict[label] = color

color_dict['TLS'] = color_dict['nearby_near_TLS_HEV']
color_dict['Stroma'] = color_dict['nearby_non_TLS_HEV']


# --- Also prepare a dict to rename the labels in the plots

rename_dict = {
    'near_TLS_HEV': 'HEV (TLS)',
    'non_TLS_HEV': 'HEV-like (stroma)',
    'nearby_spots': f'Perivascular bins',
    'nearby_near_TLS_HEV': f'Perivascular bins (TLS)',
    'nearby_non_TLS_HEV': f'Perivascular bins (stroma)',
}


# --- Functions for plotting

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

def join_lr_multi(axs, lr_pairs, offset=0.3, offset_start=0.2):
    """
    Draws brackets connecting ligands and receptors, extending the bracket
    further left if a ligand has multiple receptors or a receptor has
    multiple ligands.

    If a vertical segment would overlap with a previously drawn one, it is
    iteratively moved to the left until no overlap is found.

    Args:
        axs (list | Axes): The matplotlib Axes object(s).
        lr_pairs (list[tuple]): A list of (ligand, receptor) string pairs.
        offset (float): Horizontal distance from the label's left edge to the end of the bracket.
        offset_start (float): Horizontal distance from the longest label to the main vertical line.
    """
    if not isinstance(axs, (list, tuple, np.ndarray)):
        axs = [axs]
    # 1. Build lookup tables for label positions and properties.
    label_to_ax, label_to_idx, label_to_bbox, ax_inv = {}, {}, {}, {}
    all_labels = set(lr for pair in lr_pairs for lr in pair)
    renderer = axs[0].get_figure().canvas.get_renderer()
    for ax in axs:
        inv = ax.transData.inverted()
        ax_inv[ax] = inv
        for idx, label in enumerate(ax.get_yticklabels()):
            text = label.get_text()
            if text in all_labels:
                label_to_ax[text] = ax
                label_to_idx[text] = idx
                label_to_bbox[text] = label.get_window_extent(renderer)
    # 2. Count ligand/receptor occurrences to identify multi-connections.
    ligand_counts = {p[0]: 0 for p in lr_pairs}
    receptor_counts = {p[1]: 0 for p in lr_pairs}
    for lr1, lr2 in lr_pairs:
        ligand_counts[lr1] += 1
        receptor_counts[lr2] += 1
    # 3. Determine the main vertical alignment line based on the LEFT-MOST label edge.
    all_x0_coords = [
        ax_inv[label_to_ax[name]].transform((label_to_bbox[name].x0, 0))[0] 
        for name in all_labels
    ]
    min_x0 = min(all_x0_coords)
    x_left_align = min_x0 - offset_start - offset
    # 4. Prepare for drawing.
    ligand_seen = {lig: 0 for lig in ligand_counts}
    receptor_seen = {rec: 0 for rec in receptor_counts}
    # NEW: Track the position and span of drawn vertical lines to check for overlaps.
    drawn_verticals = []
    indent_step = offset * 2 # How far to move left to avoid an overlap.
    # 5. Main drawing loop.
    for lr1, lr2 in lr_pairs:
        # Get positional data for the current ligand-receptor pair
        ax_lig, ax_rec = label_to_ax[lr1], label_to_ax[lr2]
        idx_lig, idx_rec = label_to_idx[lr1], label_to_idx[lr2]
        bbox_lig, bbox_rec = label_to_bbox[lr1], label_to_bbox[lr2]
        inv_lig, inv_rec = ax_inv[ax_lig], ax_inv[ax_rec]
        # Calculate the end-point for each horizontal line, next to its label
        x_lig_end = inv_lig.transform((bbox_lig.x0, 0))[0] - offset
        x_rec_end = inv_rec.transform((bbox_rec.x0, 0))[0] - offset
        # Get counts for determining vertical spread
        n_total_lig, n_seen_lig = ligand_counts[lr1], ligand_seen[lr1]
        n_total_rec, n_seen_rec = receptor_counts[lr2], receptor_seen[lr2]
        # Calculate vertical positions for the horizontal lines
        y_spread_lig, y_spread_rec = 0.3, 0.3
        y_lig_data = idx_lig + 0.5 + ((n_seen_lig - (n_total_lig - 1) / 2) * y_spread_lig if n_total_lig > 1 else 0)
        y_rec_data = idx_rec + 0.5
        if n_total_rec > 1:
            order = (n_total_rec - 1) - n_seen_rec
            y_rec_data += (order - (n_total_rec - 1) / 2) * y_spread_rec
        # NEW: Iteratively find a non-overlapping x-position for the vertical segment.
        y_span = tuple(sorted((y_lig_data, y_rec_data)))
        # Determine the initial x-position to try. Multi-connections start further left.
        current_x = x_left_align
        if n_total_lig > 1 or n_total_rec > 1:
            current_x -= indent_step
        while True:
            # Check for any overlap with existing lines at the current x-position
            has_overlap = any(
                max(y_span[0], s['span'][0]) < min(y_span[1], s['span'][1])
                for s in drawn_verticals if s['x'] == current_x
            )
            if not has_overlap:
                x_left_draw = current_x  # This position is clear, lock it in.
                break
            else:
                current_x -= indent_step # Overlap found, move left and re-check.
        # Draw the bracket
        x_points = [x_lig_end, x_left_draw, x_left_draw, x_rec_end]
        y_points = [y_lig_data, y_lig_data, y_rec_data, y_rec_data]
        ax_lig.plot(x_points, y_points, color='black', linewidth=2, zorder=10, clip_on=False, solid_joinstyle='miter')
        # NEW: Record the segment that was just drawn.
        drawn_verticals.append({'x': x_left_draw, 'span': y_span})
        # Update seen counters
        ligand_seen[lr1] += 1
        receptor_seen[lr2] += 1
    # UPDATED: Calculate the furthest point reached for potential axis adjustments
    if drawn_verticals:
        furthest_x = min(d['x'] for d in drawn_verticals)
    else:
        furthest_x = x_left_align 
    return furthest_x

# Position legend function
def position_legend(axs):
    """
    Given an array (or list) of axes, find the rightmost right_edge and the highest top_edge
    in figure coordinates, for robust legend placement.
    Returns (max_right_edge, max_top_edge)
    """
    # Accept a single axis as well
    if not isinstance(axs, (list, np.ndarray)):
        axs = [axs]
    fig = axs[0].figure
    fig.canvas.draw()
    right_edges = []
    top_edges = []
    left_edges = []
    bottom_edges = []
    for ax in axs:
        bbox = ax.get_tightbbox(fig.canvas.get_renderer())
        inv_trans = fig.transFigure.inverted()
        bbox_fig = inv_trans.transform(bbox)
        # bbox_fig[0] = (x0, y0) = bottom left; bbox_fig[1] = (x1, y1) = top right
        left_edges.append(bbox_fig[0][0])
        bottom_edges.append(bbox_fig[0][1])
        right_edges.append(bbox_fig[1][0])
        top_edges.append(bbox_fig[1][1])  # y1 is the top in figure coordinates
    max_right_edge = max(right_edges)
    max_top_edge = max(top_edges)
    max_left_edge = max(left_edges)
    max_bottom_edge = max(bottom_edges)
    return max_right_edge, max_top_edge, max_left_edge, max_bottom_edge

def style_dotplot_legend(size_ax, keep_indices):
    circles_collection = size_ax.collections[0]
    original_xticklabels = [label.get_text() for label in size_ax.get_xticklabels()]
    # Get the original data
    original_offsets = circles_collection.get_offsets()
    original_sizes = circles_collection.get_sizes()
    # Filter the data to keep only the desired indices
    new_offsets = original_offsets[keep_indices]
    new_sizes = original_sizes[keep_indices]
    # Make sure the distances between circles are equal
    new_offsets[:, 0] = np.linspace(new_offsets[0, 0], new_offsets[-1, 0], len(keep_indices))
    # Update the collection with the filtered data
    circles_collection.set_offsets(new_offsets)
    circles_collection.set_sizes(new_sizes)
    # Remove xticklabels for the ones we removed
    size_ax.set_xticks(new_offsets[:, 0])  # Set ticks to match the new offsets
    # Reset the labels
    new_labels = [original_xticklabels[i] for i in keep_indices]
    size_ax.set_xticklabels(new_labels)
    # Adjust xlim
    size_ax.get_xlim()
    size_ax.set_xlim(keep_indices[0] - 0.8, keep_indices[-1] + 1.8)


# --- First make a plot of all together to get vmax & size max

genes_dict = {
    'Lymphotoxin': ['TNFSF14', 'LTB', 'LTBR', 'TNFRSF14', 'LTA', 'TNFRSF1A', 'TNFRSF1B'],
}


lr_pairs = [('LTB', 'LTBR'), ('TNFSF14', 'LTBR'), ('TNFSF14', 'TNFRSF14'), ('LTA', 'TNFRSF1A'), ('LTA', 'TNFRSF1B')]


top_genes = genes_dict['Lymphotoxin']
vmin = 0

sc.set_figure_params(fontsize=10)

# Disable grid
plt.rcParams['axes.grid'] = False

fig, ax = plt.subplots(figsize=(5, 5))

dp = sc.pl.dotplot(all_adata, var_names=top_genes, groupby=['nearby_spots_revisited'], layer=use_data, swap_axes=True,
                   ax=ax, show=False)

vmax = dp['color_legend_ax'].get_xlim()[1]
size_max = int(dp['size_legend_ax'].get_xticklabels()[-1].get_text())/100

plt.close(fig)


# --- Plot 2: Dotplot of genes corresponding to top & bottom LRs

# Split the categories into near_TLS and non_TLS for better visualization
near_cats = [cat for cat in all_adata.obs['nearby_spots_revisited'].cat.categories if 'near_TLS' in cat]
non_cats = [cat for cat in all_adata.obs['nearby_spots_revisited'].cat.categories if 'non_TLS' in cat]

all_cats = [near_cats] + [non_cats]

ncols = 2
nrows = 1

fig, axs = plt.subplots(figsize=(5, 5), ncols=ncols, nrows=nrows)

heatmap_axs = []

for n, cat_type in enumerate(all_cats):
    ax = axs[n]
    adata_sub = all_adata[all_adata.obs['nearby_spots_revisited'].isin(cat_type)]
    dp = sc.pl.dotplot(adata_sub, var_names=top_genes, groupby=['nearby_spots_revisited'], layer=use_data, swap_axes=True,
                    ax=ax, show=False, vmin=vmin, vmax=vmax, dot_max=size_max, dot_min=None, cmap='Reds')
    # Adjust the size of the axes
    auto_resize_dotplot(dp['mainplot_ax'], space_per_dot=0.04)
    main_ax = dp['mainplot_ax']
    # Adjust the ylims
    buffer = 0.25
    main_ax.set_ylim(len(main_ax.get_yticks()) + buffer, -buffer)
    # If not the last ax, remove the colorlegend & sizelegend
    if n < ncols - 1:
        dp['color_legend_ax'].remove()
        dp['size_legend_ax'].remove()
    else:
        size_ax = dp['size_legend_ax']
        cbar_ax = dp['color_legend_ax']
    heatmap_axs.append(main_ax)

# Move the plots
for n, heatmap_ax in enumerate(heatmap_axs):
    pos = heatmap_ax.get_position()
    if n == 0:
        ax_pos_prev = pos
    else:
        heatmap_ax.set_position([ax_pos_prev.x1 + 0.01, pos.y0, pos.width, pos.height])
        ax_pos_prev = heatmap_ax.get_position()
        # Remove yticks
        heatmap_ax.set_yticks([])
        heatmap_ax.set_ylabel('')

size_ax.set_title('Fraction of bins\nin group (%)', y=0.5, fontsize=8.5)

# Position the cbar
main_ax = heatmap_axs[-1]
ax_pos = main_ax.get_position()
space_per_dot = 0.04
cbar_height = space_per_dot * 0.6
cbar_width = space_per_dot * 3
cbar_ax.set_position([ax_pos.x1 + 0.05, ax_pos.y0, cbar_width, cbar_height])

# Position the size legend just above the cbar
cbar_pos = cbar_ax.get_position()
size_ax.set_position([cbar_pos.x0, cbar_pos.y1 + 0.12, cbar_width, cbar_height*3])

# Join the ligand-receptor pairs
main_ax = heatmap_axs[0]


x_left_global = join_lr_multi([main_ax], lr_pairs, offset=0.25)

# Add an ax just above for counts per group
count_axs = []

for n, main_ax in enumerate(heatmap_axs):
    ax_pos = main_ax.get_position()
    count_ax = fig.add_axes([ax_pos.x0, ax_pos.y1 + 0.01, ax_pos.width, space_per_dot*2])  # [left, bottom, width, height]
    xticklabels = [label.get_text() for label in main_ax.get_xticklabels()]
    count_df = all_adata.obs['nearby_spots_revisited'].value_counts().reindex(xticklabels, fill_value=0)
    sns.barplot(x=count_df.index, y=count_df.values, ax=count_ax, palette=color_dict, edgecolor='black', linewidth=0.5, saturation=1)
    # Set count_ax xlim
    count_ax.set_xlim(-0.5-buffer, len(xticklabels)-0.5+buffer)
    # Remove x ticks & labels
    count_ax.set_xticks([])
    count_ax.set_xlabel('')
    count_ax.set_yticks([])
    if n == 0:
        count_ax.set_ylabel('# bins', fontsize=8)
    # Remove all spines except bottom
    for spine in ['top', 'right', 'left']:
        count_ax.spines[spine].set_visible(False)
    # Annotate the bars with cell counts
    for rect, count in zip(count_ax.patches, count_df[xticklabels]):
        height = rect.get_height()
        if height > 0:
            # Rework the number format if >1000
            if count >= 1000:
                # Round to 1 decimal
                count = f'{count/1000:.1f}k'
            count_ax.text(
                rect.get_x() + rect.get_width() / 2, height + count_ax.get_ylim()[1]*0.02,
                f'{count}',
                ha='center', va='bottom',
                fontsize=6
            )
    count_axs.append(count_ax)

# Set same ylims for all count_axs
ymax = max(count_ax.get_ylim()[1] for count_ax in count_axs)
for count_ax in count_axs:
    count_ax.set_ylim(0, ymax)

# Adjust the xticklabels according to the rename_dict
for main_ax in heatmap_axs:
    new_labels = []
    for label in main_ax.get_xticklabels():
        original_text = label.get_text()
        if original_text in rename_dict:
            new_label = rename_dict[original_text]
            # Remove the bits in brackets to shorten the label
            new_label = re.sub(r'\s*\(.*?\)\s*', '', new_label)
            new_labels.append(new_label)
        else:
            new_labels.append(original_text)
    main_ax.set_xticklabels(new_labels)

# Plot rectangles just above count axes
labels = ['TLS', 'Stroma']

# Draw rectangles from top to bottom of each ax
for n, ax in enumerate(count_axs):
    # Add vertical text
    label = labels[n]
    face_color = color_dict[label]
    # Convert x_left_global from data to figure coordinates
    pos = ax.get_position()
    x0 = pos.x0
    x1 = pos.x1
    y1 = pos.y1 + 0.02
    width = abs(x1 - x0)
    height = 0.01  # Height in figure coordinates
    highlight_rect = mpl.patches.Rectangle(
        (x0, y1), width, height,
        transform=fig.transFigure,
        facecolor=face_color,
        edgecolor='black',
        linewidth=0,
        alpha=1
    )
    fig.patches.append(highlight_rect)
    # Calculate midpoint x
    x_text = (x0 + x1) / 2
    y_text = y1 + height
    # Add text
    fig.text(
        x_text, y_text,
        label,
        ha='center', va='bottom',
        color='black',
        fontsize=8.5,
    )

# Limit size_ax to 3 circles
keep_indices = [0, len(size_ax.collections[0].get_sizes()) // 2, len(size_ax.collections[0].get_sizes()) - 1]
style_dotplot_legend(size_ax, keep_indices)

# Position title
title_x = (heatmap_axs[0].get_position().x0 + heatmap_axs[-1].get_position().x1) / 2
title_y = count_axs[0].get_position().y1 + 0.1

plt.suptitle(f'All samples - Lymphotoxin system members in TLS HEVs vs stromal HEVs (bin_size={bin_size}µm)', fontsize=10, x=title_x, y=title_y, fontweight='bold')


# --- Save the figure

out_file = os.path.join(res_dir, 'figures', 'Fig_4I.pdf')
os.makedirs(os.path.dirname(out_file), exist_ok=True)
plt.savefig(out_file, bbox_inches='tight', dpi=300)

print(f"Figure saved to: {out_file}")
print("Script successfully completed")

