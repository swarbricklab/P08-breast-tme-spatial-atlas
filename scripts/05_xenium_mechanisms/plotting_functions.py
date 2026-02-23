#!/usr/bin/env python3

"""
plotting_functions.py for 05_xenium_mechanisms

Utility functions for spatial plotting and coordinate transformations.
"""

import os
import numpy as np
import pandas as pd
from shapely.ops import transform as shapely_transform

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import seaborn as sns
from matplotlib.colors import TwoSlopeNorm
from matplotlib.colors import Normalize

import scipy.cluster.hierarchy as sch


# --- Functions

def transform_geojson_to_image_space(gdf, adata, sample_id, project_dir):
    """
    Transform GeoJSON geometries from Xenium micron space to image pixel space.
    Steps:
    1. Convert from microns to pixels
    2. Apply inverse transformation matrix
    3. Scale by hires scalefactor
    4. Shift by crop coordinates
    """
    hires_factor = adata.uns['spatial'][sample_id]['scalefactors']['tissue_hires_scalef']
    microns_per_pixel = 0.2125  # From experiment.xenium file
    # Load transformation matrix
    path_to_transform_file = os.path.join(project_dir, 'resources', 'image_alignment', 'xenium', 'alignment_matrix', f'{sample_id}_matrix.csv')
    transform_matrix = pd.read_csv(path_to_transform_file, header=None).to_numpy()
    inverse_transform = np.linalg.inv(transform_matrix)
    # Load crop coordinates
    crop_coords_file = os.path.join(project_dir, 'resources', 'image_alignment', 'xenium', 'crop_coords', f'{sample_id}_crop_coordinates.csv')
    crop_coords = pd.read_csv(crop_coords_file)
    x_start = crop_coords['x_start'].values[0]
    y_start = crop_coords['y_start'].values[0]
    def transform_func(x, y, z=None):
        coords = np.column_stack([x, y])
        coords_pixels = coords / microns_per_pixel
        hom_coords = np.hstack([coords_pixels, np.ones((coords_pixels.shape[0], 1))])
        transformed = (inverse_transform @ hom_coords.T).T[:, :2]
        scaled = transformed * hires_factor
        shifted = scaled - np.array([x_start * hires_factor, y_start * hires_factor])
        return shifted[:, 0], shifted[:, 1]
    gdf['geometry'] = gdf['geometry'].apply(lambda geom: shapely_transform(transform_func, geom))
    return gdf


# Function for beautifying heatmap
def beautify_heatmap(g, col_colors=True, row_colors=True, cbar_title='Pearson correlation',
                     fontsize_title=13, label_size=12, cbar_title_pad=10, cbar_pad=6, title=None, x_rot=90,
                     cbar_height_square=5, space_per_square=0.07, title_ax_n=None, cbar_width_square=0.1,
                     width_to_height_ratio=1, colors_div=1, offset=0.01):
    # Get the current figure
    fig = plt.gcf()
    # Extract the heatmap ax
    heatmap_ax = g.ax_heatmap
    # -- Alter the main heatmap axis
    heatmap_ax.set_ylabel('')
    heatmap_ax.set_xlabel('')
    heatmap_ax.tick_params(axis='x', labelsize=label_size, rotation=x_rot)  # Set font size for xticklabels
    heatmap_ax.tick_params(axis='y', labelsize=label_size)  # Set font size for yticklabels
    heatmap_ax.set_facecolor('darkgray') # Set facecolor to gray so masked values appear in gray
    for spine in heatmap_ax.spines.values():
        spine.set_visible(True)
    heatmap_ax_pos = heatmap_ax.get_position()
    # --- Alter the top dendrogram (put the title above of the dendrogram)
    top_dendrogram_ax = g.ax_col_dendrogram
    if title_ax_n is None:
        heatmap_ax.set_title(f'{title}', fontsize=fontsize_title, fontweight='bold')
    else:
        fig.axes[title_ax_n].set_title(f'{title}', fontsize=fontsize_title, fontweight='bold')
    left_dendrogram_ax = g.ax_row_dendrogram
    # -- Adjust colorbar columns
    if col_colors:
        col_colors_ax = g.ax_col_colors
        for spine in col_colors_ax.spines.values():
            spine.set_visible(True)
        # Move column colors & dendrogram up
        top_dendrogram_ax_pos = top_dendrogram_ax.get_position()
        col_colors_ax_pos = col_colors_ax.get_position()
        col_colors_ax.set_position((heatmap_ax_pos.x0, heatmap_ax_pos.y1 + offset, col_colors_ax_pos.width, col_colors_ax_pos.height))
        top_dendrogram_ax.set_position((heatmap_ax_pos.x0, heatmap_ax_pos.y1 + col_colors_ax_pos.height + offset, top_dendrogram_ax_pos.width, top_dendrogram_ax_pos.height))
     # -- Adjust row color
    if row_colors:
        row_colors_ax = g.ax_row_colors
        for spine in row_colors_ax.spines.values():
            spine.set_visible(True)
    # Adjust heatmap ax based on number of samples
    xlabels = heatmap_ax.get_xticklabels()
    ylabels = heatmap_ax.get_yticklabels()
    new_width = space_per_square * width_to_height_ratio * len(xlabels)
    new_height = space_per_square * len(ylabels)
    heatmap_ax.set_position((heatmap_ax_pos.x0, heatmap_ax_pos.y0, new_width, new_height))
    # Adjust position other axes
    new_heatmap_pos = heatmap_ax.get_position()
    if col_colors:
        ylabels_col = col_colors_ax.get_yticklabels()
        new_height = space_per_square * len(ylabels_col) / colors_div
        col_colors_ax.set_position((new_heatmap_pos.x0, new_heatmap_pos.y1 + offset, new_heatmap_pos.width, new_height))
        new_col_pos = col_colors_ax.get_position()
        top_dendrogram_ax.set_position((new_col_pos.x0, new_col_pos.y1, new_col_pos.width, top_dendrogram_ax_pos.height))
    if row_colors:
        xlabels_row = row_colors_ax.get_xticklabels()
        new_width = space_per_square * len(xlabels_row) / colors_div
        row_colors_ax_pos = row_colors_ax.get_position()
        row_colors_ax.set_position((new_heatmap_pos.x0 - new_width - offset, new_heatmap_pos.y0, new_width,
                                    new_heatmap_pos.height))
        row_colors_ax_pos = row_colors_ax.get_position()
        left_dendrogram_ax_pos = left_dendrogram_ax.get_position()
        left_dendrogram_ax.set_position((row_colors_ax_pos.x0 - left_dendrogram_ax_pos.width, new_heatmap_pos.y0, left_dendrogram_ax_pos.width,
                                        new_heatmap_pos.height))
    else:
        left_dendrogram_ax_pos = left_dendrogram_ax.get_position()
        left_dendrogram_ax.set_position((new_heatmap_pos.x0 - left_dendrogram_ax_pos.width, new_heatmap_pos.y0, left_dendrogram_ax_pos.width,
                                        new_heatmap_pos.height))
    # -- Adjust colorbar axis
    cbar_ax = g.ax_cbar
    cbar_ax.set_ylabel(cbar_title, labelpad=cbar_title_pad, fontsize=label_size)
    cbar_ax.yaxis.set_label_position("left")
    cbar_ax.tick_params(axis='y', labelsize=label_size)  # Set font size for yticklabels
    for spine in cbar_ax.spines.values():
        spine.set_visible(True)
    # Adjust position
    cbar_ax_pos = cbar_ax.get_position()
    heatmap_ax_pos = heatmap_ax.get_position()
    cbar_ax.set_position((heatmap_ax_pos.x1 + space_per_square*cbar_pad, heatmap_ax_pos.y0, space_per_square*cbar_width_square, space_per_square*cbar_height_square))
    return fig.axes

# Function to add legend
def add_legend(ax, ordered_list, color_lut, title='Lineage', x_offset=1.5, y_offset=1.5, loc='lower left', fontsize=10, edgecolor='black'):
    # Create a list of handles and labels for the legend
    handles_subtype = []
    labels_subtype = []
    for subtype in ordered_list:
        handles_subtype.append(mpatches.Patch(facecolor=color_lut[subtype], edgecolor='black', linewidth=0.8))
        labels_subtype.append(subtype)
    # Create the legend using the Patch objects
    legend1 = ax.legend(handles=handles_subtype, labels=labels_subtype, bbox_to_anchor=(1, 1), loc=loc, fontsize=fontsize, 
                                title=title, handlelength=0.7, frameon=False)
    # Legend position
    x_offset = x_offset
    y_offset = y_offset
    legend1.set_bbox_to_anchor((x_offset, y_offset))
    legend1.get_title().set_fontweight('bold')
    # Set the alignment of the legend title to the left
    legend1._legend_box.align = "left"
    return(legend1)

# Some functions for plotting
def get_sensible_size_max(df, size_column='proportion', increment=0.1):
    # Get the maximum value from the size column
    max_value = df[size_column].max()
    # Round up the maximum value to the nearest increment
    size_max = np.ceil(max_value / increment) * increment
    return size_max

# Dotplot function
def seaborn_dotplot(df_plot, dotplot_ax, x='category', y='cell_type', hue='total_count', size='proportion',
                    title=None, vmax=None, vmin=0, title_fontsize=8, midpoint=None, cmap='Reds', max_size_px=100,
                    edgecolor='#232323', linewidth=0.5, size_max=None, size_min=0):
    # # Specify the size of the figure (width, height) in inches
    # figsize = (width, height)
    # # Create a figure and axis with the specified size
    # Make a copy of the df to avoid overwriting it when doing df.clip
    df = df_plot.copy()
    if vmax == None:
        vmax = df[hue].max()
    else:
        df[hue] = df[hue].clip(upper=vmax)
    if vmin == None:
        vmin = df[hue].min()
    # Get sensible size_max
    if size_max is None:
        size_max = get_sensible_size_max(df_plot, size_column=size)
    if midpoint is not None:
        norm = TwoSlopeNorm(vmin=vmin, vcenter=midpoint, vmax=vmax)
    else:
        norm = Normalize(vmin=vmin, vmax=vmax)    # Plot the dotplot
    sns.scatterplot(data=df, x=x, y=y, hue=hue, size=size, hue_norm=norm, size_norm=(size_min, size_max),
                    sizes=(0, max_size_px), edgecolor=edgecolor, linewidth=linewidth, palette=cmap, ax=dotplot_ax,)
    # -- Beautify
    # Add title
    dotplot_ax.set_title(title, fontweight='bold', fontsize=title_fontsize)
    # Rotate the x-axis tick labels by 90 degrees
    dotplot_ax.set_xticklabels(dotplot_ax.get_xticklabels(), rotation=45, ha='right')
    dotplot_ax.set_xlabel('')
    # Remove y-label
    dotplot_ax.set_ylabel('')
    # Return the ax
    return dotplot_ax, size_max

# Separate function for adding the dotplot legend
def add_size_legend(dotplot_ax, target_ax, size='proportion', size_title='Proportion (per layer)', limit_circles=False, labelspacing=0.5,
                    handletextpad=0.5, alignment='center', bbox_to_anchor=(4, 0.5), loc='center left'):
    # -- Add size legend
    # Move the legend to the top of the plot & remove hue variable legend
    handles, labels = dotplot_ax.get_legend_handles_labels()
    # Find the index of 'proportion'
    index = labels.index(size)
    # Get everything after 'proportion'
    size_labels = labels[index + 1:]
    size_handles = handles[index + 1:]
    # Limit to 3 circles (smallest, middle, largest)
    if limit_circles:
        if len(size_handles) > 3:
            size_indices = [0, len(size_handles) // 2, -1]  # Pick first, middle, and last
            size_handles = [size_handles[i] for i in size_indices]
            size_labels = [size_labels[i] for i in size_indices]
    # Make the spots gray
    for handle in size_handles:
        handle.set_markerfacecolor('lightgray')
        handle.set_linewidth(0.8) 
    # Remove original legend
    dotplot_ax.legend().remove()
    # Make new size legend
    legend = target_ax.legend(size_handles, size_labels, loc=loc, labelspacing=labelspacing, bbox_to_anchor=bbox_to_anchor, 
                              ncol=1, frameon=False, title=f'{size_title}', handletextpad=handletextpad, alignment=alignment)
    for txt in legend.get_texts():
        txt.set_ha("left")
    return legend

# Function to move plots so they don't overlap, moving both in x & y direction
def move_plots_x_y(axs, x_offset=0.02, y_offset=0.04, remove_ylabels=True, remove_xlabels=True, remove_title=True):
    for i, axs_row in enumerate(axs):
        for j, ax in enumerate(axs_row):
            if j > 0:
                if remove_ylabels:
                    ax.set_ylabel('')
                    ax.set_yticklabels('')
            if i != (len(axs) - 1):
                if remove_xlabels:
                    ax.set_xlabel('')
                    ax.set_xticklabels('')
            if i != 0:
                if remove_title:
                    ax.set_title('')
            # Adjust the positions of the subplots
            ax_pos = ax.get_position()
            # Shift the subplot to the desired location + set the width so it's equal for all plots
            if i == 0:
                new_pos_y = ax_pos.y0
            else:
                new_pos_y = ax_pos_row_prev.y0 - ax_pos.height - y_offset
            if j == 0:
                new_pos_x = ax_pos.x0
            else:
                new_pos_x = ax_pos_prev.x0 + ax_pos_prev.width + x_offset
            ax.set_position([new_pos_x, new_pos_y, ax_pos.width, ax_pos.height])
            ax_pos_prev = ax.get_position()
        ax_pos_row_prev = axs_row[0].get_position()

def cut_dendrogram_and_draw_lines(g, axis='row', threshold_dend=None, color='black', linestyle='--',
                                  alpha=0.7, linewidth=1, ax_main=None, line_adjust=0):
    """
    Cuts the dendrogram at a given threshold and draws separation lines at cluster boundaries.
    Parameters:
    g : seaborn.clustermap
        The clustered heatmap object from seaborn.
    axis : str, optional
        Whether to apply to the 'row' or 'col' dendrogram.
    threshold_dend : float, optional
        Custom threshold for cutting the dendrogram. If None, defaults are used.
    """
    # If no ax defined, draw lines on the heatmap ax
    if ax_main is None:
        ax_main = g.ax_heatmap
    # Select the appropriate dendrogram and axes
    if axis == 'row':
        linkage_matrix = g.dendrogram_row.linkage
        reordered_ind = g.dendrogram_row.reordered_ind
        ax_colors = g.ax_row_colors
    elif axis == 'col':
        linkage_matrix = g.dendrogram_col.linkage
        reordered_ind = g.dendrogram_col.reordered_ind
        ax_colors = g.ax_col_colors
    else:
        raise ValueError("axis must be either 'row' or 'col'.")
    # Set default threshold if not provided
    if threshold_dend is None:
        threshold_dend = 2
    # Compute clusters
    clusters = sch.fcluster(linkage_matrix, threshold_dend, criterion='distance')
    # Identify positions where clusters change
    sorted_clusters = [clusters[i] for i in reordered_ind]
    change_indices = [i for i in range(1, len(sorted_clusters)) if sorted_clusters[i] != sorted_clusters[i - 1]]
    # Draw separation lines
    for idx in change_indices:
        if axis == 'row':
            ax_main.axhline(idx+line_adjust, color=color, linestyle=linestyle, alpha=alpha, linewidth=linewidth)
            ax_colors.axhline(idx, color=color, linestyle=linestyle, alpha=alpha, linewidth=linewidth)
        else:  # axis == 'col'
            ax_main.axvline(idx+line_adjust, color=color, linestyle=linestyle, alpha=alpha, linewidth=linewidth)
            ax_colors.axvline(idx, color=color, linestyle=linestyle, alpha=alpha, linewidth=linewidth)

# Function do adjust the dendrograms
def adjust_dendrograms(g, space_per_square=0.05, num_square=4, title=''):
    heatmap_ax = g.ax_heatmap
    heatmap_pos = heatmap_ax.get_position()
    # Adjust row dendrogram if it exists
    if hasattr(g, "ax_row_dendrogram") and g.ax_row_dendrogram is not None:
        if hasattr(g, "ax_row_colors") and g.ax_row_colors is not None:
            row_colors_ax_pos = g.ax_row_colors.get_position()
        else:
            row_colors_ax_pos = heatmap_pos  # fallback if no row_colors
        
        left_dendrogram_ax = g.ax_row_dendrogram
        new_width = space_per_square * num_square
        left_dendrogram_ax.set_position((
            row_colors_ax_pos.x0 - new_width,
            row_colors_ax_pos.y0,
            new_width,
            row_colors_ax_pos.height
        ))
    # Adjust column dendrogram if it exists
    if hasattr(g, "ax_col_dendrogram") and g.ax_col_dendrogram is not None:
        if hasattr(g, "ax_col_colors") and g.ax_col_colors is not None:
            col_colors_ax_pos = g.ax_col_colors.get_position()
        else:
            col_colors_ax_pos = heatmap_pos  # fallback if no col_colors
        top_dendrogram_ax = g.ax_col_dendrogram
        new_height = space_per_square * num_square
        top_dendrogram_ax.set_position((
            col_colors_ax_pos.x0,
            col_colors_ax_pos.y1,
            col_colors_ax_pos.width,
            new_height
        ))
        # Set title on top dendrogram
        top_dendrogram_ax.set_title(title, fontweight='bold')

# Function to reformat numbers
def reformat_number(number, decimals=1):
    if number >= 1000000:
        # Divide the number by 1,000,000 and format it to one decimal place
        formatted_number = "{:.1f}m".format(number / 1000000)
    elif number >= 1000:
        # Divide the number by 1000 and format it to one decimal place
        formatted_number = "{:.1f}k".format(number / 1000)
    else:
        formatted_number = str(round(number, decimals))
    return formatted_number

# Define a function for annotation the bars
def annotate_bars(ax, df, x, y, plot_order, annotation_size=8, margin=1.05, round_int=True, offset=None):
    value_dict = df.set_index(x)[y].to_dict()
    for n, item in enumerate(plot_order):
        value = value_dict[item]
        if round_int:
            text = reformat_number(int(value))
        else:
            text = reformat_number(float(value), decimals=2)
        if offset is None:
            shift = 0
        else:
            shift = offset[n]
        ax.annotate(text, xy=(n, value+shift), ha='center', va='bottom', fontsize=annotation_size)
    # Adjust ymax to accommodate highest value
    ymax = ax.get_ylim()[1]
    ax.set_ylim(0, ymax * margin)

def beautify_plot(ax, alpha=0.7):
    # Beautify the plot
    for bar in ax.patches: # Adjust the colors of the boxplot
        # Get the RGBA values of the fill color
        fill_color = bar.get_facecolor()
        fill_color = list(fill_color)
        darkened_color = fill_color.copy()
        for i in range(len(darkened_color)-1):
            darkened_color[i] *= 0.8
        fill_color[-1] = alpha # make the fill color a bit transparent
        bar.set_color(tuple(fill_color))
        bar.set_edgecolor(tuple(darkened_color))

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
    bottom_edges = []
    left_edges = []
    for ax in axs:
        bbox = ax.get_tightbbox(fig.canvas.get_renderer())
        inv_trans = fig.transFigure.inverted()
        bbox_fig = inv_trans.transform(bbox)
        right_edges.append(bbox_fig[1][0])
        top_edges.append(bbox_fig[1][1])  # y1 is the top in figure coordinates
        bottom_edges.append(bbox_fig[0][1])  # y0 is the bottom in figure coordinates
        left_edges.append(bbox_fig[0][0])  # x0 is the left in figure coordinates
    max_right_edge = max(right_edges)
    max_top_edge = max(top_edges)
    min_bottom_edge = min(bottom_edges)
    min_left_edge = min(left_edges)
    return max_right_edge, max_top_edge, min_bottom_edge, min_left_edge

# Position the cbar so it doesn't overlap with any labels
def position_cbar(target_ax, cbar_ax, offset=0.1):
    # Draw the canvas to compute tight bbox
    fig = plt.gcf()
    fig.canvas.draw()
    # Get position of ax1 including labels (in display coords)
    bbox = target_ax.get_tightbbox(fig.canvas.get_renderer())
    # Transform display coords to figure coords
    inv_trans = fig.transFigure.inverted()
    bbox_fig = inv_trans.transform(bbox)
    # Get right edge of ax1's bbox (includes yticklabels)
    x0, y0 = bbox_fig[0]
    x1, y1 = bbox_fig[1]
    right_edge = x1
    cbar_ax_pos = cbar_ax.get_position()
    cbar_ax.set_position((right_edge + offset, cbar_ax_pos.y0, cbar_ax_pos.width, cbar_ax_pos.height))

# Function to add legend
def add_color_legend(fig, labels, color_dict, title='', x=0, y=0, offset=0.01, spacing_per_label=0.05, loc='center left'):
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
        borderpad=0,
        bbox_transform=fig.transFigure,  # <-- This is the key fix!
    )
    legend.get_title().set_fontweight('bold')
    legend._legend_box.align = "left"
    # Return updated vertical position for next stacked legend
    return legend