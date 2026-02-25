# Description: Some functions used by the scripts in this directory

# --- Load packages

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import seaborn as sns
import numpy as np
import pandas as pd


import scanpy.plotting.palettes as sc_palettes

import sys
import os
sys.path.insert(0, os.path.dirname(__file__))
from helper_functions import generate_color_variations


# --- Functions

def create_plotting_colors(category_dict, ordered_categories, available_layers, mixed=True, other=True, major_categories=False, category_combinations=None):
    # Get some colors
    if len(category_dict) < (21-2):
        colors = sc_palettes.default_20
        # Remove red (already using this for the mixed category)
        colors = colors[0:3] + colors[4:]
    else:
        colors = sc_palettes.godsnot_102[:len(category_dict.keys())]
    top_cat_colors = {}
    # If it's major categories: keep the colors consistent (epithelial blue, immune orange etc.)
    if major_categories:
        # Get the keys of the annotations of interest dict - this is the order to keep
        key_list = available_layers
        for n, cat in enumerate(category_dict):
            i = key_list.index(cat)
            top_cat_colors[cat] = colors[i]
    else:
        for n, cat in enumerate(category_dict):
            top_cat_colors[cat] = colors[n]
    # Get some color variations of the original colors (so variations of the same annotation are a similar color)
    color_variations = {}
    for n, (key, values) in enumerate(category_dict.items()):
        if major_categories:
            i = key_list.index(key)
            color = colors[i]
        else:
            color = colors[n]
        num_cols = len(values)
        if num_cols > 4:
            step_size=0.1
        else:
            step_size=0.2
        variations = generate_color_variations(color, num_cols, step_size=step_size)
        color_variations[key] = variations
    # Create a new dictionary mapping values from dict1 to values from dict2
    color_dict = {val: color_variations[key][idx] for key, values in category_dict.items() for idx, val in enumerate(values)}
    if mixed:
        # Make anything that's mixed brown
        num_cols = len(category_combinations)
        if num_cols > 4:
            step_size=0.1
        else:
            step_size=0.2
        variations = generate_color_variations('#800000', num_cols, step_size=step_size)
        for n, cat in enumerate(category_combinations):
            color_dict[cat] = variations[n]
    if other:
        # Add 'no_annotation' as gray
        color_dict['no_annotation'] = 'lightgray'
    # Extract values from the dictionary based on the keys
    plotting_colors = [color_dict[key] for key in ordered_categories]
    return color_dict, plotting_colors, top_cat_colors

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


# Define a function to add margins to the dotplot
def add_margin(ax, x=0.05,y=0.05):
    # Add a specific margin to the plot (so no dots will be cut off)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.set_xlim(xlim[0]-x,xlim[1]+x)
    ax.set_ylim(ylim[0]+y,ylim[1]-y)

# Function to prepare the axs for the dotplot
def prepare_dotplot(ordered_major_cat_dict, ncols=None, nrows=1, width=3, height=1.2, sharey=False):
    if ncols is None:
        ncols = len(ordered_major_cat_dict)
    # Calculate the width of each subtype subplot based on how many samples per subtype
    grid_spec = []
    for layer in ordered_major_cat_dict:
        grid_spec.append(len(ordered_major_cat_dict[layer]))
    fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(width, height), sharey=sharey, gridspec_kw={'width_ratios': grid_spec, 'wspace': 0.1})
    # If ncols = 1: need to reshape to make axs iterable
    if ncols == 1:
        axs = np.reshape(axs, 1)
    return axs, grid_spec

# Dotplot function
def seaborn_dotplot(df_plot, dotplot_ax, x='category', y='cell_type', hue='total_count', size='proportion', ymargin=-0.8, xmargin=0.6, cbar_height=0.8,
                    cbar_width=2, space_per_dot=0.2, title=None, vmax=None, cbar_title='Number of cells', legend=True, title_fontsize=8,
                    size_title='Proportion (per layer)', size_max=1):
    # # Specify the size of the figure (width, height) in inches
    # figsize = (width, height)
    # # Create a figure and axis with the specified size
    # Make a copy of the df to avoid overwriting it when doing df.clip
    df = df_plot.copy()
    if vmax:
        df[hue] = df[hue].clip(upper=vmax)
    # Plot the dotplot
    sns.scatterplot(data=df, x=x, y=y, hue=hue, size=size, hue_norm=(0, vmax), size_norm=(0, size_max),
                    sizes=(0, 220), edgecolor='black', linewidth=0.8, palette='Reds', ax=dotplot_ax,)
    # -- Adjust the dimensions to the dimensions of the df
    dim_x = len(df[x].unique())
    dim_y = len(df[y].unique())
    ax_pos = dotplot_ax.get_position()
    # Shift the subplot to the desired location + set the width so it's equal for all plots
    dotplot_ax.set_position([ax_pos.x0, ax_pos.y0, space_per_dot*dim_x/2, space_per_dot*dim_y])
    # Add a margin to the dotplot
    # add_margin(dotplot_ax, x=xmargin, y=ymargin)
    # Set the margins
    num_cats_y = len(df[y].unique())
    num_cats_x = len(df[x].unique())
    dotplot_ax.set_xlim(-0.75, num_cats_x-1+0.75)
    dotplot_ax.set_ylim(num_cats_y-1+0.75, -0.75)
    if legend:
        _, cax = add_dotplot_legend(dotplot_ax, df, hue=hue, size=size, vmax=vmax, space_per_dot=space_per_dot, cbar_title=cbar_title,
                           cbar_height=cbar_height, cbar_width=cbar_width, size_title=size_title)
    # -- Beautify
    # Add title
    dotplot_ax.set_title(title, fontweight='bold', fontsize=title_fontsize)
    # Rotate the x-axis tick labels by 90 degrees
    dotplot_ax.set_xticklabels(dotplot_ax.get_xticklabels(), rotation=45, ha='right')
    dotplot_ax.set_xlabel('')
    # Remove y-label
    dotplot_ax.set_ylabel('')
    # Return the ax
    if legend:
        return dotplot_ax, cax
    else:
        return dotplot_ax

# Function to move plots so they don't overlap (only move in x direction)
def move_plots(axs, offset=0.02, remove_yticks=True):
    for n, ax in enumerate(axs):
        if n > 0:
            ax.set_ylabel('')
            ax.set_yticklabels('')
            if remove_yticks:
                ax.set_yticks([])
        # Adjust the positions of the subplots
        ax_pos = axs[n].get_position()
        # Shift the subplot to the desired location + set the width so it's equal for all plots
        if n == 0:
            axs[n].set_position([ax_pos.x0, ax_pos.y0, ax_pos.width, ax_pos.height])
        else:
            new_pos = ax_pos_prev.x0 + ax_pos_prev.width + offset
            axs[n].set_position([new_pos, ax_pos.y0, ax_pos.width, ax_pos.height])
        ax_pos_prev = axs[n].get_position()

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


def move_plots_y_only(axs, y_offset=0.04):
    for n, ax in enumerate(axs):
        if n != (len(axs) - 1):
            ax.set_xlabel('')
            ax.set_xticklabels('')
            ax.set_xticks([])
        if n != 0:
            ax.set_title('')
        # Adjust the positions of the subplots
        ax_pos = ax.get_position()
        # Shift the subplot to the desired location + set the width so it's equal for all plots
        if n == 0:
            new_pos_y = ax_pos.y0
        else:
            new_pos_y = ax_pos_prev.y0 - ax_pos.height - y_offset
        ax.set_position([ax_pos.x0, new_pos_y, ax_pos.width, ax_pos.height])
        ax_pos_prev = ax.get_position()

# Separate function for adding the dotplot legend
def add_dotplot_legend(dotplot_ax, df, space_per_dot, vmin=0, vmax=None, hue='total_count', size='proportion', cbar_title='Number of cells', size_title='Distribution (per layer)', cbar_height=0.8,
                    cbar_width=2):
    fig = plt.gcf()
    # -- Add a colorbar
    if vmax == None:
        vmax = df[hue].max()
    if vmin == None:
        vmin = df[hue].min()
    norm = plt.Normalize(vmin, vmax)
    sm = plt.cm.ScalarMappable(cmap='Reds', norm=norm)
    sm.set_array([])
    cbar_width_adj = space_per_dot*cbar_width
    cbar_height_adj = space_per_dot*cbar_height
    cax = fig.add_axes([dotplot_ax.get_position().x1 + 0.1, dotplot_ax.get_position().y0, cbar_width_adj, cbar_height_adj])
    dotplot_ax.figure.colorbar(sm, cax=cax, orientation='horizontal')
    cax.set_title(cbar_title, fontsize=10)
    # -- Add size legend
    # Move the legend to the top of the plot & remove hue variable legend
    handles, labels = dotplot_ax.get_legend_handles_labels()
    # Find the index of 'proportion'
    index = labels.index(size)
    # Get everything after 'proportion'
    size_labels = labels[index + 1:]
    size_handles = handles[index + 1:]
    # Make the spots gray
    for handle in size_handles:
        handle.set_markerfacecolor('lightgray')
        handle.set_linewidth(0.8) 
    # Remove original legend
    dotplot_ax.legend().remove()
    # Make new size legend
    legend = cax.legend(size_handles, size_labels, loc=(0.2, 4), labelspacing=0, 
                            ncol=len(size_handles), frameon=False, title=f'{size_title}\n', handlelength=-2)
    # Move the legend labels under the markers
    for txt in legend.get_texts():
        txt.set_ha("center") # horizontal alignment of text item
        txt.set_x(-1) # x-position
        txt.set_y(-13) # y-position
    return dotplot_ax, cax

# Function to format the counts ax
def format_counts_ax(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_facecolor('none') 
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_yticklabels('')
    ax.set_xticklabels('')
    ax.set_xlabel('')
    ax.set_ylabel('')

# Function to format numbers
def reformat_number(number, dec=1):
    if number >= 1000000:
        # Divide the number by 1,000,000 and format it to one decimal place
        formatted_number = "{:.1f}m".format(number / 1000000)
    elif number >= 1000:
        # Divide the number by 1000 and format it to one decimal place
        formatted_number = "{:.1f}k".format(number / 1000)
    else:
        if dec == 0:
            formatted_number = str(int(round(number, dec)))
        else:
            formatted_number = str(round(number, dec))
    return formatted_number

# Function to add counts on x & y
def add_counts_bars(dotplot_ax, cax, count_df_x, count_df_y, space_per_dot=0.2, annotate=False, annotation_size=8,
                    palette_x=None, palette_y=None, alpha_x=1, alpha_y=1, counts_x=True, counts_y=True, order_ylabels=None,
                    palette_x_dict=None, order_x=None):
    # Get the current figure to add axes
    fig = plt.gcf()
    dot_ax_pos = dotplot_ax.get_position()
    new_ax_x = None
    new_ax_y = None
    # -- Counts y:
    if counts_y:
        new_ax_y = fig.add_axes([0.5, 0.5, 0.4, 0.4])  # [left, bottom, width, height]
        new_ax_y.set_position([dot_ax_pos.x1, dot_ax_pos.y0, 1.5*space_per_dot,  dot_ax_pos.height])
        # Sort count_df_according to dotplot_ax
        count_df_y = count_df_y.reindex(order_ylabels)
        # Make a horizontal barplot
        if palette_y == None:
            palette_y = 'salmon'
        count_df_y['total_count'].plot(kind='barh', ax=new_ax_y, legend=False, width=0.6, edgecolor='black', linewidth=1, color=palette_y,
                        alpha=alpha_y)
        # Set ylim to be equal to dotplotax
        dotplot_ylim = dotplot_ax.get_ylim()
        new_ax_y.set_ylim((dotplot_ylim[0], dotplot_ylim[1]))
        # Format count ax
        format_counts_ax(new_ax_y)
        # Move the legend
        cax_pos_prev = cax.get_position()
        new_ax_pos = new_ax_y.get_position()
        cax.set_position([new_ax_pos.x1 + 0.01, cax_pos_prev.y0, cax_pos_prev.width, cax_pos_prev.height])
    # -- Counts x:
    if counts_x:
        # Remove the title from the dotplot ax
        title = dotplot_ax.get_title()
        if title == 'Additional_features':
            y = 1.5
        else:
            y = None
        dotplot_ax.set_title('')
        # Add ax for the counts plot
        new_ax_x = fig.add_axes([0.5, 0.5, 0.4, 0.4])  # [left, bottom, width, height]
        new_ax_x.set_position([dot_ax_pos.x0, dot_ax_pos.y1, dot_ax_pos.width,  1.5*space_per_dot])
        # Sort according to dotplot_ax
        if order_x == None:
            order_x = [label.get_text() for label in dotplot_ax.get_xticklabels()]
        count_df_x = count_df_x.reindex(index=order_x)
        # Make the barplot (vertical)
        if palette_x == None:
            palette_x = 'salmon'
        if palette_x_dict:
            palette_x = [palette_x_dict[label] for label in order_x]
        count_df_x.plot(kind='bar', y='count', ax=new_ax_x, width=0.6, legend=False, edgecolor='black', linewidth=1, color=palette_x,
                        alpha=alpha_x)
        # Set xlim equal to dotplot
        dotplot_xlim = dotplot_ax.get_xlim()
        new_ax_x.set_xlim((dotplot_xlim[0], dotplot_xlim[1]))
        # Format count ax
        format_counts_ax(new_ax_x)
        # Set the title
        new_ax_x.set_title(title, fontweight='bold', y=y)
    # Annotate the bars
    if annotate:
        if counts_x:
            # Annotate the vertical bars
            for index, value in enumerate(count_df_x['count']):
                reformat_value = reformat_number(value)
                new_ax_x.annotate(str(reformat_value), xy=(index, value), ha='center', va='bottom', fontsize=annotation_size)
            # Increase the margin to accommodate
            margin = 1.15  # Set the margin
            ymax = new_ax_x.get_ylim()[1]
            new_ax_x.set_ylim(0, ymax * margin)
        if counts_y:
            # Annotate the horizontal bars
            for index, value in enumerate(count_df_y['total_count']):
                reformat_value = reformat_number(value, dec=0)
                new_ax_y.annotate(str(reformat_value), xy=(value, index), xytext=(2, 0),  # 2 points horizontal offset.
                                textcoords='offset points', ha='left', va='center_baseline', fontsize=annotation_size)
    return new_ax_x, new_ax_y

# Make new dfs with total counts & proportion
def calculate_stats(df, layer, cell_type):
    # Calculate total count per category
    total_counts = df.groupby(layer, observed=False)[cell_type].sum().reset_index()
    total_counts.columns = ['category', 'total_count']
    # Calculate proportion of counts in each category
    total_counts['proportion'] = total_counts['total_count'] / total_counts['total_count'].sum()
    total_counts['cell_type'] = cell_type
    return total_counts

# Make new dfs with total counts & proportion - per sample
def calculate_stats_per_sample(df, layer, cell_type, subtype_dict, sample_column='sample_id', subtype_column='subtype'):
    # Calculate total count per category
    total_counts = df.groupby([layer, sample_column], observed=False)[cell_type].sum().reset_index()
    total_counts.columns = ['category', sample_column, 'total_count']
    # Calculate proportion of counts in each category
    sum_per_sample = total_counts.groupby([sample_column], observed=False)['total_count'].sum().to_frame().reset_index()
    sum_per_sample.columns = [sample_column, 'total_per_sample']
    total_counts = pd.merge(total_counts, sum_per_sample, on=[sample_column])
    total_counts['proportion'] = total_counts['total_count'] / total_counts['total_per_sample']
    total_counts['cell_type'] = cell_type
    total_counts[subtype_column] = total_counts['sample_id'].map(subtype_dict)
    return total_counts

def add_lineage_legend(lineage_lut, ax):
    # Create a list of handles and labels for the legend
    handles_lineage = []
    labels_lineage = []
    for cat, color in lineage_lut.items():
        handles_lineage.append(mpatches.Patch(facecolor=color, edgecolor='black', linewidth=0.8))
        labels_lineage.append(cat)
    # Create the legend using the Patch objects
    legend1 = ax.legend(handles=handles_lineage, labels=labels_lineage, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10, 
                                title='Lineage', handlelength=0.7, frameon=False)
    # Legend position
    x_offset = 1.4
    y_offset = 1
    legend1.set_bbox_to_anchor((x_offset, y_offset))
    legend1.get_title().set_fontweight('bold')
    # Set the alignment of the legend title to the left
    legend1._legend_box.align = "left"

def add_subtype_legend(ax, sorted_subtypes, color_lut):
    # Create a list of handles and labels for the legend
    handles_subtype = []
    labels_subtype = []
    for subtype in sorted_subtypes:
        handles_subtype.append(mpatches.Patch(facecolor=color_lut.loc[subtype, 'color'], edgecolor='black', linewidth=0.8))
        labels_subtype.append(subtype)
    # Create the legend using the Patch objects
    legend1 = ax.legend(handles=handles_subtype, labels=labels_subtype, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10, 
                                title='Subtype', handlelength=0.7, frameon=False)
    # Legend position
    x_offset = 1.1
    y_offset = 1.5
    legend1.set_bbox_to_anchor((x_offset, y_offset))
    legend1.get_title().set_fontweight('bold')
    # Set the alignment of the legend title to the left
    legend1._legend_box.align = "left"

# Now wrap all of these into a single function to avoid repetition
def make_complete_dotplot(dot_plot_df_dict, count_df_x_dict, count_df_y_dict, available_layers, vmax, color_dict, lineage_palette, lineage_lut, all_cell_types, space_per_dot=0.2, 
                          hue='total_count', size='proportion', cbar_title='Average # cells per spot', size_title='Distribution (per layer)',
                          offset=0.02, plot_title=None, adjust_size_max=False, title_suffix=None, cbar_height=0.8):
    # Prepare the axs
    axs, grid_spec = prepare_dotplot(dot_plot_df_dict)
    # Adjust max size of dots
    if adjust_size_max:
        highest_value = 0
        for layer in available_layers:
            key = f'{layer}'
            dot_plot_df = dot_plot_df_dict[key]
            max_value = dot_plot_df[size].max()
            if max_value > highest_value:
                highest_value = max_value
        if 1 > highest_value:
            size_max = round(highest_value, 2)
        else:
            size_max = 1
    else:
        size_max = 1
    # Make the dotplot
    for n, layer in enumerate(available_layers):
        key = f'{layer}'
        dot_plot_df = dot_plot_df_dict[key]
        seaborn_dotplot(dot_plot_df, dotplot_ax=axs[n], hue=hue, size=size, space_per_dot=space_per_dot, ymargin=-0.8, xmargin=0.6,
            title=f'{layer.capitalize()}', vmax=vmax, cbar_width=1.8, legend=False, cbar_title=cbar_title, size_title=size_title,
            size_max=size_max)
    # Move plots to stop them from overlapping
    move_plots(axs, offset=offset)
    # Add a legend to the last plot
    dotplot_ax, cax = add_dotplot_legend(dotplot_ax=axs[-1], df=dot_plot_df, hue=hue, size=size, cbar_title=cbar_title,
                                          space_per_dot=space_per_dot, vmax=vmax, size_title=size_title, cbar_height=cbar_height)
    # Remove the other legends
    for ax in axs[:-1]:
        ax.get_legend().remove()
    new_ax_x_masterlist = []
    # Add the counts barplots
    for n, layer in enumerate(available_layers):
        key = f'{layer}'
        count_df_x = count_df_x_dict[key]
        count_df_y = count_df_y_dict[key]
        ymax = 0
        if n < len(available_layers)-1:
            new_ax_x, new_ax_y = add_counts_bars(axs[n], cax, count_df_x, count_df_y, counts_y = False, space_per_dot=space_per_dot, annotate=True, palette_y=lineage_palette,
                                                palette_x_dict=color_dict)
        else:
            new_ax_x, new_ax_y = add_counts_bars(axs[n], cax, count_df_x, count_df_y, space_per_dot=space_per_dot, annotate=True, palette_y=lineage_palette,
                                                palette_x_dict=color_dict, order_ylabels=all_cell_types)
            # Position the colorbar
            position_cbar(new_ax_y, cax, offset=0.05)
            add_lineage_legend(lineage_lut, new_ax_y)
        # Make sure the top counts share the same y axis
        new_ax_x_masterlist.append(new_ax_x)
        if new_ax_x.get_ylim()[1] > ymax:
            ymax = new_ax_x.get_ylim()[1]
    # Set top ax y equal
    for ax in new_ax_x_masterlist:
        ax.set_ylim(top=ymax)
    # Add labels
    new_ax_x_masterlist[0].set_ylabel('# spots', fontsize=8)
    new_ax_y.set_title('# cells', fontsize=8, y=0.99, x=0.4)
    # Add a title
    title_x = (axs[0].get_position().x0 + new_ax_y.get_position().x1)/2
    title_y = new_ax_x.get_position().y1 + 0.6
    plt.suptitle(f'{plot_title} (vmax={vmax}, size_max={size_max}, {title_suffix})', x=title_x, y=title_y, fontweight='bold')


# Define some functions to add lines between the annotation layers & cell type lineages
# Function to add lines
def add_v_lines(ordered_major_cat_dict, ax):
    where_to_draw = []
    keys = list(ordered_major_cat_dict.keys())
    for i, key in enumerate(keys):
        if i == len(keys) - 1:
            break  # Don't draw after the last group
        total = sum(where_to_draw)
        pos = total + len(ordered_major_cat_dict[key])
        where_to_draw.append(len(ordered_major_cat_dict[key]))
        ax.axvline(pos - 0.5, color='gray', linewidth=0.5, linestyle='--')

# Function to add lines
def add_h_lines(classifications_sorted, ax, decon_level):
    values = list(classifications_sorted['lineage'].unique())
    # Drop duplicates based on decon level
    classifications_sorted = classifications_sorted.drop_duplicates(subset=[f'{decon_level}'])
    value_counts = classifications_sorted['lineage'].value_counts().to_frame()
    where_to_draw = []
    for value in values:
        if pd.isna(value):
            continue  # Skip nan values
        total = sum(where_to_draw)
        pos = total + value_counts.loc[value, 'count']
        where_to_draw.append(pos)
        ax.axhline(pos - 0.5, color='gray', linewidth=0.5, linestyle='--')


# Function for making a single dotplot
def make_single_dotplot(dot_plot_df, count_df_x, count_df_y, plotting_colors, vmax, all_cell_types, lineage_lut, lineage_palette, ordered_major_cat_dict, classifications_sorted,
                        space_per_dot=0.2, hue='total_count', size='proportion', decon_level='minor',
                        cbar_title='Average # cells per spot', size_title='Distribution', offset=0.02, title='# cells per annotation', plot_title=None,
                        lines=True, adjust_size_max=False, title_suffix=None):
    # Prepare the axs
    fig, ax = plt.subplots(figsize=(3, 1.3),)
    # Adjust max size of dots
    if adjust_size_max:
        highest_value = dot_plot_df[size].max()
        if 1 > highest_value:
            size_max = round(highest_value, 2)
        else:
            size_max = 1
    else:
        size_max = 1
    # Make the dotplot
    dotplot_ax, cax = seaborn_dotplot(dot_plot_df, dotplot_ax=ax, hue=hue, cbar_title=cbar_title, size=size, space_per_dot=space_per_dot, ymargin=-0.8, xmargin=0,
        title=f'{title}', vmax=vmax, cbar_width=1.8, legend=True, size_title=size_title, size_max=size_max)
    # Add the counts barplots
    new_ax_x, new_ax_y = add_counts_bars(dotplot_ax, cax, count_df_x, count_df_y, space_per_dot=space_per_dot, annotate=True, palette_y=lineage_palette,
                                        palette_x=plotting_colors, order_ylabels=all_cell_types)
    # Position the colorbar
    position_cbar(new_ax_y, cax, offset=0.05)
    # Add lineage legend
    add_lineage_legend(lineage_lut, new_ax_y)
    # Add labels
    new_ax_x.set_ylabel('# spots', fontsize=8)
    new_ax_y.set_title('# cells', fontsize=8, y=0.99, x=0.4)
    # Add lines between layers
    if lines:
        add_v_lines(ordered_major_cat_dict, ax)
        add_h_lines(classifications_sorted, ax, decon_level)
    # Add a title
    title_x = (ax.get_position().x0 + ax.get_position().x1)/2
    title_y = new_ax_x.get_position().y1 + 0.6
    plt.suptitle(f'{plot_title} (vmax={vmax}, size_max={size_max}, {title_suffix})', x=title_x, y=title_y, fontweight='bold')
    return ax

def add_subtype_colors(axs, sorted_subtypes, color_lut, alpha=0.7):
    # Add subtype colors
    for n, subtype in enumerate(sorted_subtypes):
        # Add subtype color at top & bottom
        # Add colors for subtypes
        color = color_lut.loc[subtype, 'color']
        order_y = [label.get_text() for label in axs[n][0].get_yticklabels()]
        # Subset & sort
        xlim = axs[n][0].get_xlim()
        ylim = axs[n][0].get_ylim()
        height = ylim[0] - ylim[1]
        width = 0.3
        # Markings at left
        rect = mpatches.Rectangle((xlim[0], ylim[1]), width=width, height=height, facecolor=color, alpha=alpha, zorder=-1)
        axs[n][0].add_patch(rect)
        # Markings at right
        xlim = axs[n][-1].get_xlim()
        rect = mpatches.Rectangle((xlim[1] - width, ylim[1]), width=width, height=height, facecolor=color, alpha=alpha, zorder=-1)
        axs[n][-1].add_patch(rect)

def add_subtype_colors_y_only(axs, sorted_subtypes, color_lut, alpha=0.7):
    # Add subtype colors
    for n, subtype in enumerate(sorted_subtypes):
        # Add subtype color at top & bottom
        # Add colors for subtypes
        color = color_lut.loc[subtype, 'color']
        order_y = [label.get_text() for label in axs[n].get_yticklabels()]
        # Subset & sort
        xlim = axs[n].get_xlim()
        ylim = axs[n].get_ylim()
        height = ylim[0] - ylim[1]
        width = 0.3
        # Markings at left
        rect = mpatches.Rectangle((xlim[0], ylim[1]), width=width, height=height, facecolor=color, alpha=alpha, zorder=-1)
        axs[n].add_patch(rect)
        # Markings at right
        xlim = axs[n].get_xlim()
        rect = mpatches.Rectangle((xlim[1] - width, ylim[1]), width=width, height=height, facecolor=color, alpha=alpha, zorder=-1)
        axs[n].add_patch(rect)

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