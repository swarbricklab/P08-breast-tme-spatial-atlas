#!/usr/bin/env python3

# --- Aim
# Plot Figure 1G: cell type composition within different pathology annotation categories
# Shows deconvolution results stratified by tissue compartments (invasive, in situ, stroma, etc.)
# This creates Figure 1G

# --- Define environment
# env: spatial-pilot

# --- Load packages

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.ticker import FuncFormatter

import matplotlib.pyplot as plt
import matplotlib as mpl


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
parser.add_argument("--N_cells_per_location", help="Value for N_cells_per_location parameter cell2location")
parser.add_argument("--decon_level", help="Level of annotations used for decon")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file")

args = vars(parser.parse_args())

# --- Arguments
project_dir = args["project_dir"]
h5ad_dir = args["h5ad_dir"]
res_dir = args["res_dir"]
sample_list = args["samples"]
N_cells_per_location = args["N_cells_per_location"]
decon_level = args["decon_level"]

print(f'Samples in python: {sample_list}')

# Import function script as a module
sys.path.append(os.path.join(project_dir, 'scripts', 'figures'))
import dotplot_functions


# --- Extract some info from the sample sheet

# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_dir = os.path.join(project_dir, 'config')
    sample_sheet_file = 'sample_sheet_visium.csv'
    sample_sheet_path = os.path.join(sample_dir, sample_sheet_file)

sample_sheet = pd.read_csv(sample_sheet_path)


# --- Also get some tissue metadata

resource_dir = os.path.join(project_dir, 'resources')
metadata_file = 'tissue_metadata.csv'

tissue_metadata = pd.read_csv(os.path.join(resource_dir, metadata_file))

# Merge the two DataFrames on the common column 'tissue_id'
merged_df = pd.merge(sample_sheet, tissue_metadata, on='tissue_id')

# Create a sample - subtype dictionary
subtype_dict = dict(zip(merged_df['section_id'], merged_df['clinical_subtype']))


# --- Some annotation parameters

# Use all available layers
annotation_layers = ['epithelial', 'immune', 'stroma', 'additional_features']


# --- Order the annotations

# Get the desired order
annotation_dir_base = os.path.join(project_dir, 'resources', 'pathology_annotations')
# Load the desired categories per layer (to get the desired order)
filename = 'pathology_categories_per_layer.json'
path_to_file = os.path.join(annotation_dir_base, filename)

with open(path_to_file, 'r') as json_file:
                annotations_of_interest_dict = json.load(json_file)


# --- Load the deconvolution results

decon_level = 'major'
decon_dir = os.path.join(project_dir, 'results', '02_cell_type_mapping', 'deconvolution', 'visium', f'{decon_level}_level')

decon_dict = {}

for sample in sample_list:
    decon_filename = f'{sample}_{decon_level}_q05_cell_abundance_w_sf.csv'
    # Load decon results
    decon_results = pd.read_csv(os.path.join(decon_dir, decon_filename), index_col=0)
    decon_dict[sample] = decon_results

all_cell_types = list(decon_dict[sample].columns)


# --- Load the pathology annotations

# Define annotation path
reworked_anno_dir = os.path.join(project_dir, 'data', 'pathology_annotations', 'visium')

anno_df_dict = {}

for sample in sample_list:
    annotation_file = f'{sample}_pathology_annotations.csv'
    annotation_path = os.path.join(reworked_anno_dir, annotation_file)
    if os.path.exists(annotation_path):
        path_anno = pd.read_csv(annotation_path, index_col=0)
        anno_df_dict[sample] = path_anno
    else:
        print(f'No annotations found for {sample} - skipping')


# --- Load miniatlas qualifications
# Load miniatlas classifications

classifications_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
classifications_file = 'miniatlas_subset_classifications.csv'

classifications = pd.read_csv(os.path.join(classifications_dir, classifications_file))

# Load the colors from the miniatlas for visualization
color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
lineage_colors_file = 'miniatlas_lineage_colors.csv'
miniatlas_major_colors_file = 'miniatlas_major_colors.csv'

lineage_colors = pd.read_csv(os.path.join(color_dir, lineage_colors_file))
miniatlas_major_colors = pd.read_csv(os.path.join(color_dir, miniatlas_major_colors_file))

# Make some lookup tables
lineage_lut = dict(zip(lineage_colors['lineage'], lineage_colors['color']))
miniatlas_major_lut = dict(zip(miniatlas_major_colors['miniatlas_major'], miniatlas_major_colors['color']))

# Add nan to lut
lineage_lut[np.nan] = 'lightgray'
miniatlas_major_lut[np.nan] = 'lightgray'

# Make a dataframe containing the colors for the rows of the heatmap (lineage and miniatlas_major classifications)
row_colors_df = pd.DataFrame(zip(classifications['lineage'].map(lineage_lut), classifications['miniatlas_major'].map(miniatlas_major_lut)), 
                             index=classifications[f'miniatlas_{decon_level}'], columns=['lineage', 'miniatlas_major'])

# Drop duplicates
row_colors_df = row_colors_df[~row_colors_df.index.duplicated(keep='first')]

# Sort the decon results on lineage and major
group_columns = ['lineage', 'miniatlas_major']

# Subset classifications to cell types present
# Standardize the names in classifications
classifications = classifications[classifications[f'miniatlas_{decon_level}'].isin(all_cell_types)]

classifications_sorted = classifications.sort_values(by=group_columns)
sorted_indices = list(classifications_sorted[f'miniatlas_{decon_level}'])

# Drop duplicates
# Build an index mapping from reference
sorted_indices = list(dict.fromkeys(sorted_indices))
order = {val: i for i, val in enumerate(sorted_indices)}

# Sort all cell types according to the sorted_indices
all_cell_types = sorted(all_cell_types, key=lambda x: order.get(x, float('inf')))


# --- Merge into a single df with sample_id

# First make single df that has annotations + decon results per sample
anno_decon_dict = {}

for sample in sample_list:
    anno_decon_dict[sample] = anno_df_dict[sample]
    anno_decon_dict[sample][all_cell_types] = decon_dict[sample][all_cell_types]
    anno_decon_dict[sample]['barcode'] = anno_decon_dict[sample].index

# Concatenate DataFrames and add 'sample_id' column
anno_decon_df = pd.concat([df.assign(sample_id=key) for key, df in anno_decon_dict.items()], ignore_index=True)


# --- Order the categories

ordered_major_cat_dict = {}

for layer in annotation_layers:
    key = f'{layer}'
    categories = list(anno_decon_df[key].dropna().unique())
    # Sort according to predefined cateogires
    sorted_categories = sorted(categories, key=lambda x: annotations_of_interest_dict[layer].index(x) if x in annotations_of_interest_dict[layer] else len(annotations_of_interest_dict[layer]))
    ordered_major_cat_dict[layer] = sorted_categories
    anno_decon_df[key] = pd.Categorical(anno_decon_df[key], categories=sorted_categories, ordered=True)


# --- Get  colors per annotation

ordered_columns_major = [item for sublist in ordered_major_cat_dict.values() for item in sublist]

color_dict, plotting_colors, top_cat_colors = dotplot_functions.create_plotting_colors(ordered_major_cat_dict, ordered_columns_major, 
                                                                                       available_layers=annotation_layers, mixed=False, major_categories=True)


# --- Prepare for plotting

# Disable the grid for all plots
plt.rcParams['axes.grid'] = False


# --- Plot 1: Average number of cells per annotation (layers separate)

# -- Prepare for plotting

# Reorder the cell types in the decon results according to the order in the classifications
row_colors_df = row_colors_df[~row_colors_df.index.duplicated(keep='first')]
lineage_palette = list(row_colors_df.reindex(all_cell_types)['lineage'])

# Find missing cell types
missing_cell_types = [ct for ct, color in zip(all_cell_types, lineage_palette) if pd.isna(color)]
print("Cell types missing from lineage_palette:", missing_cell_types)

# Make a plotting df per layer
dot_plot_df_dict = {}
count_df_x_dict = {}
count_df_y_dict = {}

for n, layer in enumerate(annotation_layers):
    key = f'{layer}'
    to_select = list(all_cell_types) + [key]
    raw_df = anno_decon_df[to_select]
    # Fill the NaNs
    na_fill = f'non_{layer}'
    raw_df[key] = raw_df[key].cat.add_categories(na_fill)
    raw_df[key] = raw_df[key].fillna(f'non_{layer}')
    # Make df for dotplot
    dot_plot_df = pd.DataFrame()
    for cell_type in all_cell_types:
        df_tmp = dotplot_functions.calculate_stats(raw_df, key, cell_type)
        dot_plot_df = pd.concat([dot_plot_df, df_tmp], axis=0)
    # Make dfs for the counts plots
    count_df_y = dot_plot_df.groupby('cell_type')['total_count'].sum().to_frame()
    count_df_x = raw_df[key].value_counts().to_frame()
    count_df_x = count_df_x.reset_index()
    count_df_x.index = count_df_x[key]
    # Save in dict
    dot_plot_df_dict[key] = dot_plot_df
    count_df_x_dict[key] = count_df_x
    count_df_y_dict[key] = count_df_y

# Plotting
space_per_dot = 0.2
vmax=3
hue = 'average_count'
size = 'average_proportion'

# Prepare for calculating proportions per spot
rowsums = raw_df[all_cell_types].sum(axis=1)
proportion_per_spot = raw_df[all_cell_types].div(rowsums, axis=0)

# Prepare the df
for n, layer in enumerate(annotation_layers):
    key = f'{layer}'
    dot_plot_df = dot_plot_df_dict[key].copy()
    count_df_x = count_df_x_dict[key].copy()
    # Normalize the cell counts by number of spots per annotation
    # Reset the index to the default integer index
    count_df_x.reset_index(drop=True, inplace=True)
    count_df_x.columns = ['category', 'number_of_spots']
    dot_plot_df = pd.merge(dot_plot_df, count_df_x, on='category')
    dot_plot_df[hue] = dot_plot_df['total_count'] / dot_plot_df['number_of_spots']
    # Calculate total count per cell type
    total_count_per_cell_type = dot_plot_df.groupby('cell_type')[hue].transform('sum')
    # Also calculate the average proportion
    prop_anno_df = proportion_per_spot.copy()
    prop_anno_df[key] = anno_decon_df[key]
    # Fill the NaNs
    na_fill = f'non_{layer}'
    prop_anno_df[key] = prop_anno_df[key].cat.add_categories(na_fill)
    prop_anno_df[key] = prop_anno_df[key].fillna(na_fill)
    # Make df for dotplot
    prop_df = pd.DataFrame()
    for cell_type in all_cell_types:
        df_tmp = dotplot_functions.calculate_stats(prop_anno_df, key, cell_type)
        prop_df = pd.concat([prop_df, df_tmp], axis=0)
    # Add the number of spots
    prop_df = pd.merge(prop_df, count_df_x, on='category')
    prop_df[size] = prop_df['total_count'] / prop_df['number_of_spots']
    # Merge the dfs
    dot_plot_df = pd.merge(dot_plot_df, prop_df[['category', 'cell_type', size]], on = ['category', 'cell_type'])
    # Sort in correct order
    dot_plot_df['cell_type'] = pd.Categorical(dot_plot_df['cell_type'], categories=all_cell_types, ordered=True)
    all_cats = ordered_major_cat_dict[layer] + [f'non_{layer}']
    dot_plot_df['category'] = pd.Categorical(dot_plot_df['category'], categories=all_cats, ordered=True)
    dot_plot_df = dot_plot_df.sort_values(by=['cell_type', 'category'])
    dot_plot_df_dict[key] = dot_plot_df

lineage_palette = [
    color if pd.notna(color) else 'lightgray'
    for color in row_colors_df.reindex(all_cell_types)['lineage']
]

if decon_level == 'major':
    dot_size_max = 200
else:
    dot_size_max = 400

# Function to format counts
def format_count(value, tick_number):
    if value >= 1000000:
        return f'{value/1000000:.1f}M'
    if value >= 1000:
        return f'{value/1000:.0f}k'
    return f'{value:.0f}'

# Calculate max count across all layers for consistent y-axis in bar plots (excluding 'non_' categories)
all_counts_spots_plot1 = []
for layer in annotation_layers:
    key = f'{layer}'
    categories_in_layer = list(dot_plot_df_dict[key]['category'].unique())
    # Exclude categories that start with 'non_'
    categories_in_layer = [cat for cat in categories_in_layer if not str(cat).startswith('non_')]
    counts = [count_df_x_dict[key].loc[cat, 'count'] if cat in count_df_x_dict[key].index else 0 for cat in categories_in_layer]
    all_counts_spots_plot1.extend(counts)

max_count_spots_plot1 = max(all_counts_spots_plot1) if all_counts_spots_plot1 else 0

# Calculate max cell count for consistent x-axis in cell count bars
max_cell_count_plot1 = max([count_df_y_dict[layer]['total_count'].max() for layer in annotation_layers])

# --- Create multi-panel plot with proportion_per_layer as size
space_per_dot_cell = 0.03
space_per_dot_cat = 0.03

# Filter out 'non_' categories for Plot 1b
dot_plot_df_dict_filtered = {}
for layer in annotation_layers:
    df_layer = dot_plot_df_dict[layer].copy()
    # Exclude categories that start with 'non_'
    df_layer = df_layer[~df_layer['category'].astype(str).str.startswith('non_')]
    dot_plot_df_dict_filtered[layer] = df_layer

# Calculate global size range for consistent sizing across all plots
all_size_values_plot1 = []
for layer in annotation_layers:
    all_size_values_plot1.extend(dot_plot_df_dict_filtered[layer][size].values)

size_min_plot1 = min(all_size_values_plot1)
size_max_plot1 = max(all_size_values_plot1)


# --- Plot 1: Dotplot of z-scores of average counts per category, with size based on proportion_per_layer

# Calculate z-scores for each cell type across all categories (all layers combined)
# First, create a combined dataframe with all layers, excluding 'non_' categories
combined_df = pd.DataFrame()
for layer in annotation_layers:
    df_layer = dot_plot_df_dict[layer].copy()
    # Exclude categories that start with 'non_'
    df_layer = df_layer[~df_layer['category'].astype(str).str.startswith('non_')]
    combined_df = pd.concat([combined_df, df_layer], axis=0)

# Calculate z-score for each cell type
combined_df_zscore = combined_df.copy()
for cell_type in all_cell_types:
    mask = combined_df_zscore['cell_type'] == cell_type
    values = combined_df_zscore.loc[mask, 'average_count']
    mean_val = values.mean()
    std_val = values.std()
    if std_val > 0:
        combined_df_zscore.loc[mask, 'zscore'] = (values - mean_val) / std_val
    else:
        combined_df_zscore.loc[mask, 'zscore'] = 0

# Split back into layers
dot_plot_df_dict_zscore = {}
for layer in annotation_layers:
    df_layer_zscore = combined_df_zscore[combined_df_zscore['category'].isin(dot_plot_df_dict[layer]['category'].unique())].copy()
    df_layer_zscore['category'] = pd.Categorical(df_layer_zscore['category'], 
                                                   categories=dot_plot_df_dict[layer]['category'].cat.categories, 
                                                   ordered=True)
    df_layer_zscore['cell_type'] = pd.Categorical(df_layer_zscore['cell_type'], 
                                                    categories=all_cell_types, 
                                                    ordered=True)
    df_layer_zscore = df_layer_zscore.sort_values(['cell_type', 'category'])
    dot_plot_df_dict_zscore[layer] = df_layer_zscore

# Calculate global z-score range for color normalization
all_zscore_values_plot1 = []
for layer in annotation_layers:
    all_zscore_values_plot1.extend(dot_plot_df_dict_zscore[layer]['zscore'].values)

zscore_min_plot1 = min(all_zscore_values_plot1)
zscore_max_plot1 = max(all_zscore_values_plot1)
zscore_abs_max_plot1 = max(abs(zscore_min_plot1), abs(zscore_max_plot1))

# Create figure
fig = plt.figure(figsize=(8, 8))
axs = []
bar_axs = []

# Calculate positions for arranging horizontally
start_x = 0.08
base_y = 0.15

for n, layer in enumerate(annotation_layers):
    df_layer = dot_plot_df_dict_zscore[layer]
    # Calculate dimensions
    num_cell_types = len(all_cell_types)
    num_categories = len(df_layer['category'].unique())
    plot_width = space_per_dot_cat * num_categories
    plot_height = space_per_dot_cell * num_cell_types
    # Position for this layer's dotplot (horizontal arrangement)
    if n == 0:
        ax_x = start_x
        prev_x = ax_x + plot_width
    else:
        ax_x = prev_x + 0.007
        prev_x = ax_x + plot_width
    # Create main dotplot axis
    ax = fig.add_axes([ax_x, base_y, plot_width, plot_height])
    axs.append(ax)
    # Create scatterplot with z-score as color
    sns.scatterplot(
        data=df_layer,
        x='category',
        y='cell_type',
        hue='zscore',
        size=size,
        sizes=(0, dot_size_max),
        size_norm=(size_min_plot1, size_max_plot1),
        palette='RdBu_r',
        edgecolor='black',
        linewidth=0.5,
        ax=ax,
        hue_norm=(-zscore_abs_max_plot1, zscore_abs_max_plot1),
        legend='brief' if n == len(annotation_layers) - 1 else False
    )
    # Format axes
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title('')
    # Format x-axis
    ax.set_xticks(range(num_categories))
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='center')
    # Format y-axis - only show labels on first plot
    ax.set_yticks(range(len(all_cell_types)))
    if n == 0:
        ax.set_yticklabels(all_cell_types, rotation=0, ha='right')
    else:
        ax.set_yticklabels([])
    # Set limits
    ax.set_xlim(-0.75, num_categories - 0.25)
    ax.set_ylim(len(all_cell_types) - 0.25, -0.75)
    # Add grid
    ax.grid(True, which='both', axis='both', linestyle='-', linewidth=0.5, alpha=1)
    ax.set_axisbelow(True)
    # Remove legend
    if ax.get_legend():
        ax.get_legend().remove()
    # Create vertical barplot above
    bar_height = 0.06
    bar_ax = fig.add_axes([ax_x, base_y + plot_height + 0.007, plot_width, bar_height])
    bar_ax.set_facecolor('none')
    bar_axs.append(bar_ax)
    # Get categories and counts
    categories_in_layer = list(df_layer['category'].unique())
    count_df_x = count_df_x_dict[layer]
    counts = [count_df_x.loc[cat, 'count'] if cat in count_df_x.index else 0 for cat in categories_in_layer]
    # Create vertical barplot
    x_positions = range(len(categories_in_layer))
    colors = [color_dict.get(cat, 'lightgray') for cat in categories_in_layer]
    bar_ax.bar(x_positions, counts, color=colors, edgecolor='black', linewidth=0.8, alpha=0.9)
    # Format barplot
    bar_ax.set_xlim(-0.5, len(categories_in_layer) - 0.5)
    bar_ax.set_ylim(0, max_count_spots_plot1 * 1.1)
    bar_ax.set_xticks([])
    bar_ax.spines['top'].set_visible(False)
    bar_ax.spines['right'].set_visible(False)
    # Add text annotations above each bar
    for i, count in enumerate(counts):
        buffer = max_count_spots_plot1 * 0.1
        if count >= 1000:
            label = f'{count/1000:.1f}k'
        else:
            label = str(count)
        bar_ax.text(i, count + buffer, label, ha='center', va='bottom', fontsize=8, rotation=90)
    # If not the first barplot, remove the ytick at 0
    if n != 0:
        bar_ax.set_yticks([bar_ax.get_yticks()[1]])
    # Only show ylabel on first barplot
    if n == 0:
        bar_ax.set_ylabel('# spots', fontsize=9)
        bar_ax.tick_params(axis='y', labelsize=9)
        if max_count_spots_plot1 > 1000:
            bar_ax.yaxis.set_major_formatter(FuncFormatter(format_count))
    else:
        bar_ax.set_yticklabels([])
        bar_ax.tick_params(axis='y', labelsize=9)
    # Add layer title
    bar_ax.set_title(layer.capitalize().replace('_', ' '), fontweight='bold', pad=20)

# Format x-axis labels
caps = ['tls', 'lvi']
for ax in axs:
    xticks = ax.get_xticklabels()
    new_labels = []
    for label in xticks:
        text = label.get_text()
        if text == 'immune_other':
            text = 'lymphocytes present'
        if text in caps:
            text = text.upper()
        else:
            text = text.capitalize()
        text = text.replace('_', ' ')
        new_labels.append(text)
    ax.set_xticklabels(new_labels, rotation=90, ha='center')

# Add cell count barplots on the right side
cell_count_bar_axs = []
for n, layer in enumerate(annotation_layers):
    # Only add cell count barplot to the last plot
    if n == len(annotation_layers) - 1:
        ax = axs[n]
        ax_pos = ax.get_position()
        bar_width = 0.06
        bar_ax = fig.add_axes([ax_pos.x1 + 0.007, ax_pos.y0, bar_width, ax_pos.height])
        bar_ax.set_facecolor('none')
        cell_count_bar_axs.append(bar_ax)
        # Get cell counts
        count_df_y = count_df_y_dict[layer]
        cell_types_reversed = all_cell_types[::-1]
        counts = [count_df_y.loc[ct, 'total_count'] if ct in count_df_y.index else 0 for ct in cell_types_reversed]
        colors_reversed = lineage_palette[::-1]
        # Create horizontal barplot
        y_positions = range(len(cell_types_reversed))
        bar_ax.barh(y_positions, counts, color=colors_reversed, edgecolor='black', linewidth=0.8, alpha=0.9)
        # Add text annotations
        for i, count in enumerate(counts):
            if count >= 1000000:
                label = f'{count/1000000:.1f}M'
            elif count >= 1000:
                label = f'{count/1000:.1f}k'
            else:
                label = str(int(count))
            offset = max_cell_count_plot1 * 0.1
            bar_ax.text(count + offset, i, label, ha='left', va='center', fontsize=8)
        # Format barplot
        bar_ax.set_ylim(-0.5, len(cell_types_reversed) - 0.5)
        bar_ax.set_xlim(0, max_cell_count_plot1 * 1.15)
        bar_ax.set_yticks([])
        bar_ax.xaxis.set_major_formatter(FuncFormatter(format_count))
        bar_ax.spines['top'].set_visible(False)
        bar_ax.spines['right'].set_visible(False)
        bar_ax.set_xlabel('# cells', fontsize=9)
        bar_ax.tick_params(axis='x', labelsize=8)

# Add colorbar for z-score
if len(annotation_layers) > 0:
    last_cell_bar_ax = cell_count_bar_axs[-1]
    cbar_height = space_per_dot_cell * 4
    cax = fig.add_axes([last_cell_bar_ax.get_position().x1 + 0.03, last_cell_bar_ax.get_position().y0, 
                        0.02, cbar_height])
    norm = plt.Normalize(-zscore_abs_max_plot1, zscore_abs_max_plot1)
    sm = plt.cm.ScalarMappable(cmap='RdBu_r', norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax, orientation='vertical')
    cax.set_title('Z-score\n(avg # cells)', ha='left', x=0)
    cax.set_ylim(zscore_min_plot1, zscore_max_plot1)
    # Add size legend
    temp_ax = fig.add_axes([0, 0, 0.1, 0.1])
    temp_scatter = sns.scatterplot(
        data=dot_plot_df_dict_zscore[annotation_layers[-1]],
        x='category',
        y='cell_type',
        size=size,
        sizes=(0, dot_size_max),
        size_norm=(size_min_plot1, size_max_plot1),
        ax=temp_ax,
        legend='brief'
    )
    size_handles, size_labels = temp_ax.get_legend_handles_labels()
    if len(size_handles) > 3:
        size_indices = [0, len(size_handles) // 2, -1]
        size_handles = [size_handles[i] for i in size_indices]
        size_labels = [size_labels[i] for i in size_indices]
    for handle in size_handles:
        handle.set_markerfacecolor('lightgray')
        handle.set_linewidth(1)
        handle.set_markeredgecolor('black')
    legend = cax.legend(
        size_handles, size_labels,
        loc='lower left',
        bbox_to_anchor=(-0.5, 1.4),
        ncol=1,
        frameon=False,
        title='Avg prop\n/ spot',
        handletextpad=0.8,
        labelspacing=0.5,
        fontsize=10,
        alignment='left'
    )
    legend.get_title().set_fontsize(10)
    legend.get_title().set_fontweight('normal')
    for txt in legend.get_texts():
        txt.set_ha("left")
    temp_ax.remove()

# Add title
title_y = bar_axs[0].get_position().y1 + 0.1
plt.suptitle(f'All samples - Cell type distribution per annotation layer (z-score, horizontal layout)', 
             x=0.5, y=title_y, fontweight='bold', fontsize=13)

# Add y-axis label
axs[0].set_ylabel('Cell type', fontweight='bold', labelpad=10)


# --- Save the figure

out_dir = os.path.join(res_dir, 'figures',)
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_1G.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')

print("Script succesfully completed")
