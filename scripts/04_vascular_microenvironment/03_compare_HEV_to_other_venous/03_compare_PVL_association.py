#!/usr/bin/env python3

# --- Aim
# Compare PVL counts and associations between HEV and non-HEV vessels.
# Quantify differences in PVL distribution and vessel overlap.

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import os
import sys
from datetime import date, datetime

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import scanpy.plotting.palettes as sc_palettes
import matplotlib.colors as mcolors

# Import custom functions
from sophiesfunctions.save_multi_image import save_multi_image
from sophiesfunctions.color_functions import generate_color_variations 

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

args = vars(parser.parse_args())

# Arguments
project_dir = args["project_dir"]
h5ad_dir = args["h5ad_dir"]
res_dir = args["res_dir"]
sample_list = args["samples"]

print(f'Samples in python: {sample_list}')

# Import function script as a module
sys.path.append(os.path.join(project_dir, 'scripts', 'TLS_model', '05_HEV_analysis', '04_cellular_composition',))
import plotting_functions


# --- Extract some info from the sample sheet

sample_dir = os.path.join(project_dir, 'config')

sample_sheet_file = 'sample_sheet.csv'

sample_sheet = pd.read_csv(os.path.join(sample_dir, sample_sheet_file))
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

# Also get PVL classes
PVL_classes = list(class_df[class_df['miniatlas_major'] == 'PVL'][f'miniatlas_{decon_level}'].unique())
PVL_classes = sorted(PVL_classes)

# Define which endothelial type want to look at
endo_of_interest = 'Endothelial ACKR1'
out_name = endo_of_interest.split(' ')[-1]

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

# Find most recent results
threshold = 0.5
threshold_adj = str(threshold).split('.')[-1]

# Specify path
TLS_group_dir = os.path.join(project_dir, 'results', 'TLS_model', 'other_signatures', 'HEV_analysis', 'per_sample', 'square_016um', 
                             f'threshold_{threshold_adj}', 'other_vessels', 'vessel_HEV_overlap', 'csv_files')

most_recent_dir_path = find_most_recent(TLS_group_dir)

PVL_counts = {}
vessel_HEV_overlap = {}

for sample in sample_list:
    filename = f'{sample}_PVL_counts_per_vessel.csv'
    vessel_classes = pd.read_csv(os.path.join(most_recent_dir_path, filename),)
    PVL_counts[sample] = vessel_classes
    filename = f'{sample}_vessel_HEV_overlap.csv'
    vessel_classes = pd.read_csv(os.path.join(most_recent_dir_path, filename), index_col=0)
    vessel_HEV_overlap[sample] = vessel_classes

# Extract a single df
HEV_overlap_df = pd.concat(
    [df.assign(sample_id=key) for key, df in vessel_HEV_overlap.items()],
    ignore_index=True
)

# Redefine HEV overlap as min >0.5 of vessel
min_overlap = 0.5
HEV_overlap_df['HEV_overlap'] = HEV_overlap_df['pct_overlap'] > min_overlap

PVL_counts_df = pd.concat(
    [df.assign(sample_id=key) for key, df in PVL_counts.items()],
    ignore_index=True
)


# --- Also load the pathology annotations
# Load pathology annotations

# Specify path
TLS_group_dir = os.path.join(project_dir, 'results', 'TLS_model', 'other_signatures', 'HEV_analysis', 'per_sample', 'square_016um', 
                             f'threshold_{threshold_adj}', 'other_vessels', 'pathology_annotations', 'csv_files')

most_recent_dir_path = find_most_recent(TLS_group_dir)

path_dict = {}

for sample in sample_list:
    filename = f'{sample}_dominant_pathology_per_vessel_group.csv'
    path_dict[sample] = pd.read_csv(os.path.join(most_recent_dir_path, filename))

# Collate all pathology annotations into a single DataFrame
pathology_df = pd.concat(
    [df.assign(sample_id=key) for key, df in path_dict.items()],
    ignore_index=True
)


# --- Define some colors

# Define some colors per endothelial type
colors_file = os.path.join(project_dir, 'resources', 'plotting_colors', 'endothelial_colors.csv')
colors_df = pd.read_csv(colors_file)

# Make a color_dict
color_dict = colors_df.set_index('endothelial_type')['color'].to_dict()

# Some colors for the pathology annotations
annotation_labels = {
    'epithelial': ['invasive', 'in_situ', 'benign', 'normal'],
    'immune': ['lymphocytes'],
    'stroma': ['stroma_fibrous', 'stroma_adipose'],
    'additional_features': ['necrosis', 'blood_vessel', 'LVI']
}

# Get some colors
all_colors = sc_palettes.default_20

# Remove red to be more colorblind-friendly
del all_colors[3] 

# Extract all pathology categories
all_cats = list(pathology_df['dominant_pathology'].dropna().unique())

top_cat_colors = {}
for n, cat in enumerate(annotation_labels):
    top_cat_colors[cat] = all_colors[n]

# Get some color variations of the original colors (so variations of the same annotation are a similar color)
color_dict_path = {}

# And also define a plot order
plot_order = []

for top_cat, labels in annotation_labels.items():
    base_color = top_cat_colors[top_cat]
    # Check how many of the labels are found in the object
    num_cols = len(labels)
    # Get the matching items
    matching_items = [item for item in labels if item in set(all_cats)]
    plot_order += matching_items
    if num_cols > 4:
        step_size=0.1
    else:
        step_size=0.2
    variations = generate_color_variations(base_color, num_cols, step_size=step_size)
    for n, item in enumerate(labels):
        color_dict_path[item] = variations[n]


# --- Plot 1: Groups with HEV vs no HEV overlap

# Enable gridlines behind the plot for all plots
plt.rcParams['axes.grid'] = True
plt.rcParams['axes.axisbelow'] = True

# Plotting
nrows = 1
ncols = 4
fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(3*ncols, 2.5))

# -- First stacked barplot showing HEV overlap per category
n=0
ax = axs[n]

# Step 1: Pivot the counts to get the count of True per group
counts_true = HEV_overlap_df.groupby(['group_type', 'HEV_overlap']).size().unstack(fill_value=0)[True]

# Step 2: Pivot the counts to get the count of False per group
counts_false = HEV_overlap_df.groupby(['group_type', 'HEV_overlap']).size().unstack(fill_value=0)[False]

# Step 3: Merge the two counts into one DataFrame
final_counts = pd.concat([counts_true, counts_false], axis=1)

#  Step 1: Create an empty DataFrame with columns corresponding to each category
transformed_df = pd.DataFrame(0, index=final_counts.index, columns=endothelial)

# Step 2: Assign True counts to corresponding columns
for category in endothelial:
    transformed_df[category] = final_counts.loc[final_counts.index.str.contains(category), True]

transformed_df['False'] = final_counts[False]
transformed_df = transformed_df.fillna(0)

color_dict['False'] = 'lightgray'

# 4. Plot
transformed_df.plot(kind='bar', stacked=True, color=color_dict, ax=ax, legend=False,
                    width=0.9)

plotting_functions.beautify_plot(ax)

plotting_functions.annotate_bars(ax, final_counts.reset_index(), x='group_type', y=True, plot_order=endothelial, margin=1, annotation_size=9)
ax.grid(axis='x', visible=False)  # Turn off vertical gridlines

ax.set_xlabel('')
ax.set_ylabel('# vessels')
ax.set_title(f'All samples\nHEV overlap', fontweight='bold')

# -- Next plot proportions
n+=1
ax = axs[n]

rowsums = transformed_df.sum(axis=1)
prop_df = transformed_df.div(rowsums, axis=0)

prop_df.plot(kind='bar', stacked=True, color=color_dict, ax=ax, legend=False,
                    width=0.9)
plotting_functions.beautify_plot(ax)
ax.grid(axis='x', visible=False)  # Turn off vertical gridlines

anno_df = prop_df[endothelial].sum().reset_index()
plotting_functions.annotate_bars(ax, anno_df, x='index', y=0, plot_order=endothelial, margin=1, annotation_size=9, round_int=False,)

ax.set_xlabel('')
ax.set_ylabel('Proportion')
ax.set_title(f'All samples\nHEV overlap', fontweight='bold')

# -- Then boxplot showing overlap size

# Boxplot size
# Properties for boxplot
whisker_props = {'color': 'black'}
median_props = {'color': 'black'}

n+=1
ax = axs[n]

sns.boxplot(HEV_overlap_df, x='group_type', y='pct_overlap', order=endothelial, palette=color_dict, ax=ax, hue='group_type',
            saturation=1, whiskerprops=whisker_props, medianprops=median_props, showfliers=False)

# Make sure the dots below min_overlap are gray
# Split the DataFrame
below = HEV_overlap_df[HEV_overlap_df['pct_overlap'] <= min_overlap]
above = HEV_overlap_df[HEV_overlap_df['pct_overlap'] > min_overlap]

# Plot gray dots (below threshold)
sns.stripplot(
    data=below, x='group_type', y='pct_overlap',
    color="gray", order=endothelial, size=2, ax=ax, dodge=False
)
# Plot black dots (above threshold)
sns.stripplot(
    data=above, x='group_type', y='pct_overlap',
    color="black", order=endothelial, size=2, ax=ax, dodge=False
)

# Beautify the plot
plotting_functions.beautify_plot(ax)

# Draw a line indicating the min_overlap threshold
ax.axhline(min_overlap, color='red', linestyle='-', linewidth=1, alpha=0.7, zorder=-1)   

ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_xlabel('')
ax.set_title(f'All samples\nHEV overlap size', fontweight='bold')

# -- Then scatterplot overlap vs size colored by type

n+=1
ax = axs[n]

sns.scatterplot(HEV_overlap_df, x='count', y='pct_overlap', palette=color_dict, ax=ax, hue='group_type', hue_order=endothelial,
                s=10, alpha=0.5, edgecolor=None)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)

ax.set_title(f'All samples\nPct overlap vs vessel size', fontweight='bold')


plt.subplots_adjust(wspace=0.3)

plt.suptitle(f'All samples - Overlap vessel groups with HEVs (min_overlap={min_overlap})', 
             y=1.2, x=0.5, fontweight='bold')


# --- Plot 2: Same as above but stacked barplots with error bars

# Calculate means per group
counts_df = HEV_overlap_df.groupby(['sample_id', 'group_type', 'HEV_overlap',]).size().unstack(fill_value=0)
counts_df = counts_df.reset_index()
counts_df.columns = counts_df.columns.astype(str)

rowsums = counts_df[['True', 'False']].sum(axis=1)
prop_df = counts_df.copy()
prop_df[['True', 'False']] = counts_df[['True', 'False']].div(rowsums, axis=0)

# Step 1: Aggregate by group_type
agg = prop_df.groupby('group_type')[['True', 'False']].agg(['mean', 'std'])

# Flatten the column MultiIndex
agg.columns = ['_'.join(col).strip() for col in agg.columns.values]
agg = agg.reset_index()

# Step 2: Extract values
x = np.arange(len(agg))
true_means = agg['True_mean']
false_means = agg['False_mean']
true_stds = agg['True_std']
false_stds = agg['False_std']

# Plotting
nrows = 1
ncols = 3
fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(3*ncols, 2))

# -- First stacked barplot showing HEV overlap per category
n=0
ax = axs[n]

# Plot True (bottom)
bar1 = ax.bar(x, true_means, yerr=true_stds, label='True', color=color_list, ecolor='black',
    error_kw=dict(lw=1, capsize=3, capthick=1), width=0.9)

# Plot False (on top of True)
bar2 = ax.bar(x, false_means, bottom=true_means, label='False', color='lightgray', width=0.9)

plotting_functions.beautify_plot(ax)

plotting_functions.annotate_bars(ax, true_means.reset_index(), x='index', y='True_mean', plot_order=true_means.index,
                                 round_int=False, offset=list(true_stds), margin=1)

# Some tweaking
ax.set_xlabel('')
ax.set_ylabel('Proportion of vessels')
ax.set_title(f'All samples\nHEV overlap', fontweight='bold')
ax.grid(axis='x', visible=False)  # Turn off vertical gridlines

ax.set_xticks(x)
ax.set_xticklabels(agg['group_type'], rotation=90)


# -- Then plot vessel composition per sample
n+=1
ax = axs[n]

counts_df = HEV_overlap_df.groupby(['sample_id', 'group_type']).size().unstack(fill_value=0)

counts_df.plot(kind='bar', stacked=True, color=color_dict, ax=ax, legend=False,
                    width=0.9)

plotting_functions.beautify_plot(ax)
ax.grid(axis='x', visible=False)  # Turn off vertical gridlines

ax.set_xlabel('')
ax.set_ylabel('Number of vessels')
ax.set_title(f'All samples\nVessel type frequencies', fontweight='bold')

# -- Now proportions
n+=1
ax = axs[n]

rowsums = counts_df.sum(axis=1)
prop_df = counts_df.div(rowsums, axis=0)

prop_df.plot(kind='bar', stacked=True, color=color_dict, ax=ax, legend=False,
                    width=0.9)

plotting_functions.beautify_plot(ax)
ax.grid(axis='x', visible=False)  # Turn off vertical gridlines

ax.set_xlabel('')
ax.set_ylabel('Proportion of vessels')
ax.set_title(f'All samples\nVessel type frequencies', fontweight='bold')

ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, handlelength=0.7)

plt.subplots_adjust(wspace=0.5)

plt.suptitle(f'All samples - Overlap vessel groups with HEVs (with error bars) (min_overlap={min_overlap})', 
             y=1.2, x=0.5, fontweight='bold')


# --- Make col colors df for heatmap

# Make an index that contains all relevant info
PVL_counts_df.columns.values[0] = 'concat_val'
PVL_counts_df['index'] = PVL_counts_df.reset_index()['concat_val'] + '-' + PVL_counts_df['sample_id']

HEV_overlap_df['index'] = HEV_overlap_df['group_type'] + '-' + HEV_overlap_df['group'] + '-' + HEV_overlap_df['sample_id']

# First make a metadata df
metadata_df = HEV_overlap_df.copy()
metadata_df['endothelial_type'] = metadata_df['group_type']
metadata_df['vessel_size'] = metadata_df['count']

metadata_df = metadata_df[['sample_id', 'endothelial_type', 'HEV_overlap', 'vessel_size', 'index']]
metadata_df = metadata_df.set_index('index')

# Also add dominant pathology annotations
pathology_df['combined_id'] = pathology_df['vessel_type'] + '-' + pathology_df['group'] + '-' + pathology_df['sample_id']
pathology_df = pathology_df.set_index('combined_id')

metadata_df['dominant_pathology'] = pathology_df['dominant_pathology']

# Generate a col colors df
col_colors = metadata_df.copy()
col_colors['endothelial_type'] = col_colors['endothelial_type'].replace(color_dict)

color_dict[True] = 'indigo'
color_dict[False] = 'lightgray'

col_colors['HEV_overlap'] = col_colors['HEV_overlap'].replace(color_dict)

# And colors for sample_id
sample_color = {}

cmap = plt.get_cmap("Set2")  # or any other cmap like 'plasma', 'tab10', etc.
hex_colors = [mcolors.to_hex(cmap(i)) for i in range(cmap.N)]

for n, sample in enumerate(sorted(sample_list)):
    sample_color[sample] = hex_colors[n]

col_colors['sample_id'] = col_colors['sample_id'].replace(sample_color)

# Pathology annotations
col_colors['dominant_pathology'] = col_colors['dominant_pathology'].replace(color_dict_path)

# Generate a cmap for the vessel sizes
# Normalize the values to [0, 1]
vmin = metadata_df['vessel_size'].min() 
vmax = metadata_df['vessel_size'].quantile(0.9)  # 99th percentile
norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

# Choose a colormap (e.g., "viridis", "coolwarm", "plasma", etc.)
cmap = plt.colormaps['Reds']

# Map the values to RGBA colors
col_colors['vessel_size'] = col_colors['vessel_size'].map(lambda x: cmap(norm(x)))

# Define 2 evenly spaced values between vmin and vmax
mid_values = np.linspace(vmin, vmax, 5)[1:-1]  # This gives 2 evenly spaced values
mid_values = np.round(mid_values).astype(int)  # Round and convert to integer

# Combine min, mid-values, and max (ensure they are integers)
legend_values = [int(vmin)] + list(mid_values) + [int(vmax)]  # Ensure all are integers
legend_values = legend_values[::-1]

# Format labels and get corresponding colors
legend_labels = [f"{val}" for val in legend_values]  # No decimal places
legend_colors = [cmap(norm(val)) for val in legend_values]

# Create legend handles (colored boxes)
legend_handles = [mpatches.Patch(facecolor=color, label=label, edgecolor='black', linewidth=0.8) 
                  for color, label in zip(legend_colors, legend_labels)]


# --- Plot 2a: Heatmap showing dominant PVL

PVL_counts_df = PVL_counts_df.set_index('index')

# Disable the grid
plt.rcParams['axes.grid'] = False

# Make heatmap
g = sns.clustermap(PVL_counts_df[PVL_classes].T, row_cluster=True, col_cluster=True, cmap='viridis', linewidths=0, figsize=(4, 4), annot=False,
            colors_ratio=(0.01, 0.1), xticklabels=True, yticklabels=True, fmt='g', row_colors=None,
            col_colors=col_colors, rasterized=True,
            dendrogram_ratio=(0.2, 0.2), z_score=1)

axs = plotting_functions.beautify_heatmap(g, row_colors=False, col_colors=True, cbar_title='z score\n(column)', label_size=10, fontsize_title=11, title='',
                    cbar_pad=10.5, cbar_height_square=len(PVL_classes), cbar_width_square=0.6, space_per_square=0.05, width_to_height_ratio=0.01,
                    colors_div=1)

# Adjust dendrograms
plotting_functions.adjust_dendrograms(g, num_square=3, title=f'All samples\nEndothelial - PVL associations')

# Add some lines separating clusters
dend_thresh=1.5
plotting_functions.cut_dendrogram_and_draw_lines(g, axis='col', threshold_dend=dend_thresh, color='white', linewidth=1, linestyle='--', alpha=1)

# Remove the xticks
heatmap_ax = g.ax_heatmap
heatmap_ax.set_xticks([])

# Add the discrete color legend
legend1 = g.ax_cbar.legend(handles=legend_handles, title="vessel_size", frameon=False, alignment='left', handlelength=0.7,
                    loc='center left', bbox_to_anchor=(4, 0.5), labelspacing=0)  # adjust position
legend1.get_title().set_fontweight('bold')

# Add legends
legend2 = plotting_functions.add_legend(ax=g.ax_cbar, ordered_list=[True, False], color_lut=color_dict, title='HEV_overlap', 
                                        y_offset=1, x_offset=4, loc='lower left')

legend3 = plotting_functions.add_legend(ax=g.ax_cbar, ordered_list=endothelial, color_lut=color_dict, title='endothelial_type', 
                                        y_offset=1.8, x_offset=4, loc='lower left')

legend4 = plotting_functions.add_legend(ax=g.ax_cbar, ordered_list=sample_list, color_lut=sample_color, title='sample_id', 
                                        y_offset=3.2, x_offset=4, loc='lower left')

legend5 = plotting_functions.add_legend(ax=g.ax_cbar, ordered_list=plot_order, color_lut=color_dict_path, title='dominant_pathology', 
                                        y_offset=0, x_offset=25, loc='lower left')

g.ax_cbar.add_artist(legend1)
g.ax_cbar.add_artist(legend2)
g.ax_cbar.add_artist(legend3)
g.ax_cbar.add_artist(legend4)

# Calculate the midpoint to position the title
top_ax = g.ax_col_dendrogram
title_x = (top_ax.get_position().x0 + top_ax.get_position().x1)/2
title_y = top_ax.get_position().y1 + 0.2

plt.suptitle(f'All samples - PVL association per vessel type (decon_level={decon_level}, z_score=col, dend_thresh={dend_thresh})', y=title_y, x=title_x, fontweight='bold')


# --- Plot 2b: Same as above but with endothelial ACKR1 only

# Filter the PVL_counts_df to only include endothelial ACKR1
# Included in the index
filtered_df = PVL_counts_df[PVL_counts_df.index.str.contains(out_name)]


# Make heatmap
g = sns.clustermap(filtered_df[PVL_classes].T, row_cluster=True, col_cluster=True, cmap='viridis', linewidths=0, figsize=(4, 4), annot=False,
            colors_ratio=(0.01, 0.1), xticklabels=True, yticklabels=True, fmt='g', row_colors=None,
            col_colors=col_colors, rasterized=True,
            dendrogram_ratio=(0.2, 0.2), z_score=1)

axs = plotting_functions.beautify_heatmap(g, row_colors=False, col_colors=True, cbar_title='z score\n(column)', label_size=10, fontsize_title=11, title='',
                    cbar_pad=10.5, cbar_height_square=len(PVL_classes), cbar_width_square=0.6, space_per_square=0.05, width_to_height_ratio=0.01,
                    colors_div=1)

# Adjust dendrograms
plotting_functions.adjust_dendrograms(g, num_square=3, title=f'All samples\nEndothelial ACKR1 - PVL associations')

# Add some lines separating clusters
dend_thresh=1.5
plotting_functions.cut_dendrogram_and_draw_lines(g, axis='col', threshold_dend=dend_thresh, color='white', linewidth=1, linestyle='--', alpha=1)

# Remove the xticks
heatmap_ax = g.ax_heatmap
heatmap_ax.set_xticks([])

# Add the discrete color legend
legend1 = g.ax_cbar.legend(handles=legend_handles, title="vessel_size", frameon=False, alignment='left', handlelength=0.7,
                    loc='center left', bbox_to_anchor=(4, 0.5), labelspacing=0)  # adjust position
legend1.get_title().set_fontweight('bold')

# Add legends
legend2 = plotting_functions.add_legend(ax=g.ax_cbar, ordered_list=[True, False], color_lut=color_dict, title='HEV_overlap', 
                                        y_offset=1, x_offset=4, loc='lower left')

legend3 = plotting_functions.add_legend(ax=g.ax_cbar, ordered_list=endothelial, color_lut=color_dict, title='endothelial_type', 
                                        y_offset=1.8, x_offset=4, loc='lower left')

legend4 = plotting_functions.add_legend(ax=g.ax_cbar, ordered_list=sample_list, color_lut=sample_color, title='sample_id', 
                                        y_offset=3.2, x_offset=4, loc='lower left')

legend5 = plotting_functions.add_legend(ax=g.ax_cbar, ordered_list=plot_order, color_lut=color_dict_path, title='dominant_pathology', 
                                        y_offset=0, x_offset=25, loc='lower left')

g.ax_cbar.add_artist(legend1)
g.ax_cbar.add_artist(legend2)
g.ax_cbar.add_artist(legend3)
g.ax_cbar.add_artist(legend4)

# Calculate the midpoint to position the title
top_ax = g.ax_col_dendrogram
title_x = (top_ax.get_position().x0 + top_ax.get_position().x1)/2
title_y = top_ax.get_position().y1 + 0.2

plt.suptitle(f'All samples - PVL association per vessel type (decon_level={decon_level}, z_score=col, dend_thresh={dend_thresh})', y=title_y, x=title_x, fontweight='bold')


# --- Plot 2c: Heatmap showing dominant PVL - HEV overlapping vessels only

overlap_with_HEV = list(metadata_df[metadata_df['HEV_overlap']].index)
PVL_sub = PVL_counts_df.loc[overlap_with_HEV,:]

# Make heatmap
g = sns.clustermap(PVL_sub[PVL_classes].T, row_cluster=True, col_cluster=True, cmap='viridis', linewidths=0, figsize=(4, 4), annot=False,
            colors_ratio=(0.01, 0.1), xticklabels=True, yticklabels=True, fmt='g', row_colors=None,
            col_colors=col_colors, rasterized=True,
            dendrogram_ratio=(0.2, 0.2), z_score=1)

axs = plotting_functions.beautify_heatmap(g, row_colors=False, col_colors=True, cbar_title='z score\n(column)', label_size=10, fontsize_title=11, title='',
                    cbar_pad=10.5, cbar_height_square=len(PVL_classes), cbar_width_square=0.6, space_per_square=0.05, width_to_height_ratio=0.1,
                    colors_div=1)

# Adjust dendrograms
plotting_functions.adjust_dendrograms(g, num_square=3, title=f'All samples\nEndothelial - PVL associations\nVessels overlapping with HEVs')

# Add some lines separating clusters
dend_thresh=1.5
plotting_functions.cut_dendrogram_and_draw_lines(g, axis='col', threshold_dend=dend_thresh, color='white', linewidth=1, linestyle='--', alpha=1)

# Remove the xticks
heatmap_ax = g.ax_heatmap
heatmap_ax.set_xticks([])

# Add the discrete color legend
legend1 = g.ax_cbar.legend(handles=legend_handles, title="vessel_size", frameon=False, alignment='left', handlelength=0.7,
                    loc='center left', bbox_to_anchor=(4, 0.5), labelspacing=0)  # adjust position
legend1.get_title().set_fontweight('bold')

# Add legends
legend2 = plotting_functions.add_legend(ax=g.ax_cbar, ordered_list=[True, False], color_lut=color_dict, title='HEV_overlap', 
                                        y_offset=1, x_offset=4, loc='lower left')

legend3 = plotting_functions.add_legend(ax=g.ax_cbar, ordered_list=endothelial, color_lut=color_dict, title='endothelial_type', 
                                        y_offset=1.8, x_offset=4, loc='lower left')

legend4 = plotting_functions.add_legend(ax=g.ax_cbar, ordered_list=sample_list, color_lut=sample_color, title='sample_id', 
                                        y_offset=3.2, x_offset=4, loc='lower left')

legend5 = plotting_functions.add_legend(ax=g.ax_cbar, ordered_list=plot_order, color_lut=color_dict_path, title='dominant_pathology',
                                        y_offset=0, x_offset=25, loc='lower left')              

g.ax_cbar.add_artist(legend1)
g.ax_cbar.add_artist(legend2)
g.ax_cbar.add_artist(legend3)
g.ax_cbar.add_artist(legend4)

# Calculate the midpoint to position the title
top_ax = g.ax_col_dendrogram
title_x = (top_ax.get_position().x0 + top_ax.get_position().x1)/2
title_y = top_ax.get_position().y1 + 0.2

plt.suptitle(f'All samples - PVL association per vessel type - vessels with HEV overlap only (decon_level={decon_level}, z_score=col, dend_thresh={dend_thresh})', y=title_y, x=title_x, fontweight='bold')


# --- Save the output

qc_dir = os.path.join(project_dir, 'results', 'TLS_model')

today = date.today()
today = today.strftime("%Y%m%d")

out_dir = os.path.join(qc_dir, 'other_signatures', 'HEV_analysis', 'all_samples', 'square_016um', f'threshold_{threshold_adj}', 'other_vessels',
                       'HEV_vs_non_HEV', out_name, 'PVL_association')

os.makedirs(out_dir, exist_ok=True)

# Save plots
out_dir_plots = os.path.join(out_dir, 'plots', today)
os.makedirs(out_dir_plots, exist_ok=True)

# Save plots
filename = 'All_samples_vessel_HEV_overlap_PVL.pdf'
save_multi_image(os.path.join(out_dir_plots, filename))


# --- Finish

print("Script successfully completed")
