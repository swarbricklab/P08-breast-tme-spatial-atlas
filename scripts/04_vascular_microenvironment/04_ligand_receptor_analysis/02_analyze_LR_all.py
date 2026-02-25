#!/usr/bin/env python3

# --- Aim
# Analyze stlearn-derived ligand-receptor interaction results across all samples.
# Identify consistently detected LR pairs and their association with HEV groups.

# --- Define environment
# env: shapely (has upsetplot)

# --- Load packages
import argparse
import math
import os
import re
import sys
from collections.abc import Iterable
from datetime import date, datetime

import anndata as ad
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import squidpy as sq
from adjustText import adjust_text
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle

# Import custom functions
from sophiesfunctions.auto_crop import auto_crop_scanpy
from sophiesfunctions.save_multi_image import save_multi_image
from upsetplot import UpSet

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

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

# Import function script as a module
sys.path.append(os.path.join(project_dir, 'scripts', '04_vascular_microenvironment'))
import plotting_functions


# --- Get some info from the sample sheet
# Use --sample_sheet if provided, otherwise default to config/sample_sheet.csv
if args["sample_sheet"] is not None:
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet.csv')

sample_sheet = pd.read_csv(sample_sheet_path)
sample_lut = sample_sheet.set_index('section_id')


# --- Load the adata object with LR results

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
resolution = 'square_016um'
bin_size = 25  # um

lr_dir = os.path.join(project_dir, 'results', '04_vascular_microenvironment', '04_ligand_receptor_analysis', '01_run_LR_analysis_stlearn', f'{bin_size}um_bin', 
                        'csv_files')

most_recent_dir_path = find_most_recent(lr_dir)

adata_dict = {}

for sample in sample_list:
    # Load the anndata object
    filename = f'{sample}_HEV_LR_analysis_stlearn.h5ad'
    adata_dict[sample] = ad.read_h5ad(os.path.join(most_recent_dir_path, filename))


# --- Calculate crop coords for optimal presentation

coords = {}
for sample in adata_dict.keys():
    coords[sample] = auto_crop_scanpy(adata_dict[sample], border_size=0.1)


# --- Also load HEV groups for each sample

# Load the expanded groups
TLS_group_dir = os.path.join(res_dir, '04_vascular_microenvironment', '01_identify_HEV_niches', '02b_expand_HEV_groups', resolution, f'threshold_{threshold_adj}',)
TLS_group_dir = os.path.join(TLS_group_dir, 'csv_files')

most_recent_dir_path = find_most_recent(TLS_group_dir)

# Read in the results & add to anndata
score_of_i = 'HEVs'
group_key_source = 'expanded_groups'

for sample in sample_list:
    # Load the anndata object
    filename = f'{sample}_{score_of_i}_groups.csv'
    TLS_results = pd.read_csv(os.path.join(most_recent_dir_path, filename), index_col=0)
    # Add HEV groups to the adata object
    adata_dict[sample].obs['HEV_group'] = TLS_results[group_key_source]


# --- Extract significant LRs across samples

lr_set = set()

# Collect all LR pairs across samples
for sample in sample_list:
    lr_summary = adata_dict[sample].uns['lr_summary']
    lr_set.update(lr_summary.index)

lr_list = sorted(lr_set)

# Build DataFrame: rows=LR, columns=samples
lr_counts = pd.DataFrame(index=lr_list, columns=sample_list)

for sample in sample_list:
    lr_summary = adata_dict[sample].uns['lr_summary']
    # If 'n_spots_sig' is a column and index is LR
    lr_counts.loc[lr_summary.index, sample] = lr_summary['n_spots_sig']

# Convert to numeric
lr_counts = lr_counts.apply(pd.to_numeric, errors='coerce')

# Count in how many samples each LR is significant (n_spots_sig > 0)
lr_counts['n_samples_sig'] = (lr_counts[sample_list] > 0).sum(axis=1)

# Reorder columns
lr_counts = lr_counts[['n_samples_sig'] + sample_list]

# Optional: reset index to have LR as a column
lr_counts = lr_counts.reset_index().rename(columns={'index': 'LR'})

# Add total spots significant across samples
lr_counts['total_bins_sig'] = lr_counts[sample_list].sum(axis=1)

# Sort df by total_bins_sig
lr_counts = lr_counts.sort_values(by='total_bins_sig', ascending=False)

sample_color = {}
cmap = plt.get_cmap("Set2")  # or any other cmap like 'plasma', 'tab10', etc.


# --- Prepare a color dict for samples

sample_color = {}
cmap = plt.get_cmap("Set2")  # or any other cmap like 'plasma', 'tab10', etc.
hex_colors = [mcolors.to_hex(cmap(i)) for i in range(cmap.N)]

# Make sure the colors are assigned in the same order as the samples appear in the sample list
plot_labels = sorted(sample_list)
plot_samples = plot_labels
for n, sample in enumerate(plot_samples):
    sample_color[sample] = hex_colors[n]


# --- Also prepare a dict to rename the labels in the plots

rename_dict = {
    'near_TLS_HEV': 'HEV (TLS)',
    'non_TLS_HEV': 'HEV-like (stroma)',
    'nearby_spots': f'Perivascular bins (<{bin_size}µm)',
}


# --- Plot 1: Upsetplot showing overlap in singlificant LRs across samples

sc.set_figure_params(fontsize=10)

upsetplot_labels = sorted(sample_list)
sample_list_ordered = upsetplot_labels

# Prepare df for upsetplot
# Create a boolean DataFrame for significance, using the ordered sample list
bool_df = lr_counts[sample_list_ordered].notna() & (lr_counts[sample_list_ordered] > 0)

# Set MultiIndex of booleans for each LR, with ordered labels
multi_index = pd.MultiIndex.from_frame(bool_df, names=upsetplot_labels)

# Build the DataFrame for upsetplot: index=MultiIndex, columns=total_bins_sig (and others if needed)
upset_df = pd.DataFrame({
    'total_bins_sig': lr_counts['total_bins_sig'].values,
    # add more columns here if needed
}, index=multi_index)

# Make the upsetplot
upset = UpSet(upset_df, show_counts=True)

# Add plot showing total spots significant per LR
upset.add_catplot(value="total_bins_sig", kind="strip", color="black", size=2)

upset.plot()
plt.suptitle(f'All samples - stlearn significant LR interactions - across samples', y=1, x=0.5, fontweight='bold')


# --- Plot 2a: Barplot showing number of spots significant per LR across samples

top = 50

ncols = 1
nrows = 1

# Subset to top 50 LRs & prepare the df for plotting
plot_df = lr_counts.head(top)
plot_df.set_index('LR', inplace=True)
plot_df = plot_df[plot_samples[::-1]]

fig, ax = plt.subplots(ncols=ncols, nrows=nrows, figsize=(10*ncols, 3*nrows))

plot_df.plot(kind='bar', stacked=True, color=sample_color, width=0.8, ax=ax)

plotting_functions.beautify_plot(ax)

# Set grid behind plot
ax.set_axisbelow(True)
ax.grid(axis='x', visible=False)  # Turn off vertical gridlines

# Some tweaking
ax.set_xlabel('')
ax.set_ylabel('Number of spots')
ax.set_title(f'Significant LR interactions\nHEVs + spots <{bin_size}µm', fontweight='bold')

# Reverse the order of the legend
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], loc='best', handlelength=0.8,
              bbox_to_anchor=(1, 1), fontsize='small', title_fontsize='medium')

plt.suptitle(f'All samples - stlearn significant LR interactions top {top} by number of spots', y=1.2, x=0.5, fontweight='bold')


# --- Plot 2b: Same as above but filtering collagens for now

# Filter out collagens bc they tend to be very high in number of spots
collagen_pattern = r'COLLAG|COL1A|COL2A|COL'
plot_df_filtered = lr_counts[~lr_counts.set_index('LR').index.str.contains(collagen_pattern, regex=True)]
plot_df_filtered.set_index('LR', inplace=True)

plot_df_filtered = plot_df_filtered[plot_samples[::-1]]
plot_df_filtered = plot_df_filtered.head(top)

fig, ax = plt.subplots(ncols=ncols, nrows=nrows, figsize=(10*ncols, 3*nrows))

plot_df_filtered.plot(kind='bar', stacked=True, color=sample_color, width=0.8, ax=ax)

plotting_functions.beautify_plot(ax)

# Set grid behind plot
ax.set_axisbelow(True)
ax.grid(axis='x', visible=False)  # Turn off vertical gridlines

# Some tweaking
ax.set_xlabel('')
ax.set_ylabel('Number of spots')
ax.set_title(f'Significant LR interactions\nHEVs + spots <{bin_size}µm\n(Collagens filtered)', fontweight='bold')

# Reverse the order of the legend
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], loc='best', handlelength=0.8,
              bbox_to_anchor=(1, 1), fontsize='small', title_fontsize='medium')

plt.suptitle(f'All samples - stlearn significant LR interactions top {top} by number of spots - collagens filtered out', y=1.2, x=0.5, fontweight='bold')


# --- Plot 3: For sanity, spatial plots showing cell typing info used in LR analysis

# Get a colordict to make sure all colors align across samples
colors = adata_dict[sample_list[0]].uns['cell_type_colors']
labels = adata_dict[sample_list[0]].obs['cell_type'].cat.categories

color_dict = {}
for label, color in zip(labels, colors):
    color_dict[label] = color

# Plotting
sc.set_figure_params(fontsize=18)

if len(adata_dict.keys()) < 5:
        ncols = len(adata_dict.keys())
else:
    ncols = 5

nrows = math.ceil((len(adata_dict.keys()))/ncols)

plt.figure(figsize=(8*ncols, 7*nrows))

for n, sample in enumerate(plot_samples):
    # Only add legend to the last plot on the first row
    if n == ncols - 1:
        legend_loc = 'right margin'
    else:
        legend_loc = None
    # add a new subplot iteratively
    ax = plt.subplot(nrows, ncols, n + 1)
    display_title = sample
    sc.pl.spatial(adata_dict[sample], color='cell_type', bw=False, alpha_img=0.5,
            title=f'{display_title}\nGrouping info', crop_coord=coords[sample], ax=ax, legend_loc=legend_loc,
            palette=color_dict,)
    ax.set_title(ax.get_title(), fontweight='bold')

plt.suptitle(f'All samples - stlearn cell typing info (HEV_threshold={threshold}, bin_size={bin_size}µm)', y=1.05, x=0.5, fontweight='bold')


# --- Redefine groups also including the nearby spots
# Then afterwards will asign near TLS vs non TLS based on dominant annotation

# squidpy's spatial_connectivities just saves which spots are neighbors
# so there will be overlap between the sublists (if a spot is a neighbor of multiple other spots, it appears multiple times)
# therefore merge the sublists until all values are unique

# Define a function to do this
def extract_neighbor_groups(connected_groups, adata, obs_key='group'):
    pooled = [set(subList) for subList in connected_groups]
    merging = True
    while merging:
        merging=False
        for i,group in enumerate(pooled):
            merged = next((g for g in pooled[i+1:] if g.intersection(group)),None)
            if not merged: continue
            group.update(merged)
            pooled.remove(merged)
            merging = True
    # Convert the list of sets to a list of lists
    list_of_lists = [list(s) for s in pooled]
    # Generate group names
    groups = [f'group{i+1}' for i in range(len(list_of_lists))]
    # Flatten the list of lists and pair each item with its corresponding group
    data_list = []
    for group, sublist in zip(groups, list_of_lists):
        for item in sublist:
            data_list.append({'index': item, 'group': group})
    # Create a DataFrame
    df = pd.DataFrame(data_list)
    # Sort DataFrame based on the 'index' column
    df_sorted = df.sort_values(by='index')
    # Transfer the barcode id to the index from the anndata object
    df_sorted.index = adata.obs_names
    # Tack the group annotation onto the anndata object
    adata.obs[obs_key] = df_sorted['group']
    # Return the grouping
    return list_of_lists

# Do spatial gropuping per sample
for sample in sample_list:
    adata = adata_dict[sample]
    adata_sub = adata[adata.obs['cell_type'].notna()]
    # Get a scalefactor for converting pixel coords to um
    spot_diameter_pixel = adata.uns['spatial'][sample]['scalefactors']['spot_diameter_fullres']
    resolution_numbers = re.findall(r'\d+', resolution)[0]
    resolution_numbers = int(resolution_numbers)
    spot_diameter_um = resolution_numbers # Square size is whatever the resolution is 
    scalefactor_um = spot_diameter_pixel/spot_diameter_um
    # Vis hd spots have no distance between them, so distance from center to center is the square size
    spot_dist_um = resolution_numbers
    spot_dist_pixel = spot_dist_um * scalefactor_um
    # But Vis hd spots are squares, I also want to catch the diagonal neighbors
    # So need to adjust the radius for this (will not catch these neighbors if I use the length of the square's edge)
    # Lucky I remember how to do pythagoras
    c = pow(spot_dist_pixel/2, 2) + pow(spot_dist_pixel/2, 2)
    c = math.sqrt(c)
    radius_max = c*2
    # Round up to whole pixel
    radius_max = math.ceil(radius_max)
    # In theory this should work, but for some reason it doesn't in some areas??
    # Make it a bit larger and then it works on all spots
    radius_max = radius_max * 1.1
    # Use squidpy to define which TLS-high spots are neighbors
    sq.gr.spatial_neighbors(adata_sub, coord_type='generic', radius=(0, radius_max))
    # Extract connected spots as groups
    connected_groups = []
    for i, neighbors in enumerate(adata_sub.obsp['spatial_connectivities'].toarray()):
        group = [j for j, is_neighbor in enumerate(neighbors) if is_neighbor]
        group.append(i)  # Include the spot itself in the group
        connected_groups.append(group)
    # Extract groups
    list_of_lists = extract_neighbor_groups(connected_groups, adata=adata_sub, obs_key='group')
    # Name groups according to size
    group_counts = adata_sub.obs['group'].value_counts()
    group_counts = group_counts[group_counts > 1]  # Only keep groups with more than
    # one spot
    group_names = {group: f'Group {i+1}' for i, group in enumerate(group_counts.index)}
    # Assign group names to the adata object
    adata_sub.obs['group'] = adata_sub.obs['group'].map(group_names)
    # Add to original adata
    adata_dict[sample].obs['new_grouping'] = adata_sub.obs['group']


# --- Plot 4: Spatial plot showing the new groupings
# Check that nothing weird happened

plt.figure(figsize=(8*ncols, 7*nrows))

for n, sample in enumerate(plot_samples):
    # add a new subplot iteratively
    ax = plt.subplot(nrows, ncols, n + 1)
    display_title = sample
    sc.pl.spatial(adata_dict[sample], color='new_grouping', bw=False, alpha_img=0.5,
        title=f'{display_title}\nGrouping info', crop_coord=coords[sample], ax=ax, legend_loc=None,
        palette='Set1',)
    ax.set_title(ax.get_title(), fontweight='bold')

plt.suptitle(f'All samples - grouping of stlearn cell typing info (HEV_threshold={threshold}, bin_size={bin_size}µm)', y=1.05, x=0.5, fontweight='bold')


# --- Assign categories to groups based on dominant HEV label

# Remove the nearby_spots label
HEV_labels = [l for l in labels if l != 'nearby_spots']

# Define a function to assign categories based on dominant HEV label
def assign_HEV_category(group, adata=adata_sub, obs_key='cell_type'):
    # Get the HEV labels for the spots in the group
    adata_sub = adata[adata.obs['new_grouping'] == group]
    # Count the occurrences of each HEV label in the group
    hev_counts = adata_sub.obs[obs_key].value_counts()
    # Figure out which of the HEV labels is dominant
    hev_counts = hev_counts[hev_counts.index.isin(HEV_labels)]
    # If there are no HEV labels in the group, return 'unknown'
    if hev_counts.empty:
        dominant_label = 'unknown'
    else:
        dominant_label = hev_counts.idxmax()
    return dominant_label

# Apply the function to each group in the adata object
for sample in sample_list:
    adata_sub = adata_dict[sample]
    # Assign dominant HEV label to each group
    adata_sub.obs['dominant_HEV_label'] = adata_sub.obs['new_grouping'].apply(assign_HEV_category, adata=adata_sub, obs_key='cell_type')
    # Add the dominant HEV label to the adata object
    adata_dict[sample].obs['dominant_HEV_label'] = adata_sub.obs['dominant_HEV_label']


# --- Plot 5: Spatial plot showing the dominant HEV label for each group

plt.figure(figsize=(8*ncols, 7*nrows))

for n, sample in enumerate(plot_samples):
    # Only add legend to the last plot on the first row
    if n == ncols - 1:
        legend_loc = 'right margin'
    else:
        legend_loc = None
    # add a new subplot iteratively
    ax = plt.subplot(nrows, ncols, n + 1)
    display_title = sample
    sc.pl.spatial(adata_dict[sample], color='dominant_HEV_label', bw=False, alpha_img=0.5,
            title=f'{display_title}\nDominant label', crop_coord=coords[sample], ax=ax, legend_loc=legend_loc,
            palette=color_dict,)
    ax.set_title(ax.get_title(), fontweight='bold')

plt.suptitle(f'All samples - grouping of stlearn cell typing info: dominant HEV label (HEV_threshold={threshold}, bin_size={bin_size}µm)', y=1.05, x=0.5, fontweight='bold')


# --- Now for each shared LR across samples, figure out how many groups have that LR

# Extract all significant lrs as list
lr_list = list(lr_counts['LR'])
len(lr_list)

# Initialize a list to collect results
spot_count_results = []
group_count_results = []

for lr in lr_list:
    print(f'Analyzing LR: {lr}')
    for sample in sample_list:
        adata = adata_dict[sample]
        # Check if the LR is in the adata object
        if 'lr_summary' not in adata.uns or lr not in adata.uns['lr_summary'].index:
            print(f'LR {lr} not found in {sample}. Skipping...')
            continue
        # Get the list of all LR pair names
        lr_pairs = adata.uns['lr_summary']
        # Find the index of the LR pair of interest
        if isinstance(lr_pairs, pd.DataFrame):
            lr_index = lr_pairs.index.get_loc(lr)
        else:
            if isinstance(lr_pairs, Iterable):
                lr_index = lr_pairs.index(lr)
            else:
                raise ValueError("lr_pairs should be either a DataFrame or a list of LR pairs.")
        # Get adjusted p-values and scores
        p_adj_matrix = adata.obsm['p_adjs'][:, lr_index]
        sig_mask = p_adj_matrix < 0.05
        barcodes_sig = adata.obs_names[sig_mask]
        if len(barcodes_sig) == 0:
            continue  # No significant spots found
        # Subset to significant barcodes
        adata_sub = adata[adata.obs_names.isin(barcodes_sig)]
        # Count dominant_HEV_label values
        label_counts = (
            adata_sub.obs['dominant_HEV_label']
            .value_counts()
            .rename_axis('dominant_HEV_label')
            .reset_index(name='count')
        )
        # Add identifying columns
        label_counts['lr_pair'] = lr
        label_counts['sample'] = sample
        # Store for concatenation
        spot_count_results.append(label_counts)
        # Do the same on group level
        group_counts = (
            adata_sub.obs[['new_grouping', 'dominant_HEV_label']]
            .value_counts()
            .reset_index(name='count')
            ['dominant_HEV_label'].value_counts()
            .reset_index(name='count')
        )
        # Add identifying columns
        group_counts['lr_pair'] = lr
        group_counts['sample'] = sample
        # Store for concatenation
        group_count_results.append(group_counts)

# Combine all results into one DataFrame
result_df_group = pd.concat(group_count_results)
result_df_spot = pd.concat(spot_count_results)

# Group by LR pair and dominant HEV label, summing counts
totals_group = result_df_group.groupby(['lr_pair', 'dominant_HEV_label'])['count'].sum().reset_index()
totals_spot = result_df_spot.groupby(['lr_pair', 'dominant_HEV_label'])['count'].sum().reset_index()


# --- Plot 6: Barplot of top LRs showing dominant HEV labels across samples

# Make font smaller
sc.set_figure_params(fontsize=10)

ncols = 1
nrows = 1

# Subset to top 50 LRs & prepare the df for plotting
top = 50
plot_df = lr_counts.head(top)
plot_df.set_index('LR', inplace=True)
plot_df = plot_df[sample_list]
lrs = plot_df.index.tolist()

# Subset the pivot_df to the top LRs
totals_top = totals_spot.set_index('lr_pair').loc[lrs]

# Format for plotting
plot_df = totals_top.pivot_table(
    index='lr_pair',
    columns='dominant_HEV_label',
    values='count',
    fill_value=0
)

plot_df = plot_df.loc[lrs]

# Plotting
fig, ax = plt.subplots(ncols=ncols, nrows=nrows, figsize=(10*ncols, 3*nrows))

plot_df.plot(kind='bar', stacked=True, color=color_dict, width=0.8, ax=ax)

plotting_functions.beautify_plot(ax)

# Set grid behind plot
ax.set_axisbelow(True)
ax.grid(axis='x', visible=False)  # Turn off vertical gridlines

# Some tweaking
ax.set_xlabel('')
ax.set_ylabel('Number of spots')
ax.set_title(f'Significant LR interactions\nHEVs + spots <{bin_size}µm', fontweight='bold')

plt.suptitle(f'All samples - stlearn significant LR interactions top {top} by number of spots', y=1.2, x=0.5, fontweight='bold')


# --- Calculate LFCs between groups

plot_df = totals_spot.pivot_table(
    index='lr_pair',
    columns='dominant_HEV_label',
    values='count',
    fill_value=0
)

plot_df['LFC'] = plot_df.apply(
    lambda row: np.log2(row.get(HEV_labels[0], 1) / row.get(HEV_labels[1], 1)),
    axis=1
)   


# --- Plot 7a: Barplot of top LRs showing dominant HEV labels across samples, sorted by LFC

# Subset to only LRs that are shared across samples
shared_lrs = lr_counts[lr_counts['n_samples_sig'] == len(sample_list)]['LR'].tolist()
plot_df = plot_df.loc[shared_lrs]

nrows = 2

fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(10*ncols, 2*nrows))

# Sort by LFC (first descending)
ax = axs[0]
plot_df = plot_df.sort_values(by='LFC', ascending=False)
plot_df_sub = plot_df.head(top)

plot_df_sub[HEV_labels].plot(kind='bar', stacked=True, color=color_dict, width=0.8, ax=ax)

# Then plot in ascending order
ax = axs[1]
plot_df = plot_df.sort_values(by='LFC', ascending=True)
plot_df_sub = plot_df.head(top)

plot_df_sub[HEV_labels].plot(kind='bar', stacked=True, color=color_dict, width=0.8, ax=ax)

# Adjustments
for n, ax in enumerate(axs):
    plotting_functions.beautify_plot(ax)
    # Set grid behind plot
    ax.set_axisbelow(True)
    ax.grid(axis='x', visible=False)  # Turn off vertical gridlines
    # Some tweaking
    ax.set_xlabel('')
    ax.set_ylabel('Number of spots')
    ax.set_title(f'Significant LR interactions\nEnriched in {HEV_labels[n]}', fontweight='bold')
    # Move the legend
    ax.legend(handlelength=0.8, loc='center left', bbox_to_anchor=(1, 0.5), alignment='left', frameon=False, markerscale=1.5)


# Move the second plot down to accommodate labels
ax_pos = ax.get_position()
ax.set_position([ax_pos.x0, ax_pos.y0 - 0.45, ax_pos.width, ax_pos.height])

plt.suptitle(f'All samples - stlearn significant LR interactions top {top} by LFC HEV label (only LRs shared across samples)', y=1.1, x=0.5, fontweight='bold')


# --- Plot 7b: Scatterplot showing number of spots

def adjust_scatterplot(ax, plot_df, max_lfc, labels='Number of spots'):
    ax.set_aspect('equal', adjustable='box')
    # Draw a diagonal line
    max_value = max(plot_df[HEV_labels[0]].max(), plot_df[HEV_labels[1]].max())
    max_value = max_value * 1.05  # Add a bit of padding
    ax.plot([0, max_value], [0, max_value], color='gray', linestyle='--', zorder=1, linewidth=1)
    # Reset the axis limits
    ax.set_xlim(0, max_value)
    ax.set_ylim(0, max_value)
    # Set grid behind plot
    ax.set_axisbelow(True)
    # Add colorbar with LFC values
    # Add a new axes for the colorbar — [left, bottom, width, height]
    # Tweak values as needed to fit your layout
    norm = plt.Normalize(-max_lfc, max_lfc)
    sm = plt.cm.ScalarMappable(cmap='coolwarm', norm=norm)
    sm.set_array([])
    cbar_ax = fig.add_axes([0.91, 0.15, 0.02, 0.7])  # Right side of the plot
    # Disable grid for the colorbar axis
    cbar_ax.grid(False)
    cbar = plt.colorbar(sm, cax=cbar_ax)
    cbar.set_label('Log2 Fold Change', rotation=270, labelpad=15)
    # Adjust labels
    ax.set_xlabel(f'{labels}\n{HEV_labels[1]}')
    ax.set_ylabel(f'{labels}\n{HEV_labels[0]}')
    cbar_ax.set_ylim(plot_df['LFC'].min(), plot_df['LFC'].max())
    # Set position
    main_ax_pos = ax.get_position()
    cbar_ax.set_position([main_ax_pos.x1 + main_ax_pos.width*0.02, main_ax_pos.y0, main_ax_pos.width*0.05, main_ax_pos.height])
    # Add title
    ax.set_title(f'Significant LR interactions\nHEVs + spots <{bin_size}µm', fontweight='bold')

def adjust_text_func(ax, sorted_df, min_distance=20, labels_to_annotate=None, lw=0.2, color='gray', fontsize=5):
    texts = []
    annotated_coords = []
    for i, row in sorted_df.iterrows():
        label = row.name
        # If a list of labels is provided, only annotate those
        if labels_to_annotate is not None and label not in labels_to_annotate:
            continue
        x = row[HEV_labels[1]]
        y = row[HEV_labels[0]]
        # Check distance to all previous annotated points
        too_close = False
        for x0, y0 in annotated_coords:
            dist = np.sqrt((x - x0)**2 + (y - y0)**2)
            if dist < min_distance:
                too_close = True
                break
        if not too_close:
            text = ax.text(x, y, label, fontsize=fontsize, ha='right', va='bottom')
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
            arrowprops=dict(arrowstyle="-", color=color, lw=lw, shrinkA=5, shrinkB=0),
            expand_points=(1.2, 1.4),
            expand_text=(1.2, 1.4),
            force_text=(0.5, 0.5),
            force_points=(0, 0),
            lim=1000,
            only_move={'points':'', 'text':'xy', 'objects':''}
        )

def draw_rectangles(ax, xlim, ylim):
    # Color the area below 0 lightgray
    # 1. x < 0, all y
    left_rect = Rectangle(
        (xlim[0], ylim[0]),
        width=0 - xlim[0],
        height=ylim[1] - ylim[0],
        facecolor='lightgray',
        alpha=0.5,
        zorder=0
    )
    # 2. x ≥ 0, y < 0
    bottom_right_rect = Rectangle(
        (0, ylim[0]),
        width=xlim[1] - 0,
        height=0 - ylim[0],
        facecolor='lightgray',
        alpha=0.5,
        zorder=0
    )
    # 3. x < 0, y ≥ 0 — already covered by left_rect (overlaps), so skip
    # Add patches to the axis
    ax.add_patch(left_rect)
    ax.add_patch(bottom_right_rect)

# Plotting
nrows = 1
ncols = 3

# Sort by absolute LFC, descending → prioritize extreme values
sorted_df = plot_df.reindex(plot_df['LFC'].abs().sort_values(ascending=False).index)

# Find the most extreme LFC (positive or negative)
max_lfc = max(abs(plot_df['LFC'].min()), abs(plot_df['LFC'].max()))

# Define a midpoint for the TwoSlopeNorm
midpoint = 0
norm = TwoSlopeNorm(vmin=-max_lfc, vcenter=midpoint, vmax=max_lfc)

fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(3*ncols, 3*nrows))

# Increase the distance between the plots to avoid overlap
spacing = 0.6  # Fraction of width to use as spacing between plots

for i, ax in enumerate(axs):
    pos = ax.get_position()
    # Shift each axis to the right by (i * spacing * width)
    new_x0 = pos.x0 + i * pos.width * spacing
    ax.set_position([new_x0, pos.y0, pos.width, pos.height])

# First plot
ax = axs[0]

# Make a scatter plot of number of spots per HEV label for each LR
sns.scatterplot(data=sorted_df.iloc[::-1], x=HEV_labels[1], y=HEV_labels[0], hue='LFC', palette='coolwarm', ax=ax, legend=False,
                 s=20, linewidth=0.2, hue_norm=norm, edgecolor='gray')

adjust_scatterplot(ax, sorted_df, max_lfc, labels='Number of spots')

# Add a second plot with same data but smaller spots & annotations
ax = axs[1]

sns.scatterplot(data=sorted_df.iloc[::-1], x=HEV_labels[1], y=HEV_labels[0], hue='LFC', palette='coolwarm', ax=ax, legend=False,
                 s=10, linewidth=0.2, hue_norm=norm, edgecolor='gray')

adjust_scatterplot(ax, plot_df, max_lfc, labels='Number of spots')

# Adjust the plot limits to give extra space for annotations
ax.set_xlim(-ax.get_xlim()[1]*0.1, ax.get_xlim()[1])
ax.set_ylim(-ax.get_xlim()[1]*0.1, ax.get_ylim()[1])

# Get current axis limits
xlim = ax.get_xlim()
ylim = ax.get_ylim()

draw_rectangles(ax, xlim, ylim)

# Annotate the points with LR names
min_distance = 20  # Adjust based on your scale
if bin_size == 25:
    min_distance = 100  # Adjust based on your scale

adjust_text_func(ax, sorted_df, min_distance=min_distance)

# Add a third plot annotating just some interactions of interest
ax = axs[2]

lr_of_interest = ['CCL19_CCR7', 'CCL21_CCR7', 'SELPLG_SELL', 'SPP1_CD44', 'CD34_SELL', 'CFH_SELL', 'CXCL14_CXCR4', 'LTB_LTBR']

sns.scatterplot(data=sorted_df.iloc[::-1], x=HEV_labels[1], y=HEV_labels[0], hue='LFC', palette='coolwarm', ax=ax, legend=False,
                 s=10, linewidth=0.2, hue_norm=norm, edgecolor='gray')

adjust_scatterplot(ax, sorted_df, max_lfc, labels='Number of groups')


# Adjust the plot limits to give extra space for annotations
ax.set_xlim(-ax.get_xlim()[1]*0.1, ax.get_xlim()[1])
ax.set_ylim(-ax.get_xlim()[1]*0.1, ax.get_ylim()[1])

# Get current axis limits
xlim = ax.get_xlim()
ylim = ax.get_ylim()

draw_rectangles(ax, xlim, ylim)

adjust_text_func(ax, sorted_df, min_distance=0, labels_to_annotate=lr_of_interest, lw=0.5, color='black')

# Edit the labels of the x & y axis
for ax in axs:
    xlabel = ax.get_xlabel()
    ylabel = ax.get_ylabel()
    for hev_label in HEV_labels:
        if hev_label in xlabel:
            new_xlabel = xlabel.replace(hev_label, rename_dict.get(hev_label, hev_label))
            ax.set_xlabel(new_xlabel)
        if hev_label in ylabel:
            new_ylabel = ylabel.replace(hev_label, rename_dict.get(hev_label, hev_label))
            ax.set_ylabel(new_ylabel)

# Calculate the midpoint to position the title
title_x = (axs[0].get_position().x0 + axs[-1].get_position().x1)/2

plt.suptitle(f'All samples - stlearn significant LR interactions top {top} by LFC HEV label (only LRs shared across samples)', y=1.1, x=title_x, fontweight='bold')


# --- Now look on group level instead of spot level

# Calculate LFCs
plot_df = totals_group.pivot_table(
    index='lr_pair',
    columns='dominant_HEV_label',
    values='count',
    fill_value=0
)

plot_df['LFC'] = plot_df.apply(
    lambda row: np.log2(row.get(HEV_labels[0], 1) / row.get(HEV_labels[1], 1)),
    axis=1
)   


# --- Plot 8a: Barplot of top LRs showing dominant HEV labels across samples, sorted by LFC

# Subset to only LRs that are shared across samples
plot_df = plot_df.loc[shared_lrs]

nrows = 2
ncols = 1

fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(10*ncols, 2*nrows))

# Sort by LFC (first descending)
ax = axs[0]
plot_df = plot_df.sort_values(by='LFC', ascending=False)
plot_df_sub = plot_df.head(top)

plot_df_sub[HEV_labels].plot(kind='bar', stacked=True, color=color_dict, width=0.8, ax=ax)

# Then plot in ascending order
ax = axs[1]
plot_df = plot_df.sort_values(by='LFC', ascending=True)
plot_df_sub = plot_df.head(top)

plot_df_sub[HEV_labels].plot(kind='bar', stacked=True, color=color_dict, width=0.8, ax=ax)

# Adjustments
for n, ax in enumerate(axs):
    plotting_functions.beautify_plot(ax)
    # Set grid behind plot
    ax.set_axisbelow(True)
    ax.grid(axis='x', visible=False)  # Turn off vertical gridlines
    # Some tweaking
    ax.set_xlabel('')
    ax.set_ylabel('# of spatial groups')
    ax.set_title(f'Significant LR interactions\nEnriched in {HEV_labels[n]}', fontweight='bold')
    # Move the legend
    ax.legend(handlelength=0.8, loc='center left', bbox_to_anchor=(1, 0.5), alignment='left', frameon=False, markerscale=1.5)

# Move the second plot down to accommodate labels
ax_pos = ax.get_position()
ax.set_position([ax_pos.x0, ax_pos.y0 - 0.45, ax_pos.width, ax_pos.height])

plt.suptitle(f'All samples - stlearn significant LR interactions top {top} by LFC HEV label (only LRs shared across samples)', y=1.1, x=0.5, fontweight='bold')


# --- Plot 8b: Scatterplot showing number of groups

# Plotting
nrows = 1
ncols = 4

# Find the most extreme LFC (positive or negative)
max_lfc = max(abs(plot_df['LFC'].min()), abs(plot_df['LFC'].max()))

# Define a midpoint for the TwoSlopeNorm
midpoint = 0
norm = TwoSlopeNorm(vmin=-max_lfc, vcenter=midpoint, vmax=max_lfc)

# Sort by absolute LFC, descending → prioritize extreme values
sorted_df = plot_df.reindex(plot_df['LFC'].abs().sort_values(ascending=False).index)

fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(3*ncols, 3*nrows))

# Increase the distance between the plots to avoid overlap
spacing = 0.6  # Fraction of width to use as spacing between plots

for i, ax in enumerate(axs):
    pos = ax.get_position()
    # Shift each axis to the right by (i * spacing * width)
    new_x0 = pos.x0 + i * pos.width * spacing
    ax.set_position([new_x0, pos.y0, pos.width, pos.height])

# First plot
ax = axs[0]

# Make a scatter plot of number of groups per HEV label for each LR
sns.scatterplot(data=sorted_df.iloc[::-1], x=HEV_labels[1], y=HEV_labels[0], hue='LFC', palette='coolwarm', ax=ax, legend=False,
                 s=20, linewidth=0.2, hue_norm=norm, edgecolor='gray')

adjust_scatterplot(ax, sorted_df, max_lfc, labels='Number of groups')

# Add a second plot with same data but smaller spots & annotations
ax = axs[1]

sns.scatterplot(data=sorted_df.iloc[::-1], x=HEV_labels[1], y=HEV_labels[0], hue='LFC', palette='coolwarm', ax=ax, legend=False,
                 s=10, linewidth=0.2, hue_norm=norm, edgecolor='gray')

adjust_scatterplot(ax, sorted_df, max_lfc, labels='Number of groups')

# Adjust the plot limits to give extra space for annotations
ax.set_xlim(-ax.get_xlim()[1]*0.1, ax.get_xlim()[1])
ax.set_ylim(-ax.get_xlim()[1]*0.1, ax.get_ylim()[1])

# Get current axis limits
xlim = ax.get_xlim()
ylim = ax.get_ylim()

draw_rectangles(ax, xlim, ylim)

# Annotate the points with LR names
min_distance = 10  # Adjust based on your scale
if bin_size == 25:
    min_distance = 30

adjust_text_func(ax, sorted_df, min_distance=min_distance)

# Add a third plot annotating just some interactions of interest
ax = axs[2]

sns.scatterplot(data=sorted_df.iloc[::-1], x=HEV_labels[1], y=HEV_labels[0], hue='LFC', palette='coolwarm', ax=ax, legend=False,
                 s=10, linewidth=0.2, hue_norm=norm, edgecolor='gray')

adjust_scatterplot(ax, sorted_df, max_lfc, labels='Number of groups')


# Adjust the plot limits to give extra space for annotations
ax.set_xlim(-ax.get_xlim()[1]*0.1, ax.get_xlim()[1])
ax.set_ylim(-ax.get_xlim()[1]*0.1, ax.get_ylim()[1])

# Get current axis limits
xlim = ax.get_xlim()
ylim = ax.get_ylim()

draw_rectangles(ax, xlim, ylim)

adjust_text_func(ax, sorted_df, min_distance=0, labels_to_annotate=lr_of_interest, lw=0.5, color='black')

# Add a fourth plot highlighting specific interactions
ax = axs[3]

key_lrs = ['LTB_LTBR', 'CCL19_CCR7', 'CCL21_CCR7']

# Plot all points with their colors first
sns.scatterplot(data=sorted_df.iloc[::-1], x=HEV_labels[1], y=HEV_labels[0], hue='LFC', palette='coolwarm', ax=ax, legend=False,
                 s=10, linewidth=0.2, hue_norm=norm, edgecolor='gray', zorder=1)

# Plot the key LRs on top with black outline
key_df = sorted_df.loc[sorted_df.index.isin(key_lrs)]
sns.scatterplot(data=key_df, x=HEV_labels[1], y=HEV_labels[0], hue='LFC', palette='coolwarm', ax=ax, legend=False,
                 s=10, linewidth=0.5, hue_norm=norm, edgecolor='black', zorder=100)

adjust_scatterplot(ax, sorted_df, max_lfc, labels='Number of groups')

# Get current axis limits - add 5% buffer on low end for visibility
xlim = ax.get_xlim()
ylim = ax.get_ylim()
x_buffer = xlim[1] * 0
y_buffer = ylim[1] * 0
ax.set_xlim(-x_buffer, xlim[1])
ax.set_ylim(-y_buffer, ylim[1])

# Update limits after resetting
xlim = ax.get_xlim()
ylim = ax.get_ylim()

draw_rectangles(ax, xlim, ylim)

adjust_text_func(ax, sorted_df, min_distance=0, labels_to_annotate=key_lrs, lw=0.5, color='black', fontsize=8)

# Edit the labels of the x & y axis
for ax in axs:
    xlabel = ax.get_xlabel()
    ylabel = ax.get_ylabel()
    for hev_label in HEV_labels:
        if hev_label in xlabel:
            new_xlabel = xlabel.replace(hev_label, rename_dict.get(hev_label, hev_label))
            ax.set_xlabel(new_xlabel)
        if hev_label in ylabel:
            new_ylabel = ylabel.replace(hev_label, rename_dict.get(hev_label, hev_label))
            ax.set_ylabel(new_ylabel)

# Calculate the midpoint to position the title
title_x = (axs[0].get_position().x0 + axs[-1].get_position().x1)/2

plt.suptitle(f'All samples - stlearn significant LR interactions top {top} by LFC HEV label (only LRs shared across samples)', y=1.1, x=title_x, fontweight='bold')


# --- Save the output

qc_dir = os.path.join(project_dir, 'results', '04_vascular_microenvironment', '04_ligand_receptor_analysis', '02_analyze_LR_all', f'{bin_size}um_bin')

today = date.today()
today = today.strftime("%Y%m%d")

out_dir = os.path.join(qc_dir)

os.makedirs(out_dir, exist_ok=True)

# Save plots
out_dir_plots = os.path.join(out_dir, 'plots', today)
os.makedirs(out_dir_plots, exist_ok=True)

# Also save the groupings as a .csv
out_dir_csv = os.path.join(out_dir, 'csv_files', today)
os.makedirs(out_dir_csv, exist_ok=True)

# Save plots
filename = f'All_samples_HEV_LR_analysis.pdf'
save_multi_image(os.path.join(out_dir_plots, filename))

# And save the dfs
filename = f'All_samples_LR_summary.csv'
lr_counts.to_csv(os.path.join(out_dir_csv, filename), index=False)

# Prepare the group counts df for saving
plot_df = totals_group.pivot_table(
    index='lr_pair',
    columns='dominant_HEV_label',
    values='count',
    fill_value=0
)

# Calculate LFCs
plot_df['LFC'] = plot_df.apply(
    lambda row: np.log2(row.get(HEV_labels[0], 1) / row.get(HEV_labels[1], 1)),
    axis=1
)   

# Add the number of samples significant
plot_df['n_samples_sig'] = lr_counts.set_index('LR')['n_samples_sig']

filename = f'All_samples_group_counts.csv'
plot_df.to_csv(os.path.join(out_dir_csv, filename), index=True)

# Same for spot counts
# Prepare the spot counts df for saving
plot_df = totals_spot.pivot_table(
    index='lr_pair',
    columns='dominant_HEV_label',
    values='count',
    fill_value=0
)

# Calculate LFCs
plot_df['LFC'] = plot_df.apply(
    lambda row: np.log2(row.get(HEV_labels[0], 1) / row.get(HEV_labels[1], 1)),
    axis=1
)   

# Add the number of samples significant
plot_df['n_samples_sig'] = lr_counts.set_index('LR')['n_samples_sig']

filename = f'All_samples_spot_counts.csv'
plot_df.to_csv(os.path.join(out_dir_csv, filename), index=True)


# --- Finish

print("Script successfully completed")
