#!/usr/bin/env python3

# --- Aim
# Investigate associations between classified vessels, PVL subtypes, and HEV overlap.
# Analyze spatial relationships and compositional trends in the vascular microenvironment.

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import os
import sys
from collections import defaultdict

import anndata as ad
import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from matplotlib import colors

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments


parser = argparse.ArgumentParser(description="Import Visium data & do basic QC",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV", default=None)
args = vars(parser.parse_args())


# --- Arguments

res_dir = args["res_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
project_dir=args["project_dir"]
sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import auto_crop_scanpy, save_multi_image

# Import function script as a module
sys.path.append(os.path.join(project_dir, 'scripts', '04_vascular_microenvironment'))
import plotting_functions


# --- Get some info from the sample sheet
# Use --sample_sheet if provided, otherwise default to config/sample_sheet_visium_hd.csv
if args["sample_sheet"] is not None:
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_visium_hd.csv')

sample_sheet = pd.read_csv(sample_sheet_path)
tissue_id = sample_sheet[sample_sheet['section_id'] == sample_id]['tissue_id'].values[0]


# --- Also get some tissue metadata

resource_dir = os.path.join(project_dir, 'resources')
metadata_file = 'tissue_metadata.csv'

tissue_metadata = pd.read_csv(os.path.join(resource_dir, metadata_file))

# Create a sample - subtype dictionary
subtype = tissue_metadata[tissue_metadata['tissue_id'] == str(tissue_id)]['clinical_subtype'].values[0]


# --- Load the adata object

use_data = 'log_norm'
# Need both the 48 & 16 res
resolution = 'square_016um'

adata = ad.read_h5ad(os.path.join(h5ad_dir, 'visium_hd', resolution, f'{sample_id}.h5ad'))


# --- Calculate crop coords for optimal presentation

coords = auto_crop_scanpy(adata, border_size=0.1)


# --- Load decon results

decon_level = 'minor'

decon_dir = os.path.join(res_dir, '02_cell_type_mapping', 'deconvolution', 'visium_hd', f'{decon_level}_level', resolution,)

decon_file = f'{sample_id}_{decon_level}_{resolution}_q05_cell_abundance_w_sf.csv'

decon_res = pd.read_csv(os.path.join(decon_dir, decon_file), index_col=0)

cell_types = list(decon_res.columns)
adata.obs[cell_types] = decon_res


# --- Load miniatlas classifications

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


# --- Load the vessel classifications

# Specify path
TLS_group_dir = os.path.join(res_dir, '04_vascular_microenvironment', '02_identify_all_vasculature', '01_classify_vessels', 'csv_files')

most_recent_dir_path = TLS_group_dir

filename = f'{sample_id}_vessel_classification.csv'
vessel_classes = pd.read_csv(os.path.join(most_recent_dir_path, filename), index_col=0)

# Collapse to single vessel class per spot
collapsed_col = vessel_classes[endothelial].apply(lambda row: row.first_valid_index() if row.first_valid_index() is not None else None, axis=1)

adata.obs['vessel_class'] = collapsed_col

# also load the groups


# --- Also load HEV annotations

# Decide whether to use the expanded groups
expanded = True

threshold_hev = 0.5
threshold_hev_adj = str(threshold_hev).split('.')[-1]

TLS_group_dir = os.path.join(res_dir, '04_vascular_microenvironment', '01_identify_HEV_niches', '02b_expand_HEV_groups', 'square_016um', f'threshold_{threshold_hev_adj}')
TLS_group_dir = os.path.join(TLS_group_dir, 'csv_files')

most_recent_dir_path = TLS_group_dir

# Read in the results & add to anndata
score_of_i = 'HEVs'
TLS_file_name = f'{sample_id}_{score_of_i}_groups.csv'
TLS_results = pd.read_csv(os.path.join(most_recent_dir_path, TLS_file_name), index_col=0)

group_key = 'group_too_small_removed'
group_key_source = 'expanded_groups'

adata.obs['HEV_groups'] = TLS_results[group_key_source]


# --- Define some colors

# Also load the colors for the endothelial types
colors_file = os.path.join(project_dir, 'resources', 'plotting_colors', 'endothelial_colors.csv')
colors_df = pd.read_csv(colors_file)

# Make a color_dict
color_dict = colors_df.set_index('endothelial_type')['color'].to_dict()


# --- Figure out which vessels overlap with HEVs

HEV_index = adata[adata.obs['HEV_groups'].notna()].obs_names

# Find HEVs that intersect with vessel classes
subset_df = vessel_classes[endothelial].loc[vessel_classes[endothelial].index.intersection(HEV_index)]

# Also add indeces where there is partial overlap with a group
# Initialize an empty set to collect matching indices
matching_indices = set()

for col in subset_df.columns:
    # Get the unique values in this column in the subset
    subset_values = subset_df[col].dropna().unique()
    # Find indices in the full df where the value in the same column matches any from the subset
    matches = vessel_classes[vessel_classes[col].isin(subset_values)].index
    # Add those indices to the overall set
    matching_indices.update(matches)


# --- Plot 1: Spatial plots endothelial cells

# Plotting
sc.set_figure_params(fontsize=15)

ncols = 2
nrows = 2

fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(7*ncols, 5*nrows))
axs = axs.ravel()

# -- First the HEVs
n=0
ax = axs[n]

adata_sub = adata[adata.obs['HEV_groups'].notna()]

sc.pl.spatial(adata_sub, color='HEV_groups', bw=False, alpha_img=0.5, crop_coord=coords, 
            title=f'{sample_id} ({subtype})\nHEVs',
            palette=['indigo'], ax=ax, legend_loc=None)

ax.set_title(ax.get_title(), fontweight='bold')

# -- Then plot vessel classes
n+=1
ax = axs[n]

adata_sub = adata[adata.obs['vessel_class'].notna()]

sc.pl.spatial(adata_sub, color='vessel_class', bw=False, alpha_img=0.5, crop_coord=coords, 
            title=f'{sample_id} ({subtype})\nVessel class',
            palette=color_dict, ax=ax,)

ax.set_title(ax.get_title(), fontweight='bold')


# -- Then plot vessels with overlap with HEVs
n+=1
ax = axs[n]

adata_sub = adata[adata.obs_names.isin(matching_indices)]

sc.pl.spatial(adata_sub, color='vessel_class', bw=False, alpha_img=0.5, crop_coord=coords, 
            title=f'{sample_id} ({subtype})\nVessels overlapping with HEVs',
            palette=color_dict, ax=ax)

ax.set_title(ax.get_title(), fontweight='bold')


# -- Then highlight areas with overlap
n+=1
ax = axs[n]

color_dict['overlap'] = 'cyan'

adata_sub.obs['overlap'] = adata_sub.obs['vessel_class']
adata_sub.obs['overlap'] = np.where(adata_sub.obs['HEV_groups'].notna(),
                                    'overlap',
                                    adata_sub.obs['overlap'])

sc.pl.spatial(adata_sub, color='overlap', bw=False, alpha_img=0.5, crop_coord=coords, 
            title=f'{sample_id} ({subtype})\nVessels overlapping with HEVs',
            palette=color_dict, ax=ax)

ax.set_title(ax.get_title(), fontweight='bold')

plt.subplots_adjust(wspace=0.6, hspace=0.3)

plt.suptitle(f'{sample_id} ({subtype}) - HEVs + other vessels overlap', y=1.2, x=0.5, fontweight='bold')


# Extract a df for later
overlap_df = adata_sub.obs[['HEV_groups', 'vessel_class', 'overlap']]


# --- Make a df recording which vessel groups overlap with HEV

group_counts_df = vessel_classes[endothelial].melt(var_name='group_type', value_name='group') \
                    .dropna() \
                    .value_counts(subset=['group_type', 'group']) \
                    .reset_index(name='count')

group_counts_df = group_counts_df.set_index(group_counts_df['group_type'] + '-' + group_counts_df['group'])

# Record HEV overlap
subset_df_melt = subset_df.melt(var_name='group_type', value_name='group').dropna()
subset_df_melt['index'] = subset_df_melt['group_type'] + '-' + subset_df_melt['group']

group_counts_df['HEV_overlap'] = np.where(group_counts_df.index.isin(subset_df_melt['index']),
                                          True,
                                          False)

# Calculate HEV overlap size per vessel group
group_counts = defaultdict(dict)

for col in endothelial:
    # Subset to the matching indices
    subset = vessel_classes.loc[vessel_classes.index.intersection(subset_df.index), col]
    # Count group occurrences (excluding None/NaN)
    group_size = subset.dropna().value_counts()
    # Store in dictionary
    group_counts[col] = group_size.to_dict()

# Convert to DataFrame
tmp_df = pd.DataFrame(group_counts).fillna(0).astype(int)

# Collapse to long format
collapsed = tmp_df.stack().reset_index()
collapsed.columns = ['group', 'original_column', 'count']

# Create new index as original_column-group
collapsed['group_id'] = collapsed['original_column'] + '-' + collapsed['group']

# Set that as the index and keep only the count column
final_df = collapsed.set_index('group_id')[['count']]

group_counts_df['HEV_overlap_size'] = final_df['count']

# Get rid of any 0s that somehow snuck in
group_counts_df['HEV_overlap_size'] = group_counts_df['HEV_overlap_size'].replace(0, np.nan)

# Checking
group_counts_df[group_counts_df['HEV_overlap']]

# Also calculate pct overlap
group_counts_df['pct_overlap'] = group_counts_df['HEV_overlap_size'] / group_counts_df['count']


# --- Now also be more stringent: look only at vessels with at >50% of spots overlapping with HEVs

# Filter the overlap summary df to only include groups with >50% overlap
min_overlap = 0.5
filtered_overlap_df = group_counts_df[group_counts_df['pct_overlap'] > min_overlap]

# Add to adata
adata.obs[f'HEV_overlap_min_{min_overlap}'] = np.nan

for endo_type in endothelial:
    overlap_groups = filtered_overlap_df[filtered_overlap_df['group_type'] == endo_type]['group'].values
    # Set the overlap groups in the adata obs
    matching_indices = vessel_classes[endo_type].isin(overlap_groups)
    matching_barcodes = vessel_classes.index[matching_indices]
    adata.obs[f'HEV_overlap_min_{min_overlap}'].loc[matching_barcodes] = endo_type


# -- Plotting
fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(7*ncols, 5*nrows))
axs = axs.ravel()

# -- First the HEVs
n=0
ax = axs[n]

adata_sub = adata[adata.obs['HEV_groups'].notna()]

sc.pl.spatial(adata_sub, color='HEV_groups', bw=False, alpha_img=0.5, crop_coord=coords, 
            title=f'{sample_id} ({subtype})\nHEVs',
            palette=['indigo'], ax=ax, legend_loc=None)

ax.set_title(ax.get_title(), fontweight='bold')

# -- Then plot vessel classes
n+=1
ax = axs[n]

adata_sub = adata[adata.obs['vessel_class'].notna()]

sc.pl.spatial(adata_sub, color='vessel_class', bw=False, alpha_img=0.5, crop_coord=coords, 
            title=f'{sample_id} ({subtype})\nVessel class',
            palette=color_dict, ax=ax,)

ax.set_title(ax.get_title(), fontweight='bold')


# -- Then plot vessels with overlap with HEVs
n+=1
ax = axs[n]

plot_key = f'HEV_overlap_min_{min_overlap}'
adata_sub = adata[adata.obs[plot_key].notna()]

sc.pl.spatial(adata_sub, color=plot_key, bw=False, alpha_img=0.5, crop_coord=coords, 
            title=f'{sample_id} ({subtype})\nmin_overlap1{min_overlap}',
            palette=color_dict, ax=ax)

ax.set_title(ax.get_title(), fontweight='bold')


# -- Then highlight areas with overlap
n+=1
ax = axs[n]

color_dict['overlap'] = 'cyan'

adata_sub.obs['overlap'] = adata_sub.obs['vessel_class']
adata_sub.obs['overlap'] = np.where(adata_sub.obs['HEV_groups'].notna(),
                                    'overlap',
                                    adata_sub.obs['overlap'])

sc.pl.spatial(adata_sub, color='overlap', bw=False, alpha_img=0.5, crop_coord=coords, 
            title=f'{sample_id} ({subtype})\nVessels overlapping with HEVs',
            palette=color_dict, ax=ax)

ax.set_title(ax.get_title(), fontweight='bold')

plt.subplots_adjust(wspace=0.6, hspace=0.3)

plt.suptitle(f'{sample_id} ({subtype}) - HEVs + other vessels overlap (min_overlap={min_overlap})', y=1.2, x=0.5, fontweight='bold')




# --- Plot 2: Groups with HEV vs no HEV overlap

sc.set_figure_params(fontsize=10)

# Enable gridlines behind the plot for all plots
plt.rcParams['axes.axisbelow'] = True

group_counts_df['HEV_overlap']

# Plotting
nrows = 1
ncols = 4
fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(3*ncols, 2.5))

# -- First stacked barplot showing HEV overlap per category
n=0
ax = axs[n]

# Step 1: Pivot the counts to get the count of True per group
counts_true = group_counts_df.groupby(['group_type', 'HEV_overlap']).size().unstack(fill_value=0)[True]

# Step 2: Pivot the counts to get the count of False per group
counts_false = group_counts_df.groupby(['group_type', 'HEV_overlap']).size().unstack(fill_value=0)[False]

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
ax.set_title(f'{sample_id}\nHEV overlap', fontweight='bold')

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
ax.set_title(f'{sample_id}\nHEV overlap', fontweight='bold')

# -- Then boxplot showing overlap size

# Boxplot size
# Properties for boxplot
whisker_props = {'color': 'black'}
median_props = {'color': 'black'}

n+=1
ax = axs[n]

sns.boxplot(group_counts_df, x='group_type', y='pct_overlap', order=endothelial, palette=color_dict, ax=ax, hue='group_type',
            saturation=1, whiskerprops=whisker_props, medianprops=median_props, showfliers=False)

sns.stripplot(data=group_counts_df, x='group_type', y='pct_overlap', color="black", order=endothelial, size=2, ax=ax, dodge=False,)

# Beautify the plot
plotting_functions.beautify_plot(ax)

ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_xlabel('')
ax.set_title(f'{sample_id}\nHEV overlap size', fontweight='bold')

# -- Then scatterplot overlap vs size colored by type

n+=1
ax = axs[n]

sns.scatterplot(group_counts_df, x='count', y='pct_overlap', palette=color_dict, ax=ax, hue='group_type', hue_order=endothelial,
                s=10, alpha=0.5, edgecolor=None)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)

ax.set_title(f'{sample_id}\nPct overlap vs vessel size', fontweight='bold')


plt.subplots_adjust(wspace=0.3)

plt.suptitle(f'{sample_id} ({subtype}) - Overlap vessel groups with HEVs', 
             y=1.2, x=0.5, fontweight='bold')


# --- Save the output

qc_dir = os.path.join(res_dir, '04_vascular_microenvironment', '02_identify_all_vasculature', '02_investigate_vessels')

out_dir = os.path.join(qc_dir)

os.makedirs(out_dir, exist_ok=True)

# Save plots
out_dir_plots = os.path.join(out_dir, 'plots')
os.makedirs(out_dir_plots, exist_ok=True)

# Also save the groupings as a .csv
out_dir_csv = os.path.join(out_dir, 'csv_files')
os.makedirs(out_dir_csv, exist_ok=True)

# Save plots
filename = f'{sample_id}_vessel_HEV_overlap.pdf'
save_multi_image(os.path.join(out_dir_plots, filename))

# And save the dfs
filename = f'{sample_id}_vessel_HEV_overlap.csv'
group_counts_df.to_csv(os.path.join(out_dir_csv, filename), index=True)

# And save the overlap df
filename = f'{sample_id}_vessel_HEV_overlap_spot_level.csv'
overlap_df.to_csv(os.path.join(out_dir_csv, filename), index=True)


# --- Finish

print("Script successfully completed")
