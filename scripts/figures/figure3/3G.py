#!/usr/bin/env python3

# --- Aim
# Plot Figure 3G: stacked barplot showing proportion of HEV-overlapping vessels per endothelial lineage
# Shows how different endothelial lineages (arterial, capillary, venous, lymphatic) associate with HEVs based on overlap

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import os
import sys
from datetime import datetime

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Import custom functions
from sophiesfunctions.save_multi_image import save_multi_image

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


# --- Extract some info from the sample sheet

# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_dir = os.path.join(project_dir, 'config')
    sample_sheet_file = 'sample_sheet_visium_hd.csv'
    sample_sheet_path = os.path.join(sample_dir, sample_sheet_file)

sample_sheet = pd.read_csv(sample_sheet_path)
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

# Sort according to desired order
plot_order = ['CXCL12', 'RGS5', 'ACKR1', 'LYVE1']

# plot_order contains substrings, endothelial contains full names
sorted_endothelial = []
for key in plot_order:
    matches = [e for e in endothelial if key in e]
    sorted_endothelial.extend(matches)

# Optionally, check for any endothelial types not matched (to append at the end)
remaining = [e for e in endothelial if e not in sorted_endothelial]
sorted_endothelial.extend(remaining)

endothelial = sorted_endothelial

# Also get PVL classes
PVL_classes = list(class_df[class_df['miniatlas_major'] == 'PVL'][f'miniatlas_{decon_level}'].unique())
PVL_classes = sorted(PVL_classes)

# And a dict for the endothelial names
endo_name_dict = {
    'CXCL12': 'Arterial',
    'RGS5': 'Capillary',
    'ACKR1': 'Venous',
    'LYVE1': 'Lymphatic',
}


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

# Specify path
TLS_group_dir = os.path.join(project_dir, 'results', '04_vascular_microenvironment', '02_identify_all_vasculature', '02_investigate_vessels', 'csv_files')
most_recent_dir_path = find_most_recent(TLS_group_dir)

vessel_HEV_overlap = {}

for sample in sample_list:
    filename = f'{sample}_vessel_HEV_overlap.csv'
    vessel_classes = pd.read_csv(os.path.join(most_recent_dir_path, filename), index_col=0)
    vessel_HEV_overlap[sample] = vessel_classes

# Extract a single df
HEV_overlap_df = pd.concat(
    [df.assign(sample_id=key) for key, df in vessel_HEV_overlap.items()],
    ignore_index=True
)


# --- Define some colors

# Also load the colors for the endothelial types
colors_file = os.path.join(project_dir, 'resources', 'plotting_colors', 'endothelial_colors.csv')
colors_df = pd.read_csv(colors_file)

# Make a color_dict
color_dict = colors_df.set_index('endothelial_type')['color'].to_dict()


# --- Plot 1: Groups with HEV vs no HEV overlap

# Disable the grid
plt.rcParams['axes.grid'] = False

# Calculate means per group
counts_df = HEV_overlap_df.groupby(['sample_id', 'group_type', 'HEV_overlap',]).size().unstack(fill_value=0)
counts_df = counts_df.reset_index()
counts_df.columns = counts_df.columns.astype(str)

rowsums = counts_df[['True', 'False']].sum(axis=1)
prop_df = counts_df.copy()
prop_df[['True', 'False']] = counts_df[['True', 'False']].div(rowsums, axis=0)

# Set a minimum overlap threshold to consider a vessel as HEV-overlapping
min_overlap = 0.5

HEV_overlap_df[f'HEV_overlap_min_{min_overlap}'] = np.where(HEV_overlap_df['pct_overlap'] > min_overlap, True, False)

# Calculate means per group
counts_df_min = HEV_overlap_df.groupby(['sample_id', 'group_type', f'HEV_overlap_min_{min_overlap}',]).size().unstack(fill_value=0)
counts_df_min = counts_df_min.reset_index()
counts_df_min.columns = counts_df_min.columns.astype(str)

rowsums = counts_df_min[['True', 'False']].sum(axis=1)
prop_df_min = counts_df_min.copy()
prop_df_min[['True', 'False']] = counts_df_min[['True', 'False']].div(rowsums, axis=0)

# Step 1: Aggregate by group_type
agg = prop_df_min.groupby('group_type')[['True', 'False']].agg(['mean', 'std'])

# Flatten the column MultiIndex
agg.columns = ['_'.join(col).strip() for col in agg.columns.values]
agg = agg.reset_index()

# Order group_type according to endothelial
agg['group_type'] = pd.Categorical(agg['group_type'], categories=endothelial, ordered=True)
# Sort
agg = agg.sort_values('group_type')

# Step 2: Extract values
x = np.arange(len(agg))
true_means = agg['True_mean']
false_means = agg['False_mean']
true_stds = agg['True_std']
false_stds = agg['False_std']

# Plotting
nrows = 1
ncols = 1
fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(2*ncols, 2))

if ncols == 1:
    axs = np.array([axs])  # Ensure axs is always a 2D array for consistent indexing

ax = axs[0]
color_list = [color_dict.get(gt, 'lightgray') for gt in agg['group_type']]

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
ax.set_title(f'Vessels identified as HEV\nper lineage', fontweight='bold')
ax.grid(axis='x', visible=False)  # Turn off vertical gridlines

ax.set_xticks(x)

# Edit the xlabels
xlabels = plot_order.copy()

for i, label in enumerate(xlabels):
    xlabels[i] = endo_name_dict.get(label, label)

ax.set_xticklabels(xlabels, rotation=90)

plt.suptitle(f'All samples - Overlap vessel groups with HEVs (min_overlap={min_overlap})', 
             y=1.2, x=0.5, fontweight='bold')


# --- Save the output

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_3G.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')


# --- Finish

print("Script successfully completed")
