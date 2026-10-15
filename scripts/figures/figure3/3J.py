#!/usr/bin/env python3

# --- Aim
# Plot Figure 3J: barplot showing PVL subtype association with HEV-overlapping vs non-overlapping vessels

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from scipy import stats

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
sys.path.append(os.path.join(project_dir, 'scripts', '04_vascular_microenvironment',))
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

decon_level = 'minor'

# To make it easier to visualize trends by lineage
class_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
class_file = 'miniatlas_subset_classifications.csv'

class_df = pd.read_csv(os.path.join(class_dir, class_file))

# Get endothelial classes
endothelial = list(class_df[class_df['miniatlas_major'] == 'Endothelial'][f'miniatlas_{decon_level}'].unique())
endothelial = sorted(endothelial)

# Decide which decon level to use
decon_level = 'subset'

# Also get PVL classes
PVL_classes = list(class_df[class_df['miniatlas_major'] == 'PVL'][f'miniatlas_{decon_level}'].unique())
PVL_classes = sorted(PVL_classes)


# --- Load composition results

decon_resolution = 'square_016um'

# Define which endothelial type want to look at
endo_of_interest = 'Endothelial ACKR1'
out_name = endo_of_interest.split(' ')[-1]

# Specify path
TLS_group_dir = os.path.join(res_dir, '04_vascular_microenvironment', '03_compare_HEV_to_other_venous', '03_compare_PVL_association', 'csv_files')
most_recent_dir_path = TLS_group_dir

# Load per sample
res_dict = {}

for sample in sample_list:
    filename = f'{sample}_PVL_composition_vs_distance.csv'
    res_dict[sample] = pd.read_csv(os.path.join(most_recent_dir_path, filename), index_col=0)

# Extract a single df
res_df = pd.concat(
    [df.assign(sample_id=key) for key, df in res_dict.items()],
    ignore_index=True
)


# --- Prepare for plotting

plt.rcParams.update({'font.size': 10})

# Prepend dist_bin_ to endothelial keys
endothelial_keys = [f'dist_bin_{cell_type}' for cell_type in endothelial]

# Enable gridlines behind the plot for all plots
sns.set_style("ticks",{'axes.grid' : True})
plt.rcParams['axes.axisbelow'] = True

# Also make the PVL names a bit more interpretable
name_dict = {
    'PVL Differentiated s3': 'PVL s3\n(Contractile)',
    'PVL Immature s1': 'PVL s1\n(Inflammatory)',
    'PVL_Immature s2': 'PVL s2\n(Hypoxia response)',
}


# --- Also get some colors

# Define some colors per endothelial type
colors_file = os.path.join(project_dir, 'resources', 'plotting_colors', 'endothelial_colors.csv')
colors_df = pd.read_csv(colors_file)

# Make a color_dict
color_dict = colors_df.set_index('endothelial_type')['color'].to_dict()

# Also load cell type colors
color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
color_file = f'miniatlas_{decon_level}_colors.csv'

color_df = pd.read_csv(os.path.join(color_dir, color_file), index_col=0)
color_dict_cell = color_df['color'].to_dict()


# --- Prepare df for plotting

group_types = ['no_overlap', 'overlap']

# Set grid behind plot
plt.rcParams['axes.axisbelow'] = True

# Update color_dict
color_dict['no_overlap'] = color_dict[endo_of_interest]

# Get the first bin
bin_order = sorted(res_df[f'dist_bin_{group_types[0]}'].unique(), key=lambda x: int(x.split('-')[0]))
first_bin = bin_order[0]

res_dict = {}

# Calculate mean PVL values per sample for each cell type in the first bin
for cell_type in group_types:
    first_bin_df = res_df[res_df[f'dist_bin_{cell_type}'] == first_bin]
    mean_per_sample = first_bin_df.groupby('sample_id')[PVL_classes].mean()
    mean_per_sample.columns = [f'{cell_type}_{col}' for col in mean_per_sample.columns]
    res_dict[cell_type] = mean_per_sample

# Merge across cell types
mean_df = pd.concat(res_dict.values(), axis=1)
mean_df = mean_df.reset_index()

# Now calculate LFCs between the two groups
group1, group2 = group_types[::-1]  # e.g. ['A', 'B']
lfc_list = []

for pvl in PVL_classes:
    col1 = f'{group1}_{pvl}'
    col2 = f'{group2}_{pvl}'
    lfc = np.log2((mean_df[col1] + 1e-6) / (mean_df[col2] + 1e-6))  # small constant to avoid divide-by-zero
    mean_df[f'LFC_{pvl}'] = lfc
    lfc_list.append(f'LFC_{pvl}')

# Optional: melt for plotting
lfc_long = mean_df.melt(id_vars='sample_id', value_vars=lfc_list,
                        var_name='PVL_class', value_name='log2FC')

# Trim the 'LFC_' prefix from the PVL_class column
lfc_long['PVL_class'] = lfc_long['PVL_class'].str.replace('LFC_', '', regex=False)

# Melt mean values for both groups
mean_long = mean_df.melt(
    id_vars='sample_id',
    value_vars=[f'{g}_{pvl}' for g in group_types for pvl in PVL_classes],
    var_name='group_PVL',
    value_name='mean_PVL_value'
)

# Split 'group_PVL' into 'group' and 'PVL_class'
mean_long[['group', 'PVL_class']] = mean_long['group_PVL'].str.extract(r'^(no_overlap|overlap)_(.*)$')
mean_long = mean_long.drop(columns='group_PVL')

# Pivot so each row has both group means
mean_wide = mean_long.pivot_table(index=['sample_id', 'PVL_class'], columns='group', values='mean_PVL_value').reset_index()

# Merge with lfc_long
lfc_long = lfc_long.merge(mean_wide, on=['sample_id', 'PVL_class'], how='left')

lfc_long['change'] = lfc_long['overlap'] - lfc_long['no_overlap']


# --- Also calculate some stats

stat_results = []
metrics_to_test = ['log2FC', 'change']

for metric in metrics_to_test:
    for pvl in lfc_long['PVL_class'].unique():
        # Get data for this specific PVL class
        data = lfc_long[lfc_long['PVL_class'] == pvl][metric].dropna()
        
        # Check if we actually have data
        if len(data) < 2:
            p_val = np.nan
        else:
            # FORCE T-TEST because n=4
            t_stat, p_val = stats.ttest_1samp(data, 0)
        
        stat_results.append({
            'PVL_class': pvl,
            'metric': metric,
            'n_samples': len(data),
            'mean': data.mean(),
            'p_value': p_val,
            # Add a 'consistency' metric: how many are positive/negative?
            'consistency': f"{sum(data > 0)}/{len(data)} pos" 
        })

stats_df = pd.DataFrame(stat_results)

# Display results sorted by significance
print(stats_df.sort_values(['metric', 'p_value']))


# --- Plot 1: Plot LFCs for the first bin

# Properties for boxplot
whisker_props = {'color': 'black'}
median_props = {'color': 'black'}

plot_key = 'log2FC'

fig, ax = plt.subplots(figsize=(2, 2))

sns.boxplot(lfc_long, x=plot_key, y='PVL_class', showfliers=False, linewidth=1, whiskerprops=whisker_props, 
                medianprops=median_props, width=0.8, palette=color_dict_cell, ax=ax, hue='PVL_class', saturation=1,
                orient='h')

# Plot the individual dots on the boxplot
sns.stripplot(lfc_long, x=plot_key, y='PVL_class', color='black', jitter=True, size=3,
                marker='o', ax=ax)

# -- Add Statistical Annotations
# Filter for the metric we are currently plotting
current_stats = stats_df[stats_df['metric'] == plot_key].set_index('PVL_class')

# Get y-tick locations and labels
y_locs = ax.get_yticks()
y_labels = [label.get_text() for label in ax.get_yticklabels()]

# Iterate through each row on the plot
for loc, label in zip(y_locs, y_labels):
    if label in current_stats.index:
        p_val = current_stats.loc[label, 'p_value']
        
        # Define significance stars
        if p_val < 0.001: star = '***'
        elif p_val < 0.01: star = '**'
        elif p_val < 0.05: star = '*'
        else: star = 'ns' # or '' if you prefer to hide non-sig
        
        # Calculate where to put the star (Right of the max value + buffer)
        # Find the max value for this specific row to place the star correctly
        row_data = lfc_long[lfc_long['PVL_class'] == label]
        max_val = row_data[plot_key].max()
        
        # Add the text
        # x position is max_val + a small buffer (e.g., 10% of x-axis range)
        x_buffer = (lfc_long[plot_key].max() - lfc_long[plot_key].min()) * 0.05
        if p_val < 0.05:
            y_offset = 0.08
        else:
            y_offset = 0
        
        # Color the significant ones black, non-sig gray
        text_color = 'black'
        font_weight = 'bold' if p_val < 0.05 else 'normal'
        
        ax.text(max_val + x_buffer, loc + y_offset, star, va='center', ha='left', 
                fontsize=10, color=text_color, fontweight=font_weight)

# Increase xmax to make space for the significance stars
xmax = ax.get_xlim()[1]
ax.set_xlim(ax.get_xlim()[0], xmax * 1.3)

# Beautify plot
plotting_functions.beautify_plot(ax)

# Draw a line at y=0
ax.axvline(0, color='lightgray', linestyle='-', linewidth=3, zorder=-1)

# Update the yticklabels
yticklabels = ax.get_yticklabels()
new_labels = [name_dict.get(label.get_text(), label.get_text()) for label in yticklabels]

# Remove the newline characters for better display
new_labels = [label.replace('\n', ' ') for label in new_labels]

ax.set_yticklabels(new_labels)

# Adjust title & labels
ax.set_xlabel(f'LFC(mean proportion)')
ax.set_ylabel('')
ax.set_title(f'PVL association\n({first_bin}µm)', fontweight='bold', pad=18, fontsize=12)

# Add arrows just below the text
arrowprops = dict(arrowstyle='-|>', color='black', lw=1, shrinkA=0, shrinkB=0, joinstyle='miter')
arrow_y = 1.03
ax.annotate('', xy=(1, arrow_y), xytext=(0.7, arrow_y), xycoords='axes fraction', textcoords='axes fraction', arrowprops=arrowprops, clip_on=False)
ax.annotate('', xy=(0, arrow_y), xytext=(0.3, arrow_y), xycoords='axes fraction', textcoords='axes fraction', arrowprops=arrowprops, clip_on=False)

# Add text just above and centered the arrows
midpoint = 0.3/2
label_size = 10
buffer = 0.03
ax.text(1-midpoint, 1.04, 'HEV', fontsize=label_size, ha='center', va='bottom', color='black', clip_on=False, transform=ax.transAxes)    
ax.text(0+midpoint, 1.04, 'Venous', fontsize=label_size, ha='center', va='bottom', color='black', clip_on=False, transform=ax.transAxes)

plt.suptitle(f'All samples - PVL composition across first bin (decon_res={decon_resolution}, decon_level={decon_level}, plot_key={plot_key})',
             y=1.5, fontweight='bold', x=0.5)


# --- Save rest of the plots

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_3J.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')


# --- Finish

print("Script successfully completed")
