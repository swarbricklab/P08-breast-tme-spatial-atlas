#!/usr/bin/env python3

# --- Aim
# Plot Fig. 3k: relationship between number of HEVs and TLS size
# Shows correlation between number of HEVs proximal to a TLS and the size of the TLS (number of 48 µm bins)

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy import stats

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments
parser = argparse.ArgumentParser(description="Fig. 3k: number of HEVs vs TLS size",
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

sample_list = sorted(sample_list)

# --- Load HEV classifications

threshold = 0.5
threshold_adj = str(threshold).split('.')[-1]

HEV_dir = os.path.join(res_dir, '04_vascular_microenvironment', '01_identify_HEV_niches', '03b_calculate_HEV_TLS_distance', 'square_016um', f'threshold_{threshold_adj}', 'csv_files')
most_recent_dir_path = HEV_dir

# Read in the results
HEV_dict = {}
TLS_dict = {}

for sample in sample_list:
    TLS_filename = f'{sample}_HEV_per_TLS_group.csv'
    TLS_dict[sample] = pd.read_csv(os.path.join(most_recent_dir_path, TLS_filename))

# Convert to single df for all samples
TLS_res = pd.concat(TLS_dict, names=['sample_id']).reset_index(level=0).rename(columns={'sample_id': 'sample_id'})


# --- Plot 1: HEVs vs TLS size

# Enable gridlines behind the plot for all plots
sns.set_style("ticks",{'axes.grid' : True})
plt.rcParams['axes.axisbelow'] = True

ncols = 1
nrows = 1

xval = 'HEV_count'
yval = 'TLS_size'

# Plotting
fig, axs = plt.subplots(nrows, ncols, figsize=(2.5*ncols, 2.5*nrows))

if ncols == 1:
    axs = [axs]

# -- Plot Relationship TLS size vs number of HEVs

n = 0
ax = axs[n]

sns.scatterplot(TLS_res, x=xval, y=yval, ax=ax, linewidth=0, alpha=0.8, hue='sample_id', s=8)

# Calculate coordinates for the text (pearson correlation value and p-value)
x_coordinate = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1] - ax.get_xlim()[0])
y_coordinate = ax.get_ylim()[1] - 0.05 * (ax.get_ylim()[1] - ax.get_ylim()[0])
df_clean = TLS_res[[xval, yval]].dropna()
pearson_corr_val, p_value = stats.pearsonr(df_clean[xval], df_clean[yval])
p_text = 'p < 0.001' if p_value < 0.001 else f'p = {p_value:.3f}'
text = ax.text(x_coordinate, y_coordinate, f'r = {pearson_corr_val:.2f}\n{p_text}', verticalalignment='top', horizontalalignment='left', fontsize=10)

# Adjust plot labels
ax.set_ylabel('TLS size\n(# 48 µm bins)')
ax.set_xlabel('# proximal HEVs')

# Set grid behind plot
ax.set_axisbelow(True)

# Add title
ax.set_title('# HEVs vs TLS size', fontweight='bold')

# Move legend
legend = ax.legend()

legend.set_title('')
ax.legend(handlelength=0.8, loc='upper left', bbox_to_anchor=(1, 1), alignment='left', frameon=False, markerscale=2)

plt.suptitle(f'All samples - Relationship HEV size & distance to closest TLS (threshold={threshold})', y=1.05, fontweight='bold', x=0.5)


# --- Save output

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_3K.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')


# --- Finish

print("Script successfully completed")
