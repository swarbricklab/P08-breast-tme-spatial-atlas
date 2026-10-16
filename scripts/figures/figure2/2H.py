#!/usr/bin/env python3

# --- Aim
# Plot Fig. 2h: barplot showing number of germinal centers per clinical subtype
# GCs identified by signature scoring

# --- Define environment
# env: spatial-pilot

# --- Load packages
import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments
parser = argparse.ArgumentParser(description="Fig. 2h: number of germinal centers per clinical subtype",
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
    sample_sheet_file = 'sample_sheet_visium.csv'
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

# And make a sample patient LUT
sample_donor_lut = {}

for sample in sample_list:
    sample_donor_lut[sample] = sample_lut.loc[sample, 'donor_id']


# --- Load the GC numbers

GC_group_dir = os.path.join(res_dir, '03_TLS_niche_discovery', '04_map_niche_signatures', '02_group_high_score_areas', 'GC_score',)

most_recent_dir_path = GC_group_dir

GC_number_dict = {}

for sample in sample_list:
    sample_file = os.path.join(most_recent_dir_path, 'csv_files', f'{sample}_GC_groups.csv')
    if os.path.isfile(sample_file):
        GC_number_dict[sample] = pd.read_csv(sample_file, index_col=0)
    else:
        print(f"Warning: {sample_file} not found, skipping.")

group_key = 'group_too_small_removed'

# Count number of unique groups per sample
GC_numbers = pd.DataFrame([
    {'sample': sample, 'number': df[group_key].nunique()}
    for sample, df in GC_number_dict.items()
])

# Add subtype information
GC_numbers['subtype'] = GC_numbers['sample'].map(subtype_dict)


# --- Normalize numbers to 5000 spots, to account for differences in total number of spots per sample

TLS_res_dir = os.path.join(res_dir, '03_TLS_niche_discovery', '01_identify_TLS', '01_run_TLS_model')

most_recent_dir_path = TLS_res_dir

# Read in the results & add to anndata
for sample in sample_list:
    TLS_file_name = f'{sample}_TLS_model.csv'
    TLS_results = pd.read_csv(os.path.join(most_recent_dir_path, TLS_file_name), index_col=0)
    total_spots = len(TLS_results)
    GC_numbers.loc[GC_numbers['sample'] == sample, 'number_norm'] = GC_numbers.loc[GC_numbers['sample'] == sample, 'number'] / total_spots * 5000


# --- Read in colors per subtype

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
color_file = 'SCSubtype_colors.csv'

color_lut = pd.read_csv(os.path.join(color_dir, color_file))
color_lut = color_lut.set_index('subtype')

alpha = 0.7

subtype_color_dict = color_lut['color'].to_dict()


# --- Prepare for plotting the result per subtype: define a plotting order based on subtype

# Extract all subtypes
all_subtypes = GC_numbers['subtype'].unique()

# Sort according to preferred plotting order
subtype_order = ['Luminal A', 'Luminal B', 'Luminal B (HER2+)', 'HER2+', 'TNBC']
sorted_subtypes = sorted(all_subtypes, key=lambda x: subtype_order.index(x) if x in subtype_order else float('inf'))


# --- Plot 1: Barplot by subtype

# Set fontsize
plt.rcParams.update({'font.size': 14})

# Disable the grid for all plots
plt.rcParams['axes.grid'] = False
err_kws={"color": 'black', 'linewidth': 0.8}

# Plotting
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3, 0.6*len(sorted_subtypes)))

# Make a horizontal barplot
sns.barplot(GC_numbers, y='subtype', x='number_norm', hue='subtype', palette=subtype_color_dict, order=sorted_subtypes,
            linewidth=0.8, ax=ax, saturation=1, err_kws=err_kws, edgecolor='black', capsize=.4, )

# Adjust alpha for facecolors
for bar in ax.patches:
    # Get the current facecolor
    r, g, b, _ = bar.get_facecolor()
    # Set new facecolor with alpha
    bar.set_facecolor((r, g, b, alpha))  # Set alpha to 0.5
    # Optionally, ensure the edgecolor remains opaque
    bar.set_edgecolor('black')

# Plot the individual dots on the boxplot
sns.stripplot(GC_numbers, y='subtype', x='number_norm', color='black', jitter=True, size=3,
                marker='o', ax=ax)

# Set xlabel
ax.set_xlabel('# GCs / 5000 spots')
ax.set_ylabel('Subtype')

ax.set_title('GC frequency', fontweight='bold')

# Add an ax on the right for number of samples per subtype
count_ax = fig.add_axes([ax.get_position().x1 + 0.02, ax.get_position().y0, ax.get_position().width/3, ax.get_position().height])

# Calculate number of samples per subtype
subtype_counts = GC_numbers['subtype'].value_counts().reindex(sorted_subtypes)
sns.barplot(x=subtype_counts.values, y=subtype_counts.index, ax=count_ax, palette=subtype_color_dict, order=sorted_subtypes,
            saturation=1, edgecolor='black', linewidth=0.8)

# Make the same adjustments as before
for bar in count_ax.patches:
    # Get the current facecolor
    r, g, b, _ = bar.get_facecolor()
    # Set new facecolor with alpha
    bar.set_facecolor((r, g, b, alpha))  # Set alpha to 0.5
    # Optionally, ensure the edgecolor remains opaque
    bar.set_edgecolor('black')

# Set top & right spines invisible
count_ax.spines['top'].set_visible(False)
count_ax.spines['right'].set_visible(False)

# Remove yticks & ylabel
count_ax.set_yticks([])
count_ax.set_ylabel('')

# Set xlabel
count_ax.set_xlabel('# samples')

# Annotate bars with counts
for i, bar in enumerate(count_ax.patches):
    width = bar.get_width()
    count_ax.text(width + count_ax.get_xlim()[1]*0.1, bar.get_y() + bar.get_height() / 2,
                  f'{int(width)}', ha='left', va='center', fontsize=12)

plt.suptitle(f'Number of spatially separate GCs normalized by total spots', y=1.1, fontweight='bold')


# --- Save the output

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_2H.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')

print("Script succesfully completed")
