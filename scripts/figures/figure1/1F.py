#!/usr/bin/env python3

# --- Aim
# Plot Figure 1F: scatterplot of T cell vs B cell proportions across samples
# Shows relationship between T cell and B cell abundance across breast cancer subtypes

# --- Define environment
# env: spatial-pilot

# --- Load packages

import argparse
import os

import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

import matplotlib as mpl
import matplotlib.colors as mcolors

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium data & do basic QC",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--data_dir", help="Path to raw data")
parser.add_argument("-s", "--samples",
    nargs="*",  # 0 or more values expected => creates a list
    type=str,
    help="List of samples")
parser.add_argument("-r", "--res_dir",
                    help="Path to results directory (optional)")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files") # Not necessary for this script
parser.add_argument("--gene_file", help="Path to a csv file with genes, genes stored in first column") # Not necessary for this script
parser.add_argument("--use_data", help="Which data to use (SCT or log_norm)") # Not necessary for this script
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file")

args = vars(parser.parse_args())

# Arguments
data_dir = args["data_dir"]
h5ad_dir = args["h5ad_dir"]
sample_list = args["samples"]
project_dir = args['project_dir']
res_dir = args["res_dir"] or os.path.join(project_dir, 'results')  # default: <project_dir>/results


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


# --- Load the decon results

level = 'major'
decon_dir = os.path.join(res_dir, '02_cell_type_mapping', 'deconvolution', 'visium', f'{level}_level')

decon_dict = {}

for sample in sample_list:
    decon_filename = f'{sample}_{level}_q05_cell_abundance_w_sf.csv'
    # Load decon results
    decon_results = pd.read_csv(os.path.join(decon_dir, decon_filename), index_col=0)
    decon_dict[sample] = decon_results

all_cell_types = list(decon_dict[sample].columns)


# --- Sort the samples per subtype

# Dictionary to store the lists of samples per type
samples_by_subtype = {}

# Group sample IDs by their type
for sample in sample_list:
    subtype = subtype_dict[sample]
    if subtype not in samples_by_subtype:
        samples_by_subtype[subtype] = []
    samples_by_subtype[subtype].append(sample)

# Sort the lists of sample IDs per type
for subtype in samples_by_subtype:
    samples_by_subtype[subtype].sort()

print(samples_by_subtype)

# Define a plotting order for subtypes
subtype_order = ['Luminal A', 'Luminal B', 'Luminal B (HER2+)', 'HER2+', 'TNBC']
sorted_subtypes = sorted(samples_by_subtype.keys(), key=lambda x: subtype_order.index(x) if x in subtype_order else float('inf'))


# --- Get colors per subtype

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
color_file = 'SCSubtype_colors.csv'

color_lut = pd.read_csv(os.path.join(color_dir, color_file))
color_lut = color_lut.set_index('subtype')

# Extract a color dict
color_dict = color_lut['color'].to_dict()


# --- Load miniatlas classifications

# Load miniatlas classifications
classifications_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
classifications_file = 'miniatlas_subset_classifications.csv'

classifications = pd.read_csv(os.path.join(classifications_dir, classifications_file))
classifications.index = classifications[f'miniatlas_{level}']

# Filter for cell types that are actually present
classifications = classifications[classifications[f'miniatlas_{level}'].isin(all_cell_types)]

# Figure out which cells are non-epithelial
non_epithelial = classifications[classifications['lineage'] != 'Epithelial'][f'miniatlas_{level}'].unique().tolist()

# Filter the present cell types
non_epithelial = [cat for cat in all_cell_types if cat in non_epithelial]


# --- Prepare for plotting

sns.set_style("ticks")


# --- Calculate the normalized number of cells per cell type for each sample

# Create a dictionary to store the final values (row sum / number of rows)
final_values = {}

# Loop through each sample and perform the operations
for sample, df in decon_dict.items():
    # Subset to group of interest
    df = df[non_epithelial]
    # Calculate total number of cells per sample
    total_cells = df.sum(axis=0)
    # Calculate proportions
    total_cells = total_cells / total_cells.sum()
    # Store the value in the dictionary
    final_values[sample] = total_cells

# Convert the dictionary into a DataFrame with samples as rows
result_df = pd.DataFrame.from_dict(final_values, orient='index')

print(result_df)


# --- Add the subtype to this df

result_df['subtype'] = result_df.index.map(subtype_dict)


# --- Prepare for Plotting

alpha = 0.7

# Enable gridlines behind the plot for all plots
sns.set_style("ticks",{'axes.grid' : True})
plt.rcParams['axes.axisbelow'] = True

# --- Plot 1: Scatterplot of B & T cells

# Check if B cells and T cells are present in the data
if 'B-cells' in result_df.columns and 'T-cells' in result_df.columns:
    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    # Create scatterplot
    for subtype in sorted_subtypes:
        subtype_data = result_df[result_df['subtype'] == subtype]
        fill_color = mcolors.to_rgba(color_dict[subtype], alpha=0.6)
        edge_color = color_dict[subtype]
        ax.scatter(subtype_data['T-cells'], subtype_data['B-cells'], 
                c=[fill_color], label=subtype, s=30,
                edgecolors=edge_color, linewidth=1)
    # Add labels and title
    ax.set_xlabel('Proportion T cells\n(of non-epithelial)')
    ax.set_ylabel('Proportion B cells\n(of non-epithelial)')
    ax.set_title(f'T cells vs B cells', fontweight='bold')
    # Add legend
    ax.legend(title='', loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, fontsize=8, handletextpad=0.1,
                borderaxespad=0, borderpad=0)
    # Add grid
    ax.grid(True, alpha=0.5)
    # Set limits to start at 0
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    # Add a suptitle
    plt.suptitle(f'All samples - cell2location: T cells vs B cells (annotation_level={level})', y=1.05, fontweight='bold',
            fontsize=12)


# --- Save the output

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_1F.pdf'

# Save the plots
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')


# --- Finish

print("Script successfully completed")