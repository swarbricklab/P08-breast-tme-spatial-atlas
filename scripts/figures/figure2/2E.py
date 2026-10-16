#!/usr/bin/env python3

# --- Aim
# Plot Figure 2E: deconvolution results per TLS cluster
# Shows cell type enrichment of selected cell types in each TLS cluster


# --- Define environment
# env: spatial-pilot

# --- Load packages
import argparse
import os
import math

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap

mpl.rcParams['pdf.fonttype'] = 42  # Use Type 42 (TrueType) fonts in PDF files


# --- Argparse arguments
parser = argparse.ArgumentParser(description="Fig. 2e: deconvolution results per TLS cluster",
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
res_dir = args["res_dir"] or os.path.join(project_dir, 'results')  # default: <project_dir>/results
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


# --- Load clustering results

# Load harmony clustering results
harmony_dir = os.path.join(res_dir, '03_TLS_niche_discovery', '02_identify_TLS_niches', '01_cluster_TLS_high')
most_recent = harmony_dir

filename = 'All_samples_leiden_clusters.csv'
cluster_res = pd.read_csv(os.path.join(most_recent, 'csv_files', filename), index_col=0)


# --- Load deconvolution results

level = 'subset'
decon_dir = os.path.join(res_dir, '02_cell_type_mapping', 'deconvolution', 'visium', f'{level}_level')

decon_dict = {}

for sample in sample_list:
    decon_filename = f'{sample}_{level}_q05_cell_abundance_w_sf.csv'
    # Load decon results
    decon_results = pd.read_csv(os.path.join(decon_dir, decon_filename), index_col=0)
    # Collapse cycling cell types bc these don't deconvolve very well
    # Identify cycling cell types
    cycling_cell_types = [ct for ct in decon_results.columns if 'cycling' in ct.lower() or 'ki67' in ct.lower()]
    if cycling_cell_types:
        # Sum cycling cell types into a new column
        decon_results['Cycling'] = decon_results[cycling_cell_types].sum(axis=1)
        # Drop the original cycling cell type columns
        decon_results = decon_results.drop(columns=cycling_cell_types)
    decon_dict[sample] = decon_results

all_cell_types = list(decon_dict[sample].columns)

# Load miniatlas classifications
classifications_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
classifications_file = 'miniatlas_subset_classifications.csv'

classifications = pd.read_csv(os.path.join(classifications_dir, classifications_file))

# Also make a dict with shortened names to update the xticklabels
subset_short_lut = dict(zip(classifications['miniatlas_subset'], classifications['miniatlas_subset_short']))


# --- Add clustering info to the decon info

for sample in sample_list:
    decon_results = decon_dict[sample]
    clustering_results = cluster_res[cluster_res['sample_id'] == sample]
    decon_results[clustering_results.columns] = clustering_results.astype(str)
    # Drop the non-TLS high spots (these don't have a cluster identity)
    decon_results = decon_results.dropna(subset='leiden')
    decon_dict[sample] = decon_results


# --- Also load top genes

filename = 'All_samples_leiden_clusters_DEGS.csv'
top_genes_table = pd.read_csv(os.path.join(most_recent, 'csv_files', filename))

# Cluster numbers may shift a bit but should be pretty clear from gene expression which is which
top_marker_dict = {
    'T cell zone': ['CCL19', 'LAMP3'],
    'B cell zone': ['CXCL13', 'MS4A1'],
    'HEVs': ['CCL21', 'ACKR1'],
    'Germinal center': ['RGS13', 'MEF2B'],
}

# Get all unique cluster labels as strings
all_clusters = sorted(top_genes_table['group'].astype(str).unique(), key=lambda x: int(x) if x.isdigit() else x)

clusters_of_interest = {}

for region, markers in top_marker_dict.items():
    for cluster in all_clusters:
        # Get the top 20 genes for this cluster
        cluster_genes = top_genes_table[top_genes_table['group'].astype(str) == cluster]['genes'].values[:20]
        if all(marker in cluster_genes for marker in markers):
            clusters_of_interest[cluster] = region
            break

print("Clusters of interest based on marker gene expression:")
for cluster, region in clusters_of_interest.items():
    print(f"{region}: Cluster {cluster}")


# --- Prepare for plotting

def plot_custom_violin_seaborn(df, plot_key, ax, groupby='leiden', hue=None, 
                               hue_palette=None, cut=2, palette=None, plot_height=2, order=None):
    # If no order is provided, use sorted unique values
    if order is None:
        order = sorted(df[groupby].unique())
    # Make the violin
    sns.violinplot(
        data=df, 
        x=groupby, y=plot_key, hue=hue, legend=False, ax=ax, linewidth=0.8, palette=palette, 
        inner=None, density_norm='width', cut=cut, saturation=1, order=order
    )
    # Calculate medians in the same order
    medians = [df[df[groupby] == cat][plot_key].median() for cat in order]
    medians = np.array(medians)
    min_median, max_median = medians.min(), medians.max()
    normalized_alphas = 0.05 + 0.85 * (medians - min_median) / (max_median - min_median)
    for collection, alpha in zip(ax.collections, normalized_alphas):  # Adjust the colors of the violin
        fill_color = collection.get_facecolor()
        darkened_color = fill_color.copy()
        darkened_color[:,:-1] *= 0.8
        collection.set_edgecolor(darkened_color) 
        fill_color[0][-1] = alpha
        collection.set_facecolor(fill_color)
    ax.set_xlabel('')
    # plot median dots
    pos = range(len(medians))
    for tick, median in zip(pos, medians):
        ax.scatter(tick, median, marker='D', color='black', s=10, zorder=3)
    # Adjust the positions of the subplots
    xtick_labels = [label.get_text() for label in ax.get_xticklabels()]
    ax_pos = ax.get_position()
    ax.set_position([ax_pos.x0, ax_pos.y0, 0.08*len(xtick_labels), plot_height])
    margin = 0.75
    ax.set_xlim([0-margin, len(xtick_labels)-1+margin])
    ax.get_xlim()
    # Add colorbar for medians
    # Create a colormap reflecting the color variations with alpha
    n_colors = 100  # Number of gradations
    alpha_values = np.linspace(0.1, 0.9, n_colors)  # Alphas from 0.1 to 0.9
    fill_color = fill_color[0]
    colors_with_alpha = [(fill_color[0], fill_color[1], fill_color[2], alpha) for alpha in alpha_values]
    # Define custom colormap
    custom_cmap = LinearSegmentedColormap.from_list("custom_alpha_cmap", colors_with_alpha, N=n_colors)
    # Create a normalization based on the median values
    norm = mpl.colors.Normalize(vmin=min_median, vmax=max_median)
    # Create an inset axis for the colorbar
    ax_pos = ax.get_position()
    fig = plt.gcf()
    cbar_ax = fig.add_axes([ax_pos.x1+0.02, ax_pos.y0, 0.04, ax_pos.height])
    # Create the colorbar using the same normalization as the violin colors
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=custom_cmap), 
                        cax=cbar_ax, orientation='vertical')
    # Only do this if there are at least 8 xticklabels
    if len(xtick_labels) >= 8:
        cbar.ax.set_title('Median', fontsize=10)
    return cbar_ax

# Function to move plots so they don't overlap, moving both in x & y direction
def move_plots_x_y(axs, x_offset=0.02, y_offset=0.04):
    for i, axs_row in enumerate(axs):
        for j, ax in enumerate(axs_row):
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


# --- Also make a colormap per zone for the final violins

# Define some colors per score
color_dict = {
    'GC': '#10e0d9',
    'B_cell_zone': mpl.colors.cnames["royalblue"],
    'T_cell_zone': mpl.colors.cnames["crimson"],
    'HEVs': '#785cae',
}


# --- Plot 1: Violin plots of the proportion of expected cell types per cluster

# Extract a single df
df_combined = pd.concat(decon_dict, ignore_index=False)
df_combined.reset_index(inplace=True, names=['sample', 'barcode'])
df_combined = df_combined.drop('sample', axis=1)

# Calculate proportions per spot
row_sums = df_combined[all_cell_types].sum(axis=1)
proportions = df_combined[all_cell_types].div(row_sums, axis=0)

# Re-add clustering info
proportions['leiden'] = df_combined['leiden'].astype(int)

# Overwrite cluster names according to clusters of interest
proportions['leiden'] = proportions['leiden'].astype(str)  # Ensure leiden is a string for filtering
proportions.loc[~proportions['leiden'].isin(clusters_of_interest.keys()), 'leiden'] = 'Other'

# Also add spot barcodes & sample ids for future reference
proportions[['barcode', 'sample_id']] = df_combined[['barcode', 'sample_id']]

cell_types_of_int = ['T_cells_c0_CD4+_CCR7', 'B cells Naive', 'Endothelial ACKR1', 'Myeloid_c0_DC_LAMP3', 'Cycling', 'PVL Immature s1',]

# Color by expected zone
zone_dict = {
    'Myeloid_c0_DC_LAMP3': 'T_cell_zone',
    'B cells Naive': 'B_cell_zone',
    'Endothelial ACKR1': 'HEVs',
    'T_cells_c0_CD4+_CCR7': 'T_cell_zone',
    'Cycling': 'GC',
    'PVL Immature s1': 'HEVs',
}

# Check if these cell types are in the current annotation level
cell_types_of_int_present = [ct for ct in cell_types_of_int if ct in all_cell_types]

# Plotting order
proportions['leiden'] = proportions['leiden'].astype(str)
order = [str(k) for k in sorted(clusters_of_interest.keys(), key=int)] + ['Other']

if len(cell_types_of_int_present) > 0:
    # Make violin plot
    ncols = 3
    nrows = math.ceil(len(cell_types_of_int_present)/ncols)  # Calculate number of columns based on number of clusters of interest
    axs = []
    cbar_axs = []
    fig = plt.figure(figsize=(3, 1.3))
    for n, cell_type in enumerate(cell_types_of_int_present):
        # Add a subplot for each cell type
        ax = plt.subplot(nrows, ncols, n + 1)
        axs.append(ax)
        cbar_ax = plot_custom_violin_seaborn(proportions, plot_key=cell_type, ax=ax, cut=0, plot_height=1, order=order,
                                                palette=[color_dict[zone_dict[cell_type]]],)
        # Adjust cell type name if in subset_short_lut
        if cell_type in subset_short_lut:
            cell_type_label = subset_short_lut[cell_type]
        else:
            cell_type_label = cell_type
        ax.set_title(cell_type_label, fontweight='bold')
        ax.set_ylabel('Proportion per spot')
        ax.set_xlabel('Leiden clusters')
        cbar_axs.append(cbar_ax)
        # Update the xlabels
        xlabels = ax.get_xticklabels()
        for label in xlabels:
            # Update the label text to include the cluster of interest
            if label.get_text() in clusters_of_interest:
                label.set_text(f'c{label.get_text()}: {clusters_of_interest[label.get_text()]}')
            else:
                label.set_text(label.get_text())  # Keep 'Other' as is
        ax.set_xticklabels([label.get_text() for label in xlabels], rotation=90)
        ax.set_xlabel('')
    # Calculate how many axes per row (may vary if not a complete grid)
    axes_per_row = ncols
    grouped_axs = [axs[i:i+axes_per_row] for i in range(0, len(axs), axes_per_row)]
    x_offset = 0.15
    move_plots_x_y(grouped_axs, x_offset=x_offset, y_offset=0.25, )
    # Keep xticklabels only on bottom row
    for row in grouped_axs[:-1]:
        for ax in row:
            ax.set_xticklabels([])
    # Remove ylabels and set a shared ylabel for all plots
    for col in range(axes_per_row):
        for row in grouped_axs:
            if col < len(row):
                ax = row[col]
                ax.set_ylabel('')
        # Set shared ylabel for the column
        # Only on first column
        if col == 0:
            shared_ylabel = 'Proportion per spot'
            mid_y = (grouped_axs[0][col].get_position().y0 + grouped_axs[-1][col].get_position().y1) / 2
            x_pos = grouped_axs[0][col].get_position().x0 - 0.2
            fig.text(x_pos, mid_y, shared_ylabel, va='center', rotation='vertical', fontsize=10, fontweight='bold')
    # Now move the colorbar axes
    for i, cbar_ax in enumerate(cbar_axs):
        # plot_ax = axs[i]
        # plot_ax.get_position()  # Get the position of the plot ax
        # cbar_ax.set_position([plot_ax.get_position().x1 + 0.01, plot_ax.get_position().y0, 0.04, plot_ax.get_position().height])
        # Set cbar_ax invisible to save space
        cbar_ax.set_visible(False)
    title_x = (grouped_axs[0][0].get_position().x0 + grouped_axs[0][-1].get_position().x1)/2
    title_y = grouped_axs[0][0].get_position().y1 + 0.5
    plt.suptitle(f'All samples - Leiden clustering of TLS-high spots (cellular composition (annotation level={level})\nProportion per spot', y=title_y, x=title_x, fontweight='bold',
                fontsize=12)


# --- Save output

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_2E.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')

# --- Finish

print("Script successfully completed")