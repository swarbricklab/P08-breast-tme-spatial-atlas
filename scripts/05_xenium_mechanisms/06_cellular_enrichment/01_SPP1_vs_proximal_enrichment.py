#!/usr/bin/env python3

# --- Aim
# Quantify the enrichment of cell types within a defined distance of HEV niches in Xenium data.
# Correlate the proportion of proximal cells with SPP1 expression in the HEV niche.
# Reveal spatial relationships between SPP1 signaling and cellular microenvironment.

# --- Define environment
# env: xen-hev

# --- Load packages
import argparse
import os
import sys

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize

from scipy.stats import pearsonr
from scipy.spatial import cKDTree

import pandas as pd
import numpy as np
import scanpy as sc

import seaborn as sns

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium data & do basic QC",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required
parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_id", help="Sample ID")

parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)

args = vars(parser.parse_args())


# --- Arguments

res_dir = args["res_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
project_dir=args["project_dir"]

# Import function script as a module
sys.path.append(os.path.join(project_dir, 'scripts', '05_xenium_mechanisms'))
import plotting_functions

sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import auto_crop_xenium, save_multi_image


# --- Get some info from the sample sheet

# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_dir = os.path.join(project_dir, 'config')
    sample_sheet_file = 'sample_sheet_xenium.csv'
    sample_sheet_path = os.path.join(sample_dir, sample_sheet_file)

sample_sheet = pd.read_csv(sample_sheet_path)
tissue_id = sample_sheet[sample_sheet['section_id'] == sample_id]['tissue_id'].values[0]


# --- Also get some tissue metadata

resource_dir = os.path.join(project_dir, 'resources')
metadata_file = 'tissue_metadata.csv'

tissue_metadata = pd.read_csv(os.path.join(resource_dir, metadata_file))

# Create a sample - subtype dictionary
subtype = tissue_metadata[tissue_metadata['tissue_id'] == str(tissue_id)]['clinical_subtype'].values[0]


# --- Load data

adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample_id}.h5ad'))


# --- Calculate crop coordinates for optimal presentation

crop_coord = auto_crop_xenium(adata)


# --- Load the HEV niches

hev_dir = os.path.join(res_dir, '05_xenium_mechanisms', '02_identify_HEV_regions',)
most_recent = hev_dir

hev_file = f'{sample_id}_HEV_niches_cell_to_niche.csv'
tmp_df = pd.read_csv(os.path.join(most_recent, 'csv_files', hev_file), index_col=0)
adata.obs['niche_id'] = tmp_df['niche_id']
adata.obs['HEV_score'] = tmp_df['HEV_score']

# Also load info whether HEV is near TLS or not
filename = f'{sample_id}_HEV_niches_info.csv'
tmp_df = pd.read_csv(os.path.join(most_recent, 'csv_files', filename),)

# Make a dict close to TLS vs not
near_tls_dict = dict(zip(tmp_df['niche_id'], tmp_df['close_to_TLS']))

# Map to adata
adata.obs['HEV_in_TLS'] = adata.obs['niche_id'].map(near_tls_dict)


# --- Also load the TLS annotations

tls_dir = os.path.join(res_dir, '05_xenium_mechanisms', '01_identify_TLS_regions',)

tls_file = f'{sample_id}_cell_TLS_annotations.csv'

tmp_df = pd.read_csv(os.path.join(tls_dir, 'csv_files', tls_file), index_col=0)
adata.obs['in_TLS'] = tmp_df['in_TLS']


# --- Load the cell type annotations

# Decide which annotation levels to use
level = 'subset'

annotation_dir = os.path.join(res_dir, '02_cell_type_mapping', 'label_transfer', 'xenium', f'{level}_level')
annotation_file = f'{sample_id}_{level}_pred_cell_types.csv'

annotations = pd.read_csv(os.path.join(annotation_dir, annotation_file), index_col=0)

# Add to adata
adata.obs[level] = annotations['pred_cell_type']


# --- Load plotting colors

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')

# Use the major level for colors
color_level = 'major'
color_file = f'miniatlas_{color_level}_colors.csv'

color_dict = pd.read_csv(os.path.join(color_dir, color_file), index_col=0).to_dict()['color']


# --- Load the miniatlas classifications

# These are needed to check the alignment between major and minor levels
miniatlas_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
miniatlas_file = 'miniatlas_subset_classifications.csv'

miniatlas_df = pd.read_csv(os.path.join(miniatlas_dir, miniatlas_file),)

# For each cell type in the miniatlas at this level, check if it is in the color_dict
# If not, get the corresponding major cell type and use its color
for cell_type in miniatlas_df[f'miniatlas_{level}'].unique():
    if cell_type not in color_dict.keys():
        # Check that the value is not nan
        if pd.isna(cell_type):
            continue
        print(f'Adding color for {cell_type}')
    # Get the corresponding major cell type
    major_type = miniatlas_df[miniatlas_df[f'miniatlas_{level}'] == cell_type][f'miniatlas_major'].values
    # Find which value is in the color_dict
    major_type_in_dict = [m for m in major_type if m in color_dict.keys()]
    if len(major_type_in_dict) > 0:
        color_dict[cell_type] = color_dict[major_type_in_dict[0]]
    else:
        color_dict[cell_type] = 'lightgrey'

# Also define a dict with shorter names for the subset level cell types
miniatlas_sub_dict = miniatlas_df.set_index('miniatlas_subset')[['miniatlas_subset_short']].to_dict()


# --- Sort by lineage, then by major

lineage_order = ['Epithelial', 'Immune', 'Stromal']

major_cell_types = miniatlas_df[f'miniatlas_major'].dropna().unique().tolist()
minor_cell_types = miniatlas_df[f'miniatlas_minor'].dropna().unique().tolist()
subset_cell_types = miniatlas_df[f'miniatlas_subset'].dropna().unique().tolist()

# Sort based on lineage first (recorded in the miniatlas df)
def get_lineage(cell_type, level='major'):
    if cell_type in miniatlas_df[f'miniatlas_{level}'].values:
        lineage = miniatlas_df[miniatlas_df[f'miniatlas_{level}'] == cell_type]['lineage'].values[0]
    else:
        lineage = 'other'
    return lineage

major_cell_types_sorted = sorted(major_cell_types, key=lambda x: (lineage_order.index(get_lineage(x)) if get_lineage(x) in lineage_order else len(lineage_order), x))


# --- For each HEV group, find barcodes of all cells within x um

# Define max distance
max_dists = [200]

# Pick a gene of interest
gene_of_interest = 'SPP1'

for max_dist in max_dists:

    # --- Plot 1: Plot the identified HEV niches

    # Set fontsize
    sc.set_figure_params(fontsize=20)

    ncols = 3
    fig, axs = plt.subplots(nrows=1, ncols=ncols, figsize=(9*ncols, 7))

    # -- First plot just the H&E
    n = 0
    ax = axs[n]

    # Create an emtpy adata to speed things up
    adata_sub = adata[[]].copy()
    sc.pl.spatial(adata_sub, color=None, size=0.5, title=f'{sample_id} ({subtype})',
                crop_coord=crop_coord, ax=ax, na_in_legend=False, show=False)

    ax.set_title(ax.get_title(), fontweight='bold')

    # -- Then plot the HEV niches on top

    n = 1
    ax = axs[n]

    # Subset adata to only spots in niches
    adata_sub = adata[adata.obs['niche_id'].notna()]
    sc.pl.spatial(adata_sub, color='niche_id', size=0.5, title=f'HEV niches', palette=['blue'],
                crop_coord=crop_coord, ax=ax, na_in_legend=False, legend_loc=None, show=False,
                alpha_img=0.5)

    ax.set_title(ax.get_title(), fontweight='bold')

    # -- Then plot HEVs + TLS niches
    n = 2
    ax = axs[n]

    # Subset adata to only TLS spots
    adata_sub = adata[adata.obs['in_TLS'] == True]
    adata_sub.obs['in_TLS'] = adata_sub.obs['in_TLS'].astype('category')
    sc.pl.spatial(adata_sub, color='in_TLS', size=0.5, title=f'HEV niches in TLS', palette=['cornflowerblue'],
                crop_coord=crop_coord, ax=ax, na_in_legend=False, legend_loc=None, show=False,
                alpha_img=0.7)

    adata_sub = adata[adata.obs['niche_id'].notna()]
    sc.pl.spatial(adata_sub, color='niche_id', size=0.5, title=f'HEV + TLS niches', palette=['blue'],
                crop_coord=crop_coord, ax=ax, na_in_legend=False, legend_loc=None, show=False,
                alpha_img=0)       

    ax.set_title(ax.get_title(), fontweight='bold')

    # Add legend with TLS & HEV
    # Make the markers circles
    color_dict_TLS = {
        'TLS': 'cornflowerblue',
        'HEV niche': 'blue',
    }

    # Create custom legend handles with circles
    handles = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=label)
        for label, color in color_dict_TLS.items()
    ]

    # Remove all x & ylabels
    for ax in axs:
        ax.set_xlabel('')
        ax.set_ylabel('')

    # Adjust the positioning of the axes
    for n, ax in enumerate(axs):
        ax_pos = ax.get_position()
        ax_width = ax_pos.width
        if n == 0:
            ax_pos_prev = ax_pos
        else:
            new_x = ax_pos_prev.x1 + ax_width*0.1
            ax.set_position([new_x, ax_pos.y0, ax_pos.width, ax_pos.height])
            ax_pos_prev = ax.get_position()

    plt.legend(handles=handles, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=16, frameon=False)

    # Add a title
    plt.suptitle(f'{sample_id} ({subtype}) - HEV niches (n={adata_sub.obs["niche_id"].nunique()}) identified based on HEV scoring', y=1.05, fontweight='bold')


    # --- Extract the coordinates in um

    coords = adata.obs[['x_centroid', 'y_centroid']].copy()
    # Rename columns
    coords.columns = ['x', 'y']
    groups = adata.obs['niche_id']

    # Convert to DataFrame for easier manipulation
    df = pd.DataFrame(coords, columns=["x", "y"])
    df["HEV_group"] = groups.values

    # Add whether cell is inside TLS
    df['in_TLS'] = adata.obs['in_TLS']

    len(df['HEV_group'].unique())

    def find_cells_within_distance_min(df, group_key="HEV_group", x_col="x", y_col="y", max_distance=200):
        """
        For each HEV group, find all cell barcodes within max_distance (µm) of any cell assigned to that group.
        The distance recorded is the minimum distance from any cell in the HEV niche.
        Returns a DataFrame with columns: cell_barcode, HEV_group, distance
        """
        results = []
        # Get all cells assigned to a HEV group (the niche members)
        group_df = df.dropna(subset=[group_key])
        # For each HEV group
        for group in group_df[group_key].unique():
            niche_cells = group_df[group_df[group_key] == group]
            niche_coords = niche_cells[[x_col, y_col]].values
            # Build KDTree for the niche
            tree = cKDTree(niche_coords)
            # For all cells, query the minimum distance to any cell in the niche
            all_coords = df[[x_col, y_col]].values
            min_distances, _ = tree.query(all_coords)
            # Find all cells within max_distance
            within_mask = min_distances <= max_distance
            for idx in np.where(within_mask)[0]:
                cell_barcode = df.index[idx]
                results.append({
                    "cell_barcode": cell_barcode,
                    "HEV_group": group,
                    "distance": min_distances[idx]
                })
        return pd.DataFrame(results)

    df_results = find_cells_within_distance_min(df, group_key="HEV_group", x_col="x", y_col="y", max_distance=max_dist)


    # --- Plot 2: Spatial overview of 5 HEV groups as a test

    num_groups = 5

    fig, axs = plt.subplots(nrows=1, ncols=5, figsize=(8*num_groups, 8))

    hev_groups_to_plot = df_results['HEV_group'].unique()[:5]

    for n, hev_group in enumerate(hev_groups_to_plot):
        ax = axs[n]
        close_barcodes = df_results[df_results['HEV_group'] == hev_group]['cell_barcode']
        adata_sub = adata[close_barcodes]
        adata_sub.obs['cell_id'] = adata_sub.obs['cell_id'].astype('category')
        sc.pl.spatial(adata_sub, color='cell_id', size=0.5, title=f'HEV niches in TLS', palette=['cornflowerblue'],
                    crop_coord=crop_coord, ax=ax, na_in_legend=False, legend_loc=None, show=False,
                    alpha_img=0.7)
        adata_sub = adata[adata.obs['niche_id'] == hev_group]
        sc.pl.spatial(adata_sub, color='niche_id', size=0.5, title=f'HEV + TLS niches', palette=['blue'],
                    crop_coord=crop_coord, ax=ax, na_in_legend=False, legend_loc=None, show=False,
                    alpha_img=0)
        ax.set_title(f'{hev_group}\nCells within {max_dist} µm', fontweight='bold')

    # Remove all x & ylabels
    for ax in axs:
        ax.set_xlabel('')
        ax.set_ylabel('')

    plt.suptitle(f'{sample_id} ({subtype}) - Cells within {max_dist} µm of HEV niches', y=1.05, fontweight='bold')


    # --- Calculate total number of cells within 100 um for each HEV group

    mapping_dict = annotations['pred_cell_type'].to_dict()

    df_results['cell_type'] = df_results['cell_barcode'].map(mapping_dict)

    # Exclude any barcodes that are part of the HEV niche itself
    df_results = df_results[~df_results['cell_barcode'].isin(adata[adata.obs['niche_id'].notna()].obs_names)].copy()

    # Make a table with counts per HEV group
    count_table = pd.crosstab(df_results['HEV_group'], df_results['cell_type'])
    all_cell_types = count_table.columns.tolist()

    # Also calculate as proportions
    prop_table = count_table.div(count_table.sum(axis=1), axis=0)


    # --- Now calculate average SPP1 expression for endothelial cells per HEV group

    # Subset adata object to ACKR1+ endothelial cells within HEVs
    adata_sub = adata[adata.obs[level] == 'Endothelial ACKR1'].copy()
    adata_sub = adata_sub[adata_sub.obs['niche_id'].notna()].copy()

    # Calculate average SPP1 expression per HEV group
    spp1_avg = adata_sub.to_df()[[gene_of_interest]].groupby(adata_sub.obs['niche_id']).mean()
    spp1_avg.columns = [f'avg_{gene_of_interest}_expression_ACKR1']

    spp1_dict = spp1_avg[f'avg_{gene_of_interest}_expression_ACKR1'].to_dict()

    # Add to adata for plotting
    adata.obs[f'avg_{gene_of_interest}_expression_ACKR1'] = adata.obs['niche_id'].map(spp1_dict)

    # Also calculate average expression across all cell types (not just endothelial ACKR1 in case imperfect annotation)
    adata_sub = adata[adata.obs['niche_id'].notna()].copy()
    spp1_avg_all = adata_sub.to_df()[[gene_of_interest]].groupby(adata_sub.obs['niche_id']).mean()
    spp1_avg_all.columns = [f'avg_{gene_of_interest}_expression_all_cells']

    # Merge with spp1_avg
    spp1_avg = spp1_avg.merge(spp1_avg_all, left_index=True, right_index=True)


    # --- Load inside vs outside TLS colors

    color_file = 'TLS_HEV_colors.csv'
    color_df = pd.read_csv(os.path.join(resource_dir, 'plotting_colors', color_file), index_col=0)

    color_dict_HEV = color_df['color'].to_dict()
    color_dict_HEV = {True: color_dict_HEV['near_TLS_HEV'], False: color_dict_HEV['non_TLS_HEV']}


    # --- Plot 3: Plot SPP1 expression for HEVs inside vs outside TLS

    # Add info whether HEV is close to TLS or not
    spp1_avg['HEV_in_TLS'] = spp1_avg.index.map(near_tls_dict)

    ncols = 2

    # Properties for boxplot
    whisker_props = {'color': 'black'}
    median_props = {'color': 'black'}

    # Set fontsize
    sc.set_figure_params(fontsize=10)

    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(3, 2))

    # First plot ACKR1 only
    sns.stripplot(data=spp1_avg, x='HEV_in_TLS', y=f'avg_{gene_of_interest}_expression_ACKR1', color='black', jitter=True, size=2,
                    marker='o', ax=axs[0],)

    sns.boxplot(data=spp1_avg, x='HEV_in_TLS', y=f'avg_{gene_of_interest}_expression_ACKR1', ax=axs[0], showfliers=False, linewidth=1, linecolor='black',
                whiskerprops=whisker_props, medianprops=median_props, saturation=1, zorder=10, palette=color_dict_HEV, hue='HEV_in_TLS', dodge=False,
                legend=False)

    axs[0].set_title(f'HEV {gene_of_interest}\n(ACKR1+ endothelial)', fontweight='bold')

    # Then plot all cells
    sns.stripplot(data=spp1_avg, x='HEV_in_TLS', y=f'avg_{gene_of_interest}_expression_all_cells', color='black', jitter=True, size=2,
                    marker='o', ax=axs[1],)

    sns.boxplot(data=spp1_avg, x='HEV_in_TLS', y=f'avg_{gene_of_interest}_expression_all_cells', ax=axs[1], showfliers=False, linewidth=1, linecolor='black',
                    whiskerprops=whisker_props, medianprops=median_props, saturation=1, zorder=10, palette=color_dict_HEV, hue='HEV_in_TLS', dodge=False,
                    legend=False)

    axs[1].set_title(f'HEV {gene_of_interest}\n(all cells)', fontweight='bold')

    ymax = 0

    for ax in axs:
        plotting_functions.beautify_plot(ax)
        ax.set_ylabel(f'{gene_of_interest} expression\n(Average)')
        ymax = max(ymax, ax.get_ylim()[1])

    # Set same ylims
    for ax in axs:
        ax.set_ylim(0, ymax)

    plt.suptitle(f'{sample_id} ({subtype}) - Average {gene_of_interest} Expression by HEV Proximity to TLS', fontweight='bold', y=1.2)

    # Increase the distance between the two plots
    plt.subplots_adjust(wspace=1.5)


    # --- Plot 4: Plot SPP1 expression vs number of LAMP3+ DCs within max dist

    # Add cell numbers to spp1_avg
    for cell_type in all_cell_types:
        spp1_avg[f'prop_{cell_type}_within_{max_dist}um'] = prop_table[cell_type]

    # Calculate pearson correlation for each cell type
    corr = {}
    pval = {}

    # Calculate pearson correlation for each cell type
    for cell_type in all_cell_types:
        corr[cell_type], pval[cell_type] = pearsonr(spp1_avg[f'avg_{gene_of_interest}_expression_ACKR1'], spp1_avg[f'prop_{cell_type}_within_{max_dist}um'].fillna(0))

    # Set grid behind
    plt.rcParams['axes.axisbelow'] = True

    ncols = 5
    nrows = int(np.ceil(len(all_cell_types) / ncols))

    plt.figure(figsize=(5, 5))

    # add a new subplot iteratively
    for n, cell_type in enumerate(all_cell_types):
        ax = plt.subplot(nrows, ncols, n + 1)
        sns.scatterplot(data=spp1_avg, x=f'avg_{gene_of_interest}_expression_ACKR1', y=f'prop_{cell_type}_within_{max_dist}um', s=20, ax=ax, alpha=0.7,
                        edgecolor=None)
        ax.set_xlabel(f'HEV {gene_of_interest} expression\n(Average)')
        ax.set_ylabel(f'Prop. of cells\n(<{max_dist} µm from HEV)')
        cell_type_short = miniatlas_sub_dict['miniatlas_subset_short'].get(cell_type, cell_type)
        ax.set_title(cell_type_short, fontweight='bold')
        # Calculate pearson correlation
        ax.text(0.05, 0.95, f'Pearson r={corr[cell_type]:.2f}', transform=ax.transAxes,
                verticalalignment='top', fontsize=10, bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

    # Resize the axs
    axs = plt.gcf().get_axes()
    for ax in axs:
        ax_pos = ax.get_position()
        ax.set_position([ax_pos.x0, ax_pos.y0, 0.3, 0.3])

    # Move plots
    # Pad the axs list with None to make it evenly divisible by ncols
    padded_axs = axs + [None] * (nrows * ncols - len(axs))

    # Reshape the padded list
    axs_grid = np.array(padded_axs).reshape(nrows, ncols)
    # list of rows (the last row can be shorter, so not a 2-D array)
    axs_grid_cleaned = [[ax for ax in row if ax is not None] for row in axs_grid if any(ax is not None for ax in row)]

    plotting_functions.move_plots_x_y(axs_grid_cleaned, x_offset=0.2, y_offset=0.2, remove_ylabels=False, remove_xlabels=False, remove_title=False)

    # Add a main title
    title_x = axs_grid_cleaned[0][0].get_position().x0 + (axs_grid_cleaned[0][-1].get_position().x1 - axs_grid_cleaned[0][0].get_position().x0) / 2
    title_y = axs_grid_cleaned[0][0].get_position().y1 + 0.1
    plt.suptitle(f'{sample_id} ({subtype}) - Correlation between HEV {gene_of_interest} expression (ACKR1+ endothelial) and cell type abundance (within {max_dist} µm of HEV)', x=title_x, y=title_y, fontweight='bold')


    # --- Save the output

    out_dir = os.path.join(res_dir, '05_xenium_mechanisms', '06_cellular_enrichment', f'{gene_of_interest}_vs_cell_abundance')
    out_dir = os.path.join(out_dir, f'maxdist_{max_dist}um')

    out_dir_plots = os.path.join(out_dir, 'plots')
    os.makedirs(out_dir_plots, exist_ok=True)

    out_dir_csv = os.path.join(out_dir, 'csv_files')
    os.makedirs(out_dir_csv, exist_ok=True)

    # Save plots
    filename = f'{sample_id}_HEV_{gene_of_interest}_vs_cell_abundance.pdf'

    save_multi_image(os.path.join(out_dir_plots, filename))
    plt.close('all')

    # And save the dfs
    filename = f'{sample_id}_HEV_{gene_of_interest}_and_cell_proportions.csv'
    spp1_avg.to_csv(os.path.join(out_dir_csv, filename), index=True)

    # Also the correlation df
    filename = f'{sample_id}_HEV_{gene_of_interest}_and_cell_proportions_corr.csv'

    # Sort by correlation
    sorted_cell_types = sorted(corr.keys(), key=lambda x: corr[x], reverse=True)

    # Also calculate the average proportion of each cell type across all HEV niches
    avg_props = prop_table.mean().to_dict()

    # Create a DataFrame for plotting
    corr_df = pd.DataFrame({
        'cell_type': sorted_cell_types,
        'pearson_r': [corr[ct] for ct in sorted_cell_types],
        'p_value': [pval[ct] for ct in sorted_cell_types],
        'avg_prop': [avg_props[ct] for ct in sorted_cell_types],
    })

    # Add a column for significance
    corr_df['p<0.05'] = corr_df['p_value'] < 0.05

    # Add a column for -log10(p-value)
    corr_df['-log10(p_value)'] = -np.log10(corr_df['p_value'] + 1e-10)  # add small value to avoid log(0)

    corr_df.to_csv(os.path.join(out_dir_csv, filename), index=False)


# --- Finish

print("Script successfully completed")