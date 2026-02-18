#!/usr/bin/env python3

# --- Aim
# Score spatial transcriptomics spots for germinal center and other immune niche gene signatures.
# Integrate TLS model and spatial groupings, test multiple gene cutoffs, and visualize scores across spatial groups.
# Save scoring results, marker gene dotplots, and gene usage for downstream analysis.

# --- Define environment
# env: spatial-pilot

# --- Load packages
import argparse
import scanpy as sc
import anndata as ad
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import date
from matplotlib import patheffects

from datetime import datetime

from matplotlib.colors import LinearSegmentedColormap
from mycolorpy import colorlist as mcp

# Import custom functions
from sophiesfunctions.auto_crop import auto_crop_scanpy  # import my custom functions that are packaged locally as sophiesfunctions
from sophiesfunctions.save_multi_image import save_multi_image

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium data & do basic QC",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required

parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_id", help="Sample ID")

parser.add_argument("--min_count", help="Min count for QC", type=int)
parser.add_argument("--min_gene", help="Min gene for QC", type=int)
parser.add_argument("--min_spots", help="Min number of spots each gene is expressed in", type=int)

args = vars(parser.parse_args())


# --- Arguments

res_dir = args["res_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
project_dir=args["project_dir"]


# --- Make compatible with lymph node samples

run_lymph_node = False


# --- Get some info from the sample sheet

# Extract the external_id
sample_dir = os.path.join(project_dir, 'config')
if run_lymph_node:
    sample_sheet_file = 'sample_sheet_lymph_nodes.csv'
else:
    sample_sheet_file = 'sample_sheet.csv'

sample_sheet = pd.read_csv(os.path.join(sample_dir, sample_sheet_file))
external_id = sample_sheet[sample_sheet['section_id'] == sample_id]['external_id'].values[0]
tissue_id = sample_sheet[sample_sheet['section_id'] == sample_id]['tissue_id'].values[0]
slide_size = sample_sheet[sample_sheet['section_id'] == sample_id]['slide_size'].values[0]

joined_id = f'{external_id}-{sample_id}'

# --- Also get some tissue metadata

resource_dir = os.path.join(project_dir, 'resources')
metadata_file = 'tissue_metadata.csv'

tissue_metadata = pd.read_csv(os.path.join(resource_dir, metadata_file))

# Extract the subtype
if run_lymph_node:
    subtype = 'lymph_node'
else:
    subtype = tissue_metadata[tissue_metadata['tissue_id'] == str(tissue_id)]['clinical_subtype'].values[0]

# --- Load the adata object

use_data = 'log_norm'

if run_lymph_node:
    adata_file = f'{sample_id}_{use_data}.h5ad'
else:
    adata_file = f'{joined_id}_{use_data}.h5ad'

adata = ad.read_h5ad(os.path.join(h5ad_dir, use_data, adata_file))

# Specify which layer to use
use_layer = 'log_norm'


# --- Calculate crop coords for optimal presentation
coords = auto_crop_scanpy(adata, border_size=0.1)


# --- Set spot size for large & small capture area
# Plotting the spots larger to make it easier to see

if slide_size == '11 mm':
    spot_size = 1.8
else:
    spot_size = 1.3


# --- Load the TLS model results

TLS_res_dir = os.path.join(res_dir, 'TLS_model')

if run_lymph_node:
    TLS_res_dir = os.path.join(TLS_res_dir, 'lymph_nodes')

# Find most recent results
subdirs = [d for d in os.listdir(TLS_res_dir) if os.path.isdir(os.path.join(TLS_res_dir, d))]

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
    most_recent_dir_path = os.path.join(TLS_res_dir, most_recent_dir)
    print("Most recent results:", most_recent_dir_path)
else:
    print("No valid directories found.")

# Read in the results & add to anndata

TLS_file_name = f'{sample_id}_TLS_model.csv'
TLS_results = pd.read_csv(os.path.join(most_recent_dir_path, TLS_file_name), index_col=0)
adata.obs['TLS_model'] = TLS_results


# --- Also load the group results

threshold = 0.01
print(f'Using threshold {threshold}')

threshold_adj = str(threshold).split('.')[-1]

TLS_group_dir = os.path.join(TLS_res_dir, 'TLS_plots_per_sample', f'threshold_{threshold_adj}', 'csv_files')

# Find most recent results
subdirs = [d for d in os.listdir(TLS_group_dir) if os.path.isdir(os.path.join(TLS_group_dir, d))]

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
    most_recent_dir_path = os.path.join(TLS_group_dir, most_recent_dir)
    print("Most recent results:", most_recent_dir_path)
else:
    print("No valid directories found.")

# Read in the results & add to anndata
TLS_file_name = f'{sample_id}_TLS_groups.csv'
TLS_results = pd.read_csv(os.path.join(most_recent_dir_path, TLS_file_name), index_col=0)
adata.obs[TLS_results.columns] = TLS_results


# --- Define some genes

# -- Load the genes from a file

gene_dir = os.path.join(project_dir, 'resources', 'gene_lists')

# --- Assess multiple scores

clusters_of_i = {
    'GC': 'germinal_center_scoring/GC_score_filtered_genes/GC_score_filtered_genes.csv',
    'B_cell_zone': 'other_scores/B_cell_zone_score.csv',
    'T_cell_zone': 'other_scores/T_cell_zone_score.csv',
}

for score_name, filename in clusters_of_i.items():

    gene_table = pd.read_csv(os.path.join(gene_dir, filename))

    # Extract the genes
    if 'gene' in gene_table.columns:
        gene_list = list(gene_table.gene)
    elif 'genes' in gene_table.columns:
        gene_list = list(gene_table.genes)
    else:
        gene_list = list(gene_table.iloc[:, 0]) # Assuming the genes are in the first column of this table

    # Remove any potential nan values
    gene_list = [x for x in gene_list if isinstance(x, str) or not np.isnan(x)]

    # Select only the genes that are present in the anndata objects
    to_plot = np.intersect1d(gene_list , adata.var.index).tolist()

    missing = list(set(gene_list) - set(to_plot))

    print(f'Genes not found: {missing}') 

    # # Divide genes into groups based on cell_type
    # gene_dict = {}

    # group_key = 'source'

    # gene_dict['IHC markers'] = list(gene_table[gene_table[group_key] == 'IHC marker'].gene.values)

    # gene_dict['Potential GC'] = list(gene_table[gene_table[group_key] != 'IHC marker'].gene.values)

    # # Also remove any missing genes
    # gene_dict = {key: [value for value in values if value not in missing] for key, values in gene_dict.items()}

    gene_dict = {f'score_genes': gene_list}


    # --- Test different cutoffs

    # How many genes is optimal?
    cutoffs = [10, 15, 20, 25, 30, 50, 75, 100, 200, 300]

    all_genes = [value for sublist in gene_dict.values() for value in sublist]

    # Set seed
    seed = 42

    # Extract the layer you want to score on

    score_obj = ad.AnnData(adata.layers[use_layer])
    score_obj.obs_names = adata.obs_names
    score_obj.var_names = adata.var_names

    for cutoff in cutoffs:
        # Get top x genes
        to_score = list(all_genes[0:cutoff])
        print(f'Top {cutoff} genes: {to_score}')
        # Save per cutoff
        score_key = f'{score_name}_score_{cutoff}'
        sc.tl.score_genes(score_obj, to_score, random_state=seed, score_name=score_key)
        # Put the score back in the original adata
        adata.obs[score_key] = score_obj.obs[score_key]


    # --- Prepare for plotting

    # Subset anndata to TLS score above threshold
    adata_sub = adata[adata.obs['TLS_model'] > threshold]

    adata.obs['group_singleton_removed'] = adata_sub.obs['group_singleton_removed']

    all_groups = adata_sub.obs['group_singleton_removed'].astype('category').unique()
    all_groups = all_groups.categories

    # Drop 'singleton' from the list
    all_groups = [x for x in all_groups if x != 'singleton']

    # Sort the groups by number
    sorted_groups = sorted(all_groups, key=lambda x: int(x[5:]))  # Sort by the number part of the group name
    print(sorted_groups)

    adata_sub.obs['group_singleton_removed'] = adata_sub.obs['group_singleton_removed'].astype('category').cat.reorder_categories(sorted_groups, ordered=True)


    # --- Plot 1: TLS score & groups

    # Define custom colormap for plotting TLS model results
    # Set everything under certain threshold as the lowest color of a color map

    cmap = "YlOrRd"
    color_list = mcp.gen_color(cmap=cmap, n=11)
    color_list_extracted = ['lightgray'] + color_list[2:]

    n_bins = 100
    cmap_name = 'my_list'
    cmap_adj = LinearSegmentedColormap.from_list(cmap_name, color_list_extracted, N=n_bins)

    cmap_adj.set_under(color='lightgray', alpha=0.4)

    # Plotting
    sc.set_figure_params(fontsize=16)

    nrows = 1
    ncols = 2

    fig, axs = plt.subplots(nrows, ncols, figsize=(8*ncols, 6.5*nrows))

    # TLS model results
    sc.pl.spatial(adata, color='TLS_model', size=spot_size, cmap=cmap_adj, bw=False, alpha_img=1, crop_coord=coords, title=f'TLS score', vmin=0, ax=axs[0])

    # Spots above threshold
    sc.pl.spatial(adata_sub, color='group_singleton_removed', size=spot_size, bw=False, alpha_img=1, crop_coord=coords, title=f'Spots >{threshold}', vmin=0, 
                ax=axs[1])

    plt.suptitle(f'{sample_id} ({subtype}) - TLS score + TLS-high spots grouped on spatial continuity', y=1.05, fontweight='bold', x=0.5)


    # --- Plot 2: Plot the score

    nrows = 1
    ncols = 3

    for cutoff in cutoffs:

        score_key = f'{score_name}_score_{cutoff}'

        fig, axs = plt.subplots(nrows, ncols, figsize=(8*ncols, 6.5*nrows))

        # Scoring results
        sc.pl.spatial(adata, color=score_key, size=spot_size, cmap=cmap_adj, bw=False, alpha_img=1, crop_coord=coords, title=f'{score_key}', vmin=0, ax=axs[0])

        # Scoring results w/ vmax
        vmax=0.8
        sc.pl.spatial(adata_sub, color=score_key, size=spot_size, bw=False, alpha_img=1, crop_coord=coords, title=f'{score_key} in TLS>{threshold} (vmax={vmax})', vmin=0, 
                    vmax=vmax, ax=axs[1], cmap='viridis')


        # Spots above threshold
        sc.pl.spatial(adata_sub, color='group_singleton_removed', size=spot_size, bw=False, alpha_img=0.5, crop_coord=coords, title=f'Spots >{threshold}', vmin=0, 
                    ax=axs[2],)

        plt.suptitle(f'{sample_id} ({subtype}) - {score_name} score + TLS-high spots grouped on spatial continuity', y=1.05, fontweight='bold', x=0.5)

        # Annotate the spatial groups plot
        scalef = adata.uns['spatial'][sample_id]['scalefactors']['tissue_hires_scalef']

        # Find the midpoint for each group
        group_midpoint = {}

        for group in all_groups:
            adata_tmp = adata[adata.obs['group_singleton_removed'] == group]
            coords_tmp = auto_crop_scanpy(adata_tmp)
            # Need to muliply by hires scalefactor to get the plot coordinates
            coords_adj = tuple(x * scalef for x in coords_tmp)
            x1, x2, y1, y2 = coords_adj
            midpoint_x = (x1 + x2) / 2
            midpoint_y = (y1 + y2) / 2
            group_midpoint[group] = (midpoint_x, midpoint_y)

        # Annotate the midpoint
        ax = axs[2]

        for group in all_groups:
            group_adj = group[5:]
            ax.annotate(group_adj, xy=group_midpoint[group],ha='center', va='center', fontsize=8,
                            path_effects=[patheffects.withStroke(linewidth=2, foreground='white')])


    # --- Plot 3: Genes per TLS group dotplot

    sc.set_figure_params(fontsize=10)

    # Make the dotplot
    dp = sc.pl.dotplot(adata_sub, var_names=gene_dict, groupby=['group_singleton_removed'], layer=use_layer, 
                        swap_axes=True, return_fig=True)

    dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5, cmap='Reds').show()

    axs = plt.gcf().get_axes()
    # Get the legend ax & adjust the legend
    dots_ax = axs[4]
    dots_ax.set_title('Fraction of spots\nin group (%)', y=0.5, fontsize=8.5)

    # Rotate xlabels
    main_ax = axs[0]
    xlabels = main_ax.get_xticklabels()
    main_ax.set_xticklabels(xlabels, rotation=45, ha='right')

    # Get the top_ax to add a title
    top_ax = axs[1]
    top_ax.set_title(f'TLS score >{threshold}', fontweight='bold', fontsize=9.5)

    plt.suptitle(f'{sample_id} ({subtype}) - TLS groups (threshold={threshold}) marker genes (use_data={use_data})', y=1, fontweight='bold')


    # --- Plot 4: Score per TLS group dotplot

    # Make a fake anndata from the score

    cat_key = 'group_singleton_removed'

    df = adata.obs[[score_key, cat_key]].dropna()

    fake_adata = ad.AnnData(df[[score_key]])
    fake_adata.obs_names = df.index
    fake_adata.var_names = [score_key]

    fake_adata.obs[cat_key] = df[cat_key]

    # Only do the plot below if there's more than 1 group otherwise it breaks stuff

    if len(all_groups) > 1:

        # Use fake adata to make a dotplot

        # Make the dotplot
        cutoff_dotplot = 0.4
        dp = sc.pl.dotplot(fake_adata, var_names=[score_key], groupby=[cat_key], swap_axes=True, return_fig=True, expression_cutoff=cutoff_dotplot, mean_only_expressed=True, vmin=cutoff_dotplot,
                        vmax=0.8)

        dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5, cmap='Reds').show()

        axs = plt.gcf().get_axes()
        # Get the legend ax & adjust the legend
        dots_ax = axs[3]
        dots_ax.set_title('Fraction of spots\nin group (%)', y=0.5, fontsize=8.5)

        # Rotate xlabels
        main_ax = axs[0]
        xlabels = main_ax.get_xticklabels()
        main_ax.set_xticklabels(xlabels, rotation=45, ha='right')

        # Get the top_ax to add a title
        top_ax = axs[1]
        top_ax.set_title(f'GC score (cutoff={cutoff_dotplot})', fontweight='bold', fontsize=9.5)

        plt.suptitle(f'{sample_id} ({subtype}) - {score_key} in TLS groups (threshold={threshold}, use_data={use_data}, cutoff={cutoff_dotplot})', y=1, fontweight='bold')


    # --- Save output

    if score_name == 'GC':
        out_name = 'GC_score_filtered'
    else:
        out_name = f'{score_name}_score'

    # Put in a timestamped folder to avoid overwriting older plots
    qc_dir = os.path.join(res_dir, 'TLS_model', 'gene_expression_plots', 'potential_GC_genes', 'scoring')

    if run_lymph_node:
        qc_dir = os.path.join(qc_dir, 'lymph_nodes')

    os.makedirs(qc_dir, exist_ok=True)

    filename = f"{sample_id}_{out_name}.pdf"
    filename = filename.replace('/', '_')  # Replace / with underscores for file naming

    today = date.today()
    today = today.strftime("%Y%m%d")

    out_dir = os.path.join(qc_dir, out_name, 'per_sample', 'plots', today)
    os.makedirs(out_dir, exist_ok=True)

    print(f'Saving plot as {filename} in {out_dir}')
    save_multi_image(os.path.join(out_dir, filename))
    plt.close('all')

    # Also save the scores + which genes were used
    out_dir = os.path.join(qc_dir, out_name, 'per_sample', 'csv_files', today)
    os.makedirs(out_dir, exist_ok=True)

    score_keys = [f'{score_name}_score_' + str(item) for item in cutoffs]
    save_score = adata.obs[score_keys]
    filename = f'{sample_id}_{out_name}.csv'
    filename = filename.replace('/', '_')  # Replace / with underscores for file naming

    save_score.to_csv(os.path.join(out_dir, filename))

    # Save genes used
    # Also save the pandas dataframe as a csv file (to be able to check which genes were used in each sample)
    genes_used = pd.DataFrame.from_dict(gene_dict, orient='index').transpose()

    genes_used.to_csv(os.path.join(out_dir, f'{sample_id}_genes_used.csv'), index=False)


print("Script succesfully completed")
