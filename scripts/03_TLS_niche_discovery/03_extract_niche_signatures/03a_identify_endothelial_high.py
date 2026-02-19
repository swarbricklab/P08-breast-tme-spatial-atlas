#!/usr/bin/env python3

# --- Aim
# Refine HEV (high endothelial venule) gene signature by filtering out genes broadly expressed in endothelial-dominant spots.
# Integrate clustering and deconvolution results, identify endothelial-high regions, and visualize spatial distributions.
# Support downstream HEV-specific marker discovery.

# --- Define environment
# env: spatial-pilot

# --- Load packages
import argparse
import scanpy as sc
import anndata as ad
import os
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from datetime import date
from datetime import datetime


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
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)

parser.add_argument("--min_count", help="Min count for QC", type=int)
parser.add_argument("--min_gene", help="Min gene for QC", type=int)
parser.add_argument("--min_spots", help="Minimal spots per gene", type=int)

args = vars(parser.parse_args())


# --- Arguments

res_dir = args["res_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
project_dir=args["project_dir"]


# --- Get some info from the sample sheet
# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet.csv')

sample_sheet = pd.read_csv(sample_sheet_path)
external_id = sample_sheet[sample_sheet['section_id'] == sample_id]['external_id'].values[0]
tissue_id = sample_sheet[sample_sheet['section_id'] == sample_id]['tissue_id'].values[0]
slide_size = sample_sheet[sample_sheet['section_id'] == sample_id]['slide_size'].values[0]


# --- Also get some tissue metadata

resource_dir = os.path.join(project_dir, 'resources')
metadata_file = 'tissue_metadata.csv'

tissue_metadata = pd.read_csv(os.path.join(resource_dir, metadata_file))

subtype = tissue_metadata[tissue_metadata['tissue_id'] == str(tissue_id)]['clinical_subtype'].values[0]


# --- Load the adata object


adata = ad.read_h5ad(os.path.join(h5ad_dir, 'visium', f'{sample_id}.h5ad'))


# --- Calculate crop coords for optimal presentation
coords = auto_crop_scanpy(adata, border_size=0.1)


# --- Set spot size for large & small capture area
# Plotting the spots larger to make it easier to see

if slide_size == '11 mm':
    spot_size = 1.8
else:
    spot_size = 1.3


# --- Load clustering results

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


# Load harmony clustering results
threshold = 0.01
threshold_adj = str(threshold).split('.')[-1]

harmony_dir = os.path.join(project_dir, 'results', 'TLS_model', 'TLS_plots_all_samples', 'gene_expression', f'threshold_{threshold_adj}', 'umap', 'harmony_integration')
most_recent = find_most_recent(harmony_dir)

filename = 'All_samples_leiden_clusters.csv'
cluster_res = pd.read_csv(os.path.join(most_recent, 'csv_files', filename), index_col=0)

# Add cluster results to anndata 
results_sub = cluster_res[cluster_res['sample_id'] == sample_id]
adata.obs['leiden'] = results_sub['leiden'].astype('category')


# --- Load the decon results

# Decide which level to plot
level = 'major'

if level == 'minor':
    out_name = '202408_benchmarking'
else:
    out_name = f'{level}_level'

# Decide which n_cells_per_location
N_cells_per_location = 20

decon_dir = os.path.join(project_dir, 'resources', 'deconvolution', 'cell2location', out_name)

# Read in the results & add to anndata
decon_filename = 'q05_cell_abundance_w_sf.csv'

n_cells_key = f'n_cells_{N_cells_per_location}'
most_recent_dir_path = os.path.join(decon_dir, sample_id, n_cells_key)

# Load decon results
decon_results = pd.read_csv(os.path.join(most_recent_dir_path, 'model_output', decon_filename), index_col=0)

cell_types = list(decon_results.columns)

adata.obs[cell_types] = decon_results


# --- Prepare for saving

qc_dir = os.path.join(res_dir, 'TLS_model', 'gene_expression_plots', 'potential_GC_genes', 'scoring', 'HEV_score_v2')

os.makedirs(qc_dir, exist_ok=True)

today = date.today()
today = today.strftime("%Y%m%d")


out_dir = os.path.join(qc_dir, 'try_to_improve_HEV_score', 'per_sample', today)
os.makedirs(out_dir, exist_ok=True)


# --- Define some colors

# Define some colors per score
color_dict = {
    'GC': mpl.colors.cnames["cyan"],
    'intermediate_TLS': mpl.colors.cnames["cornflowerblue"],
    'B_cell_zone': mpl.colors.cnames["cornflowerblue"],
    'T_cell_zone': mpl.colors.cnames["crimson"],
    'mixed': mpl.colors.cnames["indigo"]
}


# --- Plot 1: Plot decon results & clustering results

HEV_cluster = 4
endothelial_key = 'Endothelial'

sc.set_figure_params(scanpy=True, fontsize=18)

nrows = 1
ncols = 3

# Make figure
fig, axs = plt.subplots(nrows, ncols, figsize=(7.5*ncols, 6*nrows))

# First plot H&E
n = 0
ax = axs[n]

sc.pl.spatial(adata, color=None, size=spot_size, ax=ax, bw=False, alpha_img=1,
    crop_coord=coords, title=f'{sample_id} ({subtype})\nH&E',)

ax.set_title(ax.get_title(), fontweight='bold')


# Then plot clustering results
n += 1
ax = axs[n]

sc.pl.spatial(adata, color='leiden', groups=[HEV_cluster], size=spot_size, ax=ax, bw=False, alpha_img=0.7,
    crop_coord=coords, title=f'{sample_id} ({subtype})\nHEV cluster', palette=['#1f77b4'])

ax.set_title(ax.get_title(), fontweight='bold')

# Then plot score classifications
n += 1
ax = axs[n]

sc.pl.spatial(adata, color=endothelial_key, size=spot_size, ax=ax, bw=False, alpha_img=0.7, cmap='magma_r',
    crop_coord=coords, title=f'{sample_id} ({subtype})\nEndothelial cells', vmax='p99.2')

ax.set_title(ax.get_title(), fontweight='bold')

plt.suptitle(f'{sample_id} ({subtype}) - clustering results & decon results (threshold={threshold})', y=1.1, fontweight='bold', x=0.5)


# --- Find endothelial-high spots

total_per_spot = decon_results.sum(axis=1)
decon_prop = decon_results.div(total_per_spot, axis=0)

# Make a new df with just the endothelial cell info
endo_df = pd.DataFrame()
endo_df['endo_cells'] = decon_results[endothelial_key]
endo_df['endo_prop'] = decon_prop[endothelial_key]

# Set some thresholds
quantile = 0.8
cells_thresh = decon_results[endothelial_key].quantile(quantile)
prop_thresh = 0.1

endo_df['pass'] = np.where((endo_df['endo_cells'] > cells_thresh) & (endo_df['endo_prop'] > prop_thresh),
                           'pass', None)

adata.obs['endo_high'] = endo_df['pass']

ncols = 5

# Make figure
fig, axs = plt.subplots(nrows, ncols, figsize=(7.5*ncols, 6*nrows))

# First plot H&E
n = 0
ax = axs[n]

sc.pl.spatial(adata, color=None, size=spot_size, ax=ax, bw=False, alpha_img=1,
    crop_coord=coords, title=f'{sample_id} ({subtype})\nH&E',)

ax.set_title(ax.get_title(), fontweight='bold')

# Then plot clustering results
n += 1
ax = axs[n]

sc.pl.spatial(adata, color='leiden', groups=[HEV_cluster], size=spot_size, ax=ax, bw=False, alpha_img=0.7,
    crop_coord=coords, title=f'{sample_id} ({subtype})\nHEV cluster', palette=['#1f77b4'])

ax.set_title(ax.get_title(), fontweight='bold')

# Then plot endothelial high
n += 1
ax = axs[n]

if len(adata.obs['endo_high'].dropna()) > 0:
    plot_key = 'endo_high'
else:
    plot_key = None

sc.pl.spatial(adata, color=plot_key, size=spot_size, ax=ax, bw=False, alpha_img=0.7,
    crop_coord=coords, title=f'{sample_id} ({subtype})\nEndothelial-high', palette=['#1f77b4'])

ax.set_title(ax.get_title(), fontweight='bold')

# Then plot decon results
n += 1
ax = axs[n]

sc.pl.spatial(adata, color=endothelial_key, size=spot_size, ax=ax, bw=False, alpha_img=0.7, cmap='magma_r',
    crop_coord=coords, title=f'{sample_id} ({subtype})\nEndothelial cells', vmax='p99.2')

ax.set_title(ax.get_title(), fontweight='bold')


# Then plot decon proportion
adata.obs['endothelial_prop'] = endo_df['endo_prop']

n += 1
ax = axs[n]

sc.pl.spatial(adata, color='endothelial_prop', size=spot_size, ax=ax, bw=False, alpha_img=0.7, cmap='magma_r',
    crop_coord=coords, title=f'{sample_id} ({subtype})\nEndothelial cells (proportion)', vmax=None)

ax.set_title(ax.get_title(), fontweight='bold')

plt.suptitle(f'{sample_id} ({subtype}) - clustering results & decon results (threshold={threshold}, endothelial-high: min_prop={prop_thresh}, min_cells={cells_thresh.round(2)} (quantile {quantile}))', y=1.1, fontweight='bold', x=0.5)


# --- Only run the rest if there are endothelial-high spots / HEV spots

adata.obs['endo_sub'] = np.where(adata.obs['endo_high'].notna(), 'Endothelial-high', None)
adata.obs['endo_sub'] = np.where(adata.obs['leiden'] == HEV_cluster, 'HEV', adata.obs['endo_sub'])

if sum(adata.obs['leiden'] == HEV_cluster) > 0:
    if len(adata.obs['endo_high'].dropna()) > 0:

        # --- Plot some marker genes

        marker_dict = {
            'endothelial': ['PECAM1', 'VWF', 'AQP1', 'COL4A1'],
            'HEV': ['CCL21', 'ACKR1', 'TSPAN7', 'MEOX2', 'RAMP3', 'CCL14']
        }

        sc.set_figure_params(scanpy=True, fontsize=10)

        dp = sc.pl.dotplot(adata, var_names=marker_dict, groupby=['endo_sub'], layer=use_data, swap_axes=False, 
                        return_fig=True, )
        dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5, cmap='Reds',).show()

        plt.suptitle(f'{sample_id} ({subtype}) - marker genes endothelial-high & HEV cluster spots (threshold={threshold}, endothelial-high: min_prop={prop_thresh}, min_cells={cells_thresh.round(2)} (quantile {quantile}))', y=1.1, fontweight='bold', x=0.5)


        # --- Plot the HEV signature

        gene_dir = os.path.join(project_dir, 'resources', 'gene_lists')

        filename = 'other_scores/HEV_score.csv'

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

        dp = sc.pl.dotplot(adata, var_names=gene_list, groupby=['endo_sub'], layer=use_data, swap_axes=False, 
                        return_fig=True, )
        dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5, cmap='Reds',).show()

        plt.suptitle(f'{sample_id} ({subtype}) - marker genes endothelial-high & HEV cluster spots (threshold={threshold}, endothelial-high: min_prop={prop_thresh}, min_cells={cells_thresh})', y=1.1, fontweight='bold', x=0.5)


        # --- Try to calculate LFCs to filter

        # Subset the data for each group
        adata_A = adata[adata.obs['endo_sub'] == 'Endothelial-high']
        adata_B = adata[adata.obs['endo_sub'] == 'HEV']

        # Compute mean expression per group
        mean_A = np.mean(adata_A[:, gene_list].to_df(layer=use_data), axis=0)
        mean_B = np.mean(adata_B[:, gene_list].to_df(layer=use_data), axis=0)

        # Compute log2 fold change
        lfc = np.log2((mean_B) / (mean_A))

        pos_lfc_genes = list(lfc[lfc > 0].index)
        neg_lfc_genes = list(lfc[lfc < 0].index)

        gene_dict = {'Positive': pos_lfc_genes,
                    'Negative': neg_lfc_genes}


        dp = sc.pl.dotplot(adata, var_names=gene_dict, groupby=['endo_sub'], layer=use_data, swap_axes=False, 
                        return_fig=True, )
        dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5, cmap='Reds',).show()

        plt.suptitle(f'{sample_id} ({subtype}) - marker genes endothelial-high & HEV cluster spots (threshold={threshold}, endothelial-high: min_prop={prop_thresh}, min_cells={cells_thresh})', y=1.1, fontweight='bold', x=0.5)


        # --- Attempt to score with just the positive LFC genes to see how that compares

        score_obj = ad.AnnData(adata.layers[use_data])
        score_obj.obs_names = adata.obs_names
        score_obj.var_names = adata.var_names

        seed = 42

        # All genes
        score_key = f'HEV_all_genes'
        sc.tl.score_genes(score_obj, gene_list, random_state=seed, score_name=score_key)
        # Put the score back in the original adata
        adata.obs[score_key] = score_obj.obs[score_key]


        # All genes
        score_key = f'HEV_pos_LFC'
        sc.tl.score_genes(score_obj, pos_lfc_genes, random_state=seed, score_name=score_key)
        # Put the score back in the original adata
        adata.obs[score_key] = score_obj.obs[score_key]


        # --- Plot 5: Plot score results

        sc.set_figure_params(scanpy=True, fontsize=18)

        nrows = 1
        ncols = 3

        # Make figure
        fig, axs = plt.subplots(nrows, ncols, figsize=(7.5*ncols, 6*nrows))

        # First plot H&E
        n = 0
        ax = axs[n]

        sc.pl.spatial(adata, color=None, size=spot_size, ax=ax, bw=False, alpha_img=1,
            crop_coord=coords, title=f'{sample_id} ({subtype})\nH&E',)

        ax.set_title(ax.get_title(), fontweight='bold')


        # Then plot clustering results
        n += 1
        ax = axs[n]

        sc.pl.spatial(adata, color='HEV_all_genes', groups=[HEV_cluster], size=spot_size, ax=ax, bw=False, alpha_img=0.7,
            crop_coord=coords, title=f'{sample_id} ({subtype})\nHEV score all genes', vmin=0)

        ax.set_title(ax.get_title(), fontweight='bold')

        # Then plot clustering results
        n += 1
        ax = axs[n]

        sc.pl.spatial(adata, color='HEV_pos_LFC', groups=[HEV_cluster], size=spot_size, ax=ax, bw=False, alpha_img=0.7,
            crop_coord=coords, title=f'{sample_id} ({subtype})\nHEV score pos LFC only', vmin=0, vmax=0.8)

        ax.set_title(ax.get_title(), fontweight='bold')


# --- Save results

# Save pots
out_dir_plots = os.path.join(out_dir, 'plots')
os.makedirs(out_dir_plots, exist_ok=True)

filename = f'{sample_id}_improve_HEV_score.pdf'
save_multi_image(os.path.join(out_dir_plots, filename))

# Save csv files
out_dir_csv = os.path.join(out_dir, 'csv_files')
os.makedirs(out_dir_csv, exist_ok=True)

if sum(adata.obs['leiden'] == HEV_cluster) > 0:
    if len(adata.obs['endo_high'].dropna()) > 0:

        # Save genes & whether they had pos or neg LFC
        lfc = pd.DataFrame(lfc)
        lfc.columns = ['LFC_HEV_vs_endothelial']
        lfc['direction'] = np.where(lfc['LFC_HEV_vs_endothelial'] > 0, 'positive', 'negative')

        # Add the original score info
        merged_df = gene_table.merge(lfc, right_index=True, left_on='genes')

        filename = f'{sample_id}_HEV_vs_endothelial_genes.csv'
        merged_df.to_csv(os.path.join(out_dir_csv, filename), index=False)

        # Also save high scoring spots (for checking with loupe)
        filename = f'{sample_id}_HEV_high_spots_loupe.csv'

        adata.obs['high_score'] = np.where(adata.obs[score_key] > 0.5, f'HEV score > 0.5', None)
        adata.obs['high_score'].to_csv(os.path.join(out_dir_csv, filename))


# Also save endothelial spots
filename = f'{sample_id}_endothelial_high_spots.csv'
to_save = adata.obs[['leiden', 'endo_high', 'endo_sub']]

to_save.to_csv(os.path.join(out_dir_csv, filename))


# --- Finish

print("Script successfully completed")
