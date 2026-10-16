#!/usr/bin/env python3

# --- Aim
# Compare gene expression profiles of vessels overlapping with HEVs versus non-HEV vessels.
# Analyze differences by vessel type across all samples.

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import math
import os
import sys

import anndata as ad
import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import scanpy.plotting.palettes as sc_palettes
from matplotlib.backends.backend_pdf import PdfPages

sc.set_figure_params(facecolor="white", figsize=(8, 8), fontsize=20)
sc.settings.verbosity = 3

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments
parser = argparse.ArgumentParser(description="Differential expression between HEV and non-HEV venous vessel bins",
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
sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import auto_crop_scanpy, save_multi_image

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


# --- Load miniatlas classifications

decon_level = 'minor'

# To make it easier to visualize trends by lineage
class_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
class_file = 'miniatlas_subset_classifications.csv'

class_df = pd.read_csv(os.path.join(class_dir, class_file))

# Get endothelial classes
endothelial = list(class_df[class_df['miniatlas_major'] == 'Endothelial'][f'miniatlas_{decon_level}'].unique())
endothelial = sorted(endothelial)


# -- Load data

use_data = 'log_norm'

resolution = 'square_016um'

# Multiple samples as a dictionary seems like an okay idea
adata_dict = {}

for sample in sample_list:
    adata_dict[sample] = ad.read_h5ad(os.path.join(h5ad_dir, 'visium_hd', resolution, f'{sample}.h5ad'))

# Add some info for plotting
for sample in sample_list:
    adata_dict[sample].obs['sample_id'] = sample
    adata_dict[sample].obs['subtype'] = subtype_dict[sample]
    adata_dict[sample].obs['donor_id'] = sample_lut.loc[sample, 'donor_id']


# --- Calculate crop coords for optimal presentation

coords_dict = {}

for sample in sample_list:
    coords_dict[sample] = auto_crop_scanpy(adata_dict[sample], border_size=0.1)


# --- Load the vessel groups

# Define which endothelial type want to look at
endo_of_interest = 'Endothelial ACKR1'
out_name = endo_of_interest.split(' ')[-1]

# Define threshold
threshold = 0.5
threshold_adj = str(threshold).split('.')[-1]

# Specify path
TLS_group_dir = os.path.join(res_dir, '04_vascular_microenvironment', '02_identify_all_vasculature', '01_classify_vessels', 'csv_files')
most_recent_dir_path = TLS_group_dir

for sample in sample_list:
    filename = f'{sample}_vessel_classification.csv'
    vessel_classes = pd.read_csv(os.path.join(most_recent_dir_path, filename), index_col=0)
    adata_dict[sample].obs[vessel_classes.columns] = vessel_classes


# --- Load HEV overlap information

TLS_group_dir = os.path.join(res_dir, '04_vascular_microenvironment', '02_identify_all_vasculature', '02_investigate_vessels', 'csv_files')

most_recent_dir_path = TLS_group_dir

for sample in sample_list:
    filename = f'{sample}_vessel_HEV_overlap_spot_level.csv'
    overlap_df = pd.read_csv(os.path.join(most_recent_dir_path, filename), index_col=0) 
    adata_dict[sample].obs['overlap'] = overlap_df['overlap']
    # Also load the file with proportion of overlap
    filename = f'{sample}_vessel_HEV_overlap.csv'
    proportion_df = pd.read_csv(os.path.join(most_recent_dir_path, filename), index_col=0)
    # Identify groups with at least 50% overlap
    min_overlap = 0.5
    # Filter for the endothelial type of interest
    proportion_df_endo = proportion_df[proportion_df['group_type'] == endo_of_interest].copy()
    # Separate groups with overlap above vs below threshold (or NaN)
    proportion_df_high = proportion_df_endo[proportion_df_endo['pct_overlap'] > min_overlap]
    # Include both low overlap and NaN (vessels with 0 or no overlap data)
    proportion_df_low = proportion_df_endo[(proportion_df_endo['pct_overlap'] <= min_overlap) | (proportion_df_endo['pct_overlap'].isna())]
    high_overlap_groups = proportion_df_high['group'].to_list()
    low_overlap_groups = proportion_df_low['group'].to_list()
    # Initialize the column
    adata_dict[sample].obs['majority_overlap'] = np.nan
    # Set 'overlap' for vessels with high overlap
    barcode_groups_high = adata_dict[sample].obs[adata_dict[sample].obs[endo_of_interest].isin(high_overlap_groups)].index
    adata_dict[sample].obs.loc[barcode_groups_high, 'majority_overlap'] = 'overlap'
    # Set vessel type name for vessels with low/no overlap
    barcode_groups_low = adata_dict[sample].obs[adata_dict[sample].obs[endo_of_interest].isin(low_overlap_groups)].index
    adata_dict[sample].obs.loc[barcode_groups_low, 'majority_overlap'] = endo_of_interest


# --- Also set colors for endothelial types

# Define some colors per endothelial type
colors_file = os.path.join(project_dir, 'resources', 'plotting_colors', 'endothelial_colors.csv')
colors_df = pd.read_csv(colors_file)

# Make a color_dict
color_dict = colors_df.set_index('endothelial_type')['color'].to_dict()

# Add a color for overlap
color_dict['overlap'] = 'indigo'
color_dict['overlap HEV'] = color_dict['overlap']


# --- Sort by subtype

# Extract all subtypes
all_subtypes = {subtype_dict[sample] for sample in sample_list if sample in subtype_dict}

# Sort according to preferred plotting order
subtype_order = ['Luminal A', 'Luminal B', 'Luminal B (HER2+)', 'HER2+', 'TNBC']
sorted_subtypes = sorted(all_subtypes, key=lambda x: subtype_order.index(x) if x in subtype_order else float('inf'))

# Plotting order per subtype
list_all_samples = sorted(sample_list, key=lambda s: (sorted_subtypes.index(subtype_dict[s]), s))


# --- Prepare for saving

qc_dir = os.path.join(res_dir, '04_vascular_microenvironment', '03_compare_HEV_to_other_venous', '02_compare_gene_expression')
os.makedirs(qc_dir, exist_ok=True)

out_dir = os.path.join(qc_dir,)
os.makedirs(out_dir, exist_ok=True)

# Save plots
out_dir_plots = os.path.join(out_dir, 'plots')
os.makedirs(out_dir_plots, exist_ok=True)


# --- Plot 1: Spatial overview of endothelial + HEV overlap

if len(adata_dict.keys()) < 5:
        ncols = len(adata_dict.keys())
else:
    ncols = 5

nrows = math.ceil((len(adata_dict.keys()))/ncols)

if nrows == 1:
    plt.figure(figsize=(8*ncols, 6*nrows))
    plt_y = 1.1
else:
    plt.figure(figsize=(8*ncols, 7*nrows))
    plt_y = 1

classification_key = 'majority_overlap'

for n, sample in enumerate(list_all_samples):
    if n == (ncols-1):
        legend_loc = 'right margin'
    else:
        legend_loc = None
    # add a new subplot iteratively
    ax = plt.subplot(nrows, ncols, n + 1)
    # filter df and plot ticker on the new subplot axis
    subtype = subtype_dict[sample]
    adata_sub = adata_dict[sample][adata_dict[sample].obs[classification_key].notna()]
    sc.pl.spatial(adata_sub, color=classification_key, bw=False, alpha_img=0.7, palette=color_dict,
                  crop_coord=coords_dict[sample], title=f'{sample} ({subtype})\n{endo_of_interest} HEV overlap', ax=ax, legend_loc=legend_loc)
    ax.set_title(ax.get_title(), fontweight='bold')

plt.suptitle(f'All samples - {endo_of_interest} HEV overlap (threshold={threshold})', y=plt_y, fontweight='bold')


# --- Save the spatial plots

filename = f'All_samples_{out_name}_HEV_overlap_spatial.pdf'
# save_multi_image(os.path.join(out_dir_plots, filename))

# Save with high resolution
filename_dir=os.path.join(out_dir_plots, filename)

pp = PdfPages(filename_dir)
fig_nums = plt.get_fignums()

figs = [plt.figure(n) for n in fig_nums]

for fig in figs:
    fig.savefig(pp, format='pdf', bbox_inches='tight', dpi=300)

pp.close()

plt.close('all')

# Reset fontsize
sc.set_figure_params(fontsize=10)


# --- Extract a single adata and subset to just the vessel spots

classification_key = 'majority_overlap'

adata_sub_dict = {}

for sample in sample_list:
    adata_sub_dict[sample] = adata_dict[sample][adata_dict[sample].obs[classification_key].notna()]

# Convert to adata list & concatenate
adata_list = list(adata_sub_dict.values())
all_adata = ad.concat(adata_list, merge='same')


# --- Now find marker genes for HEV classes

all_keys = []
rank_genes_key = 'rank_genes_HEV_overlap'
all_keys.append(rank_genes_key)

sc.tl.rank_genes_groups(all_adata, classification_key, use_raw = False, layer=use_data, method='wilcoxon', key_added=rank_genes_key,
                        pts=True)


# --- Plot 2: Plot top genes

sc.pl.rank_genes_groups_dotplot(all_adata, n_genes=50, groupby=classification_key, key=rank_genes_key, layer=use_data)

plt.suptitle(f'All samples - {endo_of_interest} HEV overlap DEGs, top 50 genes per cluster (use_data={use_data})', y=1.1, fontweight='bold')


# --- Plot 3: 1 vs rest plots

sc.pl.rank_genes_groups(all_adata, n_genes=25, groupby=classification_key, key=rank_genes_key, layer=use_data)

plt.suptitle(f'All samples - {endo_of_interest} HEV overlap DEGs 1 vs rest plots, 25 genes (use_data={use_data})', y=1.1, fontweight='bold')


# --- Save plots so far

filename = f'All_samples_{out_name}_HEV_overlap_DEGs.pdf'
save_multi_image(os.path.join(out_dir_plots, filename))

plt.close('all')


# --- Save output

out_dir_csv = os.path.join(out_dir, 'csv_files')
os.makedirs(out_dir_csv, exist_ok=True)

#  Save the differential genes
for rank_genes_key in all_keys:

    dfs = []

    result = all_adata.uns[rank_genes_key]
    groups = result['names'].dtype.names
    reference = result['params']['reference']

    # Iterate through the groups
    for group in groups:
        # Create a DataFrame for the current group
        df = pd.DataFrame({
            'group': group,
            'reference': reference,
            'genes': result['names'][group],
            'pvals': result['pvals'][group],
            'logfoldchanges': result['logfoldchanges'][group],
            'pts': result['pts'][group],
            'pvals_adj': result['pvals_adj'][group],
            'z-score': result['scores'][group]
        })
        # Append the DataFrame to the list
        dfs.append(df)

    # Concatenate the DataFrames from all groups into one DataFrame
    final_df = pd.concat(dfs, ignore_index=True)

    filename = f'All_samples_{rank_genes_key}.csv'
    final_df.to_csv(os.path.join(out_dir_csv, filename), index=False)

# Save the whole object as well for future reference
# Convert all object or category columns in .obs to string
for col in all_adata.obs.columns:
    if all_adata.obs[col].dtype.name in ['object', 'category']:
        all_adata.obs[col] = all_adata.obs[col].astype(str)

all_adata.write(os.path.join(out_dir_csv, f"All_samples_HEV_{out_name}.h5ad"))


# --- Finish

print("Script successfully completed")