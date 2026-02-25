#!/usr/bin/env python3

# --- Aim
# Plot dotplot showing top LR pairs enriched in TLS-HEVs (Figure 4K)
# Shows the source→target cell type combinations for each LR pair

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import os
import sys
from datetime import datetime

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments
parser = argparse.ArgumentParser(description="Plot top LR pairs dotplot for Figure 4K",
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


# --- Import custom plotting functions

# Import function script as a module
sys.path.append(os.path.join(project_dir, 'scripts', '04_vascular_microenvironment'))
import plotting_functions


# --- Extract some info from the sample sheet
# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_visium_hd.csv')

sample_sheet = pd.read_csv(sample_sheet_path)
sample_lut = sample_sheet.set_index('section_id')


# --- Parameters

threshold = 0.5
threshold_adj = str(threshold).split('.')[-1]
bin_size = 25  # um
resolution = 'square_016um'


# --- Function to find most recent results

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

# Find most recent results
resolution = 'square_016um'
threshold = 0.5
threshold_adj = str(threshold).split('.')[-1]
bin_size = 25  # um

lr_dir = os.path.join(res_dir, '04_vascular_microenvironment', '04_ligand_receptor_analysis', '01_run_LR_analysis_stlearn', f'{bin_size}um_bin', 'csv_files')

most_recent_dir_path = find_most_recent(lr_dir)

adata_dict_stlearn = {}

for sample in sample_list:
    # Load the anndata object
    filename = f'{sample}_HEV_LR_analysis_stlearn.h5ad'
    adata_dict_stlearn[sample] = ad.read_h5ad(os.path.join(most_recent_dir_path, filename))


# --- Load and add the decon results

decon_resolution = resolution

resolution_numbers = re.findall(r'\d+', decon_resolution)[0]
resolution_numbers = int(resolution_numbers)

decon_level = 'subset'
decon_dir = os.path.join(res_dir, '02_cell_type_mapping', 'deconvolution', 'visium_hd', f'{decon_level}_level', decon_resolution,)

for sample in sample_list:
    decon_file = f'{sample}_{decon_level}_{resolution}_q05_cell_abundance_w_sf.csv'
    decon_res = pd.read_csv(os.path.join(decon_dir, decon_file), index_col=0)
    # Convert to proportions within each spot
    decon_res = decon_res.div(decon_res.sum(axis=1), axis=0)
    cell_types = list(decon_res.columns)
    adata_dict_stlearn[sample].obs[cell_types] = decon_res


# --- Also load the revisited spot labels & add to adata
# These distinguish whether a nearby spot is nearby to TLS_HEV or non_TLS_HEV

lr_dir = os.path.join(res_dir, '04_vascular_microenvironment', '04_ligand_receptor_analysis', '03_plot_LR_per_sample', f'{bin_size}um_bin', 'csv_files')

most_recent_dir_path = find_most_recent(lr_dir)

for sample in sample_list:
    spot_label_file = f'{sample}_HEV_spot_grouping.csv'
    label_df = pd.read_csv(os.path.join(most_recent_dir_path, spot_label_file), index_col=0)
    adata_dict_stlearn[sample].obs['nearby_spots_revisited'] = label_df['nearby_spots_revisited']


# --- Load the miniatlas classifications

# These are needed to check the alignment between major and minor levels
miniatlas_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
miniatlas_file = 'miniatlas_subset_classifications.csv'

miniatlas_df = pd.read_csv(os.path.join(miniatlas_dir, miniatlas_file),)

# Also define a dict with shorter names for the subset level cell types
miniatlas_sub_dict = miniatlas_df.set_index('miniatlas_subset')[['miniatlas_subset_short']].to_dict()

# Order miniatlas_major by lineage
lineage_order = ['Epithelial', 'Immune', 'Stromal']

major_cell_types = miniatlas_df[f'miniatlas_major'].dropna().unique().tolist()

# Sort based on lineage first (recorded in the miniatlas df)
def get_lineage(cell_type, level='major'):
    if cell_type in miniatlas_df[f'miniatlas_{level}'].values:
        lineage = miniatlas_df[miniatlas_df[f'miniatlas_{level}'] == cell_type]['lineage'].values[0]
    else:
        lineage = 'other'
    return lineage

major_cell_types_sorted = sorted(major_cell_types, key=lambda x: (lineage_order.index(get_lineage(x)) if get_lineage(x) in lineage_order else len(lineage_order), x))


# --- Load plotting colors

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')

# Use the major level for colors
color_level = 'major'
color_file = f'miniatlas_{color_level}_colors.csv'

color_dict = pd.read_csv(os.path.join(color_dir, color_file), index_col=0).to_dict()['color']

# Add the other levels to the color_dict
for cell_type in miniatlas_df[f'miniatlas_{decon_level}'].unique():
    if cell_type not in color_dict.keys():
        # Check that the value is not nan
        if pd.isna(cell_type):
            continue
        print(f'Adding color for {cell_type}')
    # Get the corresponding major cell type
    major_type = miniatlas_df[miniatlas_df[f'miniatlas_{decon_level}'] == cell_type][f'miniatlas_major'].values
    # Find which value is in the color_dict
    major_type_in_dict = [m for m in major_type if m in color_dict.keys()]
    if len(major_type_in_dict) > 0:
        color_dict[cell_type] = color_dict[major_type_in_dict[0]]
    else:
        color_dict[cell_type] = 'lightgrey'


# --- Find the spots that send LTB signals

adata = adata_dict_stlearn[sample_list[0]]

# Extract all significant lrs as list
lr = 'LTB_LTBR'

print(f'Analyzing LR: {lr}')

avg_comp_dict = {}

for sample in sample_list:
    adata = adata_dict_stlearn[sample]
    # Check if the LR is in the adata object
    if 'lr_summary' not in adata.uns or lr not in adata.uns['lr_summary'].index:
        print(f'LR {lr} not found in {sample}. Skipping...')
        continue
    # Get the list of all LR pair names
    lr_pairs = adata.uns['lr_summary']
    # Find the index of the LR pair of interest
    if isinstance(lr_pairs, pd.DataFrame):
        lr_index = lr_pairs.index.get_loc(lr)
    else:
        if isinstance(lr_pairs, Iterable):
            lr_index = lr_pairs.index(lr)
        else:
            raise ValueError("lr_pairs should be either a DataFrame or a list of LR pairs.")
    # Get adjusted p-values and scores
    p_adj_matrix = adata.obsm['p_adjs'][:, lr_index]
    sig_mask = p_adj_matrix < 0.05
    barcodes_sig = adata.obs_names[sig_mask]
    if len(barcodes_sig) == 0:
        continue  # No significant spots found
    # Subset to significant barcodes
    adata_sub = adata[adata.obs_names.isin(barcodes_sig)]
    # Subset to nearby spots near TLS-HEV as these appear to be main senders
    # Going to exclude the non-TLS HEV nearby spots bc these will have very different composition
    adata_sub = adata_sub[adata_sub.obs['nearby_spots_revisited'] == 'nearby_near_TLS_HEV']
    # Calculate average composition of these spots
    avg_comp = adata_sub.obs[cell_types].mean()
    # Also calculate the average comp of all nearby spots for comparison
    adata_nearby = adata[adata.obs['nearby_spots_revisited'] == 'nearby_near_TLS_HEV']
    # Remove the spots sending LTB signals
    adata_nearby = adata_nearby[~adata_nearby.obs_names.isin(barcodes_sig)]
    avg_comp_all = adata_nearby.obs[cell_types].mean()
    # Calculate fold change
    fold_change = avg_comp / avg_comp_all
    # Make log2 fold change
    fold_change = np.log2(fold_change.replace(0, np.nan))  # Replace 0 with NaN to avoid -inf
    # Also record the absolute change
    abs_change = avg_comp - avg_comp_all
    # Make a df with all info
    avg_comp_df = pd.DataFrame({
        'cell_type': cell_types,
        'avg_composition_sending_nearby': avg_comp.values,
        'avg_composition_all_nearby': avg_comp_all.values,
        'LFC': fold_change.values,
        'abs_change': abs_change.values,
        'sample_id': sample,
    })
    avg_comp_dict[sample] = avg_comp_df

# Extract a single df for all samples
all_avg_comp_df = pd.concat(avg_comp_dict.values(), ignore_index=True)

# Sort by median LFC
all_avg_comp_df['cell_type'] = pd.Categorical(all_avg_comp_df['cell_type'], 
                                              categories=all_avg_comp_df.groupby('cell_type')['LFC'].median().sort_values(ascending=False).index,
                                              ordered=True)

sorted_cell_types = all_avg_comp_df['cell_type'].cat.categories.tolist()


# --- Statistical testing for enrichment/depletion per cell type

def analyze_cell_type_enrichment(lfc_values, cell_name):
    """
    Performs a one-sample t-test on LFC values to check for significant
    enrichment (mean > 0) or depletion (mean < 0).
    """
    print(f"--- Analyzing {cell_name} ---")
    n = len(lfc_values)
    if n < 3: 
        print(f"Cannot perform t-test: only {n} samples. Need at least 3.")
        return None, None, None
    # Perform the one-sample t-test
    t_statistic, p_value = stats.ttest_1samp(lfc_values, popmean=0)
    # Get the mean LFC
    mean_lfc = np.mean(lfc_values)
    print("\n")
    return t_statistic, p_value, mean_lfc

stat_df = pd.DataFrame(columns=['cell_type', 't_statistic', 'p_value', 'mean_correlation'])

for cell_type in sorted_cell_types:
    LFC_df = all_avg_comp_df[all_avg_comp_df['cell_type'] == cell_type]['LFC'].values
    t_statistic, p_value, mean_correlation = analyze_cell_type_enrichment(LFC_df, cell_type)
    stat_df = stat_df._append({'cell_type': cell_type, 't_statistic': t_statistic, 'p_value': p_value, 'mean_correlation': mean_correlation}, ignore_index=True)

alpha = 0.05
stat_df['significant'] = stat_df['p_value'] < alpha


# --- Get top 10 cell types for plotting

median_LFC = all_avg_comp_df.groupby('cell_type')['LFC'].median()
top10_cell_types = median_LFC[median_LFC > 0].sort_values(ascending=False).head(10).index.tolist()

# Filter df
all_avg_comp_df_top10 = all_avg_comp_df[all_avg_comp_df['cell_type'].isin(top10_cell_types)]

# Reorder the cell types based on median LFC
all_avg_comp_df_top10['cell_type'] = pd.Categorical(all_avg_comp_df_top10['cell_type'],
                                                categories=top10_cell_types,
                                                ordered=True)


# --- Plot: Top 10 enriched cell types (Figure 4K)

sc.set_figure_params(fontsize=10)

# Properties for boxplot
whisker_props = {'color': 'black'}
median_props = {'color': 'black'}

fig, ax = plt.subplots(figsize=(5, 3))

sns.boxplot(data=all_avg_comp_df_top10, x='cell_type', y='LFC', showfliers=False, linewidth=1, whiskerprops=whisker_props, 
                        medianprops=median_props, width=0.8, palette=color_dict, ax=ax, saturation=1)

sns.stripplot(data=all_avg_comp_df_top10, x='cell_type', y='LFC', color='black', jitter=True, size=3,
                        marker='o', ax=ax)

# Beautify the plot
for bar in ax.patches: # Adjust the colors of the boxplot
    # Get the RGBA values of the fill color
    fill_color = bar.get_facecolor()
    fill_color = list(fill_color)
    darkened_color = fill_color.copy()
    for i in range(len(darkened_color)-1):
        darkened_color[i] *= 0.8
    fill_color[-1] = 0.8 # make the fill color a bit transparent
    bar.set_color(tuple(fill_color))
    bar.set_edgecolor(tuple(darkened_color))

# Automatically resize based on number of xticks
n_ticks = len(ax.get_xticks())
ax_pos = ax.get_position()
ax.set_position([ax_pos.x0, ax_pos.y0, 0.05*n_ticks, ax_pos.height*0.5])

# Rotate xlabels
ax.set_xticklabels(ax.get_xticklabels(), rotation=90,)

# Adjust xticklabels using miniatlas_short notations
new_labels = [miniatlas_sub_dict['miniatlas_subset_short'].get(label.get_text(), label.get_text()) for label in ax.get_xticklabels()]
ax.set_xticklabels(new_labels)

# Add thick gray line at y=0
ax.axhline(0, color='lightgray', linewidth=2, zorder=0)

# Add an asterisk for significant cell types or p-value for non-significant
ylim = ax.get_ylim()
range_y = ylim[1] - ylim[0]

offset = range_y * 0.02  # 2% of the y-axis range

for i, cell_type in enumerate(top10_cell_types):
    cell_stats = stat_df[stat_df['cell_type'] == cell_type]
    if not cell_stats.empty:
        # Get the y position of the top whisker
        cell_data = all_avg_comp_df[all_avg_comp_df['cell_type'] == cell_type]['LFC']
        whisker_top = cell_data.max()
        if cell_stats['significant'].values[0]:
            # Place the asterisk slightly above the whisker
            ax.text(i, whisker_top + offset, '*', ha='center', va='bottom', color='darkred', fontsize=9, fontweight='bold')
        else:
            # Place the p-value above the whisker, rounded to 3 digits
            pval = cell_stats['p_value'].values[0]
            ax.text(i, whisker_top + offset, f"{pval:.3f}".lstrip('0'), ha='center', va='bottom', color='black', fontsize=7)

# Increase ymax a bit to make space for asterisks/p-values
ax.set_ylim(top=ax.get_ylim()[1]*1.2)

# Calculate equal positioning x & y
offset = 0.04
aspect = ax.get_position().width / ax.get_position().height
x_pos = 1 - offset / aspect
y_pos = 1 - offset

# Add some text at top right of ax explaining * p<0.05
stat_df_top10 = stat_df[stat_df['cell_type'].isin(top10_cell_types)]
if stat_df_top10['significant'].any():
    # Increase ymax a bit to make space for asterisks
    ax.set_ylim(top=ax.get_ylim()[1]*1.1)
    ax.text(x_pos, y_pos, '* p<0.05', ha='right', va='top', transform=ax.transAxes, fontsize=8, color='darkred',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

# Create legend handles (colored boxes with labels)
cell_types_plotted_top10 = all_avg_comp_df_top10['cell_type'].unique().tolist()

# Find matching major labels using miniatlas_df
miniatlas_subset_to_major = miniatlas_df.set_index('miniatlas_subset')['miniatlas_major'].to_dict()
major_cell_types_plotted_top10 = []
for ct in cell_types_plotted_top10:
    if ct in miniatlas_subset_to_major:
        major_ct = miniatlas_subset_to_major[ct]
        if major_ct not in major_cell_types_plotted_top10:
            major_cell_types_plotted_top10.append(major_ct)

labels_top10 = [label for label in major_cell_types_sorted if label in major_cell_types_plotted_top10]

handles_top10 = [
    mpl.patches.Patch(facecolor=color_dict[label], label=label, edgecolor='black', linewidth=0.8)
    for label in labels_top10 if label in color_dict
]

legend = ax.legend(handles=handles_top10, bbox_to_anchor=(1, 0.5), loc='center left', frameon=False, handlelength=0.8, fontsize=10)

# Make the legend boxes look the same as the barplot
for legpatch in legend.get_patches():
        fill_color = legpatch.get_facecolor()
        fill_color = list(fill_color)
        darkened_color = fill_color.copy()
        for i in range(len(darkened_color)-1):
            darkened_color[i] *= 0.8
        fill_color[-1] = 0.8 # make the fill color a bit transparent
        legpatch.set_color(tuple(fill_color))
        legpatch.set_edgecolor(tuple(darkened_color))
        legpatch.set_linewidth(1)

lr_adj = lr.replace('_', '→').replace('LTB', 'LTβ').replace('LTBR', 'LTβR')
ax.set_title(f'Cellular enrichment\nBins sending {lr_adj} signals to HEVs', fontweight='bold')

# Set labels
ax.set_ylabel('LFC', fontweight='bold')
ax.set_xlabel('Cell type', fontweight='bold')

# Position title
title_x = ax.get_position().x0 + ax.get_position().width / 2
title_y = ax.get_position().y1 + 0.3

plt.suptitle(f'All samples - Cell type composition of spots sending {lr_adj} signals in nearby_TLS_HEV vs all_TLS_HEV (bin_size={bin_size}µm)', fontsize=10, x=title_x, y=title_y, fontweight='bold')


# --- Save the output

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_4K.pdf'

plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')


# --- Finish

print("Script successfully completed")
print(f'Saved plot to: {os.path.join(out_dir, filename)}')
