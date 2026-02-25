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
        return most_recent_dir_path
    else:
        print("No valid directories found.")
        return None


# --- Load the adata objects with LR results

lr_sample_dir = os.path.join(res_dir, '04_vascular_microenvironment', '04_ligand_receptor_analysis',
                              '01_run_LR_analysis_stlearn', f'{bin_size}um_bin', 'csv_files')

most_recent_sample_path = find_most_recent(lr_sample_dir)

adata_dict = {}

for sample in sample_list:
    filename = f'{sample}_HEV_LR_analysis_stlearn.h5ad'
    adata_dict[sample] = ad.read_h5ad(os.path.join(most_recent_sample_path, filename))


# --- Load LR summary data

lr_dir = os.path.join(res_dir, '04_vascular_microenvironment', '04_ligand_receptor_analysis', 
                        '02_analyze_LR_all', f'{bin_size}um_bin', 'csv_files')

most_recent_path = find_most_recent(lr_dir)

# Load LR summary
lr_summary_file = 'All_samples_LR_summary.csv'
lr_summary = pd.read_csv(os.path.join(most_recent_path, lr_summary_file), index_col=0)

# Load LFC data
lr_lfc_file = 'All_samples_group_counts.csv'
lr_lfc = pd.read_csv(os.path.join(most_recent_path, lr_lfc_file), index_col=0)

# Sort the lr_list by LFC
lr_list = lr_lfc.sort_values(by='LFC', ascending=False).index.tolist()

# Extract list of lrs found in all samples
lr_shared = lr_summary[lr_summary['n_samples_sig'] == len(sample_list)].index.tolist()

# Subset to shared LRs
lr_list = [lr for lr in lr_list if lr in lr_shared]

print(f"Found {len(lr_shared)} LR pairs significant in all {len(sample_list)} samples")


# --- Build pivot_df from adata objects

# Extract cell type counts per sample
all_rows = []

for sample_id, adata in adata_dict.items():
    for lr_pair, df in adata.uns['per_lr_cci_cell_type'].items():
        for source in df.index:
            for target in df.columns:
                count = df.loc[source, target]
                if count > 0:  # skip zeroes if desired
                    all_rows.append({
                        'LR_pair': lr_pair,
                        'source_cell_type': source,
                        'target_cell_type': target,
                        'sample_id': sample_id,
                        'n_interactions': count
                    })

flat_df = pd.DataFrame(all_rows)

# Pivot to wide format
pivot_df = flat_df.pivot_table(
    index=['LR_pair', 'source_cell_type', 'target_cell_type'],
    columns='sample_id',
    values='n_interactions',
    fill_value=0
)

# Optional: flatten the multi-index columns
pivot_df.columns.name = None
pivot_df = pivot_df.reset_index()

# Sum across all sample columns
sample_cols = [col for col in pivot_df.columns if col not in ['LR_pair', 'source_cell_type', 'target_cell_type']]
pivot_df['total_n_interactions'] = pivot_df[sample_cols].sum(axis=1)


# --- Load HEV colors

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
color_file = 'TLS_HEV_colors.csv'

color_df = pd.read_csv(os.path.join(color_dir, color_file), index_col=0)

# Make a color_dict
color_dict = {}
for label, color in zip(color_df.index, color_df['color']):
    color_dict[label] = color


# --- Prepare rename dict for labels

rename_dict = {
    'near_TLS_HEV': 'HEV (TLS)',
    'non_TLS_HEV': 'HEV-like (stroma)',
    'nearby_spots': f'Perivascular bins',
}


# --- Select top LR pairs for plotting

top = 15
lr_sub = lr_list[:top]

df_sub = pivot_df[pivot_df['LR_pair'].isin(lr_sub)]

# Create the combined source→target column
df_sub['source_target'] = df_sub['source_cell_type'] + '→' + df_sub['target_cell_type']

# Get all unique sources and targets
sources = df_sub['source_cell_type'].unique()
targets = sorted(df_sub['target_cell_type'].unique())

# Build the desired order
ordered_source_targets = []
for source in sources:
    for target in targets:
        ordered_source_targets.append(f"{source}→{target}")

# Set as categorical with the desired order
# Subset to only categories that are in the df_sub
ordered_source_targets = [st for st in ordered_source_targets if st in df_sub['source_target'].unique()]
# Set the categorical order
df_sub['source_target'] = pd.Categorical(df_sub['source_target'], categories=ordered_source_targets, ordered=True)

merged = df_sub

# Reorder based on LR pair list
merged['LR_pair'] = pd.Categorical(merged['LR_pair'], categories=lr_sub, ordered=True)

size_max = merged['total_n_interactions'].max()


# --- Setup plotting parameters

size_per_spot = 0.1  # Size of each dot in the dotplot
cmap = 'magma_r'
xval = 'LR_pair'
yval = 'source_target'

# Disable grid
plt.rcParams['axes.grid'] = False


# --- Create the plot

nrows = 1
ncols = 1

fig, ax = plt.subplots(ncols=ncols, nrows=nrows, figsize=(3, 3))

# Extract number of cell types and LRs
n_cell_types = len(merged[yval].unique())
n_lrs = len(merged[xval].unique())

ax_pos = ax.get_position()
ax_pos = ax.set_position([ax_pos.x0, ax_pos.y0, top * size_per_spot, n_cell_types * size_per_spot*0.8])

dotplot_ax, _ = plotting_functions.seaborn_dotplot(
    merged,
    x=xval,
    y=yval,
    hue='total_n_interactions',
    size='total_n_interactions',
    cmap=cmap,
    dotplot_ax=ax,
    max_size_px=300,
    size_max=size_max,
)

# Adjust limits of plot
ax.set_xlim(-0.75, len(merged[xval].unique()) - 0.25)
ax.set_ylim(len(merged[yval].unique()) - 0.25, -0.75)

# Remove the legend
dotplot_ax.legend_.remove()

# Add labels
ax.set_xlabel('LR pair', fontweight='bold')
ax.set_ylabel('', fontweight='bold')

# Set xticklabels
ax.set_xticklabels(dotplot_ax.get_xticklabels(), rotation=90, ha='center')


# --- Add colorbar

norm = mpl.colors.Normalize(vmin=merged['total_n_interactions'].min(), vmax=merged['total_n_interactions'].max())
sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

cbar_ax = fig.add_axes([0.91, 0.15, 0.02, 0.7])
cbar_ax.grid(False)
cbar = plt.colorbar(sm, cax=cbar_ax)
cbar.set_label('# of interactions', rotation=270, labelpad=15)

# Set colorbar position
main_ax_pos = ax.get_position()
cbar_ax.set_position([main_ax_pos.x1 + size_per_spot*0.2, main_ax_pos.y0, size_per_spot*0.5, main_ax_pos.height])


# --- Add rectangles showing source→target relationships

# Create a new axis on the left side of the main axis
rect_ax = fig.add_axes([ax.get_position().x0 - 0.02 - size_per_spot*1, ax.get_position().y0, 
                        size_per_spot*1, ax.get_position().height])
rect_ax.set_ylim(ax.get_ylim())

# Add rectangles according to yticklabels
# Add small rectangles on the left to indicate the group labels
for i, label in enumerate(ax.get_yticklabels()):
    label_text = label.get_text()
    # Split on the arrow
    if '→' in label_text:
        source, target = label_text.split('→')
    else:
        continue
    
    # Left rectangle (source)
    if source in color_dict:
        color = color_dict[source]
        rect_ax.add_patch(
            patches.Rectangle((0, i - 0.35), 0.3, 0.75, facecolor=color, 
                            transform=rect_ax.transData, edgecolor='black', linewidth=0.8)
        )
    
    # Right rectangle (target)
    if target in color_dict:
        color = color_dict[target]
        rect_ax.add_patch(
            patches.Rectangle((0.7, i - 0.35), 0.3, 0.75, facecolor=color, 
                            transform=rect_ax.transData, edgecolor='black', linewidth=0.8)
        )
    
    # Add an arrow from left rectangle to right rectangle
    rect_ax.annotate('', 
        xy=(0.7, i), xytext=(0.35, i),
        arrowprops=dict(facecolor='black', width=0.5, headwidth=5, headlength=4),
        ha='center'
    )

# Make the spines of rect_ax invisible
rect_ax.spines['top'].set_visible(False)
rect_ax.spines['right'].set_visible(False)
rect_ax.spines['left'].set_visible(False)
rect_ax.spines['bottom'].set_visible(False)

# Remove ticks and labels from the rectangle axis
rect_ax.set_xticks([])

# Take the y-ticks from the main axis
rect_ax.set_yticks(ax.get_yticks())
rect_ax.set_yticklabels(ax.get_yticklabels())

# Remove y-tick labels from the main axis
ax.set_yticks([])


# --- Create legend

# Figure out the labels from the yticklabels
legend_labels = []

for label in rect_ax.get_yticklabels():
    label_text = label.get_text()
    # Split on the arrow
    if '→' in label_text:
        source, target = label_text.split('→')
        if source not in legend_labels:
            legend_labels.append(source)
        if target not in legend_labels:
            legend_labels.append(target)

# Make handles based on color_dict
handles = []
for label in legend_labels:
    if label in color_dict:
        color = color_dict[label]
        handles.append(patches.Patch(facecolor=color, label=label, edgecolor='black', linewidth=0.8))

# Add legend to the axis
x_offset = 1/len(ax.get_xticklabels())*3  # Adjust x_offset as needed
ax.legend(handles=handles, loc='lower left', bbox_to_anchor=(0, 1),
          frameon=False, handlelength=0.7, handletextpad=0.5, borderpad=0)

# Update the labels of the legend
for text in ax.get_legend().get_texts():
    if text.get_text() in rename_dict:
        text.set_text(rename_dict[text.get_text()])

# Update the yticklabels of the rect_ax
new_yticklabels = []
for label in rect_ax.get_yticklabels():
    label_text = label.get_text()
    # Labels to update will be substrings
    for key, value in rename_dict.items():
        if key in label_text:
            label_text = label_text.replace(key, value)
    # Also make sure to put space around the arrow
    if '→' in label_text:
        source, target = label_text.split('→')
        label_text = f'{source} → {target}'
    new_yticklabels.append(label_text)

rect_ax.set_yticklabels(new_yticklabels)

# Calculate the midpoint to position the title
title_x = (ax.get_position().x0 + ax.get_position().x1)/2
title_y = ax.get_position().y1 + 0.3

plt.suptitle(f'Top {top} LR pairs enriched in TLS-HEVs', 
             y=title_y, x=title_x, fontweight='bold')


# --- Save plot

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_4K.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')


# --- Finish

print("Script successfully completed")
print(f"Saved plot to: {os.path.join(out_dir, filename)}")

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
