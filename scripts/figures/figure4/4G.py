#!/usr/bin/env python3

# --- Aim
# Generate Figure 4G: Scatterplot showing key LR interactions between HEV groups (TLS-associated vs stromal HEV-like)

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from adjustText import adjust_text
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments
parser = argparse.ArgumentParser(description="Generate Figure 4G",
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


# --- Get some info from the sample sheet
# Use --sample_sheet if provided, otherwise default to config/sample_sheet_visium_hd.csv
if args["sample_sheet"] is not None:
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_visium_hd.csv')

sample_sheet = pd.read_csv(sample_sheet_path)
sample_lut = sample_sheet.set_index('section_id')


# --- Parameters

bin_size = 25  # um


# --- Load group counts and LR summary from CSV files

lr_summary_dir = os.path.join(res_dir, '04_vascular_microenvironment', '04_ligand_receptor_analysis', 
                               '02_analyze_LR_all', f'{bin_size}um_bin', 'csv_files')
most_recent_lr_summary = lr_summary_dir

# Load group counts CSV
group_counts_file = 'All_samples_group_counts.csv'
plot_df = pd.read_csv(os.path.join(most_recent_lr_summary, group_counts_file), index_col=0)

# Load LR summary to get shared LRs
lr_summary_file = 'All_samples_LR_summary.csv'
lr_summary = pd.read_csv(os.path.join(most_recent_lr_summary, lr_summary_file))

# Get shared LRs (LRs found in all samples)
shared_lrs = lr_summary[lr_summary['n_samples_sig'] == len(sample_list)]['LR'].tolist()

# Subset to only LRs that are shared across samples
plot_df = plot_df.loc[shared_lrs]

# Get HEV labels from the columns (excluding 'LFC')
HEV_labels = [col for col in plot_df.columns if col != 'LFC']


# --- Helper functions for plotting

def adjust_scatterplot(ax, plot_df, max_lfc, labels='Number of groups'):
    ax.set_aspect('equal', adjustable='box')
    # Draw a diagonal line
    max_value = max(plot_df[HEV_labels[0]].max(), plot_df[HEV_labels[1]].max())
    max_value = max_value * 1.05  # Add a bit of padding
    ax.plot([0, max_value], [0, max_value], color='gray', linestyle='--', zorder=1, linewidth=1)
    # Reset the axis limits
    ax.set_xlim(0, max_value)
    ax.set_ylim(0, max_value)
    # Set grid behind plot
    ax.set_axisbelow(True)
    # Add colorbar with LFC values
    norm = plt.Normalize(-max_lfc, max_lfc)
    sm = plt.cm.ScalarMappable(cmap='coolwarm', norm=norm)
    sm.set_array([])
    cbar_ax = fig.add_axes([0.91, 0.15, 0.02, 0.7])  # Right side of the plot
    # Disable grid for the colorbar axis
    cbar_ax.grid(False)
    cbar = plt.colorbar(sm, cax=cbar_ax)
    cbar.set_label('Log2 Fold Change', rotation=270, labelpad=15)
    # Adjust labels
    ax.set_xlabel(f'{labels}\n{HEV_labels[1]}')
    ax.set_ylabel(f'{labels}\n{HEV_labels[0]}')
    cbar_ax.set_ylim(plot_df['LFC'].min(), plot_df['LFC'].max())
    # Set position
    main_ax_pos = ax.get_position()
    cbar_ax.set_position([main_ax_pos.x1 + main_ax_pos.width*0.02, main_ax_pos.y0, main_ax_pos.width*0.05, main_ax_pos.height])
    # Add title
    ax.set_title(f'Significant LR interactions\nHEVs + spots <{bin_size}µm', fontweight='bold')

def adjust_text_func(ax, sorted_df, min_distance=20, labels_to_annotate=None, lw=0.2, color='gray', fontsize=5):
    texts = []
    annotated_coords = []
    for i, row in sorted_df.iterrows():
        label = row.name
        # If a list of labels is provided, only annotate those
        if labels_to_annotate is not None and label not in labels_to_annotate:
            continue
        x = row[HEV_labels[1]]
        y = row[HEV_labels[0]]
        # Check distance to all previous annotated points
        too_close = False
        for x0, y0 in annotated_coords:
            dist = np.sqrt((x - x0)**2 + (y - y0)**2)
            if dist < min_distance:
                too_close = True
                break
        if not too_close:
            text = ax.text(x, y, label, fontsize=fontsize, ha='right', va='bottom')
            texts.append(text)
            annotated_coords.append((x, y))
    # Now bring text in front of the arrows
    for text in texts:
        text.set_zorder(200)  # Bring text to front
    # Adjust remaining text to avoid overlap
    if texts:
        adjust_text(
            texts,
            ax=ax,
            arrowprops=dict(arrowstyle="-", color=color, lw=lw, shrinkA=5, shrinkB=0),
            expand_points=(1.2, 1.4),
            expand_text=(1.2, 1.4),
            force_text=(0.5, 0.5),
            force_points=(0, 0),
            lim=1000,
            only_move={'points':'', 'text':'xy', 'objects':''}
        )

def draw_rectangles(ax, xlim, ylim):
    # Color the area below 0 lightgray
    # 1. x < 0, all y
    left_rect = Rectangle(
        (xlim[0], ylim[0]),
        width=0 - xlim[0],
        height=ylim[1] - ylim[0],
        facecolor='lightgray',
        alpha=0.5,
        zorder=0
    )
    # 2. x ≥ 0, y < 0
    bottom_right_rect = Rectangle(
        (0, ylim[0]),
        width=xlim[1] - 0,
        height=0 - ylim[0],
        facecolor='lightgray',
        alpha=0.5,
        zorder=0
    )
    # Add patches to the axis
    ax.add_patch(left_rect)
    ax.add_patch(bottom_right_rect)


# --- Prepare for plotting

# Sort by absolute LFC, descending → prioritize extreme values
sorted_df = plot_df.reindex(plot_df['LFC'].abs().sort_values(ascending=False).index)

# Find the most extreme LFC (positive or negative)
max_lfc = max(abs(plot_df['LFC'].min()), abs(plot_df['LFC'].max()))

# Define a midpoint for the TwoSlopeNorm
midpoint = 0
norm = TwoSlopeNorm(vmin=-max_lfc, vcenter=midpoint, vmax=max_lfc)

# Define rename dict for labels
rename_dict = {
    'near_TLS_HEV': 'HEV (TLS)',
    'non_TLS_HEV': 'HEV-like (stroma)',
}


# --- Generate Figure 4G

fig, ax = plt.subplots(figsize=(3, 3))

key_lrs = ['LTB_LTBR', 'CCL19_CCR7', 'CCL21_CCR7']

# Plot all points with their colors first
sns.scatterplot(data=sorted_df.iloc[::-1], x=HEV_labels[1], y=HEV_labels[0], hue='LFC', palette='coolwarm', ax=ax, legend=False,
                 s=10, linewidth=0.2, hue_norm=norm, edgecolor='gray', zorder=1)

# Plot the key LRs on top with black outline
key_df = sorted_df.loc[sorted_df.index.isin(key_lrs)]
sns.scatterplot(data=key_df, x=HEV_labels[1], y=HEV_labels[0], hue='LFC', palette='coolwarm', ax=ax, legend=False,
                 s=10, linewidth=0.5, hue_norm=norm, edgecolor='black', zorder=100)

adjust_scatterplot(ax, sorted_df, max_lfc, labels='Number of groups')

# Get current axis limits
xlim = ax.get_xlim()
ylim = ax.get_ylim()
x_buffer = xlim[1] * 0
y_buffer = ylim[1] * 0
ax.set_xlim(-x_buffer, xlim[1])
ax.set_ylim(-y_buffer, ylim[1])

# Update limits after resetting
xlim = ax.get_xlim()
ylim = ax.get_ylim()

draw_rectangles(ax, xlim, ylim)

adjust_text_func(ax, sorted_df, min_distance=0, labels_to_annotate=key_lrs, lw=0.5, color='black', fontsize=8)

# Edit the labels of the x & y axis
xlabel = ax.get_xlabel()
ylabel = ax.get_ylabel()
for hev_label in HEV_labels:
    if hev_label in xlabel:
        new_xlabel = xlabel.replace(hev_label, rename_dict.get(hev_label, hev_label))
        ax.set_xlabel(new_xlabel)
    if hev_label in ylabel:
        new_ylabel = ylabel.replace(hev_label, rename_dict.get(hev_label, hev_label))
        ax.set_ylabel(new_ylabel)


# --- Save the output

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_4G.pdf'

plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')


# --- Finish

print("Script successfully completed")
print(f'Saved plot to: {os.path.join(out_dir, filename)}')
