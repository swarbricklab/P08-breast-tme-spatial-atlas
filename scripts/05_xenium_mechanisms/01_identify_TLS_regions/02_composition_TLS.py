#!/usr/bin/env python3

# --- Aim
# Summarise the cellular composition of the identified Xenium TLS regions across the cohort.
# Produces the proportion of each cell type located inside vs outside TLS regions (S7B) and
# the average relative abundance of each cell type among cells within TLS regions (S7C).

# --- Define environment
# env: xen-hev

# --- Load packages

import argparse
import os
import sys

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Cellular composition of Xenium TLS regions (S7B, S7C)",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-s", "--samples",
    nargs="*",
    type=str,
    help="List of Xenium sample IDs (section_id)")
parser.add_argument("-r", "--res_dir", help="Path to results directory")
parser.add_argument("--h5ad_dir", help="Path to directory with .h5ad objects")
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)
parser.add_argument("--sc_ref", help="Single-cell reference used for cell-type classifications", default='miniatlas')
parser.add_argument("--annotation_levels", nargs="*", default=['major', 'minor', 'subset'],
                    help="Annotation levels to summarise")

args = vars(parser.parse_args())

project_dir = args["project_dir"]
res_dir = args["res_dir"]
sample_list = args["samples"]
sc_ref = args["sc_ref"]
levels = args["annotation_levels"]


# --- Load the cell-type classifications
# Used to align annotation levels, derive lineage ordering and short subset names.

classif_dir = os.path.join(project_dir, 'resources', 'cell_types', sc_ref)
classif_file = f'{sc_ref}_subset_classifications.csv'
classif_df = pd.read_csv(os.path.join(classif_dir, classif_file))


# --- Load plotting colours

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
color_dict = pd.read_csv(os.path.join(color_dir, f'{sc_ref}_major_colors.csv'),
                         index_col=0).to_dict()['color']

# Propagate the major-level colour to each minor / subset cell type
for level in levels:
    for cell_type in classif_df[f'{sc_ref}_{level}'].unique():
        if pd.isna(cell_type):
            continue
        if cell_type in color_dict:
            continue
        major_type = classif_df[classif_df[f'{sc_ref}_{level}'] == cell_type][f'{sc_ref}_major'].values
        major_type_in_dict = [m for m in major_type if m in color_dict]
        if major_type_in_dict:
            color_dict[cell_type] = color_dict[major_type_in_dict[0]]
        else:
            color_dict[cell_type] = 'lightgrey'

# Short names for subset-level cell types
sub_short_dict = classif_df.set_index(f'{sc_ref}_subset')[[f'{sc_ref}_subset_short']].to_dict()


# --- Load the cell-type labels and TLS annotations per sample

def load_label_transfer(sample, level):
    annot_path = os.path.join(res_dir, '02_cell_type_mapping', 'label_transfer', 'xenium',
                              f'{level}_level', f'{sample}_{level}_pred_cell_types.csv')
    return pd.read_csv(annot_path, index_col=0)['pred_cell_type']

tls_csv_dir = os.path.join(res_dir, '05_xenium_mechanisms', '01_identify_TLS_regions', 'csv_files')

tls_df = pd.DataFrame()

for sample in sample_list:
    tmp_df = pd.read_csv(os.path.join(tls_csv_dir, f'{sample}_cell_TLS_annotations.csv'), index_col=0)
    # Prepend sample id to keep cell indices unique across samples
    tmp_df.index = [f'{sample}_{idx}' for idx in tmp_df.index]
    tmp_df['sample_id'] = sample
    for level in levels:
        annotations = load_label_transfer(sample, level)
        annotations.index = [f'{sample}_{idx}' for idx in annotations.index]
        tmp_df[level] = annotations.reindex(tmp_df.index)
    tls_df = pd.concat([tls_df, tmp_df])

# Normalise the in_TLS flag to boolean (stored as text in the annotation CSV)
tls_df['in_TLS'] = tls_df['in_TLS'].astype(str).map({'True': True, 'False': False}).fillna(False).astype(bool)


# --- Sort cell types by lineage, then by major type within lineage

lineage_order = ['Epithelial', 'Immune', 'Stromal']

major_cell_types = classif_df[f'{sc_ref}_major'].dropna().unique().tolist()
minor_cell_types = classif_df[f'{sc_ref}_minor'].dropna().unique().tolist()
subset_cell_types = classif_df[f'{sc_ref}_subset'].dropna().unique().tolist()


def get_lineage(cell_type, level='major'):
    if cell_type in classif_df[f'{sc_ref}_{level}'].values:
        return classif_df[classif_df[f'{sc_ref}_{level}'] == cell_type]['lineage'].values[0]
    return 'other'


def get_major(cell_type, level, upper_level='major'):
    if cell_type in classif_df[f'{sc_ref}_{level}'].values:
        return classif_df[classif_df[f'{sc_ref}_{level}'] == cell_type][f'{sc_ref}_{upper_level}'].values[0]
    return 'other'


major_cell_types_sorted = sorted(
    major_cell_types,
    key=lambda x: (lineage_order.index(get_lineage(x)) if get_lineage(x) in lineage_order else len(lineage_order), x))
minor_cell_types_sorted = sorted(
    minor_cell_types,
    key=lambda x: (major_cell_types_sorted.index(get_major(x, level='minor'))
                   if get_major(x, level='minor') in major_cell_types_sorted else len(major_cell_types_sorted), x))
subset_cell_types_sorted = sorted(
    subset_cell_types,
    key=lambda x: (major_cell_types_sorted.index(get_major(x, level='subset'))
                   if get_major(x, level='subset') in major_cell_types_sorted else len(major_cell_types_sorted), x))


# --- Plot 1 (S7B): proportion of each cell type located inside vs outside TLS regions

def inside_outside_proportions(tls_df, annotation_key):
    """Per-sample proportion of each cell type's cells that fall inside TLS regions."""
    sample_prop_all = pd.DataFrame()
    for sample_id in sample_list:
        sample_tls_df = tls_df[tls_df['sample_id'] == sample_id]
        sample_cell_counts = sample_tls_df.groupby([annotation_key, 'in_TLS']).size().unstack(fill_value=0)
        sample_cell_counts = sample_cell_counts.reindex(columns=[True, False], fill_value=0)
        sample_cell_props = sample_cell_counts.div(sample_cell_counts.sum(axis=1), axis=0)
        sample_cell_props.columns = sample_cell_props.columns.astype(str)
        sample_cell_props = sample_cell_props[['True', 'False']]
        sample_cell_props['sample_id'] = sample_id
        sample_cell_props = sample_cell_props.reset_index()
        sample_prop_all = pd.concat([sample_prop_all, sample_cell_props])
    mean_true = sample_prop_all.groupby(annotation_key)['True'].mean()
    std_true = sample_prop_all.groupby(annotation_key)['True'].std()
    mean_false = 1 - mean_true
    return mean_true, std_true, mean_false


ncols = len(levels)
nrows = 1

plt.rcParams['axes.axisbelow'] = True

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4, 2))
if ncols == 1:
    axs = [axs]

for n, level in enumerate(levels):
    ax = axs[n]
    mean_true, std_true, mean_false = inside_outside_proportions(tls_df, level)
    sorted_cells = mean_true.sort_values(ascending=False).index
    mean_true = mean_true.loc[sorted_cells]
    std_true = std_true.loc[sorted_cells]
    mean_false = mean_false.loc[sorted_cells]
    colors = [color_dict.get(cell_type, 'lightgrey') for cell_type in sorted_cells]
    ax.bar(sorted_cells, mean_true, color=colors, yerr=std_true, ecolor='black',
           error_kw=dict(lw=0.8, capsize=3, capthick=1), width=0.9)
    ax.bar(sorted_cells, mean_false, bottom=mean_true, label='Outside TLS', color='lightgray', width=0.9)
    ax.set_title(f'Cells in TLS regions\n(Annotations: {level})', fontweight='bold')
    ax.set_ylabel('Proportion of cells')

for ax in axs:
    ax.set_xlim(-0.7, len(ax.get_xticklabels()) - 0.3)
    ax.set_ylim(0, 1.02)

for ax in axs:
    num_cells = len(ax.get_xticklabels())
    ax_pos = ax.get_position()
    ax.set_position([ax_pos.x0, ax_pos.y0, 0.05 * num_cells, ax_pos.height])

for n, ax in enumerate(axs):
    if n > 0:
        ax_pos = ax.get_position()
        prev_ax_pos = axs[n - 1].get_position()
        new_x0 = prev_ax_pos.x1 + 0.2
        ax.set_position([new_x0, ax_pos.y0, ax_pos.width, ax_pos.height])

for ax in axs:
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

# Replace subset-level tick labels with short names on the last panel
if 'subset' in levels:
    last_ax = axs[levels.index('subset')]
    new_labels = [sub_short_dict[f'{sc_ref}_subset_short'].get(label.get_text(), label.get_text())
                  for label in last_ax.get_xticklabels()]
    last_ax.set_xticklabels(new_labels)

labels = major_cell_types_sorted
handles = [plt.Rectangle((0, 0), 0.5, 0.5, color=color_dict[label]) for label in labels if label in color_dict]
axs[-1].legend(handles, labels, title='', bbox_to_anchor=(1, 0.5), loc='center left',
               ncol=1, frameon=False, handlelength=0.8, fontsize=9)

title_x = (axs[0].get_position().x0 + axs[-1].get_position().x1) / 2
title_y = axs[0].get_position().y1 + 0.4
plt.suptitle('All samples - Cells inside vs outside of TLS regions', x=title_x, y=title_y, fontweight='bold')


# --- Plot 2 (S7C): average relative abundance of each cell type within TLS regions

def tls_composition(tls_df, annotation_key):
    """Per-sample relative abundance of each cell type among cells inside TLS regions."""
    avg_tls_compos = pd.DataFrame()
    for sample in sample_list:
        sample_tls_df = tls_df[tls_df['sample_id'] == sample]
        sample_cell_counts = sample_tls_df[sample_tls_df['in_TLS']].groupby([annotation_key]).size()
        sample_cell_props = sample_cell_counts / sample_cell_counts.sum()
        sample_cell_props.name = sample
        avg_tls_compos = pd.concat([avg_tls_compos, sample_cell_props], axis=1)
    return avg_tls_compos


plt.rcParams['axes.axisbelow'] = True

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4, 2))
if ncols == 1:
    axs = [axs]

for n, level in enumerate(levels):
    ax = axs[n]
    avg_tls_compos = tls_composition(tls_df, level)
    mean_true = avg_tls_compos.mean(axis=1)
    std_true = avg_tls_compos.std(axis=1)
    sorted_cells = mean_true.sort_values(ascending=False).index
    mean_true = mean_true.loc[sorted_cells]
    std_true = std_true.loc[sorted_cells]
    colors = [color_dict.get(cell_type, 'lightgrey') for cell_type in sorted_cells]
    ax.bar(sorted_cells, mean_true, color=colors, yerr=std_true, ecolor='black',
           error_kw=dict(lw=0.8, capsize=3, capthick=1), width=0.9)
    ax.set_title(f'TLS composition\n(Annotations: {level})', fontweight='bold', fontsize=10)
    ax.set_ylabel('Proportion of cells')

for ax in axs:
    ax.set_xlim(-0.7, len(ax.get_xticklabels()) - 0.3)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.1)

for ax in axs:
    num_cells = len(ax.get_xticklabels())
    ax_pos = ax.get_position()
    ax.set_position([ax_pos.x0, ax_pos.y0, 0.05 * num_cells, ax_pos.height])
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.7)
    ax.xaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.7)

for n, ax in enumerate(axs):
    if n > 0:
        ax_pos = ax.get_position()
        prev_ax_pos = axs[n - 1].get_position()
        new_x0 = prev_ax_pos.x1 + 0.2
        ax.set_position([new_x0, ax_pos.y0, ax_pos.width, ax_pos.height])

for ax in axs:
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

if 'subset' in levels:
    last_ax = axs[levels.index('subset')]
    new_labels = [sub_short_dict[f'{sc_ref}_subset_short'].get(label.get_text(), label.get_text())
                  for label in last_ax.get_xticklabels()]
    last_ax.set_xticklabels(new_labels)

handles = [plt.Rectangle((0, 0), 0.5, 0.5, color=color_dict[label]) for label in labels if label in color_dict]
axs[-1].legend(handles, labels, title='', bbox_to_anchor=(1, 0.5), loc='center left',
               ncol=1, frameon=False, handlelength=0.8, fontsize=9)

title_x = (axs[0].get_position().x0 + axs[-1].get_position().x1) / 2
title_y = axs[0].get_position().y1 + 0.4
plt.suptitle('All samples - Average TLS composition', x=title_x, y=title_y, fontweight='bold')


# --- Save the output

sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import save_multi_image

out_dir = os.path.join(res_dir, '05_xenium_mechanisms', '01_identify_TLS_regions',
                       'all_samples', 'cellular_composition')
plots_dir = os.path.join(out_dir, 'plots')
csv_dir = os.path.join(out_dir, 'csv_files')
os.makedirs(plots_dir, exist_ok=True)
os.makedirs(csv_dir, exist_ok=True)

save_multi_image(os.path.join(plots_dir, 'All_samples_TLS_cellular_composition.pdf'))

# Save the underlying per-cell-type statistics for each annotation level
for level in levels:
    mean_true, std_true, _ = inside_outside_proportions(tls_df, level)
    inside_stats = pd.DataFrame({'mean_prop_inside_TLS': mean_true, 'std_prop_inside_TLS': std_true})
    inside_stats.to_csv(os.path.join(csv_dir, f'All_samples_{level}_inside_vs_outside_TLS.csv'))

    compos = tls_composition(tls_df, level)
    compos_stats = pd.DataFrame({'mean_abundance_in_TLS': compos.mean(axis=1),
                                 'std_abundance_in_TLS': compos.std(axis=1)})
    compos_stats.to_csv(os.path.join(csv_dir, f'All_samples_{level}_TLS_composition.csv'))


# --- Finish

print("Script successfully completed")
