#!/usr/bin/env python3

# --- Aim
# Identify and group regions with high immune niche scores (e.g., GC score) in spatial transcriptomics data.
# Visualize spatial continuity, annotate groups, and filter by group size.
# Save spatial plots, groupings, and settings for downstream characterization of immune aggregates.

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
from datetime import datetime
import squidpy as sq

from matplotlib.colors import LinearSegmentedColormap
from mycolorpy import colorlist as mcp
from matplotlib import patheffects
import json


# Import custom functions
from sophiesfunctions.auto_crop import auto_crop_scanpy  # import my custom functions that are packaged locally as sophiesfunctions
from sophiesfunctions.auto_crop import auto_crop  
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
parser.add_argument("--min_spots", help="Minimal spots per gene", type=int)

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

# Create a sample - subtype dictionary
if run_lymph_node:
    subtype = 'lymph_node'
else:
    subtype = tissue_metadata[tissue_metadata['tissue_id'] == str(tissue_id)]['clinical_subtype'].values[0]


# --- Load the adata object

use_data = 'log_norm'

if run_lymph_node:
    adata = ad.read_h5ad(os.path.join(project_dir, 'data', 'lymph_node_anndata_objects', use_data, f'{sample_id}_{use_data}.h5ad'))
else:
    adata = ad.read_h5ad(os.path.join(h5ad_dir, 'raw', f'{joined_id}_raw.h5ad'))


# --- Calculate crop coords for optimal presentation
coords = auto_crop_scanpy(adata, border_size=0.1)

# Also calculate some coords for use w/ squidpy
coord_squidpy = auto_crop(adata, border_size=0.1)

# --- Set spot size for large & small capture area
# Plotting the spots larger to make it easier to see

if slide_size == '11 mm':
    spot_size = 1.8
else:
    spot_size = 1.3


# --- Load TLS results

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

# Find most recent results
TLS_res_dir = os.path.join(res_dir, 'TLS_model')

if run_lymph_node:
    TLS_res_dir = os.path.join(TLS_res_dir, 'lymph_nodes')

most_recent_dir_path = find_most_recent(TLS_res_dir)

TLS_file_name = f'{sample_id}_TLS_model.csv'
TLS_results = pd.read_csv(os.path.join(most_recent_dir_path, TLS_file_name), index_col=0)
TLS_model_key = 'TLS_model'
adata.obs[TLS_model_key] = TLS_results


# --- Make script compatible with multiple scores

clusters_of_i = {
    'GC': 'germinal_center_scoring/GC_score_filtered_genes/GC_score_filtered_genes.csv',
}

for score_name, filename in clusters_of_i.items():

    if score_name == 'GC':
        out_name = 'GC_score_filtered'
    else:
        out_name = f'{score_name}_score'

    # --- Also load the GC scores

    if run_lymph_node:
        GC_res_dir = os.path.join(res_dir, 'TLS_model', 'gene_expression_plots', 'potential_GC_genes', 'scoring', 'lymph_nodes', out_name, 'per_sample',
                                'csv_files')
    else:
        GC_res_dir = os.path.join(res_dir, 'TLS_model', 'gene_expression_plots', 'potential_GC_genes', 'scoring', out_name, 'per_sample',
                                'csv_files')

    most_recent_dir_path = find_most_recent(GC_res_dir)

    # Read results and add to adata
    GC_file_name = f'{sample_id}_{out_name}.csv'
    GC_results = pd.read_csv(os.path.join(most_recent_dir_path, GC_file_name), index_col=0)

    top_gene = 20

    GC_key = f'{score_name}_score_{top_gene}'

    adata.obs[GC_key] = GC_results[GC_key]


    # --- Prepare for saving the plots

    # Put in a timestamped folder to avoid overwriting older plots
    qc_dir = os.path.join(res_dir, 'TLS_model')
    os.makedirs(qc_dir, exist_ok=True)

    today = date.today()
    today = today.strftime("%Y%m%d")


    # --- Plot 1: Plot the TLS model score + spots above a certain threshold

    # Define custom colormap for plotting TLS model results
    # Set everything under certain threshold as the lowest color of a color map

    cmap = "YlOrRd"
    color_list = mcp.gen_color(cmap=cmap, n=11)
    color_list_extracted = ['lightgray'] + color_list[2:]

    n_bins = 100
    cmap_name = 'my_list'
    cmap_adj = LinearSegmentedColormap.from_list(cmap_name, color_list_extracted, N=n_bins)

    cmap_adj.set_under(color='lightgray', alpha=0.4)

    # Also define a threshold for selecting GC-high spots
    threshold = 0.7
    # Subset anndata to GC score above threshold
    adata_sub = adata[adata.obs[GC_key] > threshold]

    # Plotting
    sc.set_figure_params(fontsize=16)

    nrows = 1
    ncols = 3

    fig, axs = plt.subplots(nrows, ncols, figsize=(8*ncols, 6.5*nrows))

    # TLS model results
    sc.pl.spatial(adata, color='TLS_model', size=spot_size, cmap=cmap_adj, bw=False, alpha_img=1, crop_coord=coords, title=f'TLS score', vmin=0, ax=axs[0])

    # GC score results
    vmax = None
    sc.pl.spatial(adata, color=GC_key, size=spot_size, cmap=cmap_adj, bw=False, alpha_img=1, crop_coord=coords, title=f'{score_name} score (vmax={vmax})', vmin=0, vmax=vmax, ax=axs[1])

    # Spots above threshold
    adata_sub.obs['in_tissue'] = adata_sub.obs['in_tissue'].astype('category')
    sc.pl.spatial(adata_sub, color='in_tissue', size=spot_size, palette=['yellow'], bw=False, alpha_img=1, crop_coord=coords, title=f'Spots >{threshold}', vmin=0, 
                ax=axs[2])

    # Remove the legend
    axs[2].get_legend().set_visible(False)

    plt.suptitle(f'{sample_id} ({subtype}) - TLS score + {score_name} score + {score_name}-high spots (GC_key={GC_key})', y=1.05, fontweight='bold', x=0.5)


    # --- If there are no spots over the threshold: need to stop here
    # Or also if there is just a single spot over the threshold

    threshold_adj = str(threshold).split('.')[-1]

    # Prepare for saving
    if run_lymph_node:
        out_dir = os.path.join(qc_dir, 'gene_expression_plots', 'potential_GC_genes', 'scoring', 'lymph_nodes', out_name, 'per_sample', f'threshold_{threshold_adj}')
    else:
        out_dir = os.path.join(qc_dir, 'gene_expression_plots', 'potential_GC_genes', 'scoring', out_name, 'per_sample', f'threshold_{threshold_adj}')

    os.makedirs(out_dir, exist_ok=True)

    # Save plots
    out_dir_plots = os.path.join(out_dir, 'plots', 'groups', today)
    os.makedirs(out_dir_plots, exist_ok=True)

    # Also save the groupings as a .csv
    out_dir_csv = os.path.join(out_dir, 'csv_files', today)
    os.makedirs(out_dir_csv, exist_ok=True)

    # -- Some settings for grouping
    # Decide whether want to sever groups connected by a single edge
    sever_groups = False 

    # Set the smallest allowable group size
    min_spots = 1

    # Save in a dict
    settings = {'score_threshold': threshold,
                'sever_groups': sever_groups,
                'min_spots': min_spots,
                'GC_key': GC_key,
                }

    time_to_stop = False

    if min_spots > 1:
        if len(adata_sub) == 0 or len(adata_sub) == 1:
            time_to_stop = True
    elif len(adata_sub) == 0:
        time_to_stop = True

    if time_to_stop:
        filename = f'{sample_id}_{score_name}_groups.pdf'
        save_multi_image(os.path.join(out_dir_plots, filename))
        to_save_df = adata_sub.obs[['TLS_model', GC_key]]
        to_save_df['group'] = ''
        to_save_df['group_too_small'] = ''
        to_save_df['group_too_small_removed'] = ''
        filename = f'{sample_id}_{score_name}_groups.csv'
        to_save_df.to_csv(os.path.join(out_dir_csv, filename))
        # Save the settings as a JSON file
        with open(os.path.join(out_dir_csv, f'{sample_id}_grouping_settings.json'), 'w') as json_file:
            json.dump(settings, json_file, indent=4)  # indent=4 makes the JSON file more readable
        print('Script succesfully completed')


    # --- Define groups of neighboring spots

    # Get a scalefactor for converting pixel coords to um
    spot_diameter_pixel = adata.uns['spatial'][sample_id]['scalefactors']['spot_diameter_fullres']
    spot_diameter_um = 55 # Regular Visium spots: 55 µm 
    scalefactor_um = spot_diameter_pixel/spot_diameter_um

    # Regular Visium spots are 100 um apart
    spot_dist_um = 100
    spot_dist_pixel = spot_dist_um * scalefactor_um

    # Decide which radius I want to use for grouping spots
    num_spots = 1 # How many spots away is part of same group
    radius_max = spot_dist_pixel*num_spots

    # Use squidpy to define which TLS-high spots are neighbors
    sq.gr.spatial_neighbors(adata_sub, coord_type='generic', radius=(0, radius_max))

    adata_sub.obs['in_tissue'] = adata_sub.obs['in_tissue'].astype('category')


    # --- Plot 2: Plot squidpy's neighborhood graph result

    if slide_size == '11 mm':
        squidpy_size = 4
        squidpy_edge = 0.5
    else:
        squidpy_size = 30
        squidpy_edge = 1

    sc.set_figure_params(fontsize=15, figsize=(6,6))

    sq.pl.spatial_scatter(adata_sub,color='in_tissue', connectivity_key="spatial_connectivities", edges_color="black", crop_coord=coord_squidpy,
                          palette='viridis', edges_width=squidpy_edge, title='Squidpy spatial connectivities')

    plt.suptitle(f'{sample_id} ({subtype}) - neighborhood graph of spots with {score_name} score >{threshold})', y=1.1, fontweight='bold', x=0.5)


    # --- Now annotate the spots as groups

    # Extract connected spots as groups
    connected_groups = []
    for i, neighbors in enumerate(adata_sub.obsp['spatial_connectivities'].toarray()):
        group = [j for j, is_neighbor in enumerate(neighbors) if is_neighbor]
        group.append(i)  # Include the spot itself in the group
        connected_groups.append(group)

    # squidpy's spatial_connectivities just saves which spots are neighbors
    # so there will be overlap between the sublists (if a spot is a neighbor of multiple other spots, it appears multiple times)
    # therefore merge the sublists until all values are unique

    # Define a function to do this
    def extract_neighbor_groups(connected_groups, adata=adata_sub, obs_key='group'):
        pooled = [set(subList) for subList in connected_groups]
        merging = True
        while merging:
            merging=False
            for i,group in enumerate(pooled):
                merged = next((g for g in pooled[i+1:] if g.intersection(group)),None)
                if not merged: continue
                group.update(merged)
                pooled.remove(merged)
                merging = True
        # Convert the list of sets to a list of lists
        list_of_lists = [list(s) for s in pooled]
        # Generate group names
        groups = [f'group{i+1}' for i in range(len(list_of_lists))]
        # Flatten the list of lists and pair each item with its corresponding group
        data_list = []
        for group, sublist in zip(groups, list_of_lists):
            for item in sublist:
                data_list.append({'index': item, 'group': group})
        # Create a DataFrame
        df = pd.DataFrame(data_list)
        # Sort DataFrame based on the 'index' column
        df_sorted = df.sort_values(by='index')
        # Transfer the barcode id to the index from the anndata object
        df_sorted.index = adata.obs_names
        # Tack the group annotation onto the anndata object
        adata.obs[obs_key] = df_sorted['group']
        # Return the grouping
        return list_of_lists

    # Extract groups
    list_of_lists = extract_neighbor_groups(connected_groups, obs_key='group')

    # If each group is just a single spot, need to stop here
    # Check if the length of each nested list is one
    is_length_one = all(len(lst) == 1 for lst in list_of_lists)

    if min_spots > 1:
        if is_length_one:
            filename = f'{sample_id}_{score_name}_groups.pdf'
            save_multi_image(os.path.join(out_dir_plots, filename))
            to_save_df = adata_sub.obs[['TLS_model', GC_key, 'group']]
            to_save_df['groups_severed'] = ''
            to_save_df['group_too_small'] = ''
            to_save_df['group_too_small_removed'] = ''
            filename = f'{sample_id}_{score_name}_groups.csv'
            to_save_df.to_csv(os.path.join(out_dir_csv, filename))
            print('Script succesfully completed')


    # --- Plot 2b: plot the neighborhood graph again after grouping

    sc.set_figure_params(fontsize=15, figsize=(6,6))

    group_key = 'group'

    # Need to convert to categorical otherwise it runs into an error for some reason
    adata_sub.obs[group_key] = adata_sub.obs[group_key].astype('category')

    # For some reason this doesn't work if just one spot
    if len(adata_sub) > 1:

        sq.pl.spatial_scatter(adata_sub, color=group_key, connectivity_key="spatial_connectivities", edges_color="black", crop_coord=coord_squidpy,
                            palette='viridis', edges_width=squidpy_edge, title='Squidpy spatial connectivities', img_alpha=0.5, outline=True,
                            outline_color=('black', 'black'), outline_width=(0, 0.05))

        plt.suptitle(f'{sample_id} ({subtype}) - neighborhood graph of spots with {score_name} score >{threshold})', y=1.1, fontweight='bold', x=0.5)


    # --- Plot 3: Plot the grouping of high GC spots

    sc.set_figure_params(fontsize=12, figsize=(6,6))

    sc.pl.spatial(adata_sub, color=group_key, size=spot_size, bw=False, alpha_img=1, crop_coord=coords, title=f'{sample_id} - groups of spots >{threshold}', )

    plt.suptitle(f'{sample_id} ({subtype}) - grouping of spots with {score_name} score >{threshold})', y=1.05, fontweight='bold', x=0.5)


    # --- Now sever single edges connecting spots
    # Sometimes there are 2 GCs that are so close they're touching
    # Try to deal with these by finding groups of spots that are only connected by a single edge
    # Sever this single edge so the groups appear as 2 separate groups despite touching

    # Define a function to do this

    def identify_single_edge_neighbors(list_of_lists, connected_groups):
        single_edge_neighbors = []
        for group in list_of_lists:
            if len(group) == 1:
                print('Singleton group')
            else:
                print(f'Checking group {group}')
                for spot_index in group:
                    print(f'Checking spot {spot_index}')
                    # Get the neighbors of the spot
                    neighbors = set(connected_groups[spot_index])
                    # Remove the spot itself from the nieghbors
                    neighbors = neighbors - set([spot_index])
                    # Check that there isn't just one neighbor
                    if len(neighbors) == 1:
                        print(f'Spot {spot_index} has only one neighbor (spot {neighbors})')
                    else:
                        # Loop through the neighbors and get their neighbors
                        for neighbor_index in neighbors:
                            neighbors_of_neighbor_index = set(connected_groups[neighbor_index])
                            # Remove the spot and neighbor index
                            neighbors_of_neighbor_index = neighbors_of_neighbor_index - set([spot_index, neighbor_index])
                            # Check that neighobrs_of_neighbor_index is not empty
                            if len(neighbors_of_neighbor_index) == 0:
                                print(f'Spot {neighbor_index} has no neighbors except for index spot {spot_index}')
                            else:
                                # First check if spots have neighbors in common:
                                overlap = neighbors_of_neighbor_index.intersection(neighbors - set([spot_index, neighbor_index]))
                                if not overlap:
                                    print(f'No neighbors in common between {spot_index} & {neighbor_index}')
                                    single_edge_neighbors.append([spot_index, neighbor_index])
        return single_edge_neighbors

    # Cut groups connected by single edge (optional)
    if sever_groups:

        # Find the single-edge neighbors
        single_edge_neighbors = identify_single_edge_neighbors(list_of_lists, connected_groups)

        # Now make it seem like these pairs aren't neighbors
        updated_connected_groups = connected_groups.copy()

        for pair in single_edge_neighbors:
            # For each pair of single-edge neighbors: make it seem like these are not neighbors
            index = pair[0]
            neighbor = pair[1]
            # Remove the neighbor for this index
            updated_neighbors = [num for num in connected_groups[index] if num != neighbor]
            updated_connected_groups[index] = updated_neighbors

        # Extract new groups
        group_key = 'groups_severed'
        extract_neighbor_groups(updated_connected_groups, obs_key=group_key)


        # --- Plot 4: Plot the grouping of high GC spots

        sc.set_figure_params(fontsize=12, figsize=(6,6))

        sc.pl.spatial(adata_sub, color=group_key, size=spot_size, bw=False, alpha_img=1, crop_coord=coords, title=f'{sample_id} - groups of spots >{threshold}', )

        plt.suptitle(f'{sample_id} ({subtype}) - grouping of spots with {score_name} score >{threshold}) - single-edge connections removed', y=1.05, fontweight='bold', x=0.5)


    # --- Only keep groups with a minimum number of spots 

    # Find groups that are too small
    too_small = adata_sub.obs[group_key].value_counts()[adata_sub.obs[group_key].value_counts() < min_spots].index

    # Make a temporary df
    temp_df = adata_sub.obs[[group_key]]

    # Set small groups as too_small
    too_small_value = 'too_small'
    temp_df[group_key] = np.where(temp_df[group_key].isin(too_small), too_small_value, temp_df[group_key])

    # Save singletons as separat df
    too_small = temp_df[temp_df[group_key] == too_small_value]

    # Remove the singletons for now
    temp_df = temp_df[temp_df[group_key] != too_small_value]

    # Update group names based on size
    group_counts = temp_df[group_key].value_counts()
    sorted_groups = group_counts.sort_values(ascending=False).index
    group_map = {group: f"group{i+1}" for i, group in enumerate(sorted_groups)}

    temp_df[group_key] = temp_df[group_key].map(group_map)

    # Add results to a new column
    adata_sub.obs['group_too_small_removed'] = temp_df[group_key]

    # Add groups + singletons to another column
    group_too_small = pd.concat([temp_df, too_small])

    adata_sub.obs['group_too_small'] = group_too_small


    # --- For GC score, also check that single-spot groups have an acceptable B cell zone score

    if (score_name == 'GC') & (min_spots == 1):
        # Check whether there are any groups of length 1
        if any(len(sublist) == 1 for sublist in list_of_lists):
            
            int_score_name = 'B_cell_zone'
            int_out_name = f'{int_score_name}_score'
            
            # Load the intermediate TLS score
            if run_lymph_node:
                GC_res_dir = os.path.join(res_dir, 'TLS_model', 'gene_expression_plots', 'potential_GC_genes', 'scoring', 'lymph_nodes', int_out_name, 'per_sample',
                                        'csv_files')
            else:
                GC_res_dir = os.path.join(res_dir, 'TLS_model', 'gene_expression_plots', 'potential_GC_genes', 'scoring', int_out_name, 'per_sample',
                                        'csv_files')

            most_recent_dir_path = find_most_recent(GC_res_dir)

            # Read results and add to adata
            GC_file_name = f'{sample_id}_{int_out_name}.csv'
            GC_results = pd.read_csv(os.path.join(most_recent_dir_path, GC_file_name), index_col=0)
            int_TLS_key = f'{int_score_name}_score_{top_gene}'
            adata_sub.obs[int_TLS_key] = GC_results[int_TLS_key]
            
            # Find single-spot groups & check intermediate TLS score for these
            value_counts = adata_sub.obs['group_too_small'].value_counts()
            items_only_once = value_counts[value_counts == 1].index.tolist()

            # Set intermediate TLS threshold
            int_TLS_threshold = 0.6
            filtered_items = list(adata_sub.obs[(adata_sub.obs['group_too_small'].isin(items_only_once)) & (adata_sub.obs[int_TLS_key] > int_TLS_threshold)]['group_too_small'].unique())

            # Remove groups that don't meet the TLS score threshold
            for group in items_only_once:
                if not group in filtered_items:
                    adata_sub.obs['group_too_small'] = np.where(adata_sub.obs['group_too_small'] == group, 'TLS_score_too_low', adata_sub.obs['group_too_small'])

            temp_df = adata_sub.obs[['group_too_small']]
            temp_df = temp_df[temp_df['group_too_small'] != 'TLS_score_too_low']

            # Update group names based on size
            group_counts = temp_df['group_too_small'].value_counts()
            sorted_groups = group_counts.sort_values(ascending=False).index
            group_map = {group: f"group{i+1}" for i, group in enumerate(sorted_groups)}

            temp_df['group_too_small'] = temp_df['group_too_small'].map(group_map)
            adata_sub.obs['group_too_small'] = adata_sub.obs['group_too_small'].map(group_map)
            adata_sub.obs['group_too_small'] = adata_sub.obs['group_too_small'].fillna('TLS_score_too_low')
            adata_sub.obs['group_too_small_removed'] = temp_df['group_too_small']
            
            # Record threshold in settings
            settings['int_TLS_score_threshold'] = int_TLS_threshold


    # --- Plot 5: Groups after removing small groups

    # Need to convert to categorical otherwise it runs into an error for some reason
    adata_sub.obs['group_too_small'] = adata_sub.obs['group_too_small'].astype('category')

    sc.pl.spatial(adata_sub, color='group_too_small', size=spot_size, bw=False, alpha_img=1, crop_coord=coords, title=f'{sample_id} - GC-high spots', )

    plt.suptitle(f'{sample_id} ({subtype}) - grouping of spots with {score_name} score >{threshold} - groups too small removed (min_spots={min_spots})', y=1.05, fontweight='bold', x=0.5)

    # Annotate the spatial groups plot
    scalef = adata.uns['spatial'][sample_id]['scalefactors']['tissue_hires_scalef']

    # Define a function for annotating on the spatial plot
    def annotate_spatial_groups(adata, scalef, group_key='group_too_small'):
        # Find the midpoint for each group
        group_midpoint = {}
        all_groups = list(adata_sub.obs[group_key].cat.categories)
        # Do not annotate groups that are too small or don't have a TLS score high enough
        if too_small_value in all_groups:
            all_groups.remove(too_small_value)
        if 'TLS_score_too_low' in all_groups:
            all_groups.remove('TLS_score_too_low')
        for group in all_groups:
            adata_tmp = adata_sub[adata_sub.obs[group_key] == group]
            coords_tmp = auto_crop_scanpy(adata_tmp)
            # Need to muliply by hires scalefactor to get the plot coordinates
            coords_adj = tuple(x * scalef for x in coords_tmp)
            x1, x2, y1, y2 = coords_adj
            midpoint_x = (x1 + x2) / 2
            midpoint_y = (y1 + y2) / 2
            group_midpoint[group] = (midpoint_x, midpoint_y)
        # Annotate the midpoint
        ax = plt.gca()
        for group in all_groups:
            group_adj = group[5:]
            ax.annotate(group_adj, xy=group_midpoint[group],ha='center', va='center', fontsize=8,
                            path_effects=[patheffects.withStroke(linewidth=2, foreground='white')])

    annotate_spatial_groups(adata_sub, scalef)


    # --- Save everything so far

    # Save plots
    filename = f'{sample_id}_{score_name}_groups.pdf'

    save_multi_image(os.path.join(out_dir_plots, filename))
    plt.close('all')

    # And save the groupings as a csv
    to_save_df = adata_sub.obs[['TLS_model', GC_key, 'group', 'group_too_small', 'group_too_small_removed']]

    filename = f'{sample_id}_{score_name}_groups.csv'
    to_save_df.to_csv(os.path.join(out_dir_csv, filename))

    # Save the grouping settings as a JSON file
    with open(os.path.join(out_dir_csv, f'{sample_id}_grouping_settings.json'), 'w') as json_file:
        json.dump(settings, json_file, indent=4)  # indent=4 makes the JSON file more readable


# --- Finish

print("Script successfully completed")
