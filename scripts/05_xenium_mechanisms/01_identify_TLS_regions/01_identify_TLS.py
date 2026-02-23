#!/usr/bin/env python3

# --- Aim
# Automatically identify tertiary lymphoid structure (TLS) regions in Xenium data based on dense spatial clustering of T and B cells.
# Output spatial segmentations and annotations for downstream microenvironment analysis.


# --- Define environment
# env: xen-hev

# --- Load packages

import argparse
import os
import sys

import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from scipy.ndimage import gaussian_filter, label
from skimage.measure import regionprops_table
from skimage.morphology import remove_small_objects
from shapely.geometry import Polygon
import geopandas as gpd

from sophiesfunctions.auto_crop import auto_crop_xenium

sc.set_figure_params(facecolor="white", figsize=(8, 8), fontsize=20)


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium data & save as .h5ad file",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required
parser.add_argument("--raw_data_dir", help="Path to raw data")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--project_dir", help="Path to the root of the project directory")
parser.add_argument("--min_count", help="Minimal number of counts per spot", type=int)
parser.add_argument("--min_gene", help="Minimal number of genes per spot", type=int)
parser.add_argument("--outs_path", help="Path to spaceranger 'outs'")
parser.add_argument("--use_data", help="Which data to use (log_norm or SCT)")
parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--seed", help="Seed for reproducibility", type=int)
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)


args = vars(parser.parse_args())


# --- Arguments

h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]

res_dir = args['res_dir']
project_dir = args['project_dir']

# Import function script as a module
sys.path.append(os.path.join(project_dir, 'scripts', '05_xenium_mechanisms'))
import plotting_functions


# --- Extract some info from the sample sheet
# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_dir = os.path.join(project_dir, 'config')
    sample_sheet_file = 'sample_sheet_xenium.csv'
    sample_sheet_path = os.path.join(sample_dir, sample_sheet_file)

sample_sheet = pd.read_csv(sample_sheet_path)
tissue_id = sample_sheet[sample_sheet['section_id'] == sample_id]['tissue_id'].values[0]

# Also load tissue metadata
metadata_dir = os.path.join(project_dir, 'resources',)
tissue_metadata_file = 'tissue_metadata.csv'
tissue_metadata = pd.read_csv(os.path.join(metadata_dir, tissue_metadata_file))

subtype = tissue_metadata[tissue_metadata['tissue_id'].astype(str) == str(tissue_id)]['clinical_subtype'].values[0]


# --- Load data

adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample_id}.h5ad'))


# --- Calculate crop coordinates for optimal presentation

crop_coord = auto_crop_xenium(adata)


# --- Load the tangram annotations

# Decide which annotation level to use
level = 'major'

annotation_dir = os.path.join(project_dir, 'results', '02_cell_type_mapping', 'label_transfer', 'xenium', f'{level}_level')
annotation_file = f'{sample_id}_{level}_pred_cell_types.csv'

annotations = pd.read_csv(os.path.join(annotation_dir, annotation_file), index_col=0)

# Add to anndata object
annotation_key = 'tangram_annotation'
probability_key = 'tangram_probability'

adata.obs[annotation_key] = annotations['pred_cell_type']
adata.obs[probability_key] = annotations['pred_cell_type_prob']


# --- Load plotting colors

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
color_file = f'miniatlas_{level}_colors.csv'

color_dict = pd.read_csv(os.path.join(color_dir, color_file), index_col=0).to_dict()['color']

all_cell_types = adata.obs[annotation_key].astype('category').cat.categories

# Make sure all cell types have a color
for cell_type in all_cell_types:
    if cell_type not in color_dict.keys():
        print(f'Adding color for {cell_type}')
        color_dict[cell_type] = 'lightgrey'


# --- Plot 1: Overview of T & B cell annotations

cell_types_of_interest = ['B-cells', 'T-cells']

ncols = 1
nrows = 1

plt.figure(figsize=(9*ncols, 7*nrows))

# add a new subplot iteratively
n=0
ax = plt.subplot(nrows, ncols, n + 1)

sc.pl.spatial(adata, color=[annotation_key], groups=cell_types_of_interest, title=f'{sample_id} ({subtype})',
            crop_coord=crop_coord, ax=ax, palette=color_dict, na_in_legend=False, alpha_img=0.3, size=0.5)

# Remove the labels
ax.set_xlabel('')
ax.set_ylabel('')

ax.set_title(f'{sample_id} ({subtype})\nT & B cells', fontweight='bold')

plt.suptitle(f'{sample_id} ({subtype}) - T & B cells based on Tangram label transfer (annotation_level={level})', y=1.1, fontweight='bold')


# --- Identify TLS regions based on density of T & B cells

# Extract the spatial coordinates in microns (not pixel space) for more interpretable distance calculations
xy = adata.obs[['x_centroid', 'y_centroid']].values
ct = adata.obs[annotation_key].astype(str)
t_labels = {'T-cells'}
b_labels = {'B-cells'}

grid_um = 25
sigma_um = 60.0         # smoothing (Gaussian sigma, microns); set None/0 to skip
min_each_count = 0      # at least N T and N B per passing cell
dominance_thresh = 0.1 # (T+B)/total >= this
t_min_prop = 0.05       # local T proportion >= this
b_min_prop = 0.05       # local B proportion >= this
min_region_area_mm2 = 0.01 # Minimum area of a TLS region to keep (in mm^2)
min_total_tb = 2  # Minimum number of total T and B cells in a grid cell to consider it for TLS identification
min_total_cells = 6  # Minimum number of total cells in a grid cell

# ----------------------------
# BINNING
# ----------------------------
xmin, ymin = xy.min(axis=0)
xmax, ymax = xy.max(axis=0)
nx = int(np.ceil((xmax - xmin) / grid_um))
ny = int(np.ceil((ymax - ymin) / grid_um))

# map cells to grid
ix = np.minimum(((xy[:,0] - xmin) / grid_um).astype(int), nx-1)
iy = np.minimum(((xy[:,1] - ymin) / grid_um).astype(int), ny-1)

Tmask = ct.isin(t_labels).values
Bmask = ct.isin(b_labels).values

# counts per grid cell
total = np.zeros((ny, nx), dtype=int)
tcount = np.zeros_like(total)
bcount = np.zeros_like(total)
for x_i, y_i, isT, isB in zip(ix, iy, Tmask, Bmask):
    total[y_i, x_i] += 1
    if isT: tcount[y_i, x_i] += 1
    if isB: bcount[y_i, x_i] += 1

# optional smoothing (approximate KDE)
if sigma_um and sigma_um > 0:
    sigma_px = sigma_um / grid_um
    total_s = gaussian_filter(total.astype(float), sigma=sigma_px, mode="nearest")
    t_s = gaussian_filter(tcount.astype(float), sigma=sigma_px, mode="nearest")
    b_s = gaussian_filter(bcount.astype(float), sigma=sigma_px, mode="nearest")
else:
    total_s, t_s, b_s = total, tcount, bcount

enough_total_cells = total_s >= min_total_cells

# proportions & basic presence
eps = 1e-9
t_prop = t_s / (total_s + eps)
b_prop = b_s / (total_s + eps)

# dominance & presence criteria
enough_total_tb = (t_s + b_s) >= min_total_tb
both_present = (t_s >= min_each_count) & (b_s >= min_each_count) & enough_total_tb
dominant = ((t_s + b_s) / (total_s + eps) >= dominance_thresh)
mix_ok = (t_prop >= t_min_prop) & (b_prop >= b_min_prop)

tls_mask = both_present & dominant & mix_ok & enough_total_cells

# clean small speckles
labels, nlab = label(tls_mask)
labels = remove_small_objects(labels, min_size=4, connectivity=1)  # 4 grid cells ~ (grid)^2 scale

# region stats
props = regionprops_table(labels, properties=("label", "area", "bbox"))
props = pd.DataFrame(props)

# compute area in mm^2
cell_area_um2 = grid_um * grid_um
props["area_mm2"] = props["area"] * cell_area_um2 / 1e6

# keep only sufficiently big regions
keep = props["area_mm2"] >= min_region_area_mm2
keep_labels = set(props.loc[keep, "label"].tolist())

tls_regions_mask = np.isin(labels, list(keep_labels)) & (labels > 0)

# build polygon per region (grid-cell rectangles merged)
def region_polygon(mask):
    # build rectangle per True grid cell
    polys = []
    ys, xs = np.where(mask)
    for y, x in zip(ys, xs):
        x0 = xmin + x*grid_um
        y0 = ymin + y*grid_um
        x1 = x0 + grid_um
        y1 = y0 + grid_um
        polys.append(Polygon([(x0,y0),(x1,y0),(x1,y1),(x0,y1)]))
    return gpd.GeoSeries(polys).unary_union.buffer(0)

geoms = []
for lab in keep_labels:
    geoms.append(region_polygon(labels == lab))

tls_gdf = gpd.GeoDataFrame({"tls_id": range(1, len(keep_labels)+1)}, geometry=geoms, crs="EPSG:3857")  # CRS arbitrary in microns


# --- Annotate which cells are in TLS regions

# Create a dataframe with the cell coordinates
cells_df = adata.obs[['x_centroid', 'y_centroid']].copy()
cells_df.columns = ['x', 'y']
cells_df['cell_id'] = adata.obs_names.values

cells_gdf = gpd.GeoDataFrame(cells_df, geometry=gpd.points_from_xy(cells_df.x, cells_df.y), crs="EPSG:3857")

# Spatial join to find which cells are in TLS regions
cells_in_tls = gpd.sjoin(cells_gdf, tls_gdf, predicate='intersects')

# Add to adata
adata.obs['in_TLS'] = adata.obs_names.isin(cells_in_tls['cell_id']).astype(str)
adata.obs['TLS_id'] = adata.obs_names.map(cells_in_tls.set_index('cell_id')['tls_id']).fillna(0).astype(int)

print(f'Total cells in TLS regions: {adata.obs["in_TLS"].value_counts().get("True", 0)}')


# --- Plot 2: Overview of TLS regions identified based on T & B cell densities

# Plot the TLS regions on top of the spatial plot

ncols = 3
nrows = 1

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9*ncols, 7*nrows))

# add a new subplot iteratively
n=0
ax = axs[n]
sc.pl.spatial(adata, color=[annotation_key], groups=cell_types_of_interest, title=f'{sample_id} ({subtype})',
            crop_coord=crop_coord, ax=ax, palette=color_dict, na_in_legend=False, alpha_img=0.3, size=0.5)

ax.set_title(f'{sample_id} ({subtype})\nTLS regions', fontweight='bold')

# Remove the labels
ax.set_xlabel('')
ax.set_ylabel('')

# If plotting on image-level coordinates, transform TLS polygons
if 'xenium_um_spatial' in adata.obsm and 'spatial' in adata.uns and sample_id in adata.uns['spatial']:
    try:
        print("Transforming TLS polygons to image space...")
        tls_gdf_plot = plotting_functions.transform_geojson_to_image_space(tls_gdf.copy(), adata, sample_id, project_dir)
        print("Transformation successful!")
        plot_tls_polygons = True
    except Exception as e:
        print(f"Warning: Could not transform TLS polygons to image space: {e}")
        plot_tls_polygons = False
else:
    tls_gdf_plot = tls_gdf
    plot_tls_polygons = True

# Plot TLS regions
if plot_tls_polygons:
    tls_gdf_plot.plot(ax=ax, facecolor="none", edgecolor="black", linewidth=1)

# Add another ax next to it
ax = axs[n+1]

# Plot all cells in light grey
sc.pl.spatial(adata, color=[annotation_key], title=f'{sample_id} ({subtype})', size=0.5,
            crop_coord=crop_coord, ax=ax, palette=['lightgray'], na_in_legend=False, legend_loc=None, alpha_img=0.3)

ax.set_title(f'{sample_id} ({subtype})\nTLS regions', fontweight='bold')

# Remove the labels
ax.set_xlabel('')
ax.set_ylabel('')

# Plot TLS regions
if plot_tls_polygons:
    tls_gdf_plot.plot(ax=ax, facecolor="mediumblue", edgecolor="black", linewidth=0, alpha=0.4)

# -- Plot which cells are in TLS regions
ax = axs[n+2]
sc.pl.spatial(adata, color=['in_TLS'], title=f'{sample_id} ({subtype})', size=0.5,
            crop_coord=crop_coord, ax=ax, palette=['lightgray', 'cornflowerblue'], na_in_legend=False, alpha_img=0.3)

ax.set_title(f'{sample_id} ({subtype})\nCells in TLS regions', fontweight='bold')

# Remove the labels
ax.set_xlabel('')
ax.set_ylabel('')

# Add suptitle
plt.suptitle(f'{sample_id} ({subtype}) - TLS regions identified based on T & B cell densities', y=1.07, fontweight='bold')


# --- Plot 3: Barplot showing number of cells in TLS regions

# Calculate cell counts per cell type
cell_counts = adata.obs.groupby([annotation_key, 'in_TLS']).size().unstack(fill_value=0)

# Sort by total number of cells
cell_counts['total'] = cell_counts.sum(axis=1)
cell_counts = cell_counts.sort_values(by='total', ascending=False)
cell_counts = cell_counts.drop(columns=['total'])

# Reverse the column order
cell_counts = cell_counts[['True', 'False']]

# Plot as a stacked barplot
ncols = 2
nrows = 1

# Set the grid behind plots
plt.rcParams['axes.axisbelow'] = True

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(6*ncols, 5*nrows))

# -- First do a stacked barplot of within TLS and outside TLS

ax = axs[0]
cell_counts.plot(kind='bar', stacked=True, ax=ax, color=['cornflowerblue', 'lightgray'], width=0.8)
ax.set_title(f'{sample_id} ({subtype})\nCells in TLS regions', fontweight='bold')
ax.set_ylabel('Number of cells')
ax.set_xlabel('')
ax.legend(title='In TLS', labels=['True', 'False'], loc='upper right')

# Annotate the percentage of True
for n, cell_type in enumerate(cell_counts.index):
    total_true = cell_counts.loc[cell_type, 'True']
    total = total_true + cell_counts.loc[cell_type, 'False']
    pct = total_true / total * 100  # in_tls is total_true
    ax.text(n, total_true + total*0.01, f'{pct:.1f}%', ha='center', va='bottom', fontsize=10)

# Make y axis scientific notation if needed
ax.set_ylabel('Number of cells (x10³)')
ax.set_yticklabels(ax.get_yticks() * 1e-3)

# -- Then do a barplot showing number of cells in TLS regions per cell type

# Sort by number of cells in TLS regions
cell_counts_tls = cell_counts.sort_values(by='True', ascending=False)

ax = axs[1]
# Extract colors from color_dict for the cell types in the order of cell_counts_tls
colors = [color_dict[cell_type] for cell_type in cell_counts_tls.index]
cell_counts_tls['True'].plot(kind='bar', ax=ax, width=0.8, color=colors)
ax.set_title(f'{sample_id} ({subtype})\nCells in TLS regions', fontweight='bold')
ax.set_ylabel('Number of cells')
ax.set_xlabel('')  

# Make y axis scientific notation if needed
# Make it x10^3 for interpretability if numbers are large
ax.set_ylabel('Number of cells (x10³)')
ax.set_yticklabels(ax.get_yticks() * 1e-3)

# Move the second plot a bit to the right to avoid overlap
ax_pos = ax.get_position()
ax.set_position([ax_pos.x0 + 0.05, ax_pos.y0, ax_pos.width, ax_pos.height])

# Add suptitle
plt.suptitle(f'{sample_id} ({subtype}) - Cells in TLS regions identified based on T & B cell densities', y=1.2, fontweight='bold')


# --- Save the output

out_dir = os.path.join(res_dir, '05_xenium_mechanisms', '01_identify_TLS_regions')
os.makedirs(out_dir, exist_ok=True)

# Make a dir for plots
plots_out_dir = os.path.join(out_dir, 'plots')
os.makedirs(plots_out_dir, exist_ok=True)

# Make a separate dir for csv files
csv_out_dir = os.path.join(out_dir, 'csv_files')
os.makedirs(csv_out_dir, exist_ok=True)

# First save the geojson of the TLS regions
tls_gdf.to_file(os.path.join(csv_out_dir, f'{sample_id}_TLS_regions.geojson'), driver="GeoJSON")

# Then save the plots
filename = f'{sample_id}_TLS_regions_spatial.pdf'

# Save with high resolution
filename_dir=os.path.join(plots_out_dir, filename)

pp = PdfPages(filename_dir)
fig_nums = plt.get_fignums()

figs = [plt.figure(n) for n in fig_nums]

for fig in figs:
    fig.savefig(pp, format='pdf', bbox_inches='tight', dpi=300)

pp.close()

# Save the cell-level TLS annotations
to_save = adata.obs[[annotation_key, 'in_TLS', 'TLS_id']]
to_save.to_csv(os.path.join(csv_out_dir, f'{sample_id}_cell_TLS_annotations.csv'))


# --- Finish

print("Script successfully completed")