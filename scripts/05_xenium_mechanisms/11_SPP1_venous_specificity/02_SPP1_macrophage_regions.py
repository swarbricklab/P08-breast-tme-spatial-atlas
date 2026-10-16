#!/usr/bin/env python3

# --- Aim
# Extended Data Fig. 10e: macrophage SPP1 outside TLS. For the representative sample (P29_Xe5k),
# show four tissue regions where LAM1 (Myeloid_c1_LAM1_FABP5) macrophages express SPP1,
# alongside a whole-tissue overview with the four zoom boxes marked. Each region is one
# column of a 3-row grid:
#   Row 1: aligned H&E crop.
#   Row 2: faint H&E + SPP1 transcripts + LAM1-macrophage cell boundaries.
#   Row 3: SPP1 transcripts + LAM1-macrophage cell boundaries, no H&E.
#
# The four regions are read from an areas-of-interest geojson whose polygons are stored in
# Xenium micron coordinates and mapped onto the aligned H&E with the SpatialData image_space
# affine (same convention as 07_TLS_zonation/04_zoom_zones.py). Regions are ordered top to
# bottom by their y-coordinate and numbered 1-4.

# --- Define environment
# env: spatdat

# --- Load packages

import argparse
import json
import os
import shutil

import anndata as ad
import pandas as pd
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib_scalebar.scalebar import ScaleBar

import PIL
import tifffile
import spatialdata as sd
from spatialdata.transformations import get_transformation

mpl.rcParams['pdf.fonttype'] = 42

SPP1_GENE      = 'SPP1'
MACRO_LABEL    = 'Myeloid_c1_LAM1_FABP5'
MACRO_COLOR    = 'indigo'


# --- Argparse arguments (representative single sample; Extended Data Fig. 10e uses P29_Xe5k)

parser = argparse.ArgumentParser(description="Macrophage SPP1 outside TLS region zoom (Extended Data Fig. 10e)",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad input files")
parser.add_argument("--sample_id", default='P29_Xe5k',
                    help="Sample (section_id); the representative sample with an annotated SPP1-macrophage geojson")
parser.add_argument("--project_dir", help="Path to the root of the project directory")
parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)
parser.add_argument("--annotation_level", default='subset')
parser.add_argument("--scalebar_um", type=float, default=50.0)

args = vars(parser.parse_args())
sample_id        = args["sample_id"]
h5ad_dir         = args["h5ad_dir"]
project_dir      = args["project_dir"]
res_dir          = args["res_dir"]
annotation_level = args["annotation_level"]
scalebar_um      = args["scalebar_um"]

if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_xenium.csv')


# --- Sample metadata

sample_sheet    = pd.read_csv(sample_sheet_path)
tissue_metadata = pd.read_csv(os.path.join(project_dir, 'resources', 'tissue_metadata.csv'))
tissue_id = sample_sheet.loc[sample_sheet['section_id'] == sample_id, 'tissue_id'].values[0]
subtype   = tissue_metadata.loc[tissue_metadata['tissue_id'].astype(str) == str(tissue_id),
                                'clinical_subtype'].values[0]


# --- Colors (SPP1 transcript color; LAM1 macrophage cells indigo)

tx_color_csv  = os.path.join(project_dir, 'resources', 'plotting_colors', 'transcript_colors.csv')
tx_color_dict = pd.read_csv(tx_color_csv, index_col=0)['color'].to_dict()

classif = pd.read_csv(os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas',
                                   'miniatlas_subset_classifications.csv'))
subset_to_short = dict(zip(classif['miniatlas_subset'], classif['miniatlas_subset_short']))
macro_short = subset_to_short.get(MACRO_LABEL, MACRO_LABEL)


# --- Load adata (for the aligned image micron/pixel scale + cell annotations)

adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample_id}.h5ad'))
microns_per_pixel = adata.uns['spatial'][sample_id]['scalefactors']['microns_per_pixel']


def celltype_csv(sample, level='subset'):
    p = os.path.join(res_dir, '02_cell_type_mapping', 'label_transfer', 'xenium',
                     f'{level}_level', f'{sample}_{level}_pred_cell_types.csv')
    return p


ann = pd.read_csv(celltype_csv(sample_id, annotation_level), index_col=0)['pred_cell_type'].astype(str)
adata.obs['cell_type'] = ann.reindex(adata.obs_names)


# --- Load the aligned full-res H&E (image_space)

image_path = os.path.join(project_dir, 'data', 'images', 'he', 'aligned_images', 'xenium',
                          f'{sample_id}.ome.tif')
PIL.Image.MAX_IMAGE_PIXELS = 24078974976
# Large whole-slide image; only small crops + a downsampled overview are shown, but the whole
# image is read (lazy tiled reads need zarr>=3, and spatdat ships zarr 2.x). High-memory node.
img = tifffile.imread(image_path)
print(f'Loaded aligned H&E: {img.shape}')


# --- Load spatialdata: SPP1 transcripts + cell boundaries in image_space

spat_dat_dir = os.path.join(project_dir, 'data', 'spatialdata_objects', 'xenium')
zarr_path = os.path.join(spat_dat_dir, f'{sample_id}.zarr')
if not os.path.exists(zarr_path):
    zip_path = os.path.join(spat_dat_dir, f'{sample_id}.zip')
    print(f'Unzipping {zip_path} -> {zarr_path}')
    shutil.unpack_archive(zip_path, zarr_path)

sdata = sd.read_zarr(zarr_path)
sdata.points['transcripts'] = sdata.points['transcripts'][
    sdata.points['transcripts']['feature_name'] == SPP1_GENE
]
transformed_sdata = sd.transform(sdata, to_coordinate_system='image_space')
tx_df = transformed_sdata.points['transcripts'].compute()
shapes = transformed_sdata.shapes['cell_boundaries']
print(f'Loaded {len(tx_df)} SPP1 transcripts; {len(shapes)} cell boundaries')


# --- Xenium micron -> aligned-image (image_space) pixel transform

affine = get_transformation(sdata.points['transcripts'],
                            to_coordinate_system='image_space').to_affine_matrix(
    input_axes=('x', 'y'), output_axes=('x', 'y'))


def micron_to_he_pixel(x_um, y_um):
    coords = np.column_stack([np.asarray(x_um, dtype=float), np.asarray(y_um, dtype=float)])
    hom = np.hstack([coords, np.ones((coords.shape[0], 1))])
    out = (affine @ hom.T).T
    return out[:, 0], out[:, 1]


# Cell centroids (micron) -> aligned-image pixels, for region membership tests.
adata.obs['px_x'], adata.obs['px_y'] = micron_to_he_pixel(
    adata.obs['x_centroid'].values, adata.obs['y_centroid'].values)


# --- Areas of interest: micron polygons -> aligned-image pixel bboxes, ordered top->bottom

geojson_path = os.path.join(project_dir, 'resources', 'areas_of_interest', 'xenium',
                            f'{sample_id}_SPP1_macrophages.geojson')
gj = json.load(open(geojson_path))
region_px = []
for feat in gj['features']:
    ring = feat['geometry']['coordinates'][0]
    px, py = micron_to_he_pixel([pt[0] for pt in ring], [pt[1] for pt in ring])
    region_px.append((int(px.min()), int(px.max()), int(py.min()), int(py.max())))
region_px.sort(key=lambda c: c[2])   # by top (y_min)
area_coords = {f'region_{i + 1}': c for i, c in enumerate(region_px)}
print(f'Regions (image px, top->bottom): {area_coords}')


# --- Helpers

def crop_image(img, coords):
    x_min = max(coords[0], 0)
    x_max = min(coords[1], img.shape[1])
    y_min = max(coords[2], 0)
    y_max = min(coords[3], img.shape[0])
    return img[y_min:y_max, x_min:x_max], (x_min, x_max, y_min, y_max)


def make_transparent(color, alpha=0.5):
    rgba = mcolors.to_rgba(color)
    return (rgba[0], rgba[1], rgba[2], alpha)


def move_plots_x_y(axs, x_offset=0.05, y_offset=0.05, remove_title=True):
    for i, axs_row in enumerate(axs):
        for j, ax in enumerate(axs_row):
            if i != 0 and remove_title:
                ax.set_title('')
            ax_pos = ax.get_position()
            new_pos_y = ax_pos.y0 if i == 0 else ax_pos_row_prev.y0 - ax_pos.height - y_offset
            new_pos_x = ax_pos.x0 if j == 0 else ax_pos_prev.x0 + ax_pos_prev.width + x_offset
            ax.set_position([new_pos_x, new_pos_y, ax_pos.width, ax_pos.height])
            ax_pos_prev = ax.get_position()
        ax_pos_row_prev = axs_row[0].get_position()


def make_scalebar():
    return ScaleBar(microns_per_pixel, units='µm', location='lower right', box_alpha=0.5,
                    fixed_value=scalebar_um, frameon=False, scale_loc='none', pad=0.3,
                    height_fraction=0.02, font_properties={"size": 18})


macro_ids = set(adata.obs_names[adata.obs['cell_type'] == MACRO_LABEL])


def region_macro_shapes(coords):
    """Cell boundaries of LAM1 macrophages within the region (image px), translated to
    the region-local frame."""
    buf = int(0.1 * (coords[1] - coords[0]))
    obs = adata.obs
    in_box = ((obs['px_x'] >= coords[0] - buf) & (obs['px_x'] < coords[1] + buf) &
              (obs['px_y'] >= coords[2] - buf) & (obs['px_y'] < coords[3] + buf))
    ids = set(obs.index[in_box]) & macro_ids
    sub = shapes[shapes.index.isin(ids)].copy()
    sub['geometry'] = sub.geometry.translate(xoff=-coords[0], yoff=-coords[2])
    return sub


def region_spp1(coords):
    m = ((tx_df['x'] >= coords[0]) & (tx_df['x'] < coords[1]) &
         (tx_df['y'] >= coords[2]) & (tx_df['y'] < coords[3]))
    sub = tx_df[m].sample(frac=1, random_state=42).reset_index(drop=True)
    sub = sub.assign(x_adj=sub['x'] - coords[0], y_adj=sub['y'] - coords[2])
    return sub


# --- Build the figure (3 rows x N region columns)

plt.rcParams.update({'font.size': 25})
ncols = len(area_coords)
fig, axs = plt.subplots(3, ncols, figsize=(5, 5))
spp1_color = tx_color_dict.get(SPP1_GENE, '#1ee23e')
spp1_edge  = tuple(max(c - 0.3, 0) for c in mcolors.to_rgb(spp1_color))

for n, (area, coords) in enumerate(area_coords.items()):
    crop_img, adj = crop_image(img, coords)
    width_crop = coords[1] - coords[0]
    s_dot = 10000 * (1.0 / width_crop)

    # Row 1: H&E only
    ax = axs[0, n]
    ax.imshow(crop_img)
    ax.add_artist(make_scalebar())
    ax.set_title(f'{n + 1}', fontsize=25, fontweight='bold')

    # Row 2: faint H&E + SPP1 + LAM1 macrophage boundaries
    ax = axs[1, n]
    ax.imshow(crop_img, alpha=0.5)
    tx_sub = region_spp1(coords)
    if not tx_sub.empty:
        ax.scatter(tx_sub['x_adj'], tx_sub['y_adj'], s=s_dot, c=spp1_color,
                   alpha=1, linewidth=1, edgecolors=[spp1_edge] * len(tx_sub))
    sub = region_macro_shapes(coords)
    if len(sub):
        sub.plot(ax=ax, facecolor=make_transparent(MACRO_COLOR, 0.3),
                 edgecolor=MACRO_COLOR, linewidth=0.5)
    ax.set_xlim(0, width_crop); ax.set_ylim(width_crop, 0)
    ax.grid(False)
    ax.add_artist(make_scalebar())

    # Row 3: SPP1 + LAM1 macrophage boundaries, no H&E
    ax = axs[2, n]
    if not tx_sub.empty:
        ax.scatter(tx_sub['x_adj'], tx_sub['y_adj'], s=s_dot, c=spp1_color,
                   alpha=1, linewidth=0, edgecolors=[spp1_edge] * len(tx_sub))
    if len(sub):
        sub.plot(ax=ax, facecolor=make_transparent(MACRO_COLOR, 0.0),
                 edgecolor=MACRO_COLOR, linewidth=0.5)
    ax.set_xlim(0, width_crop); ax.set_ylim(width_crop, 0)
    ax.grid(False)
    ax.add_artist(make_scalebar())

for ax in axs.flat:
    ax.set_xlabel(''); ax.set_ylabel('')
    ax.set_xticks([]); ax.set_yticks([])
for ax in axs.flat:
    ax_pos = ax.get_position()
    ax.set_position([ax_pos.x0, ax_pos.y0, 0.8, 0.8])
move_plots_x_y(axs, x_offset=0.05, y_offset=0.05, remove_title=True)

# Row labels
axs[0, 0].set_ylabel('H&E', fontweight='bold', fontsize=25)
axs[1, 0].set_ylabel('Segmentation mask\n+ transcripts', fontweight='bold', fontsize=20)
axs[2, 0].set_ylabel('Segmentation mask\n+ transcripts', fontweight='bold', fontsize=20)
for r in range(3):
    axs[r, 0].yaxis.set_visible(True)
    axs[r, 0].set_yticks([])


# --- Overview: whole tissue (downsampled aligned H&E) with numbered zoom boxes

cx = adata.obs['px_x'].values
cy = adata.obs['px_y'].values
ov_x0 = max(int(cx.min()) - 200, 0); ov_x1 = min(int(cx.max()) + 200, img.shape[1])
ov_y0 = max(int(cy.min()) - 200, 0); ov_y1 = min(int(cy.max()) + 200, img.shape[0])
ds = max(1, int(max(ov_x1 - ov_x0, ov_y1 - ov_y0) / 1800))   # downsample to ~1800 px
ov_img = img[ov_y0:ov_y1:ds, ov_x0:ov_x1:ds]

ax_top_y    = axs[0, 0].get_position().y1
ax_bottom_y = axs[-1, 0].get_position().y0
ax_height   = ax_top_y - ax_bottom_y
ax_x        = axs[0, 0].get_position().x0 - 0.05 - ax_height
overview_ax = fig.add_axes([ax_x, ax_bottom_y, ax_height, ax_height])
overview_ax.imshow(ov_img)
overview_ax.set_xticks([]); overview_ax.set_yticks([])
overview_ax.set_title('Tissue overview', fontweight='bold', fontsize=25)

for n, (area, coords) in enumerate(area_coords.items()):
    rx = (coords[0] - ov_x0) / ds
    ry = (coords[2] - ov_y0) / ds
    rw = (coords[1] - coords[0]) / ds
    rh = (coords[3] - coords[2]) / ds
    overview_ax.add_patch(mpatches.Rectangle((rx, ry), rw, rh, linewidth=3,
                                             edgecolor='black', facecolor='none'))
    overview_ax.text(rx + rw / 2, ry - 0.01 * overview_ax.get_ylim()[0], str(n + 1),
                     color='black', fontsize=20, fontweight='bold', ha='center', va='bottom',
                     backgroundcolor=(1, 1, 1, 0.7))


# --- Legends (rows 2 and 3, last column)

def macro_legend(alpha_fill):
    handles = [mpatches.Patch(facecolor=make_transparent(MACRO_COLOR, alpha_fill),
                              edgecolor=MACRO_COLOR, linewidth=2, label=macro_short),
               Line2D([0], [0], marker='o', color='w', label=SPP1_GENE,
                      markerfacecolor=spp1_color, markersize=15, linewidth=2,
                      markeredgecolor=spp1_edge)]
    return handles


axs[1, -1].legend(handles=macro_legend(0.3), loc='center left', bbox_to_anchor=(1, 0.5),
                  fontsize=20, frameon=False, handlelength=0.7)
axs[2, -1].legend(handles=macro_legend(0.0), loc='center left', bbox_to_anchor=(1, 0.5),
                  fontsize=20, frameon=False, handlelength=0.7)


# --- Suptitle + save

title_x = (axs[0, 0].get_position().x0 + axs[0, -1].get_position().x1) / 2
title_y = axs[0, 0].get_position().y1 + 0.2
plt.suptitle(f'Extended Data Fig. 10e — {sample_id} ({subtype}): macrophage SPP1 outside TLS '
             f'(scalebar = {scalebar_um:g} µm)',
             x=title_x, y=title_y, fontsize=30, fontweight='bold')

plots_dir = os.path.join(res_dir, '05_xenium_mechanisms', '11_SPP1_venous_specificity', 'plots')
os.makedirs(plots_dir, exist_ok=True)
pdf_path = os.path.join(plots_dir, f'ED10E_{sample_id}_SPP1_macrophages.pdf')
plt.savefig(pdf_path, bbox_inches='tight', dpi=500)
plt.close('all')
print(f'\nSaved: {pdf_path}')
print('Script successfully completed')
