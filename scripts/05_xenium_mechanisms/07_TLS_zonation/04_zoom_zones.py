#!/usr/bin/env python3

# --- Aim
# Representative single-transcript zoom of adjacent B cell zone and T cell zone fields in one
# sample (Extended Data Fig. 9f,g). Produces one page:
#   - Overview (Extended Data Fig. 9f): aligned H&E cropped around the two fields, with the 50 um zone bins
#     (B cell zone / T cell zone / mixed) overlaid and the two zoom fields boxed.
#   - Grid (Extended Data Fig. 9g): rows = B cell zone (top) and T cell zone (bottom); columns = H&E, B/T cell
#     markers (MS4A1, CD3E), CXCL13, CCL19, CCL21, FDC markers (CR1, CR2) and HEV markers
#     (TSPAN7, MEOX2), each an H&E crop with the transcripts of that group overlaid.
#
# The two fields come from an areas-of-interest geojson whose polygons are stored in Xenium
# micron coordinates; they are mapped onto the aligned H&E with the same spatialdata image_space
# affine used for the transcripts (identical convention to 03_zoom_GCs.py), so boxes, bins,
# transcripts and image all share one pixel frame.
#
# Inputs : areas-of-interest geojson ({sample_id}_zone_zoom.geojson, micron polygons named
#          'B cell zone' / 'T cell zone'); BT zone bin CSV from 01_classify_zones.py; spatialdata
#          zarr; aligned full-res H&E (.ome.tif); Xenium h5ad (scalefactors); transcript colours.
# Output : {res_dir}/05_xenium_mechanisms/07_TLS_zonation/04_zoom_zones/plots/{sample_id}_zone_zoom.pdf


# --- Define environment
# env: spatdat

# --- Load packages

import argparse
import json
import os
import shutil
import sys

import anndata as ad
import pandas as pd
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib_scalebar.scalebar import ScaleBar

import spatialdata as sd
from spatialdata.transformations import get_transformation

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments
# Single representative sample (Extended Data Fig. 9f,g uses P22_Xe5k, the sample whose B/T fields were
# annotated). Dispatches through submit_array_jobs_nci.sh -> run_single_sample_job.sh.

parser = argparse.ArgumentParser(description="Zone zoom (B vs T cell zone fields), Extended Data Fig. 9f,g",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad input files")
parser.add_argument("--sample_id", default='P22_Xe5k',
                    help="Sample (section_id); the representative sample with an annotated zone-zoom geojson")
parser.add_argument("--project_dir", help="Path to the root of the project directory")
parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)
parser.add_argument("--grid_um", type=float, default=50.0, help="Bin size used in step 01 (microns)")
parser.add_argument("--overview_margin_um", type=float, default=75.0,
                    help="Minimum margin (microns) between a zoom field edge and the square overview crop edge")
parser.add_argument("--scalebar_overview_um", type=float, default=100.0)
parser.add_argument("--scalebar_zoom_um", type=float, default=20.0)

args = vars(parser.parse_args())
sample_id          = args["sample_id"]
h5ad_dir           = args["h5ad_dir"]
project_dir        = args["project_dir"]
res_dir            = args["res_dir"]
grid_um            = args["grid_um"]
overview_margin_um = args["overview_margin_um"]
scalebar_val_overview = args["scalebar_overview_um"]
scalebar_val_zoom     = args["scalebar_zoom_um"]

sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import read_tiff_region

if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_xenium.csv')


# --- Figure spec (Extended Data Fig. 9g column layout: H&E, B/T cell, CXCL13, CCL19, CCL21, FDC, HEV)

area_names_zoom = ['B cell zone', 'T cell zone']   # zoom rows, top -> bottom

zoom_col_specs = [
    ('H&E',        None),
    ('B / T cell', ['MS4A1', 'CD3E']),
    ('CXCL13',     ['CXCL13']),
    ('CCL19',      ['CCL19']),
    ('CCL21',      ['CCL21']),
    ('FDC',        ['CR1', 'CR2']),
    ('HEV',        ['TSPAN7', 'MEOX2']),
]
transcripts_to_plot = sorted({g for _, gs in zoom_col_specs if gs for g in gs})

scalebar_units         = 'µm'
scalebar_height_main   = 0.02
scalebar_border_pad    = 0.3
scalebar_pad           = 0.08
scalebar_box_color     = 'white'
scalebar_box_alpha     = 0.85
transcript_image_alpha = 0.2    # faint H&E behind the transcript scatter
transcript_marker_size = 30
ax_size      = 0.55
ax_spacing   = 0.05
overview_gap = 0.15
zoom_rect_lw = 5

zone_color_dict = {
    'B_cell_zone': mcolors.cnames['cornflowerblue'],
    'T_cell_zone': mcolors.cnames['crimson'],
    'mixed':       '#921e95',
}
zone_label_dict = {'B_cell_zone': 'B cell zone', 'T_cell_zone': 'T cell zone', 'mixed': 'Mixed'}
zone_alpha = 0.7


# --- Sample metadata

sample_sheet    = pd.read_csv(sample_sheet_path)
tissue_metadata = pd.read_csv(os.path.join(project_dir, 'resources', 'tissue_metadata.csv'))
tissue_id = sample_sheet.loc[sample_sheet['section_id'] == sample_id, 'tissue_id'].values[0]
subtype   = tissue_metadata.loc[tissue_metadata['tissue_id'].astype(str) == str(tissue_id),
                                'clinical_subtype'].values[0]


# --- Transcript colours

tx_color_csv  = os.path.join(project_dir, 'resources', 'plotting_colors', 'transcript_colors.csv')
tx_color_dict = pd.read_csv(tx_color_csv, index_col=0)['color'].to_dict()


# --- Load adata (for the aligned image micron/pixel scale)

adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample_id}.h5ad'))
microns_per_pixel = adata.uns['spatial'][sample_id]['scalefactors']['microns_per_pixel']


# --- Load the aligned full-res H&E (image_space)

image_path = os.path.join(project_dir, 'data', 'images', 'he', 'aligned_images', 'xenium',
                          f'{sample_id}.ome.tif')
# The aligned H&E is a large whole-slide image; only the zoom and overview crops are read,
# tile by tile at full resolution (read_tiff_region in crop_image below).


# --- Load spatialdata + filter transcripts to the plotted genes (image_space)

spat_dat_dir = os.path.join(project_dir, 'data', 'spatialdata_objects', 'xenium')
zarr_path = os.path.join(spat_dat_dir, f'{sample_id}.zarr')
if not os.path.exists(zarr_path):
    zip_path = os.path.join(spat_dat_dir, f'{sample_id}.zip')
    print(f'Unzipping {zip_path} -> {zarr_path}')
    shutil.unpack_archive(zip_path, zarr_path)

sdata = sd.read_zarr(zarr_path)
sdata.points['transcripts'] = sdata.points['transcripts'][
    sdata.points['transcripts']['feature_name'].isin(transcripts_to_plot)
]
# Transform only the transcripts: sd.transform(sdata, ...) would also resample the
# full-resolution label rasters into image_space
tx_df = sd.transform(sdata.points['transcripts'], to_coordinate_system='image_space').compute()
print(f'Loaded {len(tx_df)} transcripts across {len(transcripts_to_plot)} genes', flush=True)


# --- Xenium micron -> aligned-image (image_space) pixel transform (as in 03_zoom_GCs.py)

affine = get_transformation(sdata.points['transcripts'],
                            to_coordinate_system='image_space').to_affine_matrix(
    input_axes=('x', 'y'), output_axes=('x', 'y'))


def micron_to_he_pixel(x_um, y_um):
    coords = np.column_stack([np.asarray(x_um, dtype=float), np.asarray(y_um, dtype=float)])
    hom = np.hstack([coords, np.ones((coords.shape[0], 1))])
    out = (affine @ hom.T).T
    return out[:, 0], out[:, 1]


# --- Areas of interest: micron-space zoom boxes -> aligned-image pixel bboxes

def load_zoom_and_overview(geojson_path):
    """Read the AOI geojson (polygons in Xenium microns), map the two field boxes to aligned-image
    pixels, and derive a square overview crop enclosing both with an even margin."""
    gj = json.load(open(geojson_path))
    zoom_coords = {}
    for feat in gj['features']:
        cls = feat.get('properties', {}).get('classification', {})
        name = cls.get('name') if cls else None
        if name not in area_names_zoom:
            continue
        ring = feat['geometry']['coordinates'][0]
        xs_um = [pt[0] for pt in ring]
        ys_um = [pt[1] for pt in ring]
        px, py = micron_to_he_pixel(xs_um, ys_um)
        zoom_coords[name] = (int(px.min()), int(px.max()), int(py.min()), int(py.max()))

    margin_px = overview_margin_um / microns_per_pixel
    x0 = min(c[0] for c in zoom_coords.values()) - margin_px
    x1 = max(c[1] for c in zoom_coords.values()) + margin_px
    y0 = min(c[2] for c in zoom_coords.values()) - margin_px
    y1 = max(c[3] for c in zoom_coords.values()) + margin_px
    side  = max(x1 - x0, y1 - y0)
    ov_cx = (x0 + x1) / 2
    ov_cy = (y0 + y1) / 2
    overview_coords = (int(ov_cx - side / 2), int(ov_cx + side / 2),
                       int(ov_cy - side / 2), int(ov_cy + side / 2))
    return zoom_coords, overview_coords


geojson_path = os.path.join(project_dir, 'resources', 'areas_of_interest', 'xenium',
                            f'{sample_id}_zone_zoom.geojson')
zoom_coords, overview_coords = load_zoom_and_overview(geojson_path)
print(f'Zoom fields (image px): {zoom_coords}')


# --- Zone bins (bin_x/bin_y in Xenium micron coords)

zone_dir  = os.path.join(res_dir, '05_xenium_mechanisms', '07_TLS_zonation',
                         '01_classify_zones', 'csv_files')
zone_bins = pd.read_csv(os.path.join(zone_dir, f'{sample_id}_BT_zone_bins.csv'))
zone_bins = zone_bins[zone_bins['TLS_zone'].isin(zone_color_dict)].copy()


# --- Helpers

def crop_image(coords):
    return read_tiff_region(image_path, *coords)


def move_plots_x_y(axs, x_offset=0.02, y_offset=0.02, remove_title=True):
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


# --- Build the figure

plots_dir = os.path.join(res_dir, '05_xenium_mechanisms', '07_TLS_zonation', '04_zoom_zones', 'plots')
os.makedirs(plots_dir, exist_ok=True)
pdf_path = os.path.join(plots_dir, f'{sample_id}_zone_zoom.pdf')

nrows = len(area_names_zoom)
ncols = len(zoom_col_specs)
fig, axs = plt.subplots(nrows, ncols, figsize=(5, 5))

he_zoom_axes = []
for i, area in enumerate(area_names_zoom):
    coords = zoom_coords[area]
    crop_img, adjusted_coords = crop_image(coords)
    width_crop = coords[1] - coords[0]

    for j, (_col_title, gene_list) in enumerate(zoom_col_specs):
        ax = axs[i, j]
        ax.imshow(crop_img, alpha=1.0 if gene_list is None else transcript_image_alpha)

        if gene_list is not None:
            mask = (
                (tx_df["x"] >= coords[0]) & (tx_df["x"] < coords[1]) &
                (tx_df["y"] >= coords[2]) & (tx_df["y"] < coords[3]) &
                (tx_df["feature_name"].isin(gene_list))
            )
            tx_sub = tx_df[mask].copy()
            tx_sub = tx_sub.sample(frac=1, random_state=42).reset_index(drop=True)
            tx_sub['x_adj'] = tx_sub['x'] - adjusted_coords[0]
            tx_sub['y_adj'] = tx_sub['y'] - adjusted_coords[2]
            if not tx_sub.empty:
                colors     = tx_sub['feature_name'].map(tx_color_dict)
                edgecolors = colors.apply(lambda c: tuple(max(comp - 0.3, 0) for comp in mcolors.to_rgb(c)))
                ax.scatter(tx_sub['x_adj'], tx_sub['y_adj'],
                           s=transcript_marker_size, c=colors, alpha=1,
                           linewidth=0.5, edgecolors=edgecolors)
        else:
            he_zoom_axes.append(ax)

        ax.set_xlim(0, width_crop)
        ax.set_ylim(width_crop, 0)
        ax.grid(False)

for ax in axs.flat:
    ax.set_xlabel(''); ax.set_ylabel('')
    ax.set_xticks([]); ax.set_yticks([])

# Equal-size square axes, then repack tightly
for ax in axs.flat:
    ax_pos = ax.get_position()
    ax.set_position([ax_pos.x0, ax_pos.y0, ax_size, ax_size])
move_plots_x_y(axs, x_offset=ax_spacing, y_offset=ax_spacing, remove_title=True)

# Row labels + column titles
for ax, area in zip(axs[:, 0], area_names_zoom):
    ax.set_ylabel(area, fontweight='bold', fontsize=25, labelpad=10)
for ax, (title, _) in zip(axs[0, :], zoom_col_specs):
    ax.set_title(title, fontweight='bold', fontsize=25)


# --- Overview ax on the left, spanning both zoom rows (Extended Data Fig. 9f)

ax_top_y    = axs[0, 0].get_position().y1
ax_bottom_y = axs[-1, 0].get_position().y0
ax_height   = ax_top_y - ax_bottom_y
ax_x        = axs[0, 0].get_position().x0 - overview_gap - ax_height
new_ax = fig.add_axes([ax_x, ax_bottom_y, ax_height, ax_height])

ov_crop_img, ov_adjusted = crop_image(overview_coords)
new_ax.imshow(ov_crop_img)
width_ov = overview_coords[1] - overview_coords[0]
new_ax.set_xlim(0, width_ov); new_ax.set_ylim(width_ov, 0)
new_ax.grid(False); new_ax.set_xticks([]); new_ax.set_yticks([])

# 50 um zone bins overlaid. Transform bin centres + two opposite corners per bin (microns ->
# pixels) so a rotated alignment still gives the right per-bin pixel size.
bin_cx_px, bin_cy_px = micron_to_he_pixel(zone_bins['bin_x'].values, zone_bins['bin_y'].values)
bin_x0_px, bin_y0_px = micron_to_he_pixel(zone_bins['bin_x'].values - grid_um / 2,
                                          zone_bins['bin_y'].values - grid_um / 2)
bin_x1_px, bin_y1_px = micron_to_he_pixel(zone_bins['bin_x'].values + grid_um / 2,
                                          zone_bins['bin_y'].values + grid_um / 2)
bin_size_px = np.hypot(bin_x1_px - bin_x0_px, bin_y1_px - bin_y0_px) / np.sqrt(2)

cx_local = bin_cx_px - ov_adjusted[0]
cy_local = bin_cy_px - ov_adjusted[2]
in_crop = ((cx_local + bin_size_px / 2 >= 0) & (cx_local - bin_size_px / 2 <= width_ov) &
           (cy_local + bin_size_px / 2 >= 0) & (cy_local - bin_size_px / 2 <= width_ov))
bins_in = zone_bins[in_crop]
if len(bins_in):
    rects  = [mpatches.Rectangle((cx - s / 2, cy - s / 2), s, s)
              for cx, cy, s in zip(cx_local[in_crop], cy_local[in_crop], bin_size_px[in_crop])]
    colors = [zone_color_dict[z] for z in bins_in['TLS_zone']]
    new_ax.add_collection(PatchCollection(rects, facecolors=colors, edgecolors='none',
                                          alpha=zone_alpha, zorder=2))

# Thin white 50 um grid over the crop
grid_x_min = zone_bins['bin_x'].min() - grid_um / 2
grid_x_max = zone_bins['bin_x'].max() + grid_um / 2
grid_y_min = zone_bins['bin_y'].min() - grid_um / 2
grid_y_max = zone_bins['bin_y'].max() + grid_um / 2
x_grid_um = np.arange(grid_x_min, grid_x_max + grid_um / 2, grid_um)
y_grid_um = np.arange(grid_y_min, grid_y_max + grid_um / 2, grid_um)
vert_starts = np.column_stack([x_grid_um, np.full_like(x_grid_um, grid_y_min)])
vert_ends   = np.column_stack([x_grid_um, np.full_like(x_grid_um, grid_y_max)])
horz_starts = np.column_stack([np.full_like(y_grid_um, grid_x_min), y_grid_um])
horz_ends   = np.column_stack([np.full_like(y_grid_um, grid_x_max), y_grid_um])
s_he_x, s_he_y = micron_to_he_pixel(np.vstack([vert_starts, horz_starts])[:, 0],
                                    np.vstack([vert_starts, horz_starts])[:, 1])
e_he_x, e_he_y = micron_to_he_pixel(np.vstack([vert_ends, horz_ends])[:, 0],
                                    np.vstack([vert_ends, horz_ends])[:, 1])
segments = np.stack([
    np.column_stack([s_he_x - ov_adjusted[0], s_he_y - ov_adjusted[2]]),
    np.column_stack([e_he_x - ov_adjusted[0], e_he_y - ov_adjusted[2]]),
], axis=1)
new_ax.add_collection(LineCollection(segments, colors='white', linewidths=1.8, alpha=0.75, zorder=2.5))

# Overview scalebar
new_ax.add_artist(ScaleBar(microns_per_pixel, units=scalebar_units, location='lower right',
                           fixed_value=scalebar_val_overview, frameon=True,
                           box_color=scalebar_box_color, box_alpha=scalebar_box_alpha,
                           scale_loc='none', pad=scalebar_pad, border_pad=scalebar_border_pad,
                           height_fraction=scalebar_height_main, font_properties={"size": 18}))

# Zoom field boxes on the overview
for area in area_names_zoom:
    x1, x2, y1, y2 = zoom_coords[area]
    new_ax.add_patch(mpatches.Rectangle((x1 - ov_adjusted[0], y1 - ov_adjusted[2]),
                                        x2 - x1, y2 - y1, linewidth=zoom_rect_lw,
                                        edgecolor='black', facecolor='none', zorder=3))

new_ax.set_title('Zone annotation (50 µm bins)', fontweight='bold', fontsize=25)

# H&E zoom scalebars (sized to match the overview visually; transcript panels stay unobstructed)
scalebar_height_zoom = scalebar_height_main * (ax_height / ax_size)
for ax in he_zoom_axes:
    ax.add_artist(ScaleBar(microns_per_pixel, units=scalebar_units, location='lower right',
                           fixed_value=scalebar_val_zoom, frameon=True,
                           box_color=scalebar_box_color, box_alpha=scalebar_box_alpha,
                           scale_loc='none', pad=scalebar_pad, border_pad=scalebar_border_pad,
                           height_fraction=scalebar_height_zoom, font_properties={"size": 18}))


# --- Legends

for j, (_col_title, gene_list) in enumerate(zoom_col_specs):
    if gene_list is None:
        continue
    handles = []
    for gene in gene_list:
        color  = tx_color_dict.get(gene, 'black')
        darker = tuple(max(c - 0.3, 0) for c in mcolors.to_rgb(color))
        handles.append(Line2D([0], [0], marker='o', color='w', label=gene,
                              markerfacecolor=color, markersize=15,
                              linewidth=2, markeredgecolor=darker))
    axs[-1, j].legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0),
                      fontsize=18, frameon=False, handlelength=0.7,
                      ncol=len(gene_list), columnspacing=0.8)

zone_handles = [Line2D([0], [0], marker='s', color='w', label=zone_label_dict[z],
                       markerfacecolor=zone_color_dict[z], markersize=15,
                       alpha=zone_alpha, linewidth=0)
                for z in ['B_cell_zone', 'T_cell_zone', 'mixed']]
new_ax.legend(handles=zone_handles, loc='upper center', bbox_to_anchor=(0.5, 0),
              fontsize=18, frameon=False, handlelength=0.9, ncol=3, columnspacing=0.8)


# --- "Transcripts" group bracket across the transcript columns

top_y1 = max(new_ax.get_position().y1, *(axs[0, c].get_position().y1 for c in range(ncols)))
trans_cols    = [j for j, (_, gs) in enumerate(zoom_col_specs) if gs is not None]
trans_x_start = axs[0, min(trans_cols)].get_position().x0
trans_x_end   = axs[0, max(trans_cols)].get_position().x1
bracket_y     = top_y1 + 0.09
fig.add_artist(Line2D([trans_x_start, trans_x_end], [bracket_y, bracket_y],
                      transform=fig.transFigure, color='black', linewidth=1.5,
                      solid_capstyle='butt'))
fig.text((trans_x_start + trans_x_end) / 2, bracket_y + 0.005, 'Transcripts',
         ha='center', va='bottom', fontsize=28, fontweight='bold')


# --- Suptitle

title_y = top_y1 + 0.30
title_x = (new_ax.get_position().x0 + axs[0, -1].get_position().x1) / 2
fig.suptitle(f'Extended Data Fig. 9f,g — {sample_id} ({subtype}): B cell zone vs T cell zone '
             f'(scalebar: {scalebar_val_overview:g} {scalebar_units} overview, '
             f'{scalebar_val_zoom:g} {scalebar_units} insets)',
             x=title_x, y=title_y, ha='center', fontweight='bold', fontsize=30)


# --- Thicker spines

for ax in list(axs.flat) + [new_ax]:
    for spine in ax.spines.values():
        spine.set_linewidth(2.0)


# --- Save

plt.savefig(pdf_path, bbox_inches='tight', dpi=500)
plt.close(fig)
print(f'\nSaved: {pdf_path}')
print('Script successfully completed')
