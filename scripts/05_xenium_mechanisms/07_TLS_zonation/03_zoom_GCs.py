#!/usr/bin/env python3

# --- Aim
# Representative zoom plots of the discrete GC regions in one sample (Extended Data Fig. 10g).
# For each discrete GC (connected component of GC-classified bins from 01_classify_zones.py)
# a page is produced with:
#   1. Full-tissue aligned H&E overview with the GC bbox marked
#   2. Zoom: H&E cropped to the GC bbox + segmented GC outline
#   3. Zoom: H&E + 50 um grid + GC bins highlighted
#   4. Zoom: H&E only
#   5. Zoom: H&E + top-5 GC-signature transcripts overlaid
#
# Inputs : GC bin CSV from 01_classify_zones.py (with the discrete_GC column);
#          spatialdata zarr; aligned full-res H&E (.ome.tif); Xenium h5ad (scalefactors)
# Output : {res_dir}/05_xenium_mechanisms/07_TLS_zonation/03_zoom_GCs/plots/{sample_id}_GC_zoom.pdf


# --- Define environment
# env: spatdat

# --- Load packages

import argparse
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
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.lines import Line2D
from matplotlib_scalebar.scalebar import ScaleBar

import spatialdata as sd
from spatialdata.transformations import get_transformation

mpl.rcParams['pdf.fonttype'] = 42

TITLE_FS  = 20
SUP_FS    = 24
LEGEND_FS = 20

# Nice round scalebar lengths in um.
NICE_SCALEBAR_UM = [10, 20, 25, 50, 100, 200, 250, 500, 1000, 2000, 2500, 5000]


def pick_scalebar_um(side_um, lo_frac=0.25, hi_frac=0.35):
    lo, hi = side_um * lo_frac, side_um * hi_frac
    in_range = [v for v in NICE_SCALEBAR_UM if lo <= v <= hi]
    if in_range:
        return max(in_range)
    target = side_um * (lo_frac + hi_frac) / 2
    return min(NICE_SCALEBAR_UM, key=lambda v: abs(v - target))


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Zoom plots for discrete GC regions (Extended Data Fig. 10g)",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Runner contract (submit_array_jobs_nci.sh -> run_single_sample_job.sh)
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad input files")
parser.add_argument("--sample_id", help="Sample (section_id) to process")
parser.add_argument("--project_dir", help="Path to the root of the project directory")
parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)

# Analysis parameters (Methods-fixed defaults; exposed but never required)
parser.add_argument("--grid_um", type=float, default=50.0, help="Bin size used in step 01 (microns)")
parser.add_argument("--top_n_genes", type=int, default=5,
                    help="Number of top GC-signature transcripts to overlay (Extended Data Fig. 10g shows top 5)")
parser.add_argument("--genes_per_ax", type=int, default=5, help="Transcripts per panel")
parser.add_argument("--pad_um", type=float, default=50.0, help="Padding around each GC bbox in microns")

args = vars(parser.parse_args())

sample_id    = args["sample_id"]
h5ad_dir     = args["h5ad_dir"]
project_dir  = args["project_dir"]
res_dir      = args["res_dir"]
grid_um      = args["grid_um"]
top_n_genes  = args["top_n_genes"]
genes_per_ax = args["genes_per_ax"]
pad_um       = args["pad_um"]


# --- Local helper import

sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import tiff_shape, read_tiff_region, read_tiff_overview


# --- Sample metadata

if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_xenium.csv')

sample_sheet    = pd.read_csv(sample_sheet_path)
tissue_metadata = pd.read_csv(os.path.join(project_dir, 'resources', 'tissue_metadata.csv'))

tissue_id = sample_sheet.loc[sample_sheet['section_id'] == sample_id, 'tissue_id'].values[0]
subtype   = tissue_metadata.loc[tissue_metadata['tissue_id'].astype(str) == str(tissue_id),
                                'clinical_subtype'].values[0]


# --- Output paths

zone_dir  = os.path.join(res_dir, '05_xenium_mechanisms', '07_TLS_zonation', '01_classify_zones', 'csv_files')
plots_dir = os.path.join(res_dir, '05_xenium_mechanisms', '07_TLS_zonation', '03_zoom_GCs', 'plots')
os.makedirs(plots_dir, exist_ok=True)
pdf_path = os.path.join(plots_dir, f'{sample_id}_GC_zoom.pdf')


# --- Load GC bins and their discrete-GC labels (computed in step 01)

bin_df = pd.read_csv(os.path.join(zone_dir, f'{sample_id}_GC_zone_bins.csv'))
gc_bins = bin_df.dropna(subset=['GC_zone']).copy() if 'GC_zone' in bin_df.columns else bin_df.iloc[0:0].copy()
if 'discrete_GC' in gc_bins.columns:
    gc_bins = gc_bins.dropna(subset=['discrete_GC']).copy()
    gc_bins['discrete_GC'] = gc_bins['discrete_GC'].astype(int)

if len(gc_bins) == 0:
    print(f'No GC bins for {sample_id} - writing placeholder.')
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.text(0.5, 0.5, f'No GC regions detected in {sample_id}\n(subtype: {subtype})',
            ha='center', va='center', fontsize=18, transform=ax.transAxes,
            bbox=dict(boxstyle='round,pad=0.6', facecolor='#f0f0f0', edgecolor='gray'))
    ax.axis('off')
    fig.suptitle(f'{sample_id} - GC zoom plots', fontweight='bold', y=1.1)
    fig.savefig(pdf_path, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved placeholder PDF: {pdf_path}')
    print("Script successfully completed")
    raise SystemExit(0)

# Order discrete GCs by size (largest first)
comp_sizes = gc_bins['discrete_GC'].value_counts()
keep_labels = comp_sizes.index.tolist()
keep_sizes  = comp_sizes.values.tolist()
print(f'{len(keep_labels)} discrete GCs to plot')


# --- Top-N GC signature genes for this sample

genes_used = pd.read_csv(os.path.join(zone_dir, f'{sample_id}_GC_genes_used.csv'))
top_genes = genes_used['genes_used'].tolist()[:top_n_genes]
gene_groups = [top_genes[i:i + genes_per_ax] for i in range(0, len(top_genes), genes_per_ax)]
print(f'Top {len(top_genes)} GC genes: {top_genes}')


# --- Transcript colors

tx_color_csv = os.path.join(project_dir, 'resources', 'plotting_colors', 'transcript_colors.csv')
tx_color_lut = pd.read_csv(tx_color_csv, index_col=0)['color'].to_dict() if os.path.exists(tx_color_csv) else {}
fallback_palette = list(plt.get_cmap('tab10').colors)
tx_colors = {}
for i, g in enumerate(top_genes):
    if g in tx_color_lut:
        tx_colors[g] = tx_color_lut[g]
    else:
        tx_colors[g] = mcolors.to_hex(fallback_palette[i % len(fallback_palette)])
        print(f'  no color in transcript_colors.csv for {g}; using {tx_colors[g]}')


# --- Load adata (for microns_per_pixel of the aligned image)

adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample_id}.h5ad'))
microns_per_pixel = adata.uns['spatial'][sample_id]['scalefactors']['microns_per_pixel']


# --- Load the aligned full-res H&E image (image_space)

image_dir  = os.path.join(project_dir, 'data', 'images', 'he', 'aligned_images', 'xenium')
image_path = os.path.join(image_dir, f'{sample_id}.ome.tif')
# The aligned H&E is a large whole-slide image, so it is never loaded whole: the GC crops are
# read tile by tile at full resolution (read_tiff_region) and the overview panel uses a
# low-resolution pyramid level drawn in full-resolution pixel coordinates.
img_h, img_w = tiff_shape(image_path)
overview_img = read_tiff_overview(image_path)
print(f'Aligned H&E: {img_h} x {img_w} px; overview level {overview_img.shape}', flush=True)


# --- Load spatialdata + filter transcripts to the top-N GC genes

spat_dat_dir = os.path.join(project_dir, 'data', 'spatialdata_objects', 'xenium')
zarr_path = os.path.join(spat_dat_dir, f'{sample_id}.zarr')
if not os.path.exists(zarr_path):
    zip_path = os.path.join(spat_dat_dir, f'{sample_id}.zip')
    print(f'Unzipping {zip_path} -> {zarr_path}')
    shutil.unpack_archive(zip_path, zarr_path)

sdata = sd.read_zarr(zarr_path)
sdata.points['transcripts'] = sdata.points['transcripts'][
    sdata.points['transcripts']['feature_name'].isin(top_genes)
]
# Transform only the transcripts: sd.transform(sdata, ...) would also resample the
# full-resolution label rasters into image_space
tx_df = sd.transform(sdata.points['transcripts'], to_coordinate_system='image_space').compute()
print(f'Loaded {len(tx_df)} transcripts for the top {len(top_genes)} GC genes', flush=True)


# --- Xenium micron -> aligned-image (image_space) pixel transform
#
# The aligned .ome.tif and the image_space-transformed transcripts share the same
# pixel coordinate system (see scripts/figures/figure5/5B.py, 5E.py). GC bins come from
# step 01 in Xenium micron coordinates (x_centroid/y_centroid frame). We map micron -> px
# with the same affine spatialdata applies to the transcripts (native -> image_space), so
# bins land on the aligned image exactly like the transcripts.
#
# Note: a per-sample alignment matrix also exists at
#   data/images/he/alignment_files/xenium/{sample_id}_matrix.csv
# but is not needed here because the .ome.tif is already registered to image_space and the
# spatialdata affine is used to map the bins. This assumes transcript-native coords share
# the bin micron frame, which holds for Xenium (cells + transcripts use one global micron
# frame).

affine = get_transformation(sdata.points['transcripts'],
                            to_coordinate_system='image_space').to_affine_matrix(
    input_axes=('x', 'y'), output_axes=('x', 'y'))


def micron_to_he_pixel(x_um, y_um):
    coords = np.column_stack([np.asarray(x_um, dtype=float), np.asarray(y_um, dtype=float)])
    hom = np.hstack([coords, np.ones((coords.shape[0], 1))])
    out = (affine @ hom.T).T
    return out[:, 0], out[:, 1]


def crop_image(x0, x1, y0, y1):
    return read_tiff_region(image_path, x0, x1, y0, y1)


# --- Plot one page per discrete GC

mpl.rcParams['axes.grid'] = False
pad_px = pad_um / microns_per_pixel

with PdfPages(pdf_path) as pdf:
    for rank, (comp_id, comp_n_bins) in enumerate(zip(keep_labels, keep_sizes), 1):
        comp_bins = gc_bins[gc_bins['discrete_GC'] == comp_id]

        # Transform the 4 corners of every bin (microns) to aligned-image pixels; the
        # axis-aligned bbox of the transformed corners is the tight zoom region.
        n_b = len(comp_bins)
        xs_um = np.repeat(comp_bins['bin_x'].values, 4)
        ys_um = np.repeat(comp_bins['bin_y'].values, 4)
        dx = np.tile([-0.5, 0.5, -0.5, 0.5], n_b) * grid_um
        dy = np.tile([-0.5, -0.5, 0.5, 0.5], n_b) * grid_um
        cx_px, cy_px = micron_to_he_pixel(xs_um + dx, ys_um + dy)

        x0_t, x1_t = cx_px.min() - pad_px, cx_px.max() + pad_px
        y0_t, y1_t = cy_px.min() - pad_px, cy_px.max() + pad_px
        side = max(x1_t - x0_t, y1_t - y0_t)
        ctr_x, ctr_y = (x0_t + x1_t) / 2, (y0_t + y1_t) / 2
        x0_px, x1_px = ctr_x - side / 2, ctr_x + side / 2
        y0_px, y1_px = ctr_y - side / 2, ctr_y + side / 2

        crop_img, (cx0, cx1, cy0, cy1) = crop_image(x0_px, x1_px, y0_px, y1_px)
        crop_w, crop_h = cx1 - cx0, cy1 - cy0

        # Bin centers + per-bin pixel size for the GC overlay
        bin_cx_px, bin_cy_px = micron_to_he_pixel(comp_bins['bin_x'].values, comp_bins['bin_y'].values)
        bin_x0_p, bin_y0_p = micron_to_he_pixel(comp_bins['bin_x'].values - grid_um / 2,
                                                comp_bins['bin_y'].values - grid_um / 2)
        bin_x1_p, bin_y1_p = micron_to_he_pixel(comp_bins['bin_x'].values + grid_um / 2,
                                                comp_bins['bin_y'].values + grid_um / 2)
        bin_size_px = np.hypot(bin_x1_p - bin_x0_p, bin_y1_p - bin_y0_p) / np.sqrt(2)

        # 50 um grid covering the crop
        anchor_x = comp_bins['bin_x'].iloc[0] - grid_um / 2
        anchor_y = comp_bins['bin_y'].iloc[0] - grid_um / 2
        buf = max(20, int(np.ceil(side * microns_per_pixel / grid_um)) + 5)
        gx_min = anchor_x - buf * grid_um
        gx_max = anchor_x + (comp_bins['bin_ix'].max() - comp_bins['bin_ix'].min() + buf) * grid_um
        gy_min = anchor_y - buf * grid_um
        gy_max = anchor_y + (comp_bins['bin_iy'].max() - comp_bins['bin_iy'].min() + buf) * grid_um
        x_grid_um = np.arange(gx_min, gx_max + grid_um / 2, grid_um)
        y_grid_um = np.arange(gy_min, gy_max + grid_um / 2, grid_um)

        vert_s = np.column_stack([x_grid_um, np.full_like(x_grid_um, gy_min)])
        vert_e = np.column_stack([x_grid_um, np.full_like(x_grid_um, gy_max)])
        horz_s = np.column_stack([np.full_like(y_grid_um, gx_min), y_grid_um])
        horz_e = np.column_stack([np.full_like(y_grid_um, gx_max), y_grid_um])
        gs_x, gs_y = micron_to_he_pixel(np.vstack([vert_s, horz_s])[:, 0], np.vstack([vert_s, horz_s])[:, 1])
        ge_x, ge_y = micron_to_he_pixel(np.vstack([vert_e, horz_e])[:, 0], np.vstack([vert_e, horz_e])[:, 1])
        grid_segments = np.stack([
            np.column_stack([gs_x - cx0, gs_y - cy0]),
            np.column_stack([ge_x - cx0, ge_y - cy0]),
        ], axis=1)

        # Outline of the GC bin union: emit each bin edge whose neighbor is not in the component
        comp_set = set(zip(comp_bins['bin_ix'].astype(int).values, comp_bins['bin_iy'].astype(int).values))
        half = grid_um / 2
        outline_a_um, outline_b_um = [], []
        for ix_, iy_, bx_, by_ in zip(comp_bins['bin_ix'].astype(int).values,
                                      comp_bins['bin_iy'].astype(int).values,
                                      comp_bins['bin_x'].values, comp_bins['bin_y'].values):
            if (ix_ - 1, iy_) not in comp_set:
                outline_a_um.append((bx_ - half, by_ - half)); outline_b_um.append((bx_ - half, by_ + half))
            if (ix_ + 1, iy_) not in comp_set:
                outline_a_um.append((bx_ + half, by_ - half)); outline_b_um.append((bx_ + half, by_ + half))
            if (ix_, iy_ - 1) not in comp_set:
                outline_a_um.append((bx_ - half, by_ - half)); outline_b_um.append((bx_ + half, by_ - half))
            if (ix_, iy_ + 1) not in comp_set:
                outline_a_um.append((bx_ - half, by_ + half)); outline_b_um.append((bx_ + half, by_ + half))
        outline_a_um = np.array(outline_a_um)
        outline_b_um = np.array(outline_b_um)
        oa_x, oa_y = micron_to_he_pixel(outline_a_um[:, 0], outline_a_um[:, 1])
        ob_x, ob_y = micron_to_he_pixel(outline_b_um[:, 0], outline_b_um[:, 1])
        outline_segments = np.stack([
            np.column_stack([oa_x - cx0, oa_y - cy0]),
            np.column_stack([ob_x - cx0, ob_y - cy0]),
        ], axis=1)

        # Transcripts inside the bbox, shuffled to avoid one gene plotting over another
        tx_in = tx_df[(tx_df['x'] >= cx0) & (tx_df['x'] < cx1)
                      & (tx_df['y'] >= cy0) & (tx_df['y'] < cy1)].copy()
        tx_in['x_adj'] = tx_in['x'] - cx0
        tx_in['y_adj'] = tx_in['y'] - cy0
        if len(tx_in):
            tx_in = tx_in.sample(frac=1, random_state=42).reset_index(drop=True)

        side_um = side * microns_per_pixel
        scalebar_um = pick_scalebar_um(side_um)

        def make_scalebar():
            return ScaleBar(microns_per_pixel, units='µm', location='lower right', frameon=False,
                            scale_loc='none', fixed_value=scalebar_um, height_fraction=0.04, border_pad=1.2)

        n_panels = 4 + len(gene_groups)
        fig, axs = plt.subplots(1, n_panels, figsize=(n_panels * 6, 6))

        # Panel 1: overview with bbox
        ax = axs[0]
        ax.imshow(overview_img, extent=(0, img_w, img_h, 0))
        ax.add_patch(mpatches.Rectangle((cx0, cy0), crop_w, crop_h, linewidth=2,
                                        edgecolor='black', facecolor='none'))
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f'{sample_id}\nfull tissue overview', fontweight='bold', fontsize=TITLE_FS)

        # Panel 2: zoom H&E + GC outline
        ax = axs[1]
        ax.imshow(crop_img)
        ax.add_collection(LineCollection(outline_segments, colors='cyan', linewidths=2.5, zorder=3))
        ax.set_xlim(0, crop_w); ax.set_ylim(crop_h, 0)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title('GC region outline', fontweight='bold', fontsize=TITLE_FS)
        ax.add_artist(make_scalebar())

        # Panel 3: zoom H&E + 50 um grid + GC bins
        ax = axs[2]
        ax.imshow(crop_img)
        ax.add_collection(LineCollection(grid_segments, colors='white', linewidths=1.5, alpha=0.75, zorder=2))
        bin_rects = [mpatches.Rectangle((bcx - cx0 - bsz / 2, bcy - cy0 - bsz / 2), bsz, bsz)
                     for bcx, bcy, bsz in zip(bin_cx_px, bin_cy_px, bin_size_px)]
        ax.add_collection(PatchCollection(bin_rects, facecolors='cyan', edgecolors='none', alpha=0.55, zorder=3))
        ax.set_xlim(0, crop_w); ax.set_ylim(crop_h, 0)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f'GC bins + {int(grid_um)} µm grid', fontweight='bold', fontsize=TITLE_FS)
        ax.add_artist(make_scalebar())

        # Panel 4: zoom H&E only
        ax = axs[3]
        ax.imshow(crop_img)
        ax.set_xlim(0, crop_w); ax.set_ylim(crop_h, 0)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f'H&E zoom (GC #{rank}, {comp_n_bins} bins)', fontweight='bold', fontsize=TITLE_FS)
        ax.add_artist(make_scalebar())

        # Transcript panels: one per ranked group of genes
        s = float(np.clip(1e-3 * crop_w * crop_h / max(1, len(tx_in)), 30, 300))
        for gi, group in enumerate(gene_groups):
            ax = axs[4 + gi]
            ax.imshow(crop_img, alpha=0.5)
            for gene in group:
                sub = tx_in[tx_in['feature_name'] == gene]
                if not sub.empty:
                    edge = tuple(max(c - 0.3, 0) for c in mcolors.to_rgb(tx_colors[gene]))
                    ax.scatter(sub['x_adj'], sub['y_adj'], s=s, c=tx_colors[gene],
                               alpha=0.9, linewidths=0.3, edgecolors=edge)
            ax.set_xlim(0, crop_w); ax.set_ylim(crop_h, 0)
            ax.set_xticks([]); ax.set_yticks([])
            g_start = gi * genes_per_ax + 1
            g_end = g_start + len(group) - 1
            ax.set_title(f'GC transcripts {g_start}–{g_end}', fontweight='bold', fontsize=TITLE_FS)
            legend_handles = [Line2D([0], [0], marker='o', color='w', label=g,
                                     markerfacecolor=tx_colors[g], markersize=14) for g in group]
            ax.legend(handles=legend_handles, loc='upper center', bbox_to_anchor=(0.5, -0.02),
                      frameon=False, ncol=3, columnspacing=1.0, handletextpad=0.2, fontsize=LEGEND_FS - 4)

        fig.suptitle(f'{sample_id} ({subtype}) - discrete GC {rank}/{len(keep_labels)} '
                     f'({comp_n_bins} bins, n_transcripts={len(tx_in)}, scalebar={scalebar_um} µm)',
                     fontweight='bold', y=1.1, fontsize=SUP_FS)
        fig.subplots_adjust(wspace=0.06)
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)


# --- Finish

print(f'\nSaved: {pdf_path}')
print("Script successfully completed")
