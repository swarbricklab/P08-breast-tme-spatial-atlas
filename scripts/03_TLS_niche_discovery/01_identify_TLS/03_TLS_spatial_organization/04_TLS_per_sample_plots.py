#!/usr/bin/env python3

# --- Aim
# Per-sample spatial overview of the discrete-TLS classification. For each Visium
# section a 3-panel figure shows: epithelial pathology annotation, TLS topology
# (intratumoral / peritumoral / stromal_distant), and DBSCAN spatial clustering
# (each cluster a distinct colour, lone TLS in grey). All sections are saved in a
# single multi-page PDF.

# --- Define environment
# env: spatial-pilot

# --- Load packages

import argparse
import os
import sys

import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches

mpl.rcParams['pdf.fonttype'] = 42
sc.set_figure_params(facecolor='white', fontsize=13)

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))
from helper_functions import auto_crop_scanpy, save_multi_image

from discrete_tls import define_discrete_tls


# --- Argparse arguments

parser = argparse.ArgumentParser(
    description="Per-sample spatial plots of TLS topology + DBSCAN clustering",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument("-s", "--samples", nargs="*", type=str, default=None,
                    help="List of section_ids (default: all in the sample sheet)")
parser.add_argument("--project_dir", required=True, help="Project directory (root)")
parser.add_argument("--h5ad_dir", required=True,
                    help="Path to directory for .h5ad files")
parser.add_argument("-r", "--res_dir", help="Path to results directory")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file")
parser.add_argument("--TLS_threshold", type=float, default=0.01)
parser.add_argument("--min_TLS_spots", type=int, default=3,
                    help="Minimum directly-adjacent spots for a discrete TLS "
                         "(Extended Data Fig. 3b stringent definition)")
parser.add_argument("--dbscan_eps_um", type=int, default=500,
                    help="DBSCAN eps (µm) whose labels are plotted; must exist as "
                         "a dbscan_label_epsX column in clustering_per_TLS.csv")

args = vars(parser.parse_args())

project_dir       = args["project_dir"]
h5ad_dir          = args["h5ad_dir"]
res_dir           = args["res_dir"] or os.path.join(project_dir, "results")
TLS_threshold     = args["TLS_threshold"]
TLS_threshold_adj = str(TLS_threshold).split('.')[-1]
min_TLS_spots     = args["min_TLS_spots"]
eps               = args["dbscan_eps_um"]


# --- Sample metadata

if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_visium.csv')

sample_sheet    = pd.read_csv(sample_sheet_path)
sample_lut      = sample_sheet.set_index('section_id')
tissue_metadata = pd.read_csv(os.path.join(project_dir, 'resources', 'tissue_metadata.csv'))
sample_meta = pd.merge(
    sample_sheet, tissue_metadata[['tissue_id', 'clinical_subtype']], on='tissue_id',
)
subtype_dict = dict(zip(sample_meta['section_id'], sample_meta['clinical_subtype']))

sample_list = args["samples"] or list(sample_sheet['section_id'])


# --- Locate the TLS_groups directory (permissive grouping input)

groups_dir = os.path.join(
    res_dir, '03_TLS_niche_discovery', '01_identify_TLS', '02_group_TLS_high',
    f'threshold_{TLS_threshold_adj}', 'csv_files',
)
print(f"Per-sample TLS_groups from: {groups_dir}")


# --- Per-TLS topology + clustering tables (produced by 01 and 02)

distances_path = os.path.join(
    res_dir, '03_TLS_niche_discovery', '01_identify_TLS',
    '03_TLS_spatial_organization', '01_TLS_distances_topology',
    f'threshold_{TLS_threshold_adj}', 'csv_files', 'distances_per_TLS.csv',
)
clustering_path = os.path.join(
    res_dir, '03_TLS_niche_discovery', '01_identify_TLS',
    '03_TLS_spatial_organization', '02_TLS_clustering',
    f'threshold_{TLS_threshold_adj}', 'csv_files', 'clustering_per_TLS.csv',
)
dist_df = (pd.read_csv(distances_path) if os.path.exists(distances_path)
           else pd.DataFrame(columns=['sample_id', 'TLS_group', 'tls_topology_2hop']))
clust_df = (pd.read_csv(clustering_path) if os.path.exists(clustering_path)
            else pd.DataFrame(columns=['sample_id', 'TLS_group', f'dbscan_label_eps{eps}']))
print(f"Topology rows: {len(dist_df)}   Clustering rows: {len(clust_df)}")


# --- Colour schemes

topo_categories = ['intratumoral', 'peritumoral', 'stromal_distant']
topo_colors = {
    'intratumoral':    '#b30000',
    'peritumoral':     '#fc8d59',
    'stromal_distant': '#74add1',
}
# Epithelial pathology labels (invasive / in_situ = cancer)
epi_categories = ['invasive', 'in_situ', 'benign', 'normal']
epi_colors = {
    'invasive': '#b30000',
    'in_situ':  '#fc8d59',
    'benign':   '#fee090',
    'normal':   '#74add1',
}
# Cluster colours: tab20 hues excluding the grey pair (indices 14, 15), which is
# reserved for lone TLS so a cluster can never be confused with a lone TLS.
_tab20 = mpl.colormaps['tab20']
CLUSTER_COLORS = [mpl.colors.to_hex(_tab20(i)) for i in range(20) if i not in (14, 15)]
LONE_COLOR   = '#4d4d4d'                # dark grey — lone (noise) TLS: neutral, reads on the
                                        # H&E, lighter than black, and no cluster uses grey

# Pretty legend labels + font sizes
topo_display = {'intratumoral': 'Intratumoral', 'peritumoral': 'Peritumoral',
                'stromal_distant': 'Stromal/distant'}
epi_display  = {'invasive': 'Invasive', 'in_situ': 'In situ',
                'benign': 'Benign', 'normal': 'Normal'}
TITLE_FS        = 18
LEGEND_FS       = 14
LEGEND_TITLE_FS = 15

cluster_label_col = f'dbscan_label_eps{eps}'


def add_legend_below(ax, labels, colors, title):
    """Add a colour-patch legend centred just below the axes (so it never
    overlaps the neighbouring panels)."""
    if not labels:
        return
    handles = [mpatches.Patch(facecolor=c, edgecolor='black', linewidth=0.4)
               for c in colors]
    leg = ax.legend(
        handles, labels, title=title, loc='upper center',
        bbox_to_anchor=(0.5, -0.02), ncol=min(len(labels), 3),
        frameon=False, fontsize=LEGEND_FS, title_fontsize=LEGEND_TITLE_FS,
        handlelength=1.1, handleheight=1.1, columnspacing=1.2, borderaxespad=0,
    )
    leg.get_title().set_fontweight('bold')


# --- Render one sample

def render_sample(sample):
    subtype    = subtype_dict.get(sample, '?')
    slide_size = sample_lut.loc[sample, 'slide_size']
    spot_size  = 1.8 if slide_size == '11 mm' else 1.3

    h5ad_path = os.path.join(h5ad_dir, 'visium', f'{sample}.h5ad')
    if not os.path.exists(h5ad_path):
        print(f'  skip {sample}: no h5ad at {h5ad_path}')
        return

    adata = ad.read_h5ad(h5ad_path)
    crop_coord = auto_crop_scanpy(adata, border_size=0.1)

    spot_diameter_pixel = adata.uns['spatial'][sample]['scalefactors']['spot_diameter_fullres']
    scalefactor_um = spot_diameter_pixel / 55  # pixels per µm
    coords_um = pd.DataFrame(
        adata.obsm['spatial'] / scalefactor_um,
        index=adata.obs_names, columns=['x', 'y'],
    )

    # Discrete-TLS labels per spot (Extended Data Fig. 3b definition)
    groups_file = os.path.join(groups_dir, f'{sample}_TLS_groups.csv')
    if os.path.exists(groups_file):
        gdf = pd.read_csv(groups_file, index_col=0)
        discrete = define_discrete_tls(gdf, coords_um, min_spots=min_TLS_spots)
        adata.obs['discrete_TLS'] = discrete.reindex(adata.obs_names)
    else:
        adata.obs['discrete_TLS'] = np.nan

    # Epithelial pathology annotation
    anno_file = os.path.join(
        project_dir, 'data', 'pathology_annotations', 'visium',
        f'{sample}_pathology_annotations.csv',
    )
    if os.path.exists(anno_file):
        adf = pd.read_csv(anno_file, index_col=0)
        adata.obs['epithelial'] = adf.reindex(adata.obs_names)['epithelial']
    else:
        adata.obs['epithelial'] = np.nan
    adata.obs['epithelial'] = pd.Categorical(
        adata.obs['epithelial'], categories=epi_categories,
    )

    # Topology per spot (mapped from per-TLS classification)
    sdist = dist_df[dist_df['sample_id'] == sample]
    map_topo = dict(zip(sdist['TLS_group'], sdist['tls_topology_2hop']))
    adata.obs['topology'] = pd.Categorical(
        adata.obs['discrete_TLS'].map(map_topo), categories=topo_categories,
    )

    # DBSCAN clustering per spot (cluster_k / lone)
    sclust = clust_df[clust_df['sample_id'] == sample]
    if cluster_label_col in sclust.columns:
        map_clu = {
            row['TLS_group']: ('lone' if int(row[cluster_label_col]) == -1
                               else f'cluster_{int(row[cluster_label_col])}')
            for _, row in sclust.iterrows() if pd.notna(row[cluster_label_col])
        }
    else:
        map_clu = {}
    clu_series = adata.obs['discrete_TLS'].map(map_clu)
    cluster_ints = sorted({
        int(c.split('_')[1]) for c in clu_series.dropna().unique() if c != 'lone'
    })
    clu_cats    = [f'cluster_{k}' for k in cluster_ints]
    clu_palette = [CLUSTER_COLORS[k % len(CLUSTER_COLORS)] for k in cluster_ints]
    if (clu_series == 'lone').any():
        clu_cats.append('lone')
        clu_palette.append(LONE_COLOR)
    adata.obs['clustering'] = pd.Categorical(clu_series, categories=clu_cats)

    fig, axs = plt.subplots(1, 3, figsize=(18, 7))

    # Panel 0: epithelial annotation
    sc.pl.spatial(
        adata, color='epithelial', size=spot_size, ax=axs[0],
        crop_coord=crop_coord, bw=False, alpha_img=0.5,
        palette=[epi_colors[c] for c in epi_categories],
        na_color=(0, 0, 0, 0), show=False, legend_loc=None,
        title='Epithelial annotation',
    )
    epi_present = [c for c in epi_categories if (adata.obs['epithelial'] == c).any()]
    add_legend_below(axs[0], [epi_display[c] for c in epi_present],
                     [epi_colors[c] for c in epi_present], 'Epithelial')

    # Panel 1: TLS topology
    if adata.obs['topology'].notna().any():
        sc.pl.spatial(
            adata, color='topology', size=spot_size, ax=axs[1],
            crop_coord=crop_coord, bw=False, alpha_img=0.5,
            palette=[topo_colors[c] for c in topo_categories],
            na_color=(0, 0, 0, 0), show=False, legend_loc=None,
            title='TLS topology (2-hop ring)',
        )
        topo_present = [c for c in topo_categories
                        if (adata.obs['topology'] == c).any()]
        add_legend_below(axs[1], [topo_display[c] for c in topo_present],
                         [topo_colors[c] for c in topo_present], 'Topology')
    else:
        sc.pl.spatial(adata, color=None, size=spot_size, ax=axs[1],
                      crop_coord=crop_coord, bw=False, alpha_img=0.5,
                      show=False, title='TLS topology (no TLS)')

    # Panel 2: DBSCAN clustering
    if clu_cats:
        sc.pl.spatial(
            adata, color='clustering', size=spot_size, ax=axs[2],
            crop_coord=crop_coord, bw=False, alpha_img=0.5,
            palette=clu_palette, na_color=(0, 0, 0, 0), show=False,
            legend_loc=None, title=f'DBSCAN clustering (eps={eps} µm)',
        )
        clu_labels = [('Lone' if c == 'lone' else f'Cluster {c.split("_")[1]}')
                      for c in clu_cats]
        add_legend_below(axs[2], clu_labels, clu_palette, 'Clustering')
    else:
        sc.pl.spatial(adata, color=None, size=spot_size, ax=axs[2],
                      crop_coord=crop_coord, bw=False, alpha_img=0.5,
                      show=False, title=f'DBSCAN clustering (no TLS)')

    # Larger, bold panel titles; drop the pixel-coordinate axis labels
    for ax in axs:
        ax.set_title(ax.get_title(), fontsize=TITLE_FS, fontweight='bold')
        ax.set_xlabel('')
        ax.set_ylabel('')

    n_tls = int(adata.obs['discrete_TLS'].nunique())
    topo_counts = sdist['tls_topology_2hop'].value_counts()
    caption = ', '.join(
        f'{c}={int(topo_counts.get(c, 0))}' for c in topo_categories
    )
    plt.suptitle(
        f'{sample} ({subtype}, {slide_size}) — {n_tls} discrete TLS\n'
        f'Topology: {caption}',
        y=1.02, fontsize=16, fontweight='bold',
    )

    del adata
    print(f'  {sample} done')


# --- Loop over samples + save one multi-page PDF

for sample in sample_list:
    render_sample(sample)

out_dir = os.path.join(
    res_dir, '03_TLS_niche_discovery', '01_identify_TLS',
    '03_TLS_spatial_organization', '04_TLS_per_sample_plots',
    f'threshold_{TLS_threshold_adj}', 'plots',
)
os.makedirs(out_dir, exist_ok=True)
pdf_path = os.path.join(out_dir, 'TLS_per_sample_topology_clustering.pdf')
save_multi_image(pdf_path)

print(f'\nSaved: {pdf_path}')
print("Script successfully completed")
