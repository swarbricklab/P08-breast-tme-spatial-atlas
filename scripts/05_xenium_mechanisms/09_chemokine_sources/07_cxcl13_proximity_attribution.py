#!/usr/bin/env python3

# --- Aim
# S9K + S9L: Tfh/FDC density across TLS zones and CXCL13 source attribution.
#
# S9K (density): Tfh and FDC cell density (cells / mm^2) in the B vs T cell zone, per-sample
# points over a box, with a two-sided paired Wilcoxon between zones. The Tfh and FDC panels
# share a common y-axis so the two cell types can be compared on the same scale. The FDC panel
# uses only samples with a recovered FDC population.
#
# S9L (proximity attribution): CXCL13 transcripts are attributed to reference-defined producers.
# Candidate producers are miniatlas reference subsets in which >= 15% of cells express CXCL13
# (Tfh CD4+, CD8+ LAG3 and cycling MKI67+ T cells), plus FDC (absent from the reference;
# recovered separately). A cell is a proximity anchor only if it is a producer type AND carries
# >= 5 CXCL13 transcripts (subset identity alone is not evidence of expression). Each CXCL13
# transcript is assigned the TLS zone of its segmented cell (nearest cell for unassigned
# transcripts). Within the B and T cell zones each transcript is attributed to its nearest
# anchor if that anchor lies within 15 um, otherwise to "other". Among the "other" transcripts,
# those whose assigned cell carries fewer than 30 CXCL13 transcripts are treated as detection
# background and dropped. The remaining transcripts are expressed as a percentage of in-zone
# transcripts per source; any source contributing < 1 % of a zone's transcripts is pooled into
# a single "Other (each <1%)" segment. The cohort mean (over samples with a recovered FDC
# population) is shown as stacked bars, B vs T zone.
#
# Samples with no FDC barcode file (or zero FDC cells) are skipped from the FDC-dependent
# outputs automatically (no sample IDs are hard-coded).

# --- Define environment
# env: spatdat

# --- Load packages

import argparse
import os
import shutil
import tempfile
import sys

import numpy as np
import pandas as pd
import anndata as ad
from scipy.spatial import cKDTree
from scipy import stats

import spatialdata as sd

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

mpl.rcParams['pdf.fonttype'] = 42
TICK_FS, SUP_FS, ANNO_FS = 8, 12, 7


# --- Argparse arguments
# Cohort: dispatches through submit_multi_sample_job_nci.sh -> run_multi_sample_job.sh, which
# passes --h5ad_dir --project_dir --res_dir --sample_sheet --samples. The Methods parameters
# (reference producers >= 15% CXCL13, anchor >= 5 CXCL13 tx, 15 um anchor radius, background
# < 30 CXCL13 tx, pool sources < 1 %) are hard-coded defaults.

parser = argparse.ArgumentParser(
    description="S9K/S9L: Tfh/FDC density + CXCL13 reference-producer proximity attribution (cohort)",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-s", "--samples", nargs='+', required=True)
parser.add_argument("--h5ad_dir", default=None, help="defaults to <project_dir>/data/anndata_objects")
parser.add_argument("--project_dir", required=True)
parser.add_argument("--res_dir", required=True)
parser.add_argument("--sample_sheet", default=None, help="Path to sample sheet CSV (else default)")
parser.add_argument("--spatdata_dir", default=None,
                    help="Path to spatialdata .zarr/.zip files (else <project_dir>/data/spatialdata_objects/xenium)")
parser.add_argument("--grid_um", type=float, default=50.0)
parser.add_argument("--min_pct_expr", type=int, default=15,
                    help="reference-producer threshold: a subset is a CXCL13 producer if >= this %% "
                         "of its reference cells express CXCL13 (reads 01_reference_producers over_<pct>)")
parser.add_argument("--min_producer_count", type=int, default=5,
                    help="a producer-type / FDC cell serves as a proximity anchor only if it carries "
                         ">= this many CXCL13 transcripts (subset identity alone is not evidence)")
parser.add_argument("--d_max", type=float, default=15.0,
                    help="max distance (um) from a CXCL13 transcript to its nearest producer anchor "
                         "to attribute it to that source; beyond this -> 'other'")
parser.add_argument("--min_cell_count_other", type=int, default=30,
                    help="'other' transcripts whose assigned cell carries fewer than this many "
                         "CXCL13 transcripts are dropped as detection background")
parser.add_argument("--min_contribution_pct", type=float, default=1.0,
                    help="sources whose cohort-mean contribution is below this %% in both zones "
                         "are pooled into a single 'Other (each <X%%)' segment")

args = vars(parser.parse_args())
project_dir          = args["project_dir"]
res_dir              = args["res_dir"]
sample_list          = args["samples"]
h5ad_dir             = args["h5ad_dir"] or os.path.join(project_dir, 'data', 'anndata_objects')
spatdata_dir         = args["spatdata_dir"] or os.path.join(project_dir, 'data',
                                                            'spatialdata_objects', 'xenium')
grid_um              = args["grid_um"]
min_pct_expr         = args["min_pct_expr"]
min_producer_count   = args["min_producer_count"]
d_max                = args["d_max"]
min_cell_count_other = args["min_cell_count_other"]
min_contribution_pct = args["min_contribution_pct"]

if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_xenium.csv')

TFH_LABEL   = 'T_cells_c3_CD4+_Tfh_CXCL13'
FDC_LABEL   = 'FDC'
TARGET_GENE = 'CXCL13'
zones_of_i  = ['B_cell_zone', 'T_cell_zone']   # zone values as written by 07_TLS_zonation/01
zones_disp  = ['B cell zone', 'T cell zone']
bin_area_mm2 = (grid_um ** 2) / 1e6

OTHER_PREFIX      = 'OTHER::'   # namespaces an 'other' cell subset before display
MERGED_POOL_LABEL = f'Other (each <{min_contribution_pct:g}%)'

zone_palette = {
    'B_cell_zone': mpl.colors.cnames['cornflowerblue'],
    'T_cell_zone': mpl.colors.cnames['crimson'],
}
# Reserved category colours; producer subsets get their miniatlas colour (below).
src_palette = {
    'FDC':             '#f5a623',
    'other':           '#cfcfcf',
    'background':      '#e8e8e8',
    'UNASSIGNED':      '#cccccc',
    MERGED_POOL_LABEL: '#cfcfcf',
}


def pretty_zone(z):
    return z.replace('_', ' ')


def format_count(value):
    v = int(value)
    if v < 1000:
        return f'{v}'
    if v < 10_000:
        return f'{v / 1000:.1f}k'
    return f'{v // 1000}k'


def stat_label(p):
    if p < 1e-4: return '****'
    if p < 1e-3: return '***'
    if p < 1e-2: return '**'
    if p < 5e-2: return '*'
    return 'ns'


def paired_wilcoxon(values_B, values_T, sample_ids_B, sample_ids_T):
    """Two-sided paired Wilcoxon between B/T zones for samples present in both."""
    common = sorted(set(sample_ids_B) & set(sample_ids_T))
    if len(common) < 3:
        return 1.0, 'ns'
    sB, sT = list(sample_ids_B), list(sample_ids_T)
    vB = np.array([values_B[sB.index(s)] for s in common])
    vT = np.array([values_T[sT.index(s)] for s in common])
    try:
        p = float(stats.wilcoxon(vB, vT, zero_method='wilcox').pvalue)
    except ValueError:
        p = 1.0
    return p, stat_label(p)


def add_stat_bracket(ax, x1, x2, y_top, label, line_h_frac=0.025, lw=0.6, fontsize=None):
    fs = fontsize or 9
    ymin, ymax = ax.get_ylim()
    line_h = (ymax - ymin) * line_h_frac
    ax.plot([x1, x1, x2, x2], [y_top, y_top + line_h, y_top + line_h, y_top],
            color='black', linewidth=lw, clip_on=False)
    fig = ax.figure
    axis_h_in = ax.get_position().height * fig.get_figheight()
    data_per_inch = (ymax - ymin) / axis_h_in if axis_h_in > 0 else 0
    glyph_drop = (fs * 0.30 / 72.0) * data_per_inch
    ax.text((x1 + x2) / 2, y_top + line_h - glyph_drop, label,
            ha='center', va='baseline', fontsize=fs, clip_on=False)


# --- Local plotting helper

sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import save_multi_image   # noqa: E402


# --- Sample metadata

sample_sheet    = pd.read_csv(sample_sheet_path)
tissue_metadata = pd.read_csv(os.path.join(project_dir, 'resources', 'tissue_metadata.csv'))
merged          = sample_sheet.merge(tissue_metadata, on='tissue_id')
subtype_dict    = dict(zip(merged['section_id'], merged['clinical_subtype']))

subtype_order = ['Luminal A', 'Luminal B', 'Luminal B (HER2+)', 'HER2+', 'TNBC']
present_subtypes = sorted({subtype_dict[s] for s in sample_list if s in subtype_dict},
                          key=lambda x: subtype_order.index(x) if x in subtype_order else float('inf'))
samples_by_subtype = {st: sorted([s for s in sample_list if subtype_dict.get(s) == st])
                      for st in present_subtypes}
sorted_samples = [s for st in present_subtypes for s in samples_by_subtype[st]]


# --- Reference cell-type classifications + colours (for the S9L producer segments)

classif = pd.read_csv(os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas',
                                   'miniatlas_subset_classifications.csv'))
classif = classif.drop_duplicates(subset='miniatlas_subset', keep='first')
subset_to_lineage = dict(zip(classif['miniatlas_subset'], classif['lineage']))
subset_to_short   = dict(zip(classif['miniatlas_subset'], classif['miniatlas_subset_short']))
lineage_color = dict(pd.read_csv(os.path.join(project_dir, 'resources', 'plotting_colors',
                                              'miniatlas_lineage_colors.csv'))[['lineage', 'color']].values)
subset_color_lookup = dict(pd.read_csv(os.path.join(project_dir, 'resources', 'plotting_colors',
                                                    'miniatlas_subset_colors.csv'))[['miniatlas_subset', 'color']].values)


# --- Reference CXCL13 producers (>= min_pct_expr% of reference cells expressing) + FDC anchor set

producer_csv = os.path.join(res_dir, '05_xenium_mechanisms', '09_chemokine_sources',
                            '01_reference_producers', 'csv_files', 'cxcl13_producers',
                            f'over_{min_pct_expr}', 'cxcl13_producers.csv')
prod_df = pd.read_csv(producer_csv)
ref_producers = prod_df.loc[prod_df['CXCL13_producer'], 'subset'].tolist()
prod_pct = dict(zip(prod_df['subset'], prod_df['CXCL13_pct']))
# Order producers by lineage then reference CXCL13 %, so bars/legends read consistently; FDC last.
lineage_rank = {'Immune': 0, 'Stromal': 1, 'Epithelial': 2}
ref_producers = sorted(ref_producers,
                       key=lambda s: (lineage_rank.get(subset_to_lineage.get(s), 9),
                                      -prod_pct.get(s, 0.0)))
producer_labels = ref_producers + [FDC_LABEL]
print(f'Reference producers (>= {min_pct_expr}% CXCL13, n={len(ref_producers)}): {ref_producers}')
print(f'+ FDC (always). Total anchor types: {len(producer_labels)}')


def _base_subset_color(cat):
    if cat in subset_color_lookup:
        return subset_color_lookup[cat]
    lineage = subset_to_lineage.get(cat, None)
    return lineage_color.get(lineage, 'lightgray')


# Palette covering every S9L category. Tfh + CD8 LAG3 reference colours are near-white, so we
# assign explicit visible colours matching the reference-producer scheme.
cat_palette = dict(src_palette)
for p in ref_producers:
    cat_palette[p] = _base_subset_color(p)
cat_palette[TFH_LABEL] = '#888EC6'
cat_palette['T_cells_c8_CD8+_LAG3'] = '#c2456e'


def display_name(label):
    if isinstance(label, str) and label.startswith(OTHER_PREFIX):
        sub = label[len(OTHER_PREFIX):]
        return subset_to_short.get(sub, sub)
    if label in (FDC_LABEL, 'other', 'background', 'UNASSIGNED', MERGED_POOL_LABEL):
        return label
    return subset_to_short.get(label, label)


def color_for_category(cat):
    if isinstance(cat, str) and cat.startswith(OTHER_PREFIX):
        cat = cat[len(OTHER_PREFIX):]
    if cat in cat_palette:
        return cat_palette[cat]
    if cat in subset_color_lookup:
        return subset_color_lookup[cat]
    lineage = subset_to_lineage.get(cat, None)
    if lineage in lineage_color:
        return lineage_color[lineage]
    return 'lightgray'


# --- Paths

zone_dir = os.path.join(res_dir, '05_xenium_mechanisms', '07_TLS_zonation',
                        '01_classify_zones', 'csv_files')
fdc_dir  = os.path.join(res_dir, '05_xenium_mechanisms', '08_FDC_recovery',
                        'per_sample', 'csv_files')


def celltype_csv(sample, level='subset'):
    """Label-transfer predicted cell types (col 'pred_cell_type'); regenerated fallback."""
    p = os.path.join(res_dir, '02_cell_type_mapping', 'label_transfer', 'xenium',
                     f'{level}_level', f'{sample}_{level}_pred_cell_types.csv')
    if not os.path.exists(p):
        p = os.path.join(res_dir, '02_cell_type_mapping_regenerated', 'label_transfer',
                         'xenium', f'{level}_level', f'{sample}_{level}_pred_cell_types.csv')
    return p


# --- Identify samples with a recovered FDC population (auto-skip; no hard-coded IDs)

samples_with_fdc = []
for sample in sorted_samples:
    fdc_csv = os.path.join(fdc_dir, f'{sample}_FDC_Tfh_barcodes.csv')
    if os.path.isfile(fdc_csv):
        ov = pd.read_csv(fdc_csv)
        if 'group' in ov.columns and (ov['group'] == 'FDC').any():
            samples_with_fdc.append(sample)
print(f'Samples with recovered FDC population: {samples_with_fdc}')
cohort_note = f'cohort means over FDC-positive samples (n={len(samples_with_fdc)})'


# --- Per-sample: cell table (zone + cell type + centroids), n_bins per zone, CXCL13 transcripts

def load_sample(sample):
    # per-cell TLS zone (index = cell_id)
    cell_csv = os.path.join(zone_dir, f'{sample}_BT_zone_classification.csv')
    cell_df  = pd.read_csv(cell_csv)
    cell_df.rename(columns={cell_df.columns[0]: 'cell_id'}, inplace=True)
    cell_df = cell_df.rename(columns={'TLS_zone': 'zone'})
    cell_df['cell_id'] = cell_df['cell_id'].astype(str)

    # subset label
    subset = pd.read_csv(celltype_csv(sample), index_col=0)
    subset = subset[['pred_cell_type']].rename(columns={'pred_cell_type': 'subset'})
    subset.index = subset.index.astype(str)
    cell_df = cell_df.merge(subset, left_on='cell_id', right_index=True, how='left')

    # FDC / Tfh membership. cell_type (Tfh/FDC/other) drives the S9K density panels; is_fdc gates
    # the FDC proximity anchor for S9L.
    fdc_csv = os.path.join(fdc_dir, f'{sample}_FDC_Tfh_barcodes.csv')
    fdc_ids = set()
    if os.path.isfile(fdc_csv):
        ov = pd.read_csv(fdc_csv)
        if 'group' in ov.columns:
            fdc_ids = set(ov.loc[ov['group'] == 'FDC', 'cell_id'].astype(str))
    cell_df['is_fdc'] = cell_df['cell_id'].isin(fdc_ids)
    cell_df['cell_type'] = 'other'
    cell_df.loc[cell_df['subset'] == TFH_LABEL, 'cell_type'] = 'Tfh'
    cell_df.loc[cell_df['is_fdc'], 'cell_type'] = 'FDC'

    # micron centroids (obsm['spatial'] is image space in this project)
    adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample}.h5ad'), backed='r')
    coords = adata.obs[['x_centroid', 'y_centroid']].copy()
    coords.columns = ['x', 'y']
    coords.index = coords.index.astype(str)
    adata.file.close()
    cell_df = cell_df.merge(coords, left_on='cell_id', right_index=True, how='left')
    cell_df = cell_df.dropna(subset=['x', 'y'])

    # n_bins per zone (occupied 50 um bins with that zone label)
    binf = pd.read_csv(os.path.join(zone_dir, f'{sample}_BT_zone_bins.csv'))
    n_bins_by_zone = {z: int((binf['TLS_zone'] == z).sum()) for z in zones_of_i}

    # CXCL13 transcript points from the SpatialData object
    zip_path     = os.path.join(spatdata_dir, f'{sample}.zip')
    zarr_persist = os.path.join(spatdata_dir, f'{sample}.zarr')
    use_persistent = (os.path.isdir(zarr_persist) and
                      sum(os.path.getsize(os.path.join(dp, f))
                          for dp, _, fns in os.walk(zarr_persist) for f in fns) > 1e9)
    tmp_root = None
    try:
        if use_persistent:
            zarr_path = zarr_persist
            print(f'  Reading existing zarr ({zarr_persist})')
        else:
            tmp_root = tempfile.mkdtemp(prefix=f'cxcl13_prox_{sample}_', dir='/scratch/a56/sl6147/tmp')
            zarr_path = os.path.join(tmp_root, f'{sample}.zarr')
            print(f'  Unzipping {zip_path} -> {zarr_path}')
            shutil.unpack_archive(zip_path, zarr_path)
        sdata = sd.read_zarr(zarr_path)
        tx = sdata.points['transcripts']
        cx = tx[tx['feature_name'] == TARGET_GENE]
        keep_cols = [c for c in ['x', 'y', 'cell_id'] if c in cx.columns]
        tx_df = cx[keep_cols].compute().reset_index(drop=True)
    finally:
        if tmp_root and os.path.isdir(tmp_root):
            shutil.rmtree(tmp_root, ignore_errors=True)
            print('  Cleaned up tmp')

    tx_df = tx_df.dropna(subset=['x', 'y'])
    if len(tx_df) == 0:
        return cell_df, None, n_bins_by_zone

    # per-cell CXCL13 totals from segmentation (UNASSIGNED -> 0). Mapped onto transcripts (to
    # filter background) and onto cells (to gate bona-fide producer anchors).
    cell_id_str = tx_df['cell_id'].astype(str)
    assigned_mask = cell_id_str != 'UNASSIGNED'
    cell_totals = tx_df.loc[assigned_mask].groupby('cell_id').size().rename('cell_total_cxcl13')
    tx_df['cell_id'] = cell_id_str
    tx_df['cell_total_cxcl13'] = tx_df['cell_id'].map(cell_totals).fillna(0).astype(int)
    cell_df['cell_cxcl13'] = cell_df['cell_id'].map(cell_totals).fillna(0).astype(int)

    # nearest bona-fide producer anchor over all producer types. A cell anchors only if it is a
    # producer subset (or FDC) AND holds >= min_producer_count CXCL13 transcripts. FDC cells
    # anchor 'FDC'; producer-subset anchors exclude FDC barcodes to avoid double-counting.
    bona_fide = cell_df['cell_cxcl13'] >= min_producer_count
    tx_xy = tx_df[['x', 'y']].values
    dist_stack = np.full((len(tx_df), len(producer_labels)), np.inf)
    n_anchor = {}
    for j, p in enumerate(producer_labels):
        if p == FDC_LABEL:
            sel = cell_df['is_fdc'] & bona_fide
        else:
            sel = (cell_df['subset'] == p) & (~cell_df['is_fdc']) & bona_fide
        pxy = cell_df.loc[sel, ['x', 'y']].values
        n_anchor[p] = len(pxy)
        if len(pxy):
            dist_stack[:, j] = cKDTree(pxy).query(tx_xy, k=1)[0]
    print('  bona fide producer anchors: ' +
          ', '.join(f'{display_name(p)}={n_anchor[p]}' for p in producer_labels if n_anchor[p]))
    nearest_idx  = np.argmin(dist_stack, axis=1)
    nearest_dist = dist_stack[np.arange(len(tx_df)), nearest_idx]
    producer_arr = np.array(producer_labels, dtype=object)
    tx_df['nearest_producer'] = producer_arr[nearest_idx]
    tx_df['nearest_dist']     = nearest_dist
    tx_df.loc[~np.isfinite(nearest_dist), 'nearest_producer'] = 'other'

    # per-transcript zone = assigned cell's zone; UNASSIGNED -> nearest cell's zone
    lookup = cell_df.set_index('cell_id')[['zone', 'subset']]
    tx_df = tx_df.merge(lookup, left_on='cell_id', right_index=True, how='left')
    no_zone = tx_df['zone'].isna()
    if no_zone.any():
        all_xy   = cell_df[['x', 'y']].values
        all_zone = cell_df['zone'].values
        all_sub  = cell_df['subset'].values
        _, near_idx = cKDTree(all_xy).query(tx_df.loc[no_zone, ['x', 'y']].values, k=1)
        tx_df.loc[no_zone, 'zone']   = all_zone[near_idx]
        tx_df.loc[no_zone, 'subset'] = all_sub[near_idx]

    tx_df = tx_df[tx_df['zone'].isin(zones_of_i)].copy()
    tx_df['sample_id'] = sample
    tx_df['subtype']   = subtype_dict.get(sample)
    return cell_df, tx_df, n_bins_by_zone


def attribute_counts(tx_in_zone):
    """Aggregate to (sample, zone, attribution, attribution_detail) counts. Transcripts within
    d_max of a producer anchor are attributed to the nearest producer; the rest are 'other',
    except those in cells with < min_cell_count_other CXCL13 transcripts, which become
    'background'. attribution_detail is the producer / 'background', or the assigned cell's
    subset for true 'other'."""
    in_range = tx_in_zone['nearest_dist'].values <= d_max
    attribution = np.where(in_range, tx_in_zone['nearest_producer'].values, 'other').astype(object)

    other_mask = attribution == 'other'
    low_count  = tx_in_zone['cell_total_cxcl13'].values < min_cell_count_other
    attribution[other_mask & low_count] = 'background'

    detail = attribution.copy()
    keep_other = attribution == 'other'
    detail[keep_other] = tx_in_zone['subset'].where(tx_in_zone['subset'].notna(),
                                                     'UNASSIGNED').values[keep_other]
    df = tx_in_zone[['sample_id', 'subtype', 'zone']].copy()
    df['attribution'] = attribution
    df['attribution_detail'] = detail
    return (df.groupby(['sample_id', 'subtype', 'zone', 'attribution', 'attribution_detail'],
                       observed=True).size().rename('n_transcripts').reset_index())


# --- Load all samples

cell_dfs, count_parts, nbins_rows = [], [], []
for s in sorted_samples:
    print(f'\n=== {s} ===')
    cell_df, tx, n_bins_by_zone = load_sample(s)
    cell_df['sample_id'] = s
    cell_df['subtype']   = subtype_dict.get(s)
    cell_dfs.append(cell_df)
    for z in zones_of_i:
        nbins_rows.append({'sample_id': s, 'zone': z, 'n_bins': n_bins_by_zone[z]})
    if tx is not None:
        count_parts.append(attribute_counts(tx))
        print(f'  {int(tx.shape[0]):,} CXCL13 transcripts in B/T zones')

all_cells   = pd.concat(cell_dfs, ignore_index=True)
counts_long = pd.concat(count_parts, ignore_index=True)
nbins_df    = pd.DataFrame(nbins_rows)


# ==========================================================================================
# S9K: Tfh and FDC cell density (cells / mm^2) per zone (shared y-axis across panels)
# ==========================================================================================

cells_in_z = all_cells[all_cells['zone'].isin(zones_of_i)]
dens = (cells_in_z.groupby(['sample_id', 'zone'], observed=True)
        .agg(n_tfh=('cell_type', lambda x: (x == 'Tfh').sum()),
             n_fdc=('cell_type', lambda x: (x == 'FDC').sum()))
        .reset_index()
        .merge(nbins_df, on=['sample_id', 'zone'], how='left'))
dens['area_mm2']    = dens['n_bins'] * bin_area_mm2
dens['tfh_per_mm2'] = dens['n_tfh'] / dens['area_mm2'].replace(0, np.nan)
dens['fdc_per_mm2'] = dens['n_fdc'] / dens['area_mm2'].replace(0, np.nan)


def boxplot_panel(ax, values_per_x, colors, x_labels, samples_per_x):
    """Boxplot with overlaid per-sample dots + gray paired-sample connecting lines."""
    bp = ax.boxplot([np.asarray(v) for v in values_per_x],
                    positions=list(range(len(x_labels))), widths=0.55,
                    patch_artist=True, showfliers=False,
                    medianprops=dict(color='black', linewidth=1.2), zorder=2)
    for patch, color in zip(bp['boxes'], colors):
        rgb = mpl.colors.to_rgb(color)
        patch.set_facecolor((*rgb, 0.55)); patch.set_edgecolor(color); patch.set_linewidth(0.8)
    for elem in (*bp['whiskers'], *bp['caps']):
        elem.set_color('black'); elem.set_linewidth(0.8)
    rng = np.random.default_rng(seed=42)
    all_samples = sorted({s for sl in samples_per_x for s in sl})
    jitter = {(s, xi): float(rng.uniform(-0.12, 0.12))
              for xi in range(len(x_labels)) for s in all_samples}
    per_x_lookup = [dict(zip(sl, vl)) for sl, vl in zip(samples_per_x, values_per_x)]
    for s in all_samples:
        vals = [d.get(s, np.nan) for d in per_x_lookup]
        if all(np.isfinite(v) for v in vals):
            xs = [xi + jitter[(s, xi)] for xi in range(len(x_labels))]
            ax.plot(xs, vals, color='gray', alpha=0.4, linewidth=0.6, zorder=1)
    for xi, (sl, vl) in enumerate(zip(samples_per_x, values_per_x)):
        xs = [xi + jitter[(s, xi)] for s in sl]
        ax.scatter(xs, vl, c='black', s=14, linewidths=0, zorder=3)
    ax.set_xticks(range(len(x_labels)))
    ax.set_xticklabels(x_labels, rotation=90, fontsize=TICK_FS)
    ax.set_xlim(-0.5, len(x_labels) - 0.5)
    ax.tick_params(axis='y', labelsize=TICK_FS); ax.grid(False)
    for sp_ in ['top', 'right']:
        ax.spines[sp_].set_visible(False)


def density_arrays(df, metric):
    vB = df.loc[df['zone'] == 'B_cell_zone', metric].values
    vT = df.loc[df['zone'] == 'T_cell_zone', metric].values
    sB = df.loc[df['zone'] == 'B_cell_zone', 'sample_id'].tolist()
    sT = df.loc[df['zone'] == 'T_cell_zone', 'sample_id'].tolist()
    return vB, vT, sB, sT


tfh_B, tfh_T, tfh_sB, tfh_sT = density_arrays(dens, 'tfh_per_mm2')
fdc_src = dens[dens['sample_id'].isin(samples_with_fdc)]
fdc_B, fdc_T, fdc_sB, fdc_sT = density_arrays(fdc_src, 'fdc_per_mm2')

panels = [
    (tfh_B, tfh_T, tfh_sB, tfh_sT, 'Tfh density'),
    (fdc_B, fdc_T, fdc_sB, fdc_sT, 'FDC density'),
]
n = len(panels)
# Tfh and FDC share a common y-axis so the two cell types are directly comparable.
shared_ymax_data = max(
    max(np.max(vB) if len(vB) else 0, np.max(vT) if len(vT) else 0)
    for (vB, vT, _sB, _sT, _t) in panels)

panel_w_in, wspace_in, margin_left, margin_right = 0.9, 0.5, 0.85, 0.55
fig_w_in = margin_left + n * panel_w_in + (n - 1) * wspace_in + margin_right
fig = plt.figure(figsize=(fig_w_in, 2.9))
panel_w = panel_w_in / fig_w_in
wspace  = wspace_in / fig_w_in
panel_h, panel_y = 0.46, 0.26
p_left = margin_left / fig_w_in
axes = []
for i, (vB, vT, sB, sT, title) in enumerate(panels):
    ax = fig.add_axes([p_left, panel_y, panel_w, panel_h])
    boxplot_panel(ax, [vB, vT], [zone_palette[z] for z in zones_of_i], zones_disp,
                  samples_per_x=[list(sB), list(sT)])
    ax.set_title(title, fontsize=11, fontweight='bold', pad=14)
    if i == 0:
        ax.set_ylabel('cells / mm²', fontsize=TICK_FS, fontweight='bold')
    sci_fmt = ScalarFormatter(useMathText=True)
    sci_fmt.set_scientific(True); sci_fmt.set_powerlimits((-2, 3))
    ax.yaxis.set_major_formatter(sci_fmt)
    ax.yaxis.get_offset_text().set_fontsize(TICK_FS)
    ymax_data = max(np.max(vB) if len(vB) else 0, np.max(vT) if len(vT) else 0)
    p_value, label = paired_wilcoxon(vB, vT, sB, sT)
    print(f'  {title}: B vs T -> p={p_value:.3g} ({label})')
    if label != 'ns':
        ax.set_ylim(0, shared_ymax_data * 1.20)
        add_stat_bracket(ax, 0, 1, ymax_data * 1.05, label, fontsize=TICK_FS + 2)
    elif shared_ymax_data > 0:
        ax.set_ylim(0, shared_ymax_data * 1.20)
    else:
        ax.set_ylim(bottom=0)
    axes.append(ax)
    p_left += panel_w + wspace
fig.canvas.draw()
title_x = (axes[0].get_position().x0 + axes[-1].get_position().x1) / 2
top_extent = max(fig.transFigure.inverted().transform(
    ax.get_tightbbox(fig.canvas.get_renderer()))[1][1] for ax in axes)
fig.suptitle('S9K — Tfh and FDC density per TLS zone',
             x=title_x, y=top_extent + 0.05, ha='center', va='bottom',
             fontsize=SUP_FS, fontweight='bold')


# ==========================================================================================
# S9L: CXCL13 reference-producer source attribution (excluding background) per zone
# ==========================================================================================

counts_no_bg = counts_long[counts_long['attribution'] != 'background'].copy()


def per_sample_proportions(counts_df, value_col, cats=None):
    """Per (sample, zone) fraction of each category over that group's total transcripts."""
    rows = []
    for (sample, zone), sub in counts_df.groupby(['sample_id', 'zone'], observed=True):
        total = int(sub['n_transcripts'].sum())
        if total == 0:
            continue
        by_cat = sub.groupby(value_col, observed=True)['n_transcripts'].sum()
        row = {'sample_id': sample, 'zone': zone, 'subtype': sub['subtype'].iloc[0],
               'n_transcripts': total}
        for c in (cats if cats is not None else by_cat.index):
            row[c] = by_cat.get(c, 0) / total
        rows.append(row)
    return pd.DataFrame(rows)


def cohort_mean(per_sample_df, cats):
    src = per_sample_df[per_sample_df['sample_id'].isin(samples_with_fdc)]
    means  = src.groupby('zone', observed=True)[cats].mean().reindex(zones_of_i)
    counts = src.groupby('zone', observed=True)['n_transcripts'].sum().reindex(zones_of_i).fillna(0)
    return means, counts


threshold_frac = min_contribution_pct / 100.0

# attribution-level table (producers + aggregate 'other') -> which producers clear the threshold
attr_cats = list(producer_labels) + ['other']
per_sample_full = per_sample_proportions(counts_no_bg, 'attribution', cats=attr_cats)
cohort_full, _  = cohort_mean(per_sample_full, attr_cats)
prod_frac = {p: (float(np.nan_to_num(cohort_full[p].max())) if p in cohort_full.columns else 0.0)
             for p in producer_labels}
kept_prod = [p for p in producer_labels if prod_frac[p] >= threshold_frac]

# expanded table: producers stay; 'other' -> OTHER::<assigned subset>
w = counts_no_bg.copy()
w['plot_cat'] = [OTHER_PREFIX + str(det) if a == 'other' else a
                 for a, det in zip(w['attribution'], w['attribution_detail'])]
exp_cats = sorted(w['plot_cat'].unique())
ps_exp = per_sample_proportions(w, 'plot_cat', cats=exp_cats)
coh_exp, _ = cohort_mean(ps_exp, exp_cats)


def _remap(c):
    if c.startswith(OTHER_PREFIX):
        frac = float(np.nan_to_num(coh_exp[c].max())) if c in coh_exp.columns else 0.0
        return c if frac >= threshold_frac else MERGED_POOL_LABEL
    return c if c in kept_prod else MERGED_POOL_LABEL


remap = {c: _remap(c) for c in exp_cats}
w['plot_cat2'] = w['plot_cat'].map(remap)
final_cats = list(dict.fromkeys(remap.values()))
per_sample_collapsed = per_sample_proportions(w, 'plot_cat2', cats=final_cats)
cohort_collapsed, _ = cohort_mean(per_sample_collapsed, final_cats)

# sort most- to least-prominent (cohort contribution over zones); pin the merged pool last
contrib = {c: float(np.nan_to_num(cohort_collapsed[c].sum())) for c in final_cats}
named   = sorted([c for c in final_cats if c != MERGED_POOL_LABEL], key=lambda c: -contrib[c])
ordered_cats = named + ([MERGED_POOL_LABEL] if MERGED_POOL_LABEL in final_cats else [])
cohort_collapsed = cohort_collapsed[ordered_cats]


def stacked_bar(fig, ax_pos, proportions_df, xtick_labels, suptitle):
    ax_s = fig.add_axes(ax_pos)
    cats = list(proportions_df.columns)
    x = np.arange(len(proportions_df))
    for xi in x:
        bot = 0.0
        for cat in cats:
            pr = float(proportions_df.iloc[xi].get(cat, 0))
            if pr > 0:
                ax_s.bar(xi, pr, bottom=bot, color=color_for_category(cat),
                         width=0.7, linewidth=0.4, edgecolor='white')
                bot += pr
    ax_s.set_xlim(-0.5, len(x) - 0.5); ax_s.set_ylim(0, 1)
    ax_s.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax_s.set_yticklabels(['0', '25', '50', '75', '100'], fontsize=TICK_FS)
    ax_s.set_xticks(x); ax_s.set_xticklabels(xtick_labels, rotation=90, fontsize=TICK_FS)
    ax_s.set_ylabel('% of transcripts', fontsize=TICK_FS, fontweight='bold')
    for sp_ in ['top', 'right']:
        ax_s.spines[sp_].set_visible(False)
    handles = [mpl.patches.Patch(facecolor=color_for_category(c), label=display_name(c), linewidth=0)
               for c in cats]
    fig.canvas.draw()
    right_edge = ax_s.get_position().x1
    top_edge   = ax_s.get_position().y1
    fig.legend(handles=handles, bbox_to_anchor=(right_edge + 0.02, top_edge),
               bbox_transform=fig.transFigure, loc='upper left', frameon=False,
               handlelength=0.8, fontsize=TICK_FS, title='Source',
               title_fontproperties={'weight': 'bold', 'size': TICK_FS}, alignment='left')
    title_x = ax_s.get_position().x0 + ax_s.get_position().width / 2
    fig.suptitle(suptitle, x=title_x, y=top_edge + 0.04,
                 ha='center', va='bottom', fontsize=SUP_FS, fontweight='bold')


# --- S9L cohort stacked bar (B vs T zone)
fig = plt.figure(figsize=(3.0, 3.2))
p_left, p_bottom = 0.28, 0.24
p_w = 0.16 * len(zones_of_i)
stacked_bar(fig, [p_left, p_bottom, p_w, 0.56],
            proportions_df=cohort_collapsed,
            xtick_labels=[pretty_zone(z) for z in cohort_collapsed.index],
            suptitle=(f'S9L — cohort {TARGET_GENE} source attribution per zone\n'
                      f'(>= {min_pct_expr}% ref producers; d_max={d_max:g} um, anchor >= '
                      f'{min_producer_count} tx, background < {min_cell_count_other} tx; {cohort_note})'))


# --- Report + save

print('\nS9L cohort source attribution (% of in-zone transcripts, FDC-positive samples):')
if not cohort_collapsed.empty:
    print((cohort_collapsed * 100).round(1).to_string())

out_base = os.path.join(res_dir, '05_xenium_mechanisms', '09_chemokine_sources',
                        '07_cxcl13_proximity_attribution')
plot_dir = os.path.join(out_base, 'plots')
csv_dir  = os.path.join(out_base, 'csv_files')
os.makedirs(plot_dir, exist_ok=True)
os.makedirs(csv_dir, exist_ok=True)

dens.to_csv(os.path.join(csv_dir, 'cxcl13_tfh_fdc_density_per_zone.csv'), index=False)
counts_long.to_csv(os.path.join(csv_dir, 'cxcl13_producer_attribution_counts_long.csv'), index=False)
per_sample_collapsed.to_csv(os.path.join(csv_dir, 'cxcl13_producer_attribution_per_sample.csv'), index=False)
cohort_collapsed.to_csv(os.path.join(csv_dir, 'cxcl13_producer_attribution_cohort.csv'))

pdf_path = os.path.join(plot_dir, 'S9K_S9L_cxcl13_proximity_attribution.pdf')
save_multi_image(pdf_path)
print(f'\nSaved: {pdf_path}')
print('Script successfully completed')
