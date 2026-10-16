# 10_HEV_GC_proximity

Stratify HEV endothelial cells by **proximity to a germinal centre (GC)** and test which HEV genes (notably SPP1) are enriched close to GCs. Produces Extended Data Fig. 10i,j.

| Script | Env | Mode | Figure |
|---|---|---|---|
| `01_classify_HEV_by_GC_dist.py` | xen-hev | `--sample_id` | — (upstream) |
| `02_DEG_close_vs_far.py` | xen-hev | `--samples` | — (DEG table) |
| `03_volcano_dotplot.py` | shapely | `--samples` | Extended Data Fig. 10i,j |

**`01_classify_HEV_by_GC_dist.py`** — for each intra-TLS HEV niche (`close_to_TLS`, ≤100 µm from a TLS), computes the distance from the niche centroid to the nearest GC bin edge (0 if inside a GC bin) and labels it close (≤500 µm) / far. Writes `{s}_intra_TLS_HEV_GC_dist.csv`. GC-negative samples yield an empty / all-far table.

**`02_DEG_close_vs_far.py`** — pools ACKR1+ HEV endothelial cells across the GC-positive samples and, for each gene detected in ≥10 % of the pooled cells, tests close (≤500 µm) vs far log-transformed expression by two-sided Mann-Whitney U + BH. Writes the per-gene table `DEG_close_vs_far.csv`.

**`03_volcano_dotplot.py`** — renders `DEG_close_vs_far.csv` as the **Extended Data Fig. 10i** volcano (mean log-transformed difference close − far vs −log10 BH-P; purple = higher close to GC, brown = lower, grey = ns; canonical HEV markers outlined + bold-labelled, all significant genes labelled via adjustText) and the **Extended Data Fig. 10j** transposed dot plot of the significant genes (up-/down-close blocks, close vs far rows, dot size = % expressing, colour = per-gene-max mean expression). Runs on the `shapely` env because it needs `adjustText`.

**Auto-skip:** only GC-containing samples contribute (4 of 6 in this cohort); samples with no GC bins are derived from the data and skipped — no sample IDs are hard-coded. GC score threshold upstream = 0.8.

**Inputs:** GC bins from `07_TLS_zonation/01_classify_zones`; HEV niches from `results/05_xenium_mechanisms/02_identify_HEV_regions/csv_files/`; `data/anndata_objects/xenium/{s}.h5ad`; cell labels `results/02_cell_type_mapping/label_transfer/xenium/subset_level/` (label `Endothelial ACKR1`).

**Outputs:** `results/05_xenium_mechanisms/10_HEV_GC_proximity/…`.

**Run order:** requires `07_TLS_zonation` and `02_identify_HEV_regions`. Then `01` (per sample) → `02` (cohort DEG table) → `03` (cohort plots; reads `02`'s `DEG_close_vs_far.csv`).
