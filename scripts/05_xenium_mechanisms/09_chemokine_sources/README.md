# 09_chemokine_sources

Identify the cellular **sources of the TLS-organizing chemokines** (CCL19, CCL21, CXCL13) and their spatial organization across the B and T cell zones. Produces most of Figure S9.

| Script | Env | Mode | Figure |
|---|---|---|---|
| `01_reference_producers.py` | xen-hev | `--samples` | S9A |
| `02_extract_transcripts.py` | spatdat | `--sample_id` | — (upstream) |
| `03_impvl_vs_icaf.py` | xen-hev | `--samples` | S9B |
| `04_zone_expression_density.py` | xen-hev | `--samples` | S9E |
| `05_chemokine_gradients.py` | spatdat | `--samples` | S9H |
| `06_reference_producer_corr.py` | xen-hev | `--samples` | S9I |
| `07_cxcl13_proximity_attribution.py` | spatdat | `--samples` | S9K, S9L |

**`01_reference_producers.py`** — defines CCL19/CCL21 and CXCL13 producer subsets from the scRNA-seq reference (both `over_10` ≥10 % and `over_15` ≥15 % expressing) and draws the reference dot plot of niche chemokines/receptors across subsets, with candidate producers (>15 %) highlighted (S9A).

**`02_extract_transcripts.py`** — extracts raw per-cell CXCL13/CCL19/CCL21/TSPAN7 counts from the SpatialData object (feeds `04` and `06`).

**`03_impvl_vs_icaf.py`** — dot plot contrasting TLS-resident imPVL s1 vs iCAF s2 (chemokines, fibroblast/CAF markers, mural/PVL markers) (S9B).

**`04_zone_expression_density.py`** — CP10K CXCL13 & CCL19 expression density per mm² across B cell zone / T cell zone / non-TLS; two-sided paired Wilcoxon (n=6) (S9E).

**`05_chemokine_gradients.py`** — spatial gradients of CCL19/CCL21 transcript density and imPVL s1 / iCAF s2 cell density vs distance from the HEV centre (within the T cell zone) and from the T-zone edge into non-TLS (S9H). `--replot_only` regenerates the figure from the per-sample distance-profile CSVs written by a previous full run, skipping the slow SpatialData extraction — use it for plot-only changes once the extraction has already been done.

**`06_reference_producer_corr.py`** — per-50 µm-bin partial Spearman correlation of candidate-producer cell density vs CCL19/CCL21 transcript density (covariate = total cells per bin), within the T cell zone; per-sample coefficients tested vs 0 by Wilcoxon (n=6) (S9I).

**`07_cxcl13_proximity_attribution.py`** — attributes each CXCL13 transcript to the nearest Tfh/FDC anchor within 15 µm (else "other"; background < 30 CXCL13/cell excluded; sources < 1 % pooled), B vs T zone, cohort-averaged over FDC-positive samples (S9L); plus Tfh & FDC cell density by zone (S9K). FDC-negative samples are auto-skipped.

**Key inputs:** `data/anndata_objects/xenium/{s}.h5ad`; SpatialData `data/spatialdata_objects/xenium/{s}.zarr` (02, 05, 07); scRNA-seq reference `data/anndata_objects/sc_ref/miniatlas_GEO.h5ad` (01); zone bins from `07_TLS_zonation/01_classify_zones`; FDC barcodes from `08_FDC_recovery/01_identify_FDC`; cell labels from `results/02_cell_type_mapping/label_transfer/xenium/subset_level/`.

**Outputs:** `results/05_xenium_mechanisms/09_chemokine_sources/<step>/`.

**Run order:** requires `07_TLS_zonation` and `08_FDC_recovery` first. Then `01` (reference producers) and `02` (per-sample transcript extraction) → `03`–`07` (cohort). `05` also needs `02_identify_HEV_regions`.
