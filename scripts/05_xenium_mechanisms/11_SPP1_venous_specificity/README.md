# 11_SPP1_venous_specificity

Endothelial and macrophage **SPP1** in the TLS/HEV context: whether endothelial SPP1 tracks HEV identity rather than TLS location (S10C), and where macrophage SPP1 sits outside TLS (S10E).

| Script | Env | Mode | Figure |
|---|---|---|---|
| `01_SPP1_across_venous_groups.py` | xen-hev | `--samples` | S10C |
| `02_SPP1_macrophage_regions.py` | spatdat | `--sample_id` | S10E |

**`01_SPP1_across_venous_groups.py`** — splits ACKR1+ venous ECs into Non-HEV venous (stroma) / HEV-like (stroma) / Non-HEV venous (in TLS) / HEV (in TLS) and draws a per-gene-max dot plot of SPP1 + venous/HEV markers (ACKR1, TSPAN7, RAMP3, MEOX2) across the four groups, with per-group cell-count bars (S10C). HEV vs non-HEV membership uses the HEV niche `close_to_TLS` flag; the in-TLS split uses per-cell `in_TLS`.

**`02_SPP1_macrophage_regions.py`** — single-sample zoom of macrophage SPP1 outside TLS in the representative sample (P29_Xe5k): four tissue regions, each a column of a 3-row grid (aligned H&E / H&E + SPP1 transcripts + LAM1-macrophage `Myeloid_c1_LAM1_FABP5` cell boundaries / boundaries + transcripts without H&E), plus a whole-tissue overview with the four zoom boxes marked (S10E). The four regions are defined in `resources/areas_of_interest/xenium/P29_Xe5k_SPP1_macrophages.geojson` (micron polygons mapped onto the aligned H&E via the SpatialData `image_space` affine, same convention as `07_TLS_zonation/04_zoom_zones`). Loads the full-res H&E, so it needs a high-memory node; only runs for samples with an annotated geojson (currently P29_Xe5k).

**Inputs:** `data/anndata_objects/xenium/{s}.h5ad`; cell labels `results/02_cell_type_mapping/label_transfer/xenium/subset_level/` (labels `Endothelial ACKR1` for `01`, `Myeloid_c1_LAM1_FABP5` for `02`); TLS annotations from `01_identify_TLS_regions` and HEV niches from `02_identify_HEV_regions` (`01`); and, for `02`, the aligned H&E (`data/images/he/aligned_images/xenium/{s}.ome.tif`), SpatialData (`data/spatialdata_objects/xenium/{s}.zarr`, transcripts + `cell_boundaries`), transcript colours, and the SPP1-macrophage geojson.

**Outputs:** `results/05_xenium_mechanisms/11_SPP1_venous_specificity/`.

**Run order:** `01` requires `01_identify_TLS_regions` and `02_identify_HEV_regions` (single cohort run); `02` is an independent single representative-sample run.
