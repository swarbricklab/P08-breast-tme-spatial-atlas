# 07_TLS_zonation

Resolve segmented TLS regions into **B cell / T cell zones** and **germinal centres (GCs)** on a 50 µm grid, by scoring per-bin pseudobulk profiles against the Visium-derived niche signatures (`resources/gene_lists/TLS_niche_signatures/niche_signatures.csv`). Upstream of the chemokine-source (09), GC-proximity (10) and SPP1 (11) analyses.

| Script | Env | Mode | Figure |
|---|---|---|---|
| `01_classify_zones.py` | xen-hev | `--sample_id` | — (upstream) |
| `02_zonation_overview.py` | xen-hev | `--samples` | S9C, S9D, S10G/H |
| `03_zoom_GCs.py` | spatdat | `--sample_id` | S10G |
| `04_zoom_zones.py` | spatdat | `--sample_id` | S9F, S9G |

**`01_classify_zones.py`** — bins cells to a 50 µm grid inside the TLS regions and scores each bin for the B cell zone / T cell zone / GC niche signatures (top-20 on-panel genes). Classifies bins as B cell zone / T cell zone / mixed (score > 0.3) and GC (score > 0.8), and labels discrete GCs (8-connectivity, merging single-bin gaps). Writes per sample: `{s}_BT_zone_bins.csv`, `{s}_BT_zone_classification.csv`, `{s}_GC_zone_bins.csv`, `{s}_GC_zone_classification.csv`.

**`02_zonation_overview.py`** — cohort overview: the B/T zone spatial map (S9C), per-bin T:B ratio distribution (S9D), GC spatial overview (S10G), and the number of discrete GCs per sample coloured by clinical subtype (S10H).

**`03_zoom_GCs.py`** — per discrete GC, an aligned-H&E panel with the GC outline plus the top-5 GC-signature transcripts overlaid (S10G). Reads only the GC crops from the full-res H&E (tile by tile) and the overview from a low-resolution pyramid level (`read_tiff_region`/`read_tiff_overview` in `scripts/helper_functions.py`), so it needs < 30 GB.

**`04_zoom_zones.py`** — single-transcript zoom of adjacent B cell zone and T cell zone fields in the representative sample (P22_Xe5k). One page: an overview (S9F) of the aligned H&E cropped around the two fields with the 50 µm zone bins overlaid and the two fields boxed, and a 2×7 grid (S9G) with rows = B cell zone / T cell zone and columns = H&E, B/T cell markers (MS4A1, CD3E), CXCL13, CCL19, CCL21, FDC markers (CR1, CR2), HEV markers (TSPAN7, MEOX2). The two fields are defined in `resources/areas_of_interest/xenium/P22_Xe5k_zone_zoom.geojson`. Reads only the overview and zoom crops from the full-res H&E (tile by tile, `read_tiff_region` in `scripts/helper_functions.py`), so it needs < 20 GB.

**Zone label values** (used downstream): `TLS_zone ∈ {B_cell_zone, T_cell_zone, mixed}`; `GC_zone = GC`.

**Inputs:** `data/anndata_objects/xenium/{s}.h5ad` (layer `filtered_counts`); TLS regions from `results/05_xenium_mechanisms/01_identify_TLS_regions/csv_files/`; niche signatures; and, for `03`/`04`, the aligned H&E (`data/images/he/aligned_images/xenium/{s}.ome.tif`) + SpatialData (`data/spatialdata_objects/xenium/{s}.zarr`). `04` also reads the zone-zoom AOI geojson (`resources/areas_of_interest/xenium/{s}_zone_zoom.geojson`) and transcript colours (`resources/plotting_colors/transcript_colors.csv`); it only runs for samples that have an annotated zone-zoom geojson (currently P22_Xe5k).

**Outputs:** `results/05_xenium_mechanisms/07_TLS_zonation/…`.

**Run order:** `01` (per sample) → `02` (cohort) / `03` (per sample) / `04` (representative sample). Run `01`–`03` before workflows 09 and 10.
