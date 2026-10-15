# 03_HEV_microenvironment

Compute each cell's **radial distance to the nearest segmented HEV niche** (overall, and split by HEVs close to / far from TLS). Upstream workflow (no figure of its own); provides the distance tables underpinning the proximity-based analyses in 04-06.

| Script | Env | Mode | Figure |
|---|---|---|---|
| `01_calculate_distance_HEV.py` | xen-hev | `--sample_id` | — (upstream) |

**`01_calculate_distance_HEV.py`** — loads the HEV niches from `02_identify_HEV_regions` and the TLS annotations, then for every cell computes the KD-tree distance to the nearest HEV niche, and separately to the nearest TLS-proximal (`close_to_TLS`) and TLS-distal HEV. Writes per-cell distance tables (`{s}_closest_HEV_distances.csv`, `{s}_HEV_distances_close_TLS.csv`, `{s}_HEV_distances_far_from_TLS.csv`) and a distance overview PDF.

**Inputs:** `data/anndata_objects/xenium/{s}.h5ad`; HEV niches from `results/05_xenium_mechanisms/02_identify_HEV_regions/csv_files/` (`{s}_HEV_niches_cell_to_niche.csv`, `{s}_HEV_niches_info.csv`); TLS annotations `results/05_xenium_mechanisms/01_identify_TLS_regions/csv_files/{s}_cell_TLS_annotations.csv`; `config/sample_sheet_xenium.csv`; `resources/tissue_metadata.csv`.

**Outputs:** `results/05_xenium_mechanisms/03_HEV_microenvironment/01_calculate_distance/…` (`csv_files/`, `plots/`).

**Run order:** requires `01_identify_TLS_regions` and `02_identify_HEV_regions`; run per sample. Upstream of the proximity-based analyses (04-06).
