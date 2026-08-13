# 02_identify_HEV_regions

Score and segment **high endothelial venule (HEV)** structures at single-cell resolution using a targeted 5-gene endothelial signature, producing per-sample HEV niches. Upstream workflow (no figure of its own); its HEV niches feed workflows 03, 04, 05, 06, 10 and 11.

| Script | Env | Mode | Figure |
|---|---|---|---|
| `01_identify_HEV.py` | xen-hev | `--sample_id` | — (upstream) |

**`01_identify_HEV.py`** — scores each cell for a 5-gene HEV signature (`TSPAN7`, `MEOX2`, `RAMP3`, `SPNS2`, `IL33`; threshold 0.6), grows spatial HEV niches from the high-scoring cores by KD-tree neighbour expansion and merging, then computes each niche's distance to the nearest TLS cell and flags `close_to_TLS` (≤100 µm). Writes per-cell niche membership (`{s}_HEV_niches_cell_to_niche.csv`), per-niche info including TLS distance (`{s}_HEV_niches_info.csv`), a parameters JSON, and an overview PDF — all under a dated `{YYYYMMDD}/` subdirectory.

**Inputs:** `data/anndata_objects/xenium/{s}.h5ad`; Tangram major- and subset-level labels `results/02_cell_type_mapping/label_transfer/xenium/{major,subset}_level/{s}_{level}_pred_cell_types.csv`; TLS annotations `results/05_xenium_mechanisms/01_identify_TLS_regions/csv_files/{s}_cell_TLS_annotations.csv`; `config/sample_sheet_xenium.csv`; `resources/tissue_metadata.csv`; colours `resources/plotting_colors/miniatlas_major_colors.csv` and `resources/cell_types/miniatlas/miniatlas_subset_classifications.csv`.

**Outputs:** `results/05_xenium_mechanisms/02_identify_HEV_regions/{YYYYMMDD}/…` (`csv_files/`, `plots/`).

**Run order:** requires `01_identify_TLS_regions`; run per sample. Downstream scripts (03-06, 10, 11) resolve the most-recent dated `02_identify_HEV_regions/{YYYYMMDD}` subdirectory.
