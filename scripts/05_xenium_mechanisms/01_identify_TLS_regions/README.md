# 01_identify_TLS_regions

Segment **tertiary lymphoid structure (TLS)** regions from dense spatial co-clustering of T and B cells. Upstream workflow (no figure of its own); its TLS regions feed workflows 03, 07, 08, 09 and 11 (and are also read directly by 02 and 06).

| Script | Env | Mode | Figure |
|---|---|---|---|
| `01_identify_TLS.py` | xen-hev | `--sample_id` | — (upstream) |

**`01_identify_TLS.py`** — bins each sample on a 25 µm grid, Gaussian-smooths the T-cell and B-cell (Tangram major-level) densities, and keeps grid regions that are T/B-dominant, mixed (both lineages present), and ≥0.01 mm² as TLS polygons; each cell is then flagged `in_TLS` and given a `TLS_id` by spatial join. Writes the TLS polygons (`{s}_TLS_regions.geojson`), the per-cell TLS annotations (`{s}_cell_TLS_annotations.csv`), and a spatial overview PDF.

**Inputs:** `data/anndata_objects/xenium/{s}.h5ad`; Tangram major-level labels `results/02_cell_type_mapping/label_transfer/xenium/major_level/{s}_major_pred_cell_types.csv` (labels `T-cells`, `B-cells`); `config/sample_sheet_xenium.csv`; `resources/tissue_metadata.csv`; colours `resources/plotting_colors/miniatlas_major_colors.csv`.

**Outputs:** `results/05_xenium_mechanisms/01_identify_TLS_regions/…` (`csv_files/`, `plots/`).

**Run order:** first step of the workflow — run per sample. The resulting `{s}_cell_TLS_annotations.csv` is consumed by `02_identify_HEV_regions`, `03_HEV_microenvironment`, `06_cellular_enrichment` and the downstream zonation/SPP1 workflows (07-09, 11).
