# 06_cellular_enrichment

Quantify **cell-type enrichment within ~200 µm of HEV niches** and correlate proximal abundance with endothelial `SPP1` expression. Generates data for Table S12.

| Script | Env | Mode | Figure |
|---|---|---|---|
| `01_SPP1_vs_proximal_enrichment.py` | xen-hev | `--sample_id` | Table S12 |

**`01_SPP1_vs_proximal_enrichment.py`** — loads the HEV niches from `02_identify_HEV_regions`, TLS annotations and Tangram subset labels; for each HEV niche finds all non-niche cells within 200 µm (KD-tree minimum distance), tabulates their cell-type proportions, computes mean ACKR1+ endothelial `SPP1` per niche, and Pearson-correlates each proximal cell-type proportion against niche SPP1. Writes per-niche SPP1 + cell proportions (`{s}_HEV_SPP1_and_cell_proportions.csv`) and the correlation table (`{s}_HEV_SPP1_and_cell_proportions_corr.csv`), plus a PDF, under `maxdist_200um/`.

**Inputs:** `data/anndata_objects/xenium/{s}.h5ad`; HEV niches from `results/05_xenium_mechanisms/02_identify_HEV_regions/csv_files/` (`{s}_HEV_niches_cell_to_niche.csv`, `{s}_HEV_niches_info.csv`); TLS annotations `results/05_xenium_mechanisms/01_identify_TLS_regions/csv_files/{s}_cell_TLS_annotations.csv`; Tangram subset labels `results/02_cell_type_mapping/label_transfer/xenium/subset_level/{s}_subset_pred_cell_types.csv`; colours `resources/plotting_colors/{miniatlas_major_colors.csv,TLS_HEV_colors.csv}` and `resources/cell_types/miniatlas/miniatlas_subset_classifications.csv`; `config/sample_sheet_xenium.csv`; `resources/tissue_metadata.csv`.

**Outputs:** `results/05_xenium_mechanisms/06_cellular_enrichment/SPP1_vs_cell_abundance/maxdist_200um/…` (`csv_files/`, `plots/`).

**Run order:** requires `01_identify_TLS_regions` and `02_identify_HEV_regions`; run per sample. Feeds Table S12.
