# 04_score_signatures

Score **STING-pathway and type-I-interferon (IFN)** gene signatures across cell types in the vascular niche. Generates data for Table S9.

| Script | Env | Mode | Figure |
|---|---|---|---|
| `01_score_STING_signatures.py` | xen-hev | `--sample_id` | Table S9 |

**`01_score_STING_signatures.py`** — loads Tangram major/subset labels and the HEV niches from `02_identify_HEV_regions`, splits ACKR1+ endothelium into Venous / HEV-far-from-TLS / HEV-close-to-TLS, `sc.tl.score_genes`-scores each STING and type-I-IFN signature, and renders per-cell-type violin plots plus dot plots (raw and 0-1-scaled mean score × % expressing) across the vessel-type groups. Writes per-cell signature scores (`{s}_STING_signature_scores.csv`), the genes used per signature (`{s}_STING_signatures_genes_used.csv`) and a PDF.

**Inputs:** `data/anndata_objects/xenium/{s}.h5ad`; Tangram major- and subset-level labels `results/02_cell_type_mapping/label_transfer/xenium/{major,subset}_level/`; HEV niches from `results/05_xenium_mechanisms/02_identify_HEV_regions/csv_files/` (`{s}_HEV_niches_cell_to_niche.csv`, `{s}_HEV_niches_info.csv`); signature gene lists `resources/gene_lists/STING/{STING_signatures.csv,type_I_interferon_signatures.csv}`; colours `resources/plotting_colors/{miniatlas_major_colors.csv,endothelial_colors.csv,TLS_HEV_colors.csv}` and `resources/cell_types/miniatlas/miniatlas_subset_classifications.csv`; `config/sample_sheet_xenium.csv`; `resources/tissue_metadata.csv`.

**Outputs:** `results/05_xenium_mechanisms/04_score_signatures/01_score_STING_signatures/…` (`plots/`, `csv_files/`).

**Run order:** requires `02_identify_HEV_regions`; run per sample. Feeds Table S9.
