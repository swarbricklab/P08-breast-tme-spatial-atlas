# 01_identify_TLS_regions

Segment **tertiary lymphoid structure (TLS)** regions from dense spatial co-clustering of T and B cells, then characterise their cellular composition and gene expression. `01_identify_TLS.py` is the upstream segmentation step (its TLS regions feed workflows 03, 07, 08, 09 and 11, and are read directly by 02 and 06); `02` and `03` are cohort summaries that reproduce **Extended Data Fig. 7b–d**.

| Script | Env | Mode | Figure |
|---|---|---|---|
| `01_identify_TLS.py` | xen-hev | `--sample_id` | — (upstream) |
| `02_composition_TLS.py` | xen-hev | `--samples` (cohort) | Extended Data Fig. 7b,c |
| `03_gene_expression_TLS.py` | xen-hev | `--samples` (cohort) | Extended Data Fig. 7d |

**`01_identify_TLS.py`** — bins each sample on a 25 µm grid, Gaussian-smooths the T-cell and B-cell (Tangram major-level) densities, and keeps grid regions that are T/B-dominant, mixed (both lineages present), and ≥0.01 mm² as TLS polygons; each cell is then flagged `in_TLS` and given a `TLS_id` by spatial join. Writes the TLS polygons (`{s}_TLS_regions.geojson`), the per-cell TLS annotations (`{s}_cell_TLS_annotations.csv`), and a spatial overview PDF.

**`02_composition_TLS.py`** — cohort cellular composition of the TLS regions at the major / minor / subset annotation levels: the proportion of each cell type located inside vs outside TLS regions (Extended Data Fig. 7b) and the average relative abundance of each cell type among cells within TLS regions (Extended Data Fig. 7c; mean ± SD across samples). Writes `All_samples_TLS_cellular_composition.pdf` and the per-level statistics CSVs.

**`03_gene_expression_TLS.py`** — cohort matrixplot of mean CD40 and CD40LG expression (min-max scaled 0–1 per gene) across the subset-level cell types comprising ≥1% of in-TLS cells, with a per-cell-type count bar on top (Extended Data Fig. 7d). Writes `All_samples_TLS_gene_expression_CD40LG_CD40.pdf`.

**Inputs:** `data/anndata_objects/xenium/{s}.h5ad`; Tangram label transfers `results/02_cell_type_mapping/label_transfer/xenium/{level}_level/{s}_{level}_pred_cell_types.csv` (col `pred_cell_type`; `01` uses the major level `T-cells`/`B-cells`, `02`/`03` use major/minor/subset); the per-cell TLS annotations from `01` (`{s}_cell_TLS_annotations.csv`).

**Outputs:** `results/05_xenium_mechanisms/01_identify_TLS_regions/…` — per-sample `csv_files/` + `plots/` from `01`; cohort summaries under `all_samples/{cellular_composition, gene_expression}/`.

**Run order:** run `01_identify_TLS.py` per sample first; `02` and `03` then run once over the cohort (`--samples`), reading `01`'s `{s}_cell_TLS_annotations.csv`. Those annotations are also consumed by `02_identify_HEV_regions`, `03_HEV_microenvironment`, `06_cellular_enrichment` and the downstream zonation/SPP1 workflows (07–09, 11).
