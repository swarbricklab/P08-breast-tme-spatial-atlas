# 05_gene_gene_correlations

Per-HEV-niche **pseudobulk correlation** of genes with endothelial `TMEM173` (STING) and `SPP1`, followed by **GSEA** on the correlated-gene lists. Together generate data for Table S10.

| Script | Env | Mode | Figure |
|---|---|---|---|
| `01_investigate_gene_correlations.py` | xen-hev | `--samples` | Table S10 |
| `02_GSEA_analysis_on_correlated_genes.py` | gseapy | `--samples` | Table S10 |

**`01_investigate_gene_correlations.py`** — pools all `--samples` HEV-niche cells (from `02_identify_HEV_regions`), builds a per-niche pseudobulk (mean-expression) matrix, and Pearson-correlates every gene against the gene of interest across niches, flagging p<0.05. Writes the ranked correlation table (`All_samples_HEV_{gene}_gene_correlations.csv`) and a rank-vs-correlation plot. The gene is set with `--gene_of_interest` (default `TMEM173`); the manuscript uses `TMEM173` (Fig. 4d) and `SPP1` (Extended Data Fig. 10d).

**`02_GSEA_analysis_on_correlated_genes.py`** — reads the `01` correlation table, filters to genes present in the panel, and runs `gseapy.prerank` GSEA against a GO gene set (`--gene_set`, default `GO_Molecular_Function_2025`; `--gene_of_interest` selects the `01` table, default `TMEM173`) using the correlation coefficient as the rank metric. Writes the full preranked GSEA results (`All_samples_HEV_{gene}_preranked_GSEA_results.csv`) and a top-terms barplot. Runs on the `gseapy` env. Manuscript runs (Supplementary Table 10): `TMEM173` × GO MF (Fig. 4d), `SPP1` × GO MF (Extended Data Fig. 10d) and `SPP1` × GO BP (Extended Data Fig. 10f).

**Inputs:** `01` — `data/anndata_objects/xenium/{s}.h5ad` for all `--samples`; HEV niches from `results/05_xenium_mechanisms/02_identify_HEV_regions/csv_files/` (`{s}_HEV_niches_cell_to_niche.csv`, `{s}_HEV_niches_info.csv`); `config/sample_sheet_xenium.csv`; `resources/tissue_metadata.csv`. `02` — the `01_investigate_gene_correlations/{gene}_correlations/csv_files/` correlation CSV; GO `.gmt` gene sets `resources/gene_lists/enrichment_analysis/`.

**Outputs:** `results/05_xenium_mechanisms/05_gene_gene_correlations/…` (`01_investigate_gene_correlations/{gene}_correlations/`, `02_GSEA_analysis/{gene}_correlations/{gene_set}/`).

**Run order:** both require `02_identify_HEV_regions` (single cohort run over all samples). `01` → `02` (which reads the `01` correlation table).
