# 05_gene_gene_correlations

Per-HEV-niche **pseudobulk correlation** of genes with endothelial `TMEM173` (STING) and `SPP1`, followed by **GSEA** on the correlated-gene lists. Together generate data for Table S10.

| Script | Env | Mode | Figure |
|---|---|---|---|
| `01_investigate_gene_correlations.py` | xen-hev | `--samples` | Table S10 |
| `02_GSEA_analysis_on_correlated_genes.py` | gseapy | `--samples` | Table S10 |

**`01_investigate_gene_correlations.py`** — pools all `--samples` HEV-niche cells (from the latest `02_identify_HEV_regions/{date}`), builds a per-niche pseudobulk (mean-expression) matrix, and Pearson-correlates every gene against the gene of interest across niches, flagging p<0.05. Writes the ranked correlation table (`All_samples_HEV_{gene}_gene_correlations.csv`) and a rank-vs-correlation plot under a dated subdirectory. The committed script sets `gene_of_i = TMEM173`; switch it to `SPP1` for the SPP1 correlations.

**`02_GSEA_analysis_on_correlated_genes.py`** — reads the latest `01` correlation table, filters to genes present in the panel, and runs `gseapy.prerank` GSEA against a GO gene set (default `GO_Molecular_Function_2025`) using the correlation coefficient as the rank metric. Writes the full preranked GSEA results (`All_samples_HEV_{gene}_preranked_GSEA_results.csv`) and a top-terms barplot under dated subdirectories. Runs on the `gseapy` env.

**Inputs:** `01` — `data/anndata_objects/xenium/{s}.h5ad` for all `--samples`; HEV niches from the latest `results/05_xenium_mechanisms/02_identify_HEV_regions/{YYYYMMDD}/csv_files/` (`{s}_HEV_niches_cell_to_niche.csv`, `{s}_HEV_niches_info.csv`); `config/sample_sheet_xenium.csv`; `resources/tissue_metadata.csv`. `02` — the latest `01_investigate_gene_correlations/{gene}_correlations/{date}/csv_files/` correlation CSV; GO `.gmt` gene sets `resources/gene_lists/enrichment_analysis/`.

**Outputs:** `results/05_xenium_mechanisms/05_gene_gene_correlations/…` (`01_investigate_gene_correlations/{gene}_correlations/{date}/`, `02_GSEA_analysis/{gene}_correlations/{gene_set}/{date}/`).

**Run order:** both require `02_identify_HEV_regions` (single cohort run over all samples). `01` → `02` (which reads the latest `01` correlation table).
