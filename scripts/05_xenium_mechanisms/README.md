# Xenium Mechanisms Analysis Workflow

This directory contains scripts for the high-resolution, single-cell spatial analysis of tertiary lymphoid structures (TLS), vascular niches, and microenvironmental signaling mechanisms using Xenium in situ profiling data.

## Directory Structure

- **01_identify_TLS_regions/** 
  Single-cell spatial identification of TLS architectures.
  - `01_identify_TLS.py`: Automatically segment TLS regions based on the density and spatial co-localization of T cells and B cells.

- **02_identify_HEV_regions/** 
  Single-cell identification of high endothelial venules.
  - `01_identify_HEV.py`: Score and segment HEV structures at single-cell resolution using a targeted 5-gene endothelial signature.

- **03_HEV_microenvironment/** 
  Spatial profiling of the HEV microenvironment.
  - `01_calculate_distance_HEV.py`: Calculate radial distances from every cell to the nearest segmented HEV for downstream proximity-based enrichment analyses.

- **04_score_signatures/** 
  Scoring of immune signaling pathways. *(Relates to Table S9)*
  - `01_score_STING_signatures.py`: Score STING-pathway and Type I Interferon (IFN) gene signatures across specific cell types and map their spatial activity in the vascular niche.

- **05_gene_gene_correlations/** 
  Gene correlation and pathway enrichment analysis. *(Generates data for Table S10)*
  - `01_investigate_gene_correlations.py`: Calculate pseudobulk expression profiles per HEV structure to investigate transcriptomic correlations with endothelial `TMEM173` (STING) and `SPP1` expression.
  - `02_GSEA_analysis_on_correlated_genes.py`: Perform Gene Set Enrichment Analysis (GSEA) on the correlated gene lists to identify functional pathways (e.g., Cell Adhesion vs. Chemokine Activity).

- **06_cellular_enrichment/** 
  Spatial recruitment and signaling analysis. *(Generates data for Table S12)*
  - `01_SPP1_vs_proximal_enrichment.py`: Quantify the spatial enrichment of specific cell types (e.g., Tfh cells, LAMP3+ DCs) within a 200 µm radius of HEVs, and correlate their abundance with endothelial `SPP1` expression.

## Workflow Overview

1. **TLS & HEV Segmentation:** 
  Identify and map the boundaries of TLS and HEV regions based on single-cell coordinates, cellular density, and transcriptomic signatures.

2. **Microenvironment Profiling:** 
  Calculate continuous spatial distances to HEVs to analyze cellular recruitment gradients and proximity effects.

3. **Signature Scoring & Correlation:** 
  Score functional pathways (STING/IFN) within specific cellular compartments and extract gene-gene correlation networks operating within the HEV niche.

4. **Pathway Enrichment:** 
  Perform GSEA on genes co-varying with key vascular drivers (`TMEM173` and `SPP1`) to reveal their distinct biological roles.

5. **Spatial Signaling Analysis:** 
  Statistically link endothelial signaling states (e.g., SPP1 production) to the local recruitment and enrichment of target immune populations.
