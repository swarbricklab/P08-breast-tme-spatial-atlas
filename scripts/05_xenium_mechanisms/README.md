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

- **07_TLS_zonation/** 
  Resolve TLS into B cell / T cell zones and germinal centres (GCs) on a 50 µm grid. *(Figures S9C–D, S9F–G, S10G–H)*
  - `01_classify_zones.py`: Score per-bin pseudobulk profiles against B/T-cell-zone and GC niche signatures; classify bins and label discrete GCs.
  - `02_zonation_overview.py`: Cohort B/T-zone map (S9C), per-bin T:B ratio distribution (S9D), GC overview (S10G) and discrete-GCs-per-sample bar (S10H).
  - `03_zoom_GCs.py`: Per-discrete-GC aligned-H&E zoom with the top GC-signature transcripts (S10G).
  - `04_zoom_zones.py`: Representative adjacent B/T-cell-zone single-transcript zoom + zone-annotation overview (S9F–G).

- **08_FDC_recovery/** 
  Recover follicular dendritic cells (FDCs) and Tfh from the B-cell / germinal-centre compartment. *(Figure S9J)*
  - `01_identify_FDC.py`: Per-sample subclustering to recover the FDC and Tfh populations.
  - `02_FDC_overview.py`: Cohort FDC / Tfh subcluster marker dot plot (S9J).

- **09_chemokine_sources/** 
  Cellular and spatial sources of the TLS-organising chemokines (CXCL13, CCL19, CCL21). *(Figures S9A–B, S9E, S9H, S9I, S9K–L)*
  - `01_reference_producers.py`: Reference (scRNA-seq) chemokine-producer subsets and niche-gene dot plot (S9A).
  - `02_extract_transcripts.py`: Per-cell chemokine transcript extraction from SpatialData (upstream input).
  - `03_impvl_vs_icaf.py`: imPVL s1 vs iCAF s2 discrimination dot plot (S9B).
  - `04_zone_expression_density.py`: CXCL13 / CCL19 expression density per zone (S9E).
  - `05_chemokine_gradients.py`: CCL19/CCL21 transcript and imPVL/iCAF cell density vs distance from HEV (S9H).
  - `06_reference_producer_corr.py`: Reference-producer abundance vs chemokine density correlation per 50 µm bin (S9I).
  - `07_cxcl13_proximity_attribution.py`: Tfh/FDC density per zone (S9K) and CXCL13 reference-producer source attribution (S9L).

- **10_HEV_GC_proximity/** 
  Stratify HEV endothelial cells by proximity to a germinal centre and test GC-proximal genes (notably `SPP1`). *(Figures S10I–J)*
  - `01_classify_HEV_by_GC_dist.py`: Label intra-TLS HEV niches close / far from the nearest GC.
  - `02_DEG_close_vs_far.py`: Cell-level close-vs-far Mann-Whitney DEG table for ACKR1+ HEV endothelial cells.
  - `03_volcano_dotplot.py`: Close-vs-far DEG volcano (S10I) and transposed dot plot of the significant genes (S10J).

- **11_SPP1_venous_specificity/** 
  Endothelial and macrophage `SPP1` in the TLS / HEV context. *(Figures S10C, S10E)*
  - `01_SPP1_across_venous_groups.py`: `SPP1` dot plot across four ACKR1+ venous EC groups (HEV / non-HEV × TLS / stroma) (S10C).
  - `02_SPP1_macrophage_regions.py`: Representative-sample region zoom of macrophage `SPP1` outside TLS (S10E).

Each workflow subdirectory has its own README with the per-script conda environment, run mode (`--sample_id` per-sample vs `--samples` cohort), inputs, and the exact supplementary figure panel(s) it produces.

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

6. **TLS Zonation & Source Attribution:** 
  Resolve TLS into B/T cell zones and germinal centres (07), recover FDC and Tfh populations (08), attribute the TLS-organising chemokines to their cellular sources and map their spatial gradients out from the HEVs (09), stratify HEV endothelial cells by germinal-centre proximity (10), and dissect endothelial versus macrophage `SPP1` across the HEV / TLS compartments (11).
