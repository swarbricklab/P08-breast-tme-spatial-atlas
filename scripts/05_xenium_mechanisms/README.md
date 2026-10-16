# Xenium Mechanisms Analysis Workflow

This directory contains scripts for the high-resolution, single-cell spatial analysis of tertiary lymphoid structures (TLS), vascular niches, and microenvironmental signaling mechanisms using Xenium in situ profiling data.

## Directory Structure

- **01_identify_TLS_regions/** 
  Single-cell spatial identification of TLS architectures, plus their cellular composition and gene expression. *(Extended Data Fig. 7b–d)*
  - `01_identify_TLS.py`: Automatically segment TLS regions based on the density and spatial co-localization of T cells and B cells.
  - `02_composition_TLS.py`: Cohort cellular composition of TLS regions — cell-type proportions inside vs outside TLS (Extended Data Fig. 7b) and average relative abundance within TLS (Extended Data Fig. 7c).
  - `03_gene_expression_TLS.py`: Cohort CD40 / CD40LG expression across the subset-level cell types within TLS regions (Extended Data Fig. 7d).

- **02_identify_HEV_regions/** 
  Single-cell identification of high endothelial venules.
  - `01_identify_HEV.py`: Score and segment HEV structures at single-cell resolution using a targeted 5-gene endothelial signature.

- **03_HEV_microenvironment/** 
  Spatial profiling of the HEV microenvironment.
  - `01_calculate_distance_HEV.py`: Calculate radial distances from every cell to the nearest segmented HEV for downstream proximity-based enrichment analyses.

- **04_score_signatures/** 
  Scoring of immune signaling pathways. *(Relates to Supplementary Table 9)*
  - `01_score_STING_signatures.py`: Score STING-pathway and Type I Interferon (IFN) gene signatures across specific cell types and map their spatial activity in the vascular niche.

- **05_gene_gene_correlations/** 
  Gene correlation and pathway enrichment analysis. *(Generates data for Supplementary Table 10)*
  - `01_investigate_gene_correlations.py`: Calculate pseudobulk expression profiles per HEV structure to investigate transcriptomic correlations with endothelial `TMEM173` (STING) and `SPP1` expression.
  - `02_GSEA_analysis_on_correlated_genes.py`: Perform Gene Set Enrichment Analysis (GSEA) on the correlated gene lists to identify functional pathways (e.g., Cell Adhesion vs. Chemokine Activity).

- **06_cellular_enrichment/** 
  Spatial recruitment and signaling analysis. *(Generates data for Supplementary Table 12)*
  - `01_SPP1_vs_proximal_enrichment.py`: Quantify the spatial enrichment of specific cell types (e.g., Tfh cells, LAMP3+ DCs) within a 200 µm radius of HEVs, and correlate their abundance with endothelial `SPP1` expression.

- **07_TLS_zonation/** 
  Resolve TLS into B cell / T cell zones and germinal centres (GCs) on a 50 µm grid. *(Extended Data Figs. 9c,d,f,g and 10g,h)*
  - `01_classify_zones.py`: Score per-bin pseudobulk profiles against B/T-cell-zone and GC niche signatures; classify bins and label discrete GCs.
  - `02_zonation_overview.py`: Cohort B/T-zone map (Extended Data Fig. 9c), per-bin T:B ratio distribution (Extended Data Fig. 9d), GC overview (Extended Data Fig. 10g) and discrete-GCs-per-sample bar (Extended Data Fig. 10h).
  - `03_zoom_GCs.py`: Per-discrete-GC aligned-H&E zoom with the top GC-signature transcripts (Extended Data Fig. 10g).
  - `04_zoom_zones.py`: Representative adjacent B/T-cell-zone single-transcript zoom + zone-annotation overview (Extended Data Fig. 9f,g).

- **08_FDC_recovery/** 
  Recover follicular dendritic cells (FDCs) and Tfh from the B-cell / germinal-centre compartment. *(Extended Data Fig. 9j)*
  - `01_identify_FDC.py`: Per-sample subclustering to recover the FDC and Tfh populations.
  - `02_FDC_overview.py`: Cohort FDC / Tfh subcluster marker dot plot (Extended Data Fig. 9j).

- **09_chemokine_sources/** 
  Cellular and spatial sources of the TLS-organising chemokines (CXCL13, CCL19, CCL21). *(Extended Data Fig. 9a,b,e,h,i,k,l)*
  - `01_reference_producers.py`: Reference (scRNA-seq) chemokine-producer subsets and niche-gene dot plot (Extended Data Fig. 9a).
  - `02_extract_transcripts.py`: Per-cell chemokine transcript extraction from SpatialData (upstream input).
  - `03_impvl_vs_icaf.py`: imPVL s1 vs iCAF s2 discrimination dot plot (Extended Data Fig. 9b).
  - `04_zone_expression_density.py`: CXCL13 / CCL19 expression density per zone (Extended Data Fig. 9e).
  - `05_chemokine_gradients.py`: CCL19/CCL21 transcript and imPVL/iCAF cell density vs distance from HEV (Extended Data Fig. 9h).
  - `06_reference_producer_corr.py`: Reference-producer abundance vs chemokine density correlation per 50 µm bin (Extended Data Fig. 9i).
  - `07_cxcl13_proximity_attribution.py`: Tfh/FDC density per zone (Extended Data Fig. 9k) and CXCL13 reference-producer source attribution (Extended Data Fig. 9l).

- **10_HEV_GC_proximity/** 
  Stratify HEV endothelial cells by proximity to a germinal centre and test GC-proximal genes (notably `SPP1`). *(Extended Data Fig. 10i,j)*
  - `01_classify_HEV_by_GC_dist.py`: Label intra-TLS HEV niches close / far from the nearest GC.
  - `02_DEG_close_vs_far.py`: Cell-level close-vs-far Mann-Whitney DEG table for ACKR1+ HEV endothelial cells.
  - `03_volcano_dotplot.py`: Close-vs-far DEG volcano (Extended Data Fig. 10i) and transposed dot plot of the significant genes (Extended Data Fig. 10j).

- **11_SPP1_venous_specificity/** 
  Endothelial and macrophage `SPP1` in the TLS / HEV context. *(Extended Data Fig. 10c,e)*
  - `01_SPP1_across_venous_groups.py`: `SPP1` dot plot across four ACKR1+ venous EC groups (HEV / non-HEV × TLS / stroma) (Extended Data Fig. 10c).
  - `02_SPP1_macrophage_regions.py`: Representative-sample region zoom of macrophage `SPP1` outside TLS (Extended Data Fig. 10e).

Each workflow subdirectory has its own README with the per-script conda environment, run mode (`--sample_id` per-sample vs `--samples` cohort), inputs, and the figure panel(s) it produces.

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
