# 05_xenium_mechanisms

Single-cell analysis of TLS, HEVs and their signalling in the Xenium Prime 5K data. Each folder has a README with its scripts, conda environments, run mode (`--sample_id` per sample or `--samples` for the cohort), inputs and outputs. Run folders 01 and 02 first: most later folders use the TLS regions and HEV niches they define. Besides the inputs listed per folder, all scripts read the sample sheet (`config/sample_sheet_xenium.csv`), `resources/tissue_metadata.csv` and the colour and cell type files in `resources/`.

| Folder | What it does | Manuscript |
|---|---|---|
| [01_identify_TLS_regions](01_identify_TLS_regions/README.md) | TLS regions from the local density of T and B cells; their composition and CD40/CD40LG expression | Extended Data Fig. 7b–d |
| [02_identify_HEV_regions](02_identify_HEV_regions/README.md) | HEV niches from a 5-gene HEV signature score | |
| [03_HEV_microenvironment](03_HEV_microenvironment/README.md) | Distance from each cell to the nearest HEV niche | |
| [04_score_signatures](04_score_signatures/README.md) | STING and type I interferon signature scores per cell type | Supplementary Table 9 |
| [05_gene_gene_correlations](05_gene_gene_correlations/README.md) | Genes correlated with endothelial TMEM173 (STING) or SPP1 across HEV niches, and GSEA | Supplementary Table 10 |
| [06_cellular_enrichment](06_cellular_enrichment/README.md) | Cell types within 200 µm of HEV niches and their correlation with endothelial SPP1 | Supplementary Table 12 |
| [07_TLS_zonation](07_TLS_zonation/README.md) | B cell zones, T cell zones and germinal centres within TLS | Extended Data Figs. 9c,d,f,g and 10g,h |
| [08_FDC_recovery](08_FDC_recovery/README.md) | Follicular dendritic cells recovered from the label-transferred Tfh cells | Extended Data Fig. 9j |
| [09_chemokine_sources](09_chemokine_sources/README.md) | Cellular sources and spatial gradients of CXCL13, CCL19 and CCL21 | Extended Data Fig. 9a,b,e,h,i,k,l |
| [10_HEV_GC_proximity](10_HEV_GC_proximity/README.md) | HEV gene expression close to vs far from germinal centres | Extended Data Fig. 10i,j |
| [11_SPP1_venous_specificity](11_SPP1_venous_specificity/README.md) | SPP1 in venous endothelial cells and macrophages, inside and outside TLS | Extended Data Fig. 10c,e |

`plotting_functions.py` holds plotting functions shared by these scripts.
