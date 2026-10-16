# 04_vascular_microenvironment

Map HEVs and the wider vasculature in the Visium HD data (16 µm bins), compare HEVs with other venous vessels, and analyse ligand-receptor interactions in the HEV niche. Scripts use the `vis-hd` environment unless noted, and run in the order listed.

- **01_identify_HEV_niches/** — segment HEVs and place them relative to TLS.
  - `01_score_signature.py`: score bins for the HEV and TLS zone signatures (published signatures by default; see [03_TLS_niche_discovery](../03_TLS_niche_discovery/README.md)).
  - `02a_group_HEV_bins.py`, `02b_expand_HEV_groups.py`: group HEV-high "seed" bins, then expand each group to neighbouring bins with lower HEV scores.
  - `03a_identify_TLS_areas.py` (`spatial-pilot`): group bins with high TLS zone scores into TLS areas.
  - `03b_calculate_HEV_TLS_distance.py`: distance from each HEV to the nearest TLS area, classifying HEVs as TLS-associated or stromal HEV-like.
  - `04_assess_HEV_pathology_annotations.py`: HEV locations against the pathology annotations.
  - `05_calculate_all_HEV_distances.py`: distance from each HEV to the surrounding bins.
- **02_identify_all_vasculature/** — classify all vessels.
  - `01_classify_vessels.py`: segment vessels and classify them by dominant deconvolved endothelial type (arterial, capillary, venous, lymphatic).
  - `02_investigate_vessels.py`: overlap of vessel types with HEVs and perivascular-like (PVL) subsets.
  - `03_calculate_distances.py`, `04_assess_PVL_association.py`: PVL composition across distance from each vessel type.
- **03_compare_HEV_to_other_venous/** — HEVs vs other venous vessels.
  - `01_compare_microenvironment.py`: surrounding tissue composition (pathology annotations).
  - `02_compare_gene_expression.py`: differential expression (Supplementary Table 8, Fig. 3i).
  - `03_compare_PVL_association.py`: PVL composition around each vessel group.
- **04_ligand_receptor_analysis/** — stLearn ligand-receptor analysis around HEVs (Supplementary Table 11). `01`–`03` take `--bin_size`: `0` for HEV bins only, `25` to add the perivascular bins within 25 µm.
  - `01_run_LR_analysis_stlearn.py`: run stLearn per sample (the longest step, up to ~4 h per sample).
  - `02_analyze_LR_all.py` (`shapely`): aggregate the interactions across samples, by direction (e.g. perivascular to HEV) and TLS context.
  - `03_plot_LR_per_sample.py`: per-sample spatial plots and dot plots by distance to HEV.
  - `04_check_composition_LTB_signaling_bins.py`: cell type composition of the bins sending LTB to HEV LTBR.

`plotting_functions.py` holds plotting functions shared by these scripts.
