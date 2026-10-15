# Rules for the figure panels of Figures 1-3 and Extended Data Figure 1 (scripts/figures/figure1-3 and
# scripts/figures/supplementary_overviews): one rule per script, named fig_<panel>.
#   Figure 1 - cohort, pathology annotations, deconvolution (published inputs only)
#   Figure 2 - TLS model and TLS niches (Visium; inputs from 03_TLS_niche_discovery)
#   Figure 3 - HEVs and vasculature (Visium HD; inputs from 04_vascular_microenvironment)
#   S1C, S1G - pathology annotation overviews (published inputs only)
#
# All figures go to {res_dir}/figures (S1C/S1G: {res_dir}/figures/supplementary_overviews).
# Inputs from outside the workflow: {H5AD}/visium, {H5AD}/visium_hd, data/anndata_objects/sc_ref,
# data/pathology_annotations, data/images/he/aligned_images, results/02_cell_type_mapping, resources/.
# Several scripts plot fixed example sections (hard-coded in the script); their sample lists below
# must match the script.
#
# No f-strings in .smk files (see common.smk); all names here start with F13_ to stay clear of the
# other rule files, which share this namespace.

F13_S = "scripts/figures"
F13_OUT = RES + "/figures"
F13_T3 = RES + "/03_TLS_niche_discovery"
F13_T4 = RES + "/04_vascular_microenvironment"

# Published inputs
F13_METADATA = "resources/tissue_metadata.csv"
F13_CLASSES = "resources/cell_types/miniatlas/miniatlas_subset_classifications.csv"
F13_PATH_CATEGORIES = "resources/pathology_annotations/pathology_categories_per_layer.json"
F13_SUBTYPE_COLORS = "resources/plotting_colors/SCSubtype_colors.csv"
F13_PATH_COLORS = "resources/plotting_colors/path_colors.csv"
F13_ENDO_COLORS = "resources/plotting_colors/endothelial_colors.csv"
F13_LINEAGE_COLORS = "resources/plotting_colors/miniatlas_lineage_colors.csv"
F13_MAJOR_COLORS = "resources/plotting_colors/miniatlas_major_colors.csv"
F13_SUBSET_COLORS = "resources/plotting_colors/miniatlas_subset_colors.csv"
F13_VIS_H5AD = H5AD + "/visium/{sample}.h5ad"
F13_HD_H5AD = H5AD + "/visium_hd/{res}/{sample}.h5ad"
F13_VIS_ANNO = "data/pathology_annotations/visium/{sample}_pathology_annotations.csv"
F13_HD_ANNO = "data/pathology_annotations/visium_hd/{sample}_pathology_annotations.csv"
F13_LN_ANNO_DIR = "data/pathology_annotations/visium_lymph_nodes"
F13_VIS_IMAGE = "data/images/he/aligned_images/visium/{sample}.ome.tif"
F13_HD_IMAGE = "data/images/he/aligned_images/visium_hd/{sample}.ome.tif"
F13_DECON = (RES + "/02_cell_type_mapping/deconvolution/visium/{level}_level/"
             "{sample}_{level}_q05_cell_abundance_w_sf.csv")

# Outputs of 03_TLS_niche_discovery
F13_TLS_MODEL = F13_T3 + "/01_identify_TLS/01_run_TLS_model/{sample}_TLS_model.csv"
F13_TLS_CLUSTER = F13_T3 + "/02_identify_TLS_niches/01_cluster_TLS_high/csv_files"
F13_SCORE = (F13_T3 + "/04_map_niche_signatures/01_score_signatures/"
             "{score}_score/csv_files/{sample}_{score}_score.csv")
F13_SCORE_LN = (F13_T3 + "/04_map_niche_signatures/01_score_signatures/lymph_nodes/"
                "{score}_score/csv_files/{sample}_{score}_score.csv")
F13_GC_GROUPS = F13_T3 + "/04_map_niche_signatures/02_group_high_score_areas/GC_score/csv_files/{sample}_GC_groups.csv"
F13_NICHE_SCORES = ["GC", "B_cell_zone", "T_cell_zone"]

# Outputs of 04_vascular_microenvironment
F13_HEV_SCORE = F13_T4 + "/01_identify_HEV_niches/01_score_signature/HEVs_score/{res}/csv_files/{sample}_HEVs_score.csv"
F13_HEV_PER_TLS = (F13_T4 + "/01_identify_HEV_niches/03b_calculate_HEV_TLS_distance/square_016um/threshold_5/"
                   "csv_files/{sample}_HEV_per_TLS_group.csv")
F13_VESSELS = F13_T4 + "/02_identify_all_vasculature/01_classify_vessels/csv_files"
F13_VESSEL_HEV_OVERLAP = F13_T4 + "/02_identify_all_vasculature/02_investigate_vessels/csv_files/{sample}_vessel_HEV_overlap.csv"
F13_PVL_ALL_VESSELS = (F13_T4 + "/02_identify_all_vasculature/04_assess_PVL_association/csv_files/"
                       "{sample}_PVL_composition_vs_distance.csv.gz")
F13_HEV_VS_VENOUS = F13_T4 + "/03_compare_HEV_to_other_venous"

# Script arguments: multi-sample scripts get the full sample list of their sheet (as
# scripts/run_multi_sample_job.sh does); single-sample scripts get --sample_id.
F13_VIS_ARGS = base_args("visium") + " --samples " + " ".join(VIS_SAMPLES)
F13_LN_ARGS = base_args("visium_lymph_nodes") + " --samples " + " ".join(LN_SAMPLES)
F13_HD_ARGS = base_args("visium_hd") + " --samples " + " ".join(HD_SAMPLES)

# Example sections drawn by the scripts (hard-coded there)
F13_1D_SAMPLES = ["P03_Vis", "P15_Vis", "P20_Vis", "P23_Vis", "P31_Vis"]
F13_2F_SAMPLES = ["P22_Vis", "P24_Vis", "P29_Vis_Rep1"]
F13_2G_SAMPLES = ["P22_Vis", "P24_Vis", "P18_Vis", "P29_Vis_Rep1"]
F13_3C_SAMPLES = ["P26_HD", "P27_HD", "P32_HD"]


# ── Figure 1 ─────────────────────────────────────────────────────────────────

rule fig_1C:
    """1C: Data modalities profiled per patient (number of sections per modality)."""
    input:
        sheets   = list(SHEETS.values()),
        metadata = F13_METADATA,
        colors   = F13_SUBTYPE_COLORS,
    output:
        pdf = F13_OUT + "/Fig_1C.pdf",
    conda:
        "spatial-pilot"
    params:
        args = F13_VIS_ARGS,
    resources:
        mem_mb  = 4000,        # peak 0.2 GB (1 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/fig_1C.log",
    shell:
        "python " + F13_S + "/figure1/1C.py {params.args} > {log} 2>&1"


rule fig_1D:
    """1D: H&E and pathology annotations of one example section per clinical subtype."""
    input:
        h5ad        = expand(F13_VIS_H5AD, sample=F13_1D_SAMPLES),
        annotations = expand(F13_VIS_ANNO, sample=F13_1D_SAMPLES),
        categories  = F13_PATH_CATEGORIES,
        metadata    = F13_METADATA,
        sheet       = SHEETS["visium"],
    output:
        pdf = F13_OUT + "/Fig_1D.pdf",
    conda:
        "spatial-pilot"
    params:
        args = F13_VIS_ARGS,
    resources:
        mem_mb  = 8000,        # peak 5.6 GB (1 jobs)
        runtime = 15,          # peak 2 min
    log:
        LOGW + "/fig_1D.log",
    shell:
        "python " + F13_S + "/figure1/1D.py {params.args} > {log} 2>&1"


rule fig_1E:
    """1E: Composition per section by pathology annotation and by deconvolution (major level)."""
    input:
        annotations    = expand(F13_VIS_ANNO, sample=VIS_SAMPLES),
        decon          = expand(F13_DECON, sample=VIS_SAMPLES, level="major"),
        categories     = F13_PATH_CATEGORIES,
        classes        = F13_CLASSES,
        colors         = F13_SUBTYPE_COLORS,
        lineage_colors = F13_LINEAGE_COLORS,
        major_colors   = F13_MAJOR_COLORS,
        metadata       = F13_METADATA,
        sheet          = SHEETS["visium"],
    output:
        pdf = F13_OUT + "/Fig_1E.pdf",
    conda:
        "spatial-pilot"
    params:
        args = F13_VIS_ARGS,
    resources:
        mem_mb  = 4000,        # peak 1.1 GB (1 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/fig_1E.log",
    shell:
        "python " + F13_S + "/figure1/1E.py {params.args} > {log} 2>&1"


rule fig_1F:
    """1F: T cell vs B cell proportions per section (deconvolution, major level)."""
    input:
        decon    = expand(F13_DECON, sample=VIS_SAMPLES, level="major"),
        classes  = F13_CLASSES,
        colors   = F13_SUBTYPE_COLORS,
        metadata = F13_METADATA,
        sheet    = SHEETS["visium"],
    output:
        pdf = F13_OUT + "/Fig_1F.pdf",
    conda:
        "spatial-pilot"
    params:
        args = F13_VIS_ARGS,
    resources:
        mem_mb  = 4000,        # peak 0.9 GB (1 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/fig_1F.log",
    shell:
        "python " + F13_S + "/figure1/1F.py {params.args} > {log} 2>&1"


rule fig_1G:
    """1G: Deconvolved cell type composition per pathology annotation category."""
    input:
        annotations    = expand(F13_VIS_ANNO, sample=VIS_SAMPLES),
        decon          = expand(F13_DECON, sample=VIS_SAMPLES, level="major"),
        categories     = F13_PATH_CATEGORIES,
        classes        = F13_CLASSES,
        lineage_colors = F13_LINEAGE_COLORS,
        major_colors   = F13_MAJOR_COLORS,
        metadata       = F13_METADATA,
        sheet          = SHEETS["visium"],
    output:
        pdf = F13_OUT + "/Fig_1G.pdf",
    conda:
        "spatial-pilot"
    params:
        args = F13_VIS_ARGS,
    resources:
        mem_mb  = 4000,        # peak 0.9 GB (1 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/fig_1G.log",
    shell:
        "python " + F13_S + "/figure1/1G.py {params.args} > {log} 2>&1"


rule fig_1H:
    """1H: B/T cell ratio in TLS/aggregate-annotated spots per clinical subtype."""
    input:
        annotations = expand(F13_VIS_ANNO, sample=VIS_SAMPLES),
        decon       = expand(F13_DECON, sample=VIS_SAMPLES, level="major"),
        colors      = F13_SUBTYPE_COLORS,
        metadata    = F13_METADATA,
        sheet       = SHEETS["visium"],
    output:
        pdf = F13_OUT + "/Fig_1H.pdf",
    conda:
        "spatial-pilot"
    params:
        args = F13_VIS_ARGS,
    resources:
        mem_mb  = 4000,        # peak 0.6 GB (1 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/fig_1H.log",
    shell:
        "python " + F13_S + "/figure1/1H.py {params.args} > {log} 2>&1"


# ── Figure 2 ─────────────────────────────────────────────────────────────────

rule fig_2A:
    """2A: TLS model score on P24_Vis with two zoomed areas (H&E, pathology annotations)."""
    input:
        h5ad        = expand(F13_VIS_H5AD, sample="P24_Vis"),
        tls_model   = expand(F13_TLS_MODEL, sample="P24_Vis"),
        image       = expand(F13_VIS_IMAGE, sample="P24_Vis"),
        annotations = expand(F13_VIS_ANNO, sample="P24_Vis"),
        categories  = F13_PATH_CATEGORIES,
        metadata    = F13_METADATA,
        sheet       = SHEETS["visium"],
    output:
        pdf = F13_OUT + "/Fig_2A.pdf",
    conda:
        "spatial-pilot"
    params:
        args = base_args("visium") + " --sample_id P24_Vis",
    resources:
        mem_mb  = 12000,       # peak 8.8 GB (1 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/fig_2A.log",
    shell:
        "python " + F13_S + "/figure2/2A.py {params.args} > {log} 2>&1"


rule fig_2B:
    """2B: TLS model score distribution per immune pathology annotation (ridge plot)."""
    input:
        annotations = expand(F13_VIS_ANNO, sample=VIS_SAMPLES),
        tls_model   = expand(F13_TLS_MODEL, sample=VIS_SAMPLES),
        metadata    = F13_METADATA,
        sheet       = SHEETS["visium"],
    output:
        pdf = F13_OUT + "/Fig_2B.pdf",
    conda:
        "spatial-pilot"
    params:
        # TLS-high threshold (0.01, as in 03_TLS_niche_discovery); the script has no default
        args = F13_VIS_ARGS + " --TLS_threshold 0.01",
    resources:
        mem_mb  = 4000,        # peak 0.9 GB (1 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/fig_2B.log",
    shell:
        "python " + F13_S + "/figure2/2B.py {params.args} > {log} 2>&1"


rule fig_2C:
    """2C: Proportion of TLS-high spots (score > 0.01) per clinical subtype."""
    input:
        tls_model = expand(F13_TLS_MODEL, sample=VIS_SAMPLES),
        colors    = F13_SUBTYPE_COLORS,
        metadata  = F13_METADATA,
        sheet     = SHEETS["visium"],
    output:
        pdf = F13_OUT + "/Fig_2C.pdf",
    conda:
        "spatial-pilot"
    params:
        args = F13_VIS_ARGS,
    resources:
        mem_mb  = 4000,        # peak 0.2 GB (1 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/fig_2C.log",
    shell:
        "python " + F13_S + "/figure2/2C.py {params.args} > {log} 2>&1"


rule fig_2D:
    """2D: Heatmap of top marker genes per TLS niche cluster."""
    input:
        h5ad     = F13_TLS_CLUSTER + "/All_samples_TLS_umap.h5ad",
        degs     = F13_TLS_CLUSTER + "/All_samples_leiden_clusters_DEGS.csv",
        metadata = F13_METADATA,
        sheet    = SHEETS["visium"],
    output:
        pdf = F13_OUT + "/Fig_2D.pdf",
    conda:
        "spatial-pilot"
    params:
        args = F13_VIS_ARGS,
    resources:
        mem_mb  = 8000,        # peak 3.7 GB (1 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/fig_2D.log",
    shell:
        "python " + F13_S + "/figure2/2D.py {params.args} > {log} 2>&1"


rule fig_2E:
    """2E: Deconvolved cell type composition (subset level) per TLS niche cluster."""
    input:
        clusters = F13_TLS_CLUSTER + "/All_samples_leiden_clusters.csv",
        degs     = F13_TLS_CLUSTER + "/All_samples_leiden_clusters_DEGS.csv",
        decon    = expand(F13_DECON, sample=VIS_SAMPLES, level="subset"),
        classes  = F13_CLASSES,
        metadata = F13_METADATA,
        sheet    = SHEETS["visium"],
    output:
        pdf = F13_OUT + "/Fig_2E.pdf",
    conda:
        "spatial-pilot"
    params:
        args = F13_VIS_ARGS,
    resources:
        mem_mb  = 4000,        # peak 0.6 GB (1 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/fig_2E.log",
    shell:
        "python " + F13_S + "/figure2/2E.py {params.args} > {log} 2>&1"


rule fig_2F:
    """2F: Three example TLSs: H&E, pathology, TLS model, niche scores, B and T cells."""
    input:
        h5ad        = expand(F13_VIS_H5AD, sample=F13_2F_SAMPLES),
        annotations = expand(F13_VIS_ANNO, sample=F13_2F_SAMPLES),
        image       = expand(F13_VIS_IMAGE, sample=F13_2F_SAMPLES),
        tls_model   = expand(F13_TLS_MODEL, sample=F13_2F_SAMPLES),
        scores      = expand(F13_SCORE, sample=F13_2F_SAMPLES, score=F13_NICHE_SCORES),
        decon       = expand(F13_DECON, sample=F13_2F_SAMPLES, level="major"),
        categories  = F13_PATH_CATEGORIES,
        sheet       = SHEETS["visium"],
    output:
        pdf = F13_OUT + "/Fig_2F.pdf",
    conda:
        "vis-hd"
    params:
        args = F13_VIS_ARGS,
    resources:
        mem_mb  = 28000,       # peak 20.9 GB (1 jobs)
        runtime = 15,          # peak 2 min
    log:
        LOGW + "/fig_2F.log",
    shell:
        "python " + F13_S + "/figure2/2F.py {params.args} > {log} 2>&1"


rule fig_2G:
    """2G: Example germinal centres (GC score > 0.7 groups) on H&E."""
    input:
        h5ad      = expand(F13_VIS_H5AD, sample=VIS_SAMPLES),
        gc_groups = expand(F13_GC_GROUPS, sample=VIS_SAMPLES),
        gc_score  = expand(F13_SCORE, sample=VIS_SAMPLES, score="GC"),
        image     = expand(F13_VIS_IMAGE, sample=F13_2G_SAMPLES),
        metadata  = F13_METADATA,
        sheet     = SHEETS["visium"],
    output:
        pdf = F13_OUT + "/Fig_2G.pdf",
    conda:
        "spatial-pilot"
    params:
        args = F13_VIS_ARGS,
    resources:
        mem_mb  = 52000,       # peak 38.0 GB (1 jobs)
        runtime = 15,          # peak 5 min
    log:
        LOGW + "/fig_2G.log",
    shell:
        "python " + F13_S + "/figure2/2G.py {params.args} > {log} 2>&1"


rule fig_2H:
    """2H: Number of germinal centres per 5,000 spots per clinical subtype."""
    input:
        gc_groups = expand(F13_GC_GROUPS, sample=VIS_SAMPLES),
        tls_model = expand(F13_TLS_MODEL, sample=VIS_SAMPLES),
        colors    = F13_SUBTYPE_COLORS,
        metadata  = F13_METADATA,
        sheet     = SHEETS["visium"],
    output:
        pdf = F13_OUT + "/Fig_2H.pdf",
    conda:
        "spatial-pilot"
    params:
        args = F13_VIS_ARGS,
    resources:
        mem_mb  = 4000,        # peak 0.5 GB (1 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/fig_2H.log",
    shell:
        "python " + F13_S + "/figure2/2H.py {params.args} > {log} 2>&1"


rule fig_2I:
    """2I: TLS niche scores and pathologist GC annotations on lymph node P12_LN_Vis."""
    input:
        h5ad        = expand(F13_VIS_H5AD, sample="P12_LN_Vis"),
        scores      = expand(F13_SCORE_LN, sample="P12_LN_Vis", score=F13_NICHE_SCORES),
        tissue_anno = F13_LN_ANNO_DIR + "/P12_LN_Vis_tissue_annotations.csv",
        gc_geojson  = F13_LN_ANNO_DIR + "/P12_LN_Vis_GC_annotations.geojson",
        image       = expand(F13_VIS_IMAGE, sample="P12_LN_Vis"),
        sheet       = SHEETS["visium_lymph_nodes"],
    output:
        pdf = F13_OUT + "/Fig_2I.pdf",
    conda:
        "shapely"
    params:
        args = base_args("visium_lymph_nodes") + " --sample_id P12_LN_Vis",
    resources:
        mem_mb  = 36000,       # peak 26.6 GB (1 jobs)
        runtime = 15,          # peak 3 min
    log:
        LOGW + "/fig_2I.log",
    shell:
        "python " + F13_S + "/figure2/2I.py {params.args} > {log} 2>&1"


rule fig_2J:
    """2J: GC score in pathologist-annotated GCs vs other spots (lymph nodes; ridge plot)."""
    input:
        gc_score = expand(F13_SCORE_LN, sample=LN_SAMPLES, score="GC"),
        gc_anno  = expand(F13_LN_ANNO_DIR + "/{sample}_GC_annotations.csv", sample=LN_SAMPLES),
        sheet    = SHEETS["visium_lymph_nodes"],
    output:
        pdf = F13_OUT + "/Fig_2J.pdf",
    conda:
        "spatial-pilot"
    params:
        args = F13_LN_ARGS,
    resources:
        mem_mb  = 4000,        # peak 0.8 GB (1 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/fig_2J.log",
    shell:
        "python " + F13_S + "/figure2/2J.py {params.args} > {log} 2>&1"


# ── Figure 3 ─────────────────────────────────────────────────────────────────
# The full-resolution Visium HD H&E images are 11-17 GB, hence the memory of 3A-3C.

F13_HD_RES = ["square_048um", "square_016um", "square_008um"]

rule fig_3A:
    """3A: HEV signature score on P29_HD at 48, 16 and 8 µm bins, with zoomed areas."""
    input:
        h5ad      = expand(F13_HD_H5AD, sample="P29_HD", res=F13_HD_RES),
        hev_score = expand(F13_HEV_SCORE, sample="P29_HD", res=F13_HD_RES),
        image     = expand(F13_HD_IMAGE, sample="P29_HD"),
        sheet     = SHEETS["visium_hd"],
    output:
        pdf = F13_OUT + "/Fig_3A.pdf",
    conda:
        "vis-hd"
    params:
        args = base_args("visium_hd") + " --sample_id P29_HD",
    resources:
        mem_mb  = 56000,       # peak 42.0 GB (1 jobs)
        runtime = 15,          # peak 4 min
    log:
        LOGW + "/fig_3A.log",
    shell:
        "python " + F13_S + "/figure3/3A.py {params.args} > {log} 2>&1"


rule fig_3B:
    """3B: Deconvolved endothelial types and pathology annotations in zoomed areas of P29_HD."""
    input:
        h5ad         = expand(F13_HD_H5AD, sample="P29_HD", res=["square_016um", "square_008um"]),
        endo_numbers = F13_VESSELS + "/P29_HD_endothelial_numbers_normalized.csv",
        annotations  = expand(F13_HD_ANNO, sample="P29_HD"),
        image        = expand(F13_HD_IMAGE, sample="P29_HD"),
        classes      = F13_CLASSES,
        colors       = F13_ENDO_COLORS,
        metadata     = F13_METADATA,
        sheet        = SHEETS["visium_hd"],
    output:
        pdf = F13_OUT + "/Fig_3B.pdf",
    conda:
        "vis-hd"
    params:
        args = base_args("visium_hd") + " --sample_id P29_HD",
    resources:
        mem_mb  = 56000,       # peak 41.1 GB (1 jobs)
        runtime = 15,          # peak 4 min
    log:
        LOGW + "/fig_3B.log",
    shell:
        "python " + F13_S + "/figure3/3B.py {params.args} > {log} 2>&1"


rule fig_3C:
    """3C: Example vessel of each endothelial type (H&E + deconvolution)."""
    input:
        h5ad         = expand(F13_HD_H5AD, sample=F13_3C_SAMPLES, res="square_016um"),
        endo_numbers = expand(F13_VESSELS + "/{sample}_endothelial_numbers_normalized.csv", sample=F13_3C_SAMPLES),
        image        = expand(F13_HD_IMAGE, sample=F13_3C_SAMPLES),
        classes      = F13_CLASSES,
        colors       = F13_ENDO_COLORS,
        metadata     = F13_METADATA,
        sheet        = SHEETS["visium_hd"],
    output:
        pdf = F13_OUT + "/Fig_3C.pdf",
    conda:
        "vis-hd"
    params:
        args = F13_HD_ARGS,
    resources:
        mem_mb  = 80000,       # peak 60.1 GB (1 jobs)
        runtime = 15,          # peak 6 min
    log:
        LOGW + "/fig_3C.log",
    shell:
        "python " + F13_S + "/figure3/3C.py {params.args} > {log} 2>&1"


rule fig_3D:
    """3D: Arterial/capillary/venous/lymphatic marker expression per vessel class (Kumar 2023)."""
    input:
        h5ad     = expand(F13_HD_H5AD, sample=HD_SAMPLES, res="square_016um"),
        vessels  = expand(F13_VESSELS + "/{sample}_vessel_classification.csv", sample=HD_SAMPLES),
        genes    = "resources/gene_lists/marker_genes/kumar_2023_endothelial_genes.csv",
        classes  = F13_CLASSES,
        metadata = F13_METADATA,
        sheet    = SHEETS["visium_hd"],
    output:
        pdf = F13_OUT + "/Fig_3D.pdf",
    conda:
        "vis-hd"
    params:
        args = F13_HD_ARGS,
    resources:
        mem_mb  = 12000,       # peak 7.9 GB (1 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/fig_3D.log",
    shell:
        "python " + F13_S + "/figure3/3D.py {params.args} > {log} 2>&1"


rule fig_3E:
    """3E: PVL subtype marker genes in the scRNA-seq reference (DEGs, matrix plot)."""
    input:
        sc_ref  = "data/anndata_objects/sc_ref/miniatlas_GEO.h5ad",
        classes = F13_CLASSES,
    output:
        pdf = F13_OUT + "/Fig_3E.pdf",
    conda:
        "vis-hd"
    params:
        args = F13_HD_ARGS,
    resources:
        mem_mb  = 12000,       # peak 7.6 GB (1 jobs)
        runtime = 15,          # peak 3 min
    log:
        LOGW + "/fig_3E.log",
    shell:
        "python " + F13_S + "/figure3/3E.py {params.args} > {log} 2>&1"


rule fig_3F:
    """3F: PVL subtype abundance vs distance from each vessel type."""
    input:
        pvl      = expand(F13_PVL_ALL_VESSELS, sample=HD_SAMPLES),
        classes  = F13_CLASSES,
        colors   = F13_ENDO_COLORS,
        metadata = F13_METADATA,
        sheet    = SHEETS["visium_hd"],
    output:
        pdf = F13_OUT + "/Fig_3F.pdf",
    conda:
        "vis-hd"
    params:
        args = F13_HD_ARGS,
    resources:
        mem_mb  = 4000,        # peak 1.4 GB (1 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/fig_3F.log",
    shell:
        "python " + F13_S + "/figure3/3F.py {params.args} > {log} 2>&1"


rule fig_3G:
    """3G: Proportion of vessels overlapping HEVs per endothelial type."""
    input:
        overlap  = expand(F13_VESSEL_HEV_OVERLAP, sample=HD_SAMPLES),
        classes  = F13_CLASSES,
        colors   = F13_ENDO_COLORS,
        metadata = F13_METADATA,
        sheet    = SHEETS["visium_hd"],
    output:
        pdf = F13_OUT + "/Fig_3G.pdf",
    conda:
        "vis-hd"
    params:
        args = F13_HD_ARGS,
    resources:
        mem_mb  = 4000,        # peak 0.3 GB (1 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/fig_3G.log",
    shell:
        "python " + F13_S + "/figure3/3G.py {params.args} > {log} 2>&1"


rule fig_3H:
    """3H: Pathology composition vs distance from HEV-overlapping vs other ACKR1+ vessels."""
    input:
        composition = expand(F13_HEV_VS_VENOUS + "/01_compare_microenvironment/csv_files/"
                             "{sample}_ACKR1_{group}_composition_vs_distance.csv",
                             sample=HD_SAMPLES, group=["no_overlap", "overlap"]),
        classes     = F13_CLASSES,
        metadata    = F13_METADATA,
        sheet       = SHEETS["visium_hd"],
    output:
        pdf = F13_OUT + "/Fig_3H.pdf",
    conda:
        "vis-hd"
    params:
        args = F13_HD_ARGS,
    resources:
        mem_mb  = 4000,        # peak 0.8 GB (1 jobs)
        runtime = 15,          # peak 2 min
    log:
        LOGW + "/fig_3H.log",
    shell:
        "python " + F13_S + "/figure3/3H.py {params.args} > {log} 2>&1"


rule fig_3I:
    """3I: Volcano plot of genes in HEV-overlapping vs other ACKR1+ vessels."""
    input:
        h5ad     = F13_HEV_VS_VENOUS + "/02_compare_gene_expression/csv_files/All_samples_HEV_ACKR1.h5ad",
        degs     = F13_HEV_VS_VENOUS + "/02_compare_gene_expression/csv_files/All_samples_rank_genes_HEV_overlap.csv",
        genes    = "resources/gene_lists/HEV_markers/HEV_sawada_2022.csv",
        classes  = F13_CLASSES,
        metadata = F13_METADATA,
        sheet    = SHEETS["visium_hd"],
    output:
        pdf = F13_OUT + "/Fig_3I.pdf",
    conda:
        "shapely"
    params:
        args = F13_HD_ARGS,
    resources:
        mem_mb  = 4000,        # peak 0.6 GB (1 jobs)
        runtime = 15,          # peak 2 min
    log:
        LOGW + "/fig_3I.log",
    shell:
        "python " + F13_S + "/figure3/3I.py {params.args} > {log} 2>&1"


rule fig_3J:
    """3J: PVL subtype enrichment around HEV-overlapping vs other ACKR1+ vessels (t-tests)."""
    input:
        pvl           = expand(F13_HEV_VS_VENOUS + "/03_compare_PVL_association/csv_files/"
                               "{sample}_PVL_composition_vs_distance.csv", sample=HD_SAMPLES),
        classes       = F13_CLASSES,
        colors        = F13_ENDO_COLORS,
        subset_colors = F13_SUBSET_COLORS,
        metadata      = F13_METADATA,
        sheet         = SHEETS["visium_hd"],
    output:
        pdf = F13_OUT + "/Fig_3J.pdf",
    conda:
        "vis-hd"
    params:
        args = F13_HD_ARGS,
    resources:
        mem_mb  = 4000,        # peak 0.8 GB (1 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/fig_3J.log",
    shell:
        "python " + F13_S + "/figure3/3J.py {params.args} > {log} 2>&1"


rule fig_3K:
    """3K: Number of HEVs vs TLS size (Pearson correlation)."""
    input:
        hev_per_tls = expand(F13_HEV_PER_TLS, sample=HD_SAMPLES),
        metadata    = F13_METADATA,
        sheet       = SHEETS["visium_hd"],
    output:
        pdf = F13_OUT + "/Fig_3K.pdf",
    conda:
        "vis-hd"
    params:
        args = F13_HD_ARGS,
    resources:
        mem_mb  = 4000,        # peak 0.7 GB (1 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/fig_3K.log",
    shell:
        "python " + F13_S + "/figure3/3K.py {params.args} > {log} 2>&1"


# ── Extended Data Figure 1 (supplementary overviews) ─────────────────────────

rule fig_S1C:
    """S1C: Presence of each pathology annotation category per section."""
    input:
        annotations = expand(F13_VIS_ANNO, sample=VIS_SAMPLES),
        categories  = F13_PATH_CATEGORIES,
        colors      = F13_SUBTYPE_COLORS,
        path_colors = F13_PATH_COLORS,
        metadata    = F13_METADATA,
        sheet       = SHEETS["visium"],
    output:
        pdf = F13_OUT + "/supplementary_overviews/Fig_S1C.pdf",
    conda:
        "spatial-pilot"
    params:
        args = F13_VIS_ARGS,
    resources:
        mem_mb  = 4000,        # peak 0.5 GB (1 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/fig_S1C.log",
    shell:
        "python " + F13_S + "/supplementary_overviews/S1C.py {params.args} > {log} 2>&1"


rule fig_S1G:
    """S1G: DEGs per pathology annotation (one vs rest, all sections) and GO BP 2025 enrichment."""
    input:
        h5ad        = expand(F13_VIS_H5AD, sample=VIS_SAMPLES),
        annotations = expand(F13_VIS_ANNO, sample=VIS_SAMPLES),
        gene_sets   = "resources/gene_lists/enrichment_analysis/GO_Biological_Process_2025.gmt",
        categories  = F13_PATH_CATEGORIES,
        path_colors = F13_PATH_COLORS,
    output:
        pdf        = F13_OUT + "/supplementary_overviews/Fig_S1G.pdf",
        markers    = F13_OUT + "/supplementary_overviews/csv_files/S1G_pathology_marker_genes.csv",
        samples    = F13_OUT + "/supplementary_overviews/csv_files/S1G_samples_included.csv",
        enrichment = F13_OUT + "/supplementary_overviews/csv_files/S1G_pathology_GO_enrichment.csv",
    conda:
        "gseapy"
    params:
        args = F13_VIS_ARGS,
    resources:
        # Concatenates all 33 sections; needs ~190 GB (script header), so hugemem with margin
        mem_mb  = 256000,
        runtime = 240,
    log:
        LOGW + "/fig_S1G.log",
    shell:
        "python " + F13_S + "/supplementary_overviews/S1G.py {params.args} > {log} 2>&1"


# ── targets ──────────────────────────────────────────────────────────────────

for _f13_rule in [rules.fig_1C, rules.fig_1D, rules.fig_1E, rules.fig_1F, rules.fig_1G, rules.fig_1H,
                  rules.fig_2A, rules.fig_2B, rules.fig_2C, rules.fig_2D, rules.fig_2E, rules.fig_2F,
                  rules.fig_2G, rules.fig_2H, rules.fig_2I, rules.fig_2J,
                  rules.fig_3A, rules.fig_3B, rules.fig_3C, rules.fig_3D, rules.fig_3E, rules.fig_3F,
                  rules.fig_3G, rules.fig_3H, rules.fig_3I, rules.fig_3J, rules.fig_3K,
                  rules.fig_S1C, rules.fig_S1G]:
    TARGETS += list(_f13_rule.output)
