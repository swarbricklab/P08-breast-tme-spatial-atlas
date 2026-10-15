# Rules for scripts/figures/figure4-5: one rule per manuscript panel (Fig_4A ... Fig_5H).
#   Xenium panels (4A-4E, 4K, 5A-5C, 5E-5H) read the Xenium objects + 05_xenium_mechanisms outputs.
#   Visium HD panels (4F-4J, 5D) read 04_vascular_microenvironment outputs, incl. the Visium HD
#   ligand-receptor analysis (up to ~4 h per sample), so their figures go into TARGETS_LONG.
#
# Multi-sample scripts get --samples (see scripts/run_multi_sample_job.sh); 4F, 5B and 5E are set up
# for a single sample (hard-coded in the script) and get --sample_id (see run_single_sample_job.sh).
# Resources: from the PBS summaries of earlier runs in logs/*.out, with margin.
# No f-strings here: snakemake 7.32 on Python 3.12 mangles them (adds spaces around {fields}).

F45_FIG = RES + "/figures"
F45_COLORS = "resources/plotting_colors"
F45_METADATA = "resources/tissue_metadata.csv"
F45_CLASSIFICATIONS = "resources/cell_types/miniatlas/miniatlas_subset_classifications.csv"

# Published inputs
F45_XEN_H5AD = expand(H5AD + "/xenium/{sample}.h5ad", sample=XEN_SAMPLES)
F45_XEN_LABELS = expand(RES + "/02_cell_type_mapping/label_transfer/xenium/subset_level/"
                         "{sample}_subset_pred_cell_types.csv", sample=XEN_SAMPLES)
F45_HD_H5AD = expand(H5AD + "/visium_hd/square_016um/{sample}.h5ad", sample=HD_SAMPLES)
F45_HD_DECON = expand(RES + "/02_cell_type_mapping/deconvolution/visium_hd/subset_level/square_016um/"
                      "{sample}_subset_square_016um_q05_cell_abundance_w_sf.csv", sample=HD_SAMPLES)

# Upstream analysis outputs (rules in 04_vascular_microenvironment.smk / 05_xenium_mechanisms.smk)
F45_XM = RES + "/05_xenium_mechanisms"
F45_HEV = expand(F45_XM + "/02_identify_HEV_regions/csv_files/{sample}_HEV_niches_{kind}.csv",
                 sample=XEN_SAMPLES, kind=["cell_to_niche", "info"])
F45_TLS = expand(F45_XM + "/01_identify_TLS_regions/csv_files/{sample}_cell_TLS_annotations.csv",
                 sample=XEN_SAMPLES)
F45_STING_CORR = (F45_XM + "/05_gene_gene_correlations/01_investigate_gene_correlations/TMEM173_correlations/"
                  "csv_files/All_samples_HEV_TMEM173_gene_correlations.csv")

F45_LR = RES + "/04_vascular_microenvironment/04_ligand_receptor_analysis"
# stLearn objects per sample; 4F-4J use the 25 um bins, 5D the 0 um bins (HEV bins only)
F45_LR_H5AD = (F45_LR + "/01_run_LR_analysis_stlearn/{bin}um_bin/csv_files/"
               "{sample}_HEV_LR_analysis_stlearn.h5ad")
F45_LR_ALL = F45_LR + "/02_analyze_LR_all/{bin}um_bin/csv_files/All_samples_{table}.csv"
F45_LR_SPOTS = expand(F45_LR + "/03_plot_LR_per_sample/25um_bin/csv_files/{sample}_HEV_spot_grouping.csv",
                      sample=HD_SAMPLES)


# ── 4A · ACKR1+ EC marker expression, HEV vs non-HEV niches (Xenium) ──────────

rule fig_4A:
    """Fig 4A: dotplot of HEV markers / LRI receptors in ACKR1+ ECs, venous vs HEV niches."""
    input:
        h5ad     = F45_XEN_H5AD,
        labels   = F45_XEN_LABELS,
        hev      = F45_HEV,
        sheet    = SHEETS["xenium"],
        metadata = F45_METADATA,
        colors   = [F45_COLORS + "/endothelial_colors.csv", F45_COLORS + "/TLS_HEV_colors.csv"],
    output:
        pdf = F45_FIG + "/Fig_4A.pdf",
    conda:
        "xen-hev"
    params:
        args    = base_args("xenium"),
        samples = " ".join(XEN_SAMPLES),
    resources:
        mem_mb  = 28000,       # peak 21.4 GB (3 jobs)
        runtime = 15,          # peak 2 min
    log:
        LOGW + "/fig_4A.log",
    shell:
        """
        python scripts/figures/figure4/4A.py {params.args} \\
            --samples {params.samples} \\
        > {log} 2>&1
        """

TARGETS += rules.fig_4A.output


# ── 4B · HEV TMEM173 (STING) expression vs distance to TLS (Xenium) ──────────

rule fig_4B:
    """Fig 4B: TMEM173 (STING) expression in HEV ECs as a function of distance to TLS."""
    input:
        h5ad     = F45_XEN_H5AD,
        labels   = F45_XEN_LABELS,
        hev      = F45_HEV,
        sheet    = SHEETS["xenium"],
        metadata = F45_METADATA,
        colors   = F45_COLORS + "/TLS_HEV_colors.csv",
    output:
        pdf = F45_FIG + "/Fig_4B.pdf",
    conda:
        "xen-hev"
    params:
        args    = base_args("xenium"),
        samples = " ".join(XEN_SAMPLES),
    resources:
        mem_mb  = 24000,       # peak 17.3 GB (3 jobs)
        runtime = 15,          # peak 2 min
    log:
        LOGW + "/fig_4B.log",
    shell:
        """
        python scripts/figures/figure4/4B.py {params.args} \\
            --samples {params.samples} \\
        > {log} 2>&1
        """

TARGETS += rules.fig_4B.output


# ── 4C · STING signature scores in venous ECs by environment (Xenium) ───────

rule fig_4C:
    """Fig 4C: STING signature scores in canonical venous, stromal HEV-like and TLS-HEV ECs."""
    input:
        h5ad       = F45_XEN_H5AD,
        labels     = F45_XEN_LABELS,
        hev        = F45_HEV,
        scores     = expand(F45_XM + "/04_score_signatures/01_score_STING_signatures/"
                            "csv_files/{sample}_STING_signature_scores.csv", sample=XEN_SAMPLES),
        signatures = ["resources/gene_lists/STING/STING_signatures.csv",
                      "resources/gene_lists/STING/type_I_interferon_signatures.csv"],
        sheet      = SHEETS["xenium"],
        metadata   = F45_METADATA,
        colors     = [F45_COLORS + "/endothelial_colors.csv", F45_COLORS + "/TLS_HEV_colors.csv"],
    output:
        pdf = F45_FIG + "/Fig_4C.pdf",
    conda:
        "xen-hev"
    params:
        args    = base_args("xenium"),
        samples = " ".join(XEN_SAMPLES),
    resources:
        mem_mb  = 36000,       # peak 25.3 GB (4 jobs)
        runtime = 15,          # peak 3 min
    log:
        LOGW + "/fig_4C.log",
    shell:
        """
        python scripts/figures/figure4/4C.py {params.args} \\
            --samples {params.samples} \\
        > {log} 2>&1
        """

TARGETS += rules.fig_4C.output


# ── 4D · genes correlated with TMEM173 per HEV + GSEA (Xenium) ───────────────

rule fig_4D_1:
    """Fig 4D (left): genes ranked by correlation with mean TMEM173 expression per HEV niche."""
    input:
        corr = F45_STING_CORR,
    output:
        pdf = F45_FIG + "/Fig_4D-1.pdf",
    conda:
        "shapely"
    params:
        args    = base_args("xenium"),
        samples = " ".join(XEN_SAMPLES),
    resources:
        mem_mb  = 4000,        # peak 0.7 GB (4 jobs)
        runtime = 15,          # peak 2 min
    log:
        LOGW + "/fig_4D_1.log",
    shell:
        """
        python scripts/figures/figure4/4D-1.py {params.args} \\
            --samples {params.samples} \\
        > {log} 2>&1
        """

TARGETS += rules.fig_4D_1.output


rule fig_4D_2:
    """Fig 4D (right): GSEA running-score plots (GO Molecular Function) for the TMEM173-ranked genes."""
    input:
        gsea     = (F45_XM + "/05_gene_gene_correlations/02_GSEA_analysis/TMEM173_correlations/"
                    "GO_Molecular_Function_2025/csv_files/All_samples_HEV_TMEM173_preranked_GSEA_results.csv"),
        corr     = F45_STING_CORR,
        gene_set = "resources/gene_lists/enrichment_analysis/GO_Molecular_Function_2025.gmt",
    output:
        pdf = F45_FIG + "/Fig_4D-2.pdf",
    conda:
        "gseapy"
    params:
        args    = base_args("xenium"),
        samples = " ".join(XEN_SAMPLES),
    resources:
        mem_mb  = 4000,        # peak 0.7 GB (4 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/fig_4D_2.log",
    shell:
        """
        python scripts/figures/figure4/4D-2.py {params.args} \\
            --samples {params.samples} \\
        > {log} 2>&1
        """

TARGETS += rules.fig_4D_2.output


# ── 4E · adhesion-receptor expression near HEVs (Xenium) ─────────────────────

rule fig_4E:
    """Fig 4E: integrin / SPN / CD38 receptor expression across HEV-enriched cell types <50 um of HEVs."""
    input:
        h5ad            = F45_XEN_H5AD,
        labels          = F45_XEN_LABELS,
        distances       = expand(F45_XM + "/03_HEV_microenvironment/01_calculate_distance/csv_files/"
                                 "{sample}_closest_HEV_distances.csv", sample=XEN_SAMPLES),
        classifications = F45_CLASSIFICATIONS,
        colors          = F45_COLORS + "/miniatlas_major_colors.csv",
    output:
        pdf = F45_FIG + "/Fig_4E.pdf",
    conda:
        "xen-hev"
    params:
        args    = base_args("xenium"),
        samples = " ".join(XEN_SAMPLES),
    resources:
        mem_mb  = 48000,       # peak 35.9 GB (4 jobs)
        runtime = 15,          # peak 2 min
    log:
        LOGW + "/fig_4E.log",
    shell:
        """
        python scripts/figures/figure4/4E.py {params.args} \\
            --samples {params.samples} \\
        > {log} 2>&1
        """

TARGETS += rules.fig_4E.output


# ── 4F · zoom TLS-HEVs vs stromal HEV-like vessels, P29_HD (Visium HD) ───────

rule fig_4F:
    """Fig 4F: H&E zoom with HEV groups used in the LR analysis (P29_HD, hard-coded in the script)."""
    input:
        lr_h5ad  = expand(F45_LR_H5AD, bin=25, sample="P29_HD"),
        image    = "data/images/he/aligned_images/visium_hd/P29_HD.ome.tif",
        area     = "resources/areas_of_interest/visium_hd/P29_HD_HEV_LR_area.geojson",
        sheet    = SHEETS["visium_hd"],
        metadata = F45_METADATA,
        colors   = F45_COLORS + "/TLS_HEV_colors.csv",
    output:
        pdf = F45_FIG + "/Fig_4F.pdf",
    conda:
        "shapely"
    params:
        args = base_args("visium_hd"),
    resources:
        mem_mb  = 44000,       # peak 31.9 GB (1 jobs)
        runtime = 15,          # peak 5 min
    log:
        LOGW + "/fig_4F.log",
    shell:
        """
        python scripts/figures/figure4/4F.py {params.args} \\
            --sample_id P29_HD \\
        > {log} 2>&1
        """

TARGETS_LONG += rules.fig_4F.output


# ── 4G · key LR interactions, TLS-HEVs vs stromal HEV-like (Visium HD) ──────

rule fig_4G:
    """Fig 4G: scatterplot of LR interactions between HEV groups (25 um bins)."""
    input:
        lr_all = expand(F45_LR_ALL, bin=25, table=["group_counts", "LR_summary"]),
        sheet  = SHEETS["visium_hd"],
    output:
        pdf = F45_FIG + "/Fig_4G.pdf",
    conda:
        "shapely"          # adjustText (only in this env)
    params:
        args    = base_args("visium_hd"),
        samples = " ".join(HD_SAMPLES),
    resources:
        mem_mb  = 4000,        # peak 0.8 GB (1 jobs)
        runtime = 15,          # peak 3 min
    log:
        LOGW + "/fig_4G.log",
    shell:
        """
        python scripts/figures/figure4/4G.py {params.args} \\
            --samples {params.samples} \\
        > {log} 2>&1
        """

TARGETS_LONG += rules.fig_4G.output


# ── 4H · top LR pairs enriched in TLS-HEVs (Visium HD) ───────────────────────

rule fig_4H:
    """Fig 4H: dotplot of top LR pairs enriched in TLS-HEVs, source-target combinations (25 um bins)."""
    input:
        lr_h5ad = expand(F45_LR_H5AD, bin=25, sample=HD_SAMPLES),
        lr_all  = expand(F45_LR_ALL, bin=25, table=["group_counts", "LR_summary"]),
        sheet   = SHEETS["visium_hd"],
        colors  = F45_COLORS + "/TLS_HEV_colors.csv",
    output:
        pdf = F45_FIG + "/Fig_4H.pdf",
    conda:
        "vis-hd"
    params:
        args    = base_args("visium_hd"),
        samples = " ".join(HD_SAMPLES),
    resources:
        mem_mb  = 4000,        # peak 2.7 GB (1 jobs)
        runtime = 15,          # peak 2 min
    log:
        LOGW + "/fig_4H.log",
    shell:
        """
        python scripts/figures/figure4/4H.py {params.args} \\
            --samples {params.samples} \\
        > {log} 2>&1
        """

TARGETS_LONG += rules.fig_4H.output


# ── 4I · lymphotoxin system, TLS-HEVs vs stromal HEVs (Visium HD) ───────────

rule fig_4I:
    """Fig 4I: expression of lymphotoxin signaling components in HEV / perivascular bins."""
    input:
        h5ad   = F45_HD_H5AD,
        spots  = F45_LR_SPOTS,
        sheet  = SHEETS["visium_hd"],
        colors = F45_COLORS + "/TLS_HEV_colors.csv",
    output:
        pdf = F45_FIG + "/Fig_4I.pdf",
    conda:
        "vis-hd"
    params:
        args    = base_args("visium_hd"),
        samples = " ".join(HD_SAMPLES),
    resources:
        mem_mb  = 20000,       # peak 12.8 GB (1 jobs)
        runtime = 15,          # peak 2 min
    log:
        LOGW + "/fig_4I.log",
    shell:
        """
        python scripts/figures/figure4/4I.py {params.args} \\
            --samples {params.samples} \\
        > {log} 2>&1
        """

TARGETS_LONG += rules.fig_4I.output


# ── 4J · cell types in spots sending LTB signals to HEVs (Visium HD) ────────

rule fig_4J:
    """Fig 4J: top 10 cell types (deconvolution) enriched in spots sending LTB signals to HEVs."""
    input:
        lr_h5ad         = expand(F45_LR_H5AD, bin=25, sample=HD_SAMPLES),
        spots           = F45_LR_SPOTS,
        decon           = F45_HD_DECON,
        sheet           = SHEETS["visium_hd"],
        metadata        = F45_METADATA,
        classifications = F45_CLASSIFICATIONS,
        colors          = F45_COLORS + "/miniatlas_major_colors.csv",
    output:
        pdf = F45_FIG + "/Fig_4J.pdf",
    conda:
        "vis-hd"
    params:
        args    = base_args("visium_hd"),
        samples = " ".join(HD_SAMPLES),
    resources:
        mem_mb  = 8000,        # peak 4.1 GB (1 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/fig_4J.log",
    shell:
        """
        python scripts/figures/figure4/4J.py {params.args} \\
            --samples {params.samples} \\
        > {log} 2>&1
        """

TARGETS_LONG += rules.fig_4J.output


# ── 4K · cell type-specific expression inside TLS (Xenium) ───────────────────

rule fig_4K:
    """Fig 4K: expression of key genes per cell type inside TLS regions."""
    input:
        h5ad            = F45_XEN_H5AD,
        labels          = F45_XEN_LABELS,
        tls             = F45_TLS,
        sheet           = SHEETS["xenium"],
        classifications = F45_CLASSIFICATIONS,
        colors          = F45_COLORS + "/miniatlas_major_colors.csv",
    output:
        pdf = F45_FIG + "/Fig_4K.pdf",
    conda:
        "xen-hev"
    params:
        args    = base_args("xenium"),
        samples = " ".join(XEN_SAMPLES),
    resources:
        mem_mb  = 36000,       # peak 27.8 GB (3 jobs)
        runtime = 15,          # peak 4 min
    log:
        LOGW + "/fig_4K.log",
    shell:
        """
        python scripts/figures/figure4/4K.py {params.args} \\
            --samples {params.samples} \\
        > {log} 2>&1
        """

TARGETS += rules.fig_4K.output


# ── 5A · HEV CCL19 / CCL21 expression vs distance to TLS (Xenium) ───────────

rule fig_5A:
    """Fig 5A: CCL19 and CCL21 expression in HEVs vs distance to TLS."""
    input:
        h5ad   = F45_XEN_H5AD,
        labels = F45_XEN_LABELS,
        hev    = F45_HEV,
        colors = F45_COLORS + "/TLS_HEV_colors.csv",
    output:
        pdf = F45_FIG + "/Fig_5A.pdf",
    conda:
        "xen-hev"
    params:
        args    = base_args("xenium"),
        samples = " ".join(XEN_SAMPLES),
    resources:
        mem_mb  = 24000,       # peak 16.0 GB (3 jobs)
        runtime = 15,          # peak 2 min
    log:
        LOGW + "/fig_5A.log",
    shell:
        """
        python scripts/figures/figure5/5A.py {params.args} \\
            --samples {params.samples} \\
        > {log} 2>&1
        """

TARGETS += rules.fig_5A.output


# ── 5B · CCL19 / CCL21 / CXCL13 / TSPAN7 transcripts, P29_Xe5k (Xenium) ──────
# 5B and 5E read data/spatialdata_objects/xenium/P29_Xe5k.zarr, unpacked once by x05_unzip_spatialdata.

rule fig_5B:
    """Fig 5B: chemokine / HEV marker transcripts on H&E in TLS & HEV regions (P29_Xe5k)."""
    input:
        h5ad     = H5AD + "/xenium/P29_Xe5k.h5ad",
        sdata    = "data/spatialdata_objects/xenium/P29_Xe5k.zarr",     # unpacked by x05_unzip_spatialdata
        image    = "data/images/he/aligned_images/xenium/P29_Xe5k.ome.tif",
        areas    = "resources/areas_of_interest/xenium/P29_Xe5k_cropped_HE_TLS_HEV_new.geojson",
        sheet    = SHEETS["xenium"],
        metadata = F45_METADATA,
        colors   = F45_COLORS + "/transcript_colors.csv",
    output:
        pdf = F45_FIG + "/Fig_5B.pdf",
    conda:
        "spatdat"
    params:
        args = base_args("xenium"),
    resources:
        mem_mb  = 56000,       # peak 40.9 GB (2 jobs)
        runtime = 15,          # peak 5 min
    log:
        LOGW + "/fig_5B.log",
    shell:
        """
        python scripts/figures/figure5/5B.py {params.args} \\
            --sample_id P29_Xe5k \\
        > {log} 2>&1
        """

TARGETS += rules.fig_5B.output


# ── 5C · chemokine / receptor expression per cell type inside TLS (Xenium) ───

rule fig_5C:
    """Fig 5C: chemokine and receptor expression per cell type inside TLS regions."""
    input:
        h5ad            = F45_XEN_H5AD,
        labels          = F45_XEN_LABELS,
        tls             = F45_TLS,
        classifications = F45_CLASSIFICATIONS,
        colors          = F45_COLORS + "/miniatlas_major_colors.csv",
    output:
        pdf = F45_FIG + "/Fig_5C.pdf",
    conda:
        "xen-hev"
    params:
        args    = base_args("xenium"),
        samples = " ".join(XEN_SAMPLES),
    resources:
        mem_mb  = 36000,       # peak 27.8 GB (3 jobs)
        runtime = 15,          # peak 3 min
    log:
        LOGW + "/fig_5C.log",
    shell:
        """
        python scripts/figures/figure5/5C.py {params.args} \\
            --samples {params.samples} \\
        > {log} 2>&1
        """

TARGETS += rules.fig_5C.output


# ── 5D · LR interactions of HEV-upregulated genes (Visium HD) ───────────────

rule fig_5D:
    """Fig 5D: LR pairs of genes upregulated in HEVs vs ACKR1+ ECs (bin_size 0: HEV bins only)."""
    input:
        lr_h5ad = expand(F45_LR_H5AD, bin=0, sample=HD_SAMPLES),
        lr_all  = expand(F45_LR_ALL, bin=0, table=["group_counts", "LR_summary"]),
        de      = (RES + "/04_vascular_microenvironment/03_compare_HEV_to_other_venous/02_compare_gene_expression/"
                   "csv_files/All_samples_rank_genes_HEV_overlap.csv"),
        colors  = F45_COLORS + "/TLS_HEV_colors.csv",
    output:
        pdf = F45_FIG + "/Fig_5D.pdf",
    conda:
        "vis-hd"
    params:
        args    = base_args("visium_hd"),
        samples = " ".join(HD_SAMPLES),
    resources:
        mem_mb  = 4000,        # peak 1.0 GB (1 jobs)
        runtime = 15,          # peak 2 min
    log:
        LOGW + "/fig_5D.log",
    shell:
        """
        python scripts/figures/figure5/5D.py {params.args} \\
            --samples {params.samples} \\
        > {log} 2>&1
        """

TARGETS_LONG += rules.fig_5D.output


# ── 5E · SPP1 transcripts in HEV regions, P29_Xe5k (Xenium) ──────────────────

rule fig_5E:
    """Fig 5E: SPP1 transcripts on H&E in TLS / HEV regions (P29_Xe5k)."""
    input:
        h5ad     = H5AD + "/xenium/P29_Xe5k.h5ad",
        sdata    = "data/spatialdata_objects/xenium/P29_Xe5k.zarr",     # unpacked by x05_unzip_spatialdata
        image    = "data/images/he/aligned_images/xenium/P29_Xe5k.ome.tif",
        areas    = "resources/areas_of_interest/xenium/P29_Xe5k_cropped_HE_TLS_SPP1.geojson",
        sheet    = SHEETS["xenium"],
        metadata = F45_METADATA,
        colors   = F45_COLORS + "/transcript_colors.csv",
    output:
        pdf = F45_FIG + "/Fig_5E.pdf",
    conda:
        "spatdat"
    params:
        args = base_args("xenium"),
    resources:
        mem_mb  = 56000,       # peak 41.4 GB (2 jobs)
        runtime = 15,          # peak 5 min
    log:
        LOGW + "/fig_5E.log",
    shell:
        """
        python scripts/figures/figure5/5E.py {params.args} \\
            --sample_id P29_Xe5k \\
        > {log} 2>&1
        """

TARGETS += rules.fig_5E.output


# ── 5F · HEV SPP1 expression vs distance to TLS (Xenium) ────────────────────

rule fig_5F:
    """Fig 5F: SPP1 expression in HEVs as a function of proximity to TLS."""
    input:
        h5ad     = F45_XEN_H5AD,
        labels   = F45_XEN_LABELS,
        hev      = F45_HEV,
        sheet    = SHEETS["xenium"],
        metadata = F45_METADATA,
        colors   = F45_COLORS + "/TLS_HEV_colors.csv",
    output:
        pdf = F45_FIG + "/Fig_5F.pdf",
    conda:
        "xen-hev"
    params:
        args    = base_args("xenium"),
        samples = " ".join(XEN_SAMPLES),
    resources:
        mem_mb  = 24000,       # peak 17.8 GB (3 jobs)
        runtime = 15,          # peak 2 min
    log:
        LOGW + "/fig_5F.log",
    shell:
        """
        python scripts/figures/figure5/5F.py {params.args} \\
            --samples {params.samples} \\
        > {log} 2>&1
        """

TARGETS += rules.fig_5F.output


# ── 5G · expression inside vs outside TLS per cell type (Xenium) ─────────────

rule fig_5G:
    """Fig 5G: cell type-specific expression in cells inside vs outside TLS regions."""
    input:
        h5ad            = F45_XEN_H5AD,
        labels          = F45_XEN_LABELS,
        tls             = F45_TLS,
        classifications = F45_CLASSIFICATIONS,
        colors          = F45_COLORS + "/miniatlas_major_colors.csv",
    output:
        pdf = F45_FIG + "/Fig_5G.pdf",
    conda:
        "xen-hev"
    params:
        args    = base_args("xenium"),
        samples = " ".join(XEN_SAMPLES),
    resources:
        mem_mb  = 76000,       # peak 56.6 GB (2 jobs)
        runtime = 15,          # peak 4 min
    log:
        LOGW + "/fig_5G.log",
    shell:
        """
        python scripts/figures/figure5/5G.py {params.args} \\
            --samples {params.samples} \\
        > {log} 2>&1
        """

TARGETS += rules.fig_5G.output


# ── 5H · HEV SPP1 vs nearby cell abundance correlations (Xenium) ─────────────

rule fig_5H:
    """Fig 5H: per-sample correlation of HEV SPP1 expression with cell type proportions <200 um."""
    input:
        corr            = expand(F45_XM + "/06_cellular_enrichment/SPP1_vs_cell_abundance/maxdist_200um/csv_files/"
                                 "{sample}_HEV_SPP1_and_cell_proportions_corr.csv", sample=XEN_SAMPLES),
        classifications = F45_CLASSIFICATIONS,
        colors          = F45_COLORS + "/miniatlas_major_colors.csv",
    output:
        pdf = F45_FIG + "/Fig_5H.pdf",
    conda:
        "xen-hev"
    params:
        args    = base_args("xenium"),
        samples = " ".join(XEN_SAMPLES),
    resources:
        mem_mb  = 4000,        # peak 0.5 GB (4 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/fig_5H.log",
    shell:
        """
        python scripts/figures/figure5/5H.py {params.args} \\
            --samples {params.samples} \\
        > {log} 2>&1
        """

TARGETS += rules.fig_5H.output
