# Rules for scripts/03_TLS_niche_discovery (Visium breast tumours + matched lymph nodes):
#   01 identify TLS          - TLS model score per spot, grouping of TLS-high spots (tumours + lymph
#                              nodes), discrete-TLS number/size/clustering/topology (Fig S3B)
#   02 identify TLS niches   - Harmony/Leiden clustering of TLS-high spots across tumour sections
#   03 extract signatures    - niche marker genes; GC filtered against invasive epithelium, HEV
#                              against endothelial-high spots; final signature files (used by 04)
#   04 map signatures        - niche scores per spot (tumours + lymph nodes), high-score areas
#
# All scripts write under --res_dir (config res_dir), so res_dir can be changed.
# No f-strings in this file: Snakemake 7.32 on Python 3.12 mangles them (pads the replacement
# fields with spaces and collapses escaped braces), so paths are built by concatenation.
# Inputs from outside the workflow: {H5AD}/visium (tumour and lymph node objects),
# data/pathology_annotations/visium, results/02_cell_type_mapping/deconvolution/visium/major_level,
# resources/ (gene lists, tissue metadata, colours, pathology categories).

T3_SCRIPTS = "scripts/03_TLS_niche_discovery"
T3 = RES + "/03_TLS_niche_discovery"
TLS_MODEL = T3 + "/01_identify_TLS/01_run_TLS_model"
TLS_GROUPS = T3 + "/01_identify_TLS/02_group_TLS_high/threshold_01"
TLS_SPATIAL = T3 + "/01_identify_TLS/03_TLS_spatial_organization"
TLS_CLUSTER = T3 + "/02_identify_TLS_niches/01_cluster_TLS_high"
NICHE_SIGS = T3 + "/03_extract_niche_signatures"
HEV_FILTER = NICHE_SIGS + "/03_filter_HEV_signature"
NICHE_MAP = T3 + "/04_map_niche_signatures"

# Niche scores computed by 04_map_niche_signatures (score_name in the scripts)
NICHE_SCORES = ["GC", "B_cell_zone", "T_cell_zone"]
# Signature files written by 04_organize_signature_files.py (HEVs is used by 04_vascular_microenvironment)
FINAL_SIGNATURES = [NICHE_SIGS + "/04_final_niche_signature_files/" + f for f in [
    "germinal_center_filtered_epithelial_genes.csv",
    "b_cell_zone_top_genes_filtered.csv",
    "t_cell_zone_top_genes_filtered.csv",
    "hevs_filtered_endothelial_genes.csv",
]]

# Tumour and lymph node sections share {H5AD}/visium and most result folders (lymph nodes go to a
# lymph_nodes/ subfolder), so per-sample rules are restricted to their own sample list.
VIS_RE = "|".join(VIS_SAMPLES)
LN_RE = "|".join(LN_SAMPLES)
VIS_SAMPLE_ARGS = " ".join(VIS_SAMPLES)

T3_METADATA = "resources/tissue_metadata.csv"
T3_SUBTYPE_COLORS = "resources/plotting_colors/SCSubtype_colors.csv"
T3_TLS_COEFS = "resources/gene_lists/andersson_2021/andersson_2021_coef_full_visium.csv"
VIS_H5AD = H5AD + "/visium/{sample}.h5ad"
VIS_ANNOTATION = "data/pathology_annotations/visium/{sample}_pathology_annotations.csv"


# ── 01 · identify TLS ────────────────────────────────────────────────────────

rule v03_run_TLS_model:
    """01: Andersson 2021 TLS model score per spot (tumour sections)."""
    input:
        h5ad  = VIS_H5AD,
        coefs = T3_TLS_COEFS,
        sheet = SHEETS["visium"],
    output:
        csv = TLS_MODEL + "/{sample}_TLS_model.csv",
    wildcard_constraints:
        sample = VIS_RE,
    conda:
        "spatial-pilot"
    params:
        args = base_args("visium"),
    resources:
        mem_mb  = 8000,        # peak 4.9 GB (33 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/v03_run_TLS_model/{sample}.log",
    shell:
        "python {T3_SCRIPTS}/01_identify_TLS/01_run_TLS_model.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


rule v03_run_TLS_model_lymph_nodes:
    """01: Andersson 2021 TLS model score per spot (lymph node sections)."""
    input:
        h5ad  = VIS_H5AD,
        coefs = T3_TLS_COEFS,
        sheet = SHEETS["visium_lymph_nodes"],
    output:
        csv = TLS_MODEL + "/lymph_nodes/{sample}_TLS_model.csv",
    wildcard_constraints:
        sample = LN_RE,
    conda:
        "spatial-pilot"
    params:
        args = base_args("visium_lymph_nodes"),
    resources:
        mem_mb  = 4000,        # peak 2.0 GB (7 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/v03_run_TLS_model_lymph_nodes/{sample}.log",
    shell:
        "python {T3_SCRIPTS}/01_identify_TLS/01_run_TLS_model.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


rule v03_group_TLS_high:
    """02: Group spatially contiguous TLS-high spots (score > 0.01; tumour sections)."""
    input:
        h5ad      = VIS_H5AD,
        tls_model = rules.v03_run_TLS_model.output.csv,
        metadata  = T3_METADATA,
        sheet     = SHEETS["visium"],
    output:
        csv = TLS_GROUPS + "/csv_files/{sample}_TLS_groups.csv",
        pdf = TLS_GROUPS + "/plots/groups/{sample}_TLS_groups.pdf",
    wildcard_constraints:
        sample = VIS_RE,
    conda:
        "spatial-pilot"
    params:
        args = base_args("visium"),
    resources:
        mem_mb  = 8000,        # peak 3.8 GB (33 jobs)
        runtime = 15,          # peak 2 min
    log:
        LOGW + "/v03_group_TLS_high/{sample}.log",
    shell:
        "python {T3_SCRIPTS}/01_identify_TLS/02_group_TLS_high.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


rule v03_group_TLS_high_lymph_nodes:
    """02: Group spatially contiguous TLS-high spots (score > 0.01; lymph node sections)."""
    input:
        h5ad      = VIS_H5AD,
        tls_model = rules.v03_run_TLS_model_lymph_nodes.output.csv,
        metadata  = T3_METADATA,
        sheet     = SHEETS["visium_lymph_nodes"],
    output:
        csv = TLS_GROUPS + "/lymph_nodes/csv_files/{sample}_TLS_groups.csv",
        pdf = TLS_GROUPS + "/lymph_nodes/plots/groups/{sample}_TLS_groups.pdf",
    wildcard_constraints:
        sample = LN_RE,
    conda:
        "spatial-pilot"
    params:
        args = base_args("visium_lymph_nodes"),
    resources:
        mem_mb  = 8000,        # peak 3.7 GB (7 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/v03_group_TLS_high_lymph_nodes/{sample}.log",
    shell:
        "python {T3_SCRIPTS}/01_identify_TLS/02_group_TLS_high.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


rule v03_TLS_distances_topology:
    """03/01: Discrete-TLS topology (intratumoral / peritumoral / stromal) and distances."""
    input:
        groups      = expand(rules.v03_group_TLS_high.output.csv, sample=VIS_SAMPLES),
        h5ad        = expand(VIS_H5AD, sample=VIS_SAMPLES),
        annotations = expand(VIS_ANNOTATION, sample=VIS_SAMPLES),
        metadata    = T3_METADATA,
        sheet       = SHEETS["visium"],
    output:
        csv = TLS_SPATIAL + "/01_TLS_distances_topology/threshold_01/csv_files/distances_per_TLS.csv",
    conda:
        "spatial-pilot"
    params:
        args    = base_args("visium"),
        samples = VIS_SAMPLE_ARGS,
    resources:
        mem_mb  = 4000,        # peak 2.7 GB (1 jobs)
        runtime = 15,          # peak 8 min
    log:
        LOGW + "/v03_TLS_distances_topology.log",
    shell:
        "python {T3_SCRIPTS}/01_identify_TLS/03_TLS_spatial_organization/01_TLS_distances_topology.py "
        "{params.args} --samples {params.samples} > {log} 2>&1"


rule v03_TLS_clustering:
    """03/02: DBSCAN (eps 500 um, edge-to-edge) clustering of discrete TLS per section."""
    input:
        groups   = expand(rules.v03_group_TLS_high.output.csv, sample=VIS_SAMPLES),
        h5ad     = expand(VIS_H5AD, sample=VIS_SAMPLES),
        metadata = T3_METADATA,
        sheet    = SHEETS["visium"],
    output:
        per_tls    = TLS_SPATIAL + "/02_TLS_clustering/threshold_01/csv_files/clustering_per_TLS.csv",
        per_sample = TLS_SPATIAL + "/02_TLS_clustering/threshold_01/csv_files/clustering_per_sample.csv",
    conda:
        "spatial-pilot"
    params:
        args    = base_args("visium"),
        samples = VIS_SAMPLE_ARGS,
    resources:
        mem_mb  = 16000,       # peak 10.3 GB (1 jobs)
        runtime = 15,          # peak 3 min
    log:
        LOGW + "/v03_TLS_clustering.log",
    shell:
        "python {T3_SCRIPTS}/01_identify_TLS/03_TLS_spatial_organization/02_TLS_clustering.py "
        "{params.args} --samples {params.samples} > {log} 2>&1"


rule v03_TLS_number_size_overview:
    """03/03: Cohort overview of TLS number, size, clustering and topology (Fig S3B)."""
    input:
        tls_model  = expand(rules.v03_run_TLS_model.output.csv, sample=VIS_SAMPLES),
        distances  = rules.v03_TLS_distances_topology.output.csv,
        clustering = rules.v03_TLS_clustering.output.per_tls,
        metadata   = T3_METADATA,
        colors     = T3_SUBTYPE_COLORS,
        sheet      = SHEETS["visium"],
    output:
        pdf     = TLS_SPATIAL + "/03_TLS_number_size_overview/threshold_01/plots/TLS_number_size_overview.pdf",
        summary = TLS_SPATIAL + "/03_TLS_number_size_overview/threshold_01/csv_files/summary_per_sample.csv",
        sizes   = TLS_SPATIAL + "/03_TLS_number_size_overview/threshold_01/csv_files/sizes_per_aggregate.csv",
    conda:
        "spatial-pilot"
    params:
        # This script takes neither --h5ad_dir nor --samples (all sections of the sample sheet)
        args = ("--project_dir " + PROJECT_DIR + " --res_dir " + PROJECT_DIR + "/" + RES
                + " --sample_sheet " + PROJECT_DIR + "/" + SHEETS["visium"]),
    resources:
        mem_mb  = 4000,        # peak 0.3 GB (1 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/v03_TLS_number_size_overview.log",
    shell:
        "python {T3_SCRIPTS}/01_identify_TLS/03_TLS_spatial_organization/03_TLS_number_size_overview.py "
        "{params.args} > {log} 2>&1"


rule v03_TLS_per_sample_plots:
    """03/04: Per-section spatial plots of epithelial annotation, TLS topology and clustering."""
    input:
        groups      = expand(rules.v03_group_TLS_high.output.csv, sample=VIS_SAMPLES),
        h5ad        = expand(VIS_H5AD, sample=VIS_SAMPLES),
        annotations = expand(VIS_ANNOTATION, sample=VIS_SAMPLES),
        distances   = rules.v03_TLS_distances_topology.output.csv,
        clustering  = rules.v03_TLS_clustering.output.per_tls,
        metadata    = T3_METADATA,
        sheet       = SHEETS["visium"],
    output:
        pdf = TLS_SPATIAL + "/04_TLS_per_sample_plots/threshold_01/plots/TLS_per_sample_topology_clustering.pdf",
    conda:
        "spatial-pilot"
    params:
        args    = base_args("visium"),
        samples = VIS_SAMPLE_ARGS,
    resources:
        mem_mb  = 24000,       # peak 16.1 GB (1 jobs)
        runtime = 15,          # peak 5 min
    log:
        LOGW + "/v03_TLS_per_sample_plots.log",
    shell:
        "python {T3_SCRIPTS}/01_identify_TLS/03_TLS_spatial_organization/04_TLS_per_sample_plots.py "
        "{params.args} --samples {params.samples} > {log} 2>&1"


# ── 02 · identify TLS niches ─────────────────────────────────────────────────

rule v03_cluster_TLS_high:
    """02/01: Harmony integration + Leiden clustering (res 1.5) of TLS-high spots, marker genes."""
    input:
        h5ad      = expand(VIS_H5AD, sample=VIS_SAMPLES),
        tls_model = expand(rules.v03_run_TLS_model.output.csv, sample=VIS_SAMPLES),
        groups    = expand(rules.v03_group_TLS_high.output.csv, sample=VIS_SAMPLES),
        metadata  = T3_METADATA,
        colors    = T3_SUBTYPE_COLORS,
        sheet     = SHEETS["visium"],
    output:
        clusters = TLS_CLUSTER + "/csv_files/All_samples_leiden_clusters.csv",
        degs     = TLS_CLUSTER + "/csv_files/All_samples_leiden_clusters_DEGS.csv",
        h5ad     = TLS_CLUSTER + "/csv_files/All_samples_TLS_umap.h5ad",
        pdf      = TLS_CLUSTER + "/plots/All_samples_TLS_groups_umap.pdf",
    conda:
        "spatial-pilot"
    params:
        args    = base_args("visium"),
        samples = VIS_SAMPLE_ARGS,
    threads: 2
    resources:
        mem_mb  = 68000,       # peak 51.1 GB (1 jobs)
        runtime = 15,          # peak 7 min
    log:
        LOGW + "/v03_cluster_TLS_high.log",
    shell:
        "python {T3_SCRIPTS}/02_identify_TLS_niches/01_cluster_TLS_high.py "
        "{params.args} --samples {params.samples} > {log} 2>&1"


# ── 03 · extract niche signatures ────────────────────────────────────────────

rule v03_extract_niche_signatures:
    """03/01: Filter cluster markers by specificity; signatures of the T cell zone, B cell zone,
    HEV and germinal centre clusters (identified by marker genes)."""
    input:
        h5ad = rules.v03_cluster_TLS_high.output.h5ad,
        degs = rules.v03_cluster_TLS_high.output.degs,
    output:
        t_cell_zone     = NICHE_SIGS + "/01_extract_niche_signatures/t_cell_zone_top_genes_filtered.csv",
        b_cell_zone     = NICHE_SIGS + "/01_extract_niche_signatures/b_cell_zone_top_genes_filtered.csv",
        hevs            = NICHE_SIGS + "/01_extract_niche_signatures/hevs_top_genes_filtered.csv",
        germinal_center = NICHE_SIGS + "/01_extract_niche_signatures/germinal_center_top_genes_filtered.csv",
        pct_expressing  = TLS_CLUSTER + "/top_genes_dotplot/post_filtering/csv_files/All_clusters_pct_expressing.csv",
        pdfs            = expand(TLS_CLUSTER + "/top_genes_dotplot/post_filtering/"
                                 "All_samples_TLS_umap_marker_genes_top_{n}_LFC_filtered.pdf",
                                 n=[20, 25, 30]),
    conda:
        "spatial-pilot"
    params:
        args    = base_args("visium"),
        samples = VIS_SAMPLE_ARGS,
    resources:
        mem_mb  = 8000,        # peak 4.3 GB (1 jobs)
        runtime = 15,          # peak 2 min
    log:
        LOGW + "/v03_extract_niche_signatures.log",
    shell:
        "python {T3_SCRIPTS}/03_extract_niche_signatures/01_extract_niche_signatures.py "
        "{params.args} --samples {params.samples} > {log} 2>&1"


rule v03_filter_GC_signature:
    """03/02: Remove GC signature genes also expressed in invasive epithelium."""
    input:
        h5ad        = expand(VIS_H5AD, sample=VIS_SAMPLES),
        annotations = expand(VIS_ANNOTATION, sample=VIS_SAMPLES),
        gc_genes    = rules.v03_extract_niche_signatures.output.germinal_center,
        categories  = "resources/pathology_annotations/pathology_categories_per_layer.json",
        metadata    = T3_METADATA,
        colors      = T3_SUBTYPE_COLORS,
        sheet       = SHEETS["visium"],
    output:
        filtered      = NICHE_SIGS + "/02_filter_GC_signature/germinal_center_filtered_epithelial_genes.csv",
        documentation = NICHE_SIGS + "/02_filter_GC_signature/germinal_center_filtering_documentation.csv",
        per_gene      = NICHE_SIGS + "/02_filter_GC_signature/csv_files/All_samples_mean_pct_expressing_per_gene.csv",
        pdf           = NICHE_SIGS + "/02_filter_GC_signature/plots/All_samples_GC_score_filtering.pdf",
    conda:
        "spatial-pilot"
    params:
        args    = base_args("visium"),
        samples = VIS_SAMPLE_ARGS,
    resources:
        mem_mb  = 124000,      # peak 96.7 GB (1 jobs)
        runtime = 30,          # peak 8 min
    log:
        LOGW + "/v03_filter_GC_signature.log",
    shell:
        "python {T3_SCRIPTS}/03_extract_niche_signatures/02_filter_GC_signature.py "
        "{params.args} --samples {params.samples} > {log} 2>&1"


rule v03_identify_endothelial_high:
    """03/03a: Endothelial-high spots per section (cell2location major level) vs the HEV cluster."""
    input:
        h5ad     = VIS_H5AD,
        clusters = rules.v03_cluster_TLS_high.output.clusters,
        hevs     = rules.v03_extract_niche_signatures.output.hevs,
        decon    = RES + "/02_cell_type_mapping/deconvolution/visium/major_level/{sample}_major_q05_cell_abundance_w_sf.csv",
        metadata = T3_METADATA,
        sheet    = SHEETS["visium"],
    output:
        csv = HEV_FILTER + "/identify_endothelial_high/csv_files/{sample}_endothelial_high_spots.csv",
        pdf = HEV_FILTER + "/identify_endothelial_high/plots/{sample}_identify_endothelial_high.pdf",
    wildcard_constraints:
        sample = VIS_RE,
    conda:
        "spatial-pilot"
    params:
        args = base_args("visium"),
    resources:
        mem_mb  = 8000,        # peak 4.4 GB (33 jobs)
        runtime = 15,          # peak 3 min
    log:
        LOGW + "/v03_identify_endothelial_high/{sample}.log",
    shell:
        "python {T3_SCRIPTS}/03_extract_niche_signatures/03a_identify_endothelial_high.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


rule v03_create_endothelial_object:
    """03/03b: Pool endothelial-high and HEV-cluster spots; LFC HEV vs endothelial-high."""
    input:
        h5ad     = expand(VIS_H5AD, sample=VIS_SAMPLES),
        endo     = expand(rules.v03_identify_endothelial_high.output.csv, sample=VIS_SAMPLES),
        degs     = rules.v03_cluster_TLS_high.output.degs,
        metadata = T3_METADATA,
        colors   = T3_SUBTYPE_COLORS,
        sheet    = SHEETS["visium"],
    output:
        h5ad         = HEV_FILTER + "/create_endothelial_object/csv_files/All_samples_endothelial_object.h5ad",
        genes        = HEV_FILTER + "/create_endothelial_object/csv_files/All_samples_HEV_vs_endothelial_genes.csv",
        positive_lfc = HEV_FILTER + "/create_endothelial_object/csv_files/All_samples_HEV_vs_endothelial_positive_LFC.csv",
        spatial_pdf  = HEV_FILTER + "/create_endothelial_object/plots/All_samples_endothelial_classification_spatial.pdf",
        genes_pdf    = HEV_FILTER + "/create_endothelial_object/plots/All_samples_endothelial_vs_HEV_genes.pdf",
    conda:
        "spatial-pilot"
    params:
        args    = base_args("visium"),
        samples = VIS_SAMPLE_ARGS,
    resources:
        mem_mb  = 100000,      # peak 77.0 GB (1 jobs)
        runtime = 15,          # peak 4 min
    log:
        LOGW + "/v03_create_endothelial_object.log",
    shell:
        "python {T3_SCRIPTS}/03_extract_niche_signatures/03b_create_endothelial_object.py "
        "{params.args} --samples {params.samples} > {log} 2>&1"


rule v03_filter_HEV_signature:
    """03/03c: Keep HEV genes expressed in HEV spots but not general in endothelial-high spots."""
    input:
        h5ad         = rules.v03_create_endothelial_object.output.h5ad,
        positive_lfc = rules.v03_create_endothelial_object.output.positive_lfc,
        metadata     = T3_METADATA,
        sheet        = SHEETS["visium"],
    output:
        filtered      = HEV_FILTER + "/filter_HEV_signature/csv_files/hevs_filtered_endothelial_genes.csv",
        documentation = HEV_FILTER + "/filter_HEV_signature/csv_files/hevs_filtering_documentation.csv",
        pdf           = HEV_FILTER + "/filter_HEV_signature/plots/HEV_signature_top_50_LFC_filtered.pdf",
    conda:
        "spatial-pilot"
    params:
        args    = base_args("visium"),
        samples = VIS_SAMPLE_ARGS,
    resources:
        mem_mb  = 8000,        # peak 4.6 GB (1 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/v03_filter_HEV_signature.log",
    shell:
        "python {T3_SCRIPTS}/03_extract_niche_signatures/03c_filter_HEV_signature.py "
        "{params.args} --samples {params.samples} > {log} 2>&1"


rule v03_organize_signature_files:
    """03/04: Top 300 genes of the final GC, B cell zone, T cell zone and HEV signatures."""
    input:
        t_cell_zone = rules.v03_extract_niche_signatures.output.t_cell_zone,
        b_cell_zone = rules.v03_extract_niche_signatures.output.b_cell_zone,
        gc          = rules.v03_filter_GC_signature.output.filtered,
        hevs        = rules.v03_filter_HEV_signature.output.filtered,
    output:
        FINAL_SIGNATURES,
    conda:
        "spatial-pilot"
    params:
        args    = base_args("visium"),
        samples = VIS_SAMPLE_ARGS,
    resources:
        mem_mb  = 4000,        # peak 0.2 GB (1 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/v03_organize_signature_files.log",
    shell:
        "python {T3_SCRIPTS}/03_extract_niche_signatures/04_organize_signature_files.py "
        "{params.args} --samples {params.samples} > {log} 2>&1"


# ── 04 · map niche signatures ────────────────────────────────────────────────

rule v03_score_signatures:
    """04/01: GC, B cell zone and T cell zone scores (top 20/50/100 genes) per spot (tumours)."""
    input:
        h5ad       = VIS_H5AD,
        tls_model  = rules.v03_run_TLS_model.output.csv,
        groups     = rules.v03_group_TLS_high.output.csv,
        signatures = "resources/gene_lists/TLS_niche_signatures/niche_signatures.csv",   # published signatures (use_resource_signatures = True)
        metadata   = T3_METADATA,
        sheet      = SHEETS["visium"],
    output:
        scores = expand(NICHE_MAP + "/01_score_signatures/{score}_score/csv_files/{sample}_{score}_score.csv",
                        score=NICHE_SCORES, allow_missing=True),
        genes  = expand(NICHE_MAP + "/01_score_signatures/{score}_score/csv_files/{sample}_genes_used.csv",
                        score=NICHE_SCORES, allow_missing=True),
        pdfs   = expand(NICHE_MAP + "/01_score_signatures/{score}_score/plots/{sample}_{score}_score.pdf",
                        score=NICHE_SCORES, allow_missing=True),
    wildcard_constraints:
        sample = VIS_RE,
    conda:
        "spatial-pilot"
    params:
        args = base_args("visium"),
    resources:
        mem_mb  = 8000,        # peak 4.4 GB (33 jobs)
        runtime = 15,          # peak 2 min
    log:
        LOGW + "/v03_score_signatures/{sample}.log",
    shell:
        "python {T3_SCRIPTS}/04_map_niche_signatures/01_score_signatures.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


rule v03_score_signatures_lymph_nodes:
    """04/01: GC, B cell zone and T cell zone scores per spot (lymph nodes; Fig 2i, 2j)."""
    input:
        h5ad       = VIS_H5AD,
        tls_model  = rules.v03_run_TLS_model_lymph_nodes.output.csv,
        groups     = rules.v03_group_TLS_high_lymph_nodes.output.csv,
        signatures = "resources/gene_lists/TLS_niche_signatures/niche_signatures.csv",   # published signatures (use_resource_signatures = True)
        metadata   = T3_METADATA,
        sheet      = SHEETS["visium_lymph_nodes"],
    output:
        scores = expand(NICHE_MAP + "/01_score_signatures/lymph_nodes/{score}_score/csv_files/{sample}_{score}_score.csv",
                        score=NICHE_SCORES, allow_missing=True),
        genes  = expand(NICHE_MAP + "/01_score_signatures/lymph_nodes/{score}_score/csv_files/{sample}_genes_used.csv",
                        score=NICHE_SCORES, allow_missing=True),
        pdfs   = expand(NICHE_MAP + "/01_score_signatures/lymph_nodes/{score}_score/plots/{sample}_{score}_score.pdf",
                        score=NICHE_SCORES, allow_missing=True),
    wildcard_constraints:
        sample = LN_RE,
    conda:
        "spatial-pilot"
    params:
        args = base_args("visium_lymph_nodes"),
    resources:
        mem_mb  = 8000,        # peak 3.2 GB (7 jobs)
        runtime = 15,          # peak 4 min
    log:
        LOGW + "/v03_score_signatures_lymph_nodes/{sample}.log",
    shell:
        "python {T3_SCRIPTS}/04_map_niche_signatures/01_score_signatures.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


rule v03_group_high_score_areas:
    """04/02: Group contiguous spots with a niche score (top 20 genes) > 0.7 (tumours; Fig 2g, 2h)."""
    input:
        h5ad      = VIS_H5AD,
        tls_model = rules.v03_run_TLS_model.output.csv,
        scores    = rules.v03_score_signatures.output.scores,
        metadata  = T3_METADATA,
        sheet     = SHEETS["visium"],
    output:
        groups   = expand(NICHE_MAP + "/02_group_high_score_areas/{score}_score/csv_files/{sample}_{score}_groups.csv",
                          score=NICHE_SCORES, allow_missing=True),
        settings = expand(NICHE_MAP + "/02_group_high_score_areas/{score}_score/csv_files/{sample}_grouping_settings.json",
                          score=NICHE_SCORES, allow_missing=True),
        pdfs     = expand(NICHE_MAP + "/02_group_high_score_areas/{score}_score/plots/groups/{sample}_{score}_groups.pdf",
                          score=NICHE_SCORES, allow_missing=True),
    wildcard_constraints:
        sample = VIS_RE,
    conda:
        "spatial-pilot"
    params:
        args = base_args("visium"),
    resources:
        mem_mb  = 8000,        # peak 4.0 GB (33 jobs)
        runtime = 15,          # peak 2 min
    log:
        LOGW + "/v03_group_high_score_areas/{sample}.log",
    shell:
        "python {T3_SCRIPTS}/04_map_niche_signatures/02_group_high_score_areas.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


# ── targets ──────────────────────────────────────────────────────────────────

def _t3_per_sample(rule_outputs, samples):
    """Expand every output of a per-sample rule over a sample list."""
    return [f for pattern in rule_outputs for f in expand(pattern, sample=samples)]

TARGETS += _t3_per_sample(rules.v03_group_TLS_high.output, VIS_SAMPLES)
TARGETS += _t3_per_sample(rules.v03_group_TLS_high_lymph_nodes.output, LN_SAMPLES)
TARGETS += list(rules.v03_TLS_number_size_overview.output)
TARGETS += list(rules.v03_TLS_per_sample_plots.output)
TARGETS += list(rules.v03_cluster_TLS_high.output)
TARGETS += list(rules.v03_extract_niche_signatures.output)
TARGETS += list(rules.v03_filter_GC_signature.output)
TARGETS += _t3_per_sample(rules.v03_identify_endothelial_high.output, VIS_SAMPLES)
TARGETS += list(rules.v03_create_endothelial_object.output)
TARGETS += list(rules.v03_filter_HEV_signature.output)
TARGETS += list(rules.v03_organize_signature_files.output)
TARGETS += _t3_per_sample(rules.v03_score_signatures.output, VIS_SAMPLES)
TARGETS += _t3_per_sample(rules.v03_score_signatures_lymph_nodes.output, LN_SAMPLES)
TARGETS += _t3_per_sample(rules.v03_group_high_score_areas.output, VIS_SAMPLES)
