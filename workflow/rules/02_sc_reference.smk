# Rule for the scRNA-seq reference (scripts/02_cell_type_mapping): builds
# {H5AD}/sc_ref/miniatlas_GEO.h5ad from the GEO archive of Wu et al. 2021 (GSE176078), which is not on
# Zenodo. Used by x05_reference_producers and fig_3E. The archive is a DVC target:
# dvc pull resources/sc_ref/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz.dvc

rule s02_build_sc_ref:
    """02: extract the GEO archive and build the reference AnnData object (raw counts in
    layers['counts'], normalized and log1p-transformed .X, Ensembl IDs as var_names)."""
    input:
        archive = "resources/sc_ref/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz",
        mapping = "resources/sc_ref/ensembl_mapping/map.csv",
        genes   = "resources/sc_ref/ensembl_mapping/genes.tsv",
    output:
        H5AD + "/sc_ref/miniatlas_GEO.h5ad",
    conda:
        "spatial-pilot"
    params:
        args = base_args("visium"),
    resources:
        mem_mb  = 20000,       # peak 13.5 GB (1 jobs)
        runtime = 15,          # peak 3 min
    log:
        LOGW + "/s02_build_sc_ref.log",
    shell:
        "tar -xzf {input.archive} -C resources/sc_ref && "
        "python scripts/02_cell_type_mapping/visium_deconvolution/01_prepare_sc_ref/01_create_GEO_miniatlas_object.py "
        "{params.args} > {log} 2>&1"
