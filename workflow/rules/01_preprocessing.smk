# Rules for scripts/01_preprocessing: build the AnnData and SpatialData objects from the 10x pipeline
# outputs deposited on GEO (GSE351850), for any object that is not already present.
#
# The analyses read the objects from {H5AD}/ and data/spatialdata_objects/, the paths `dvc pull` fills
# from Zenodo. The object rules only match the objects that are missing when the workflow starts: objects
# pulled from Zenodo (or built in an earlier run) are used as they are and never rebuilt. To rebuild an
# object from GEO, delete it first.
#
# Building needs the GEO imports (dvc pull data/geo/<platform>.dvc) and, for Xenium, the aligned H&E
# images and alignment matrices (data/images/he/aligned_images/xenium.dvc, data/images/he/alignment_files.dvc).
#   p01_unpack_geo            GEO files -> data/processed/ (Visium HD 16 and 48 um bins rebuilt)
#   Visium, Visium HD         import + QC -> log-normalization
#   Xenium                    AnnData object + H&E image; SpatialData object -> registered to the H&E

P01_S = "scripts/01_preprocessing"
P01_UNPACKED = "data/processed/.unpacked_{platform}"
P01_VIS = H5AD + "/visium/{sample}.h5ad"
P01_HD_RES = ["square_008um", "square_016um", "square_048um"]      # the bin sizes the analyses use
P01_HD = expand(H5AD + "/visium_hd/{res}/{{sample}}.h5ad", res=P01_HD_RES)
P01_XEN = H5AD + "/xenium/{sample}.h5ad"
P01_XEN_SDATA = "data/spatialdata_objects/xenium/{sample}.zip"
P01_XEN_IMAGES = ["data/images/he/aligned_images/xenium/{sample}.ome.tif",
                  "data/images/he/alignment_files/xenium/{sample}_matrix.csv"]


def p01_missing(samples, objects):
    """Wildcard constraint matching the samples with any of objects (paths with {sample}) missing."""
    missing = [s for s in samples if not all(os.path.exists(o.format(sample=s)) for o in objects)]
    return "|".join(missing) or "(?!)"     # (?!) matches nothing: every object is present


rule p01_unpack_geo:
    """00: unpack the GEO files of one platform into the Space Ranger / Xenium Onboard Analysis layout."""
    input:
        "data/geo/{platform}",
    output:
        touch(P01_UNPACKED),
    wildcard_constraints:
        platform = "visium|visium_hd|xenium",
    conda:
        "vis-hd"
    resources:
        mem_mb  = lambda wc: {"visium": 4000, "visium_hd": 32000, "xenium": 4000}[wc.platform],   # peak 20.4 GB (visium_hd)
        runtime = lambda wc: {"visium": 15, "visium_hd": 45, "xenium": 45}[wc.platform],          # peak 14 min (visium_hd), 10 min (xenium)
    log:
        LOGW + "/p01_unpack_geo/{platform}.log",
    shell:
        "python " + P01_S + "/00_unpack_geo.py --project_dir " + PROJECT_DIR
        + " --modality {wildcards.platform} > {log} 2>&1"


def p01_visium_args(wildcards):
    sheet = "visium_lymph_nodes" if wildcards.sample in LN_SAMPLES else "visium"
    return base_args(sheet) + " --sample_id " + wildcards.sample


rule p01_visium_import:
    """visium/01: import the Space Ranger outputs and run basic QC."""
    input:
        expand(P01_UNPACKED, platform="visium"),
    output:
        temp(H5AD + "/visium/raw/{sample}_raw.h5ad"),
    wildcard_constraints:
        sample = "|".join(VIS_SAMPLES + LN_SAMPLES),
    conda:
        "spatial-pilot"
    params:
        args = p01_visium_args,
    resources:
        mem_mb  = 4000,        # peak 1.0 GB (2 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/p01_visium_import/{sample}.log",
    shell:
        "python " + P01_S + "/visium/01_import_basic_qc.py {params.args} > {log} 2>&1"


rule p01_visium_object:
    """visium/02: log-normalize -> the Visium AnnData object."""
    input:
        H5AD + "/visium/raw/{sample}_raw.h5ad",
    output:
        P01_VIS,
    wildcard_constraints:
        sample = p01_missing(VIS_SAMPLES + LN_SAMPLES, [P01_VIS]),
    conda:
        "spatial-pilot"
    params:
        args = p01_visium_args,
    resources:
        mem_mb  = 8000,        # peak 2.3 GB (2 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/p01_visium_object/{sample}.log",
    shell:
        "python " + P01_S + "/visium/02_log_norm.py {params.args} > {log} 2>&1"


rule p01_visium_hd_import:
    """visium_hd/01: import the Space Ranger outputs of every bin size and run basic QC."""
    input:
        expand(P01_UNPACKED, platform="visium_hd"),
    output:
        temp(expand(H5AD + "/visium_hd/{res}/raw/{{sample}}_raw.h5ad",
                    res=["square_002um"] + P01_HD_RES)),
    wildcard_constraints:
        sample = "|".join(HD_SAMPLES),
    conda:
        "vis-hd"
    params:
        args = lambda wc: base_args("visium_hd") + " --sample_id " + wc.sample,
    resources:
        mem_mb  = 32000,       # peak 12.4 GB (1 job)
        runtime = 30,          # peak 4 min
    log:
        LOGW + "/p01_visium_hd_import/{sample}.log",
    shell:
        "python " + P01_S + "/visium_hd/01_import_basic_qc.py {params.args} > {log} 2>&1"


rule p01_visium_hd_object:
    """visium_hd/02: log-normalize -> the Visium HD AnnData objects (the 2 um one is written too)."""
    input:
        expand(H5AD + "/visium_hd/{res}/raw/{{sample}}_raw.h5ad", res=["square_002um"] + P01_HD_RES),
    output:
        P01_HD,
    wildcard_constraints:
        sample = p01_missing(HD_SAMPLES, P01_HD),
    conda:
        "vis-hd"
    params:
        args = lambda wc: base_args("visium_hd") + " --sample_id " + wc.sample,
    resources:
        mem_mb  = 48000,       # peak 19.2 GB (1 job)
        runtime = 30,          # peak 4 min
    log:
        LOGW + "/p01_visium_hd_object/{sample}.log",
    shell:
        "python " + P01_S + "/visium_hd/02_log_norm.py {params.args} > {log} 2>&1"


rule p01_xenium_object:
    """xenium/01a + 02_handle_images/01: the Xenium AnnData object, with the H&E image added."""
    input:
        expand(P01_UNPACKED, platform="xenium"),
        P01_XEN_IMAGES,
    output:
        P01_XEN,
    wildcard_constraints:
        sample = p01_missing(XEN_SAMPLES, [P01_XEN]),
    conda:
        "xen-hev"
    params:
        args = lambda wc: base_args("xenium") + " --sample_id " + wc.sample,
    resources:
        mem_mb  = 32000,       # peak 8.7 GB (1 job)
        runtime = 30,          # peak 3 min
    log:
        LOGW + "/p01_xenium_object/{sample}.log",
    shell:
        "(python " + P01_S + "/xenium/01_make_objects/01a_make_adata_objects.py {params.args} && "
        "python " + P01_S + "/xenium/02_handle_images/01_make_adata_objects_with_images.py {params.args}) > {log} 2>&1"


rule p01_xenium_spatialdata_unregistered:
    """xenium/01b: the Xenium SpatialData object before registration to the H&E image."""
    input:
        expand(P01_UNPACKED, platform="xenium"),
    output:
        temp("data/spatialdata_objects/xenium/without_image/{sample}.zip"),
    wildcard_constraints:
        sample = "|".join(XEN_SAMPLES),
    conda:
        "spatdat"
    params:
        args = lambda wc: base_args("xenium") + " --sample_id " + wc.sample,
    resources:
        mem_mb  = 64000,       # peak 24.2 GB (1 job)
        runtime = 60,          # peak 12 min
    log:
        LOGW + "/p01_xenium_spatialdata_unregistered/{sample}.log",
    shell:
        "python " + P01_S + "/xenium/01_make_objects/01b_make_spatialdata_objects.py {params.args} > {log} 2>&1"


rule p01_xenium_spatialdata:
    """xenium/02_handle_images/02: register the SpatialData object to the H&E image."""
    input:
        "data/spatialdata_objects/xenium/without_image/{sample}.zip",
        P01_XEN_IMAGES,
    output:
        P01_XEN_SDATA,
    wildcard_constraints:
        sample = p01_missing(XEN_SAMPLES, [P01_XEN_SDATA]),
    conda:
        "spatdat"
    params:
        args = lambda wc: base_args("xenium") + " --sample_id " + wc.sample,
    resources:
        mem_mb  = 64000,       # peak 28.8 GB (1 job)
        runtime = 60,          # peak 8 min
    log:
        LOGW + "/p01_xenium_spatialdata/{sample}.log",
    shell:
        "python " + P01_S + "/xenium/02_handle_images/02_transform_spatialdata_object.py {params.args} > {log} 2>&1"


# Every object the analyses read, for `snakemake objects`
P01_OBJECTS = (expand(P01_VIS, sample=VIS_SAMPLES + LN_SAMPLES)
               + expand(P01_HD, sample=HD_SAMPLES)
               + expand(P01_XEN, sample=XEN_SAMPLES)
               + expand(P01_XEN_SDATA, sample=XEN_SAMPLES))
