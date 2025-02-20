# Libraries
from os.path import join

# Includes
include: "load_samples.smk"

# Variables
QC_FILTERED_DIR = join(config["qc_dir"],"filtered")

SCT_DIR = join(config["downstream_dir"],"normalized")
MERGED_DIR = join(config["downstream_dir"],"merged")
CELLXGENE_DIR = config["cellxgene_dir"]
PLOTS_DIR = join(config["plots_dir"],"integration")

plots_list = ["umap_integrated", "umap_unknown", "umap_celltypes"]

all_output_list = [join(MERGED_DIR,"merged.qs"),
                   join(CELLXGENE_DIR,"merged.h5ad"),
                   *[join(PLOTS_DIR,"{0}.png".format(e)) for e in plots_list]
                  ]

# Rules
rule integration_all:
    input:
        expand(all_output_list,
                sample = SAMPLES_DF.index)
    default_target: True


rule integration_sct:
    conda: 
        config["conda_env"]
    input:
        os.path.join(QC_FILTERED_DIR,"{sample}.rds")
    output:
        os.path.join(SCT_DIR,"{sample}.rds")
    log:
        err="logs/{sample}/integration_sct.err",
        log="logs/{sample}/integration_sct.log",
    script:
        "../../scripts/integration_sct.R"


rule integration_merge:
    conda: 
        config["conda_env"]
    input:
        samples_rds = expand(rules.integration_sct.output,
                             sample = SAMPLES_DF.index)
    output:
        integrated = join(MERGED_DIR,"merged.qs"),
        features = join(MERGED_DIR,"merged_features.qs"),
        anchors = join(MERGED_DIR,"merged_anchors.qs")
    params:
        samples = SAMPLES_DF.index
    log:
        err="logs/integration_merge.err",
        log="logs/integration_merge.log",
    script:
        "../../scripts/integration_merge.R"


rule integration_plots:
    conda: 
        config["conda_env"]
    input:
        rules.integration_merge.output.integrated
    output:
        [join(PLOTS_DIR,"{0}.png".format(e)) for e in plots_list]
    params:
        plots_list = plots_list
    log:
        err="logs/integration_plot.err",
        log="logs/integration_plot.log"
    script:
        "../../scripts/make_plots_integration.R"

rule integration_cellxgene:
    conda: 
        config["conda_env"]
    input:
        rules.integration_merge.output.integrated
    output:
        join(CELLXGENE_DIR,"merged.h5ad")
    log:
        err="logs/integration_cellxgene.err",
        log="logs/integration_cellxgene.log"
    script:
        "../../scripts/integration_cellxgene.R"