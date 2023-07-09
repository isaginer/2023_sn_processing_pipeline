from snakemake.utils import min_version
min_version("6.0")

include: "rules/load_samples.smk"
include: "rules/00_cellranger.smk"
include: "rules/01_decontx.smk"
include: "rules/02_doublets.smk"
include: "rules/03_qc_filtering.smk"


all_input_list = [rules.cellranger_all.input,
                  rules.doublets_all.input,
                  rules.qc_all.input]

if (not(config.get("skip_soup", False))):
    all_input_list.append([rules.decontx_all.input])


rule all:
    input:
        all_input_list
    default_target: True

# rule all:
    # input:
        # os.path.join(config["results_dir"],"cellranger_stats.csv"),
        # expand(os.path.join(config["doublets_dir"],"params","{sample}","doublet_finder_pk.rds"),
        #        sample = SAMPLES_DF.index),
        # expand(os.path.join(config["after_qc_dir"],"metadata/{sample}.rds"),
        #        sample = SAMPLES_DF.index),
        # expand(os.path.join(config["downstream_dirs"]["normalization"],"{sample}.rds"),
        #        sample = SAMPLES_DF.index),
        # expand(os.path.join(config["doublets_dir_2"],f'{{sample}}.rds'),
        #        sample = SAMPLES_DF.index),
        # expand(os.path.join(config["re_doublets_dir"],f'{{sample}}.rds'),
        #        sample = SAMPLES_DF.index),
        # os.path.join(config["downstream_dirs"]["integration"],"merged.qs"),
        # os.path.join(config["downstream_dirs"]["garnett"],"merged.qs"),
        # os.path.join(config["cellxgene_dir"],"merged.h5ad")

# rule sct_normalization:
#     conda: "2023_aus_brain"
#     input:
#         os.path.join(config["after_qc_dir"],"{sample}.rds")
#     output:
#         os.path.join(config["downstream_dirs"]["normalization"],"{sample}.rds")
#     log:
#         "logs/{sample}/sct_normalization.log"
#     script:
#         "../scripts/normalization.R"

# rule integrate:
#     conda: "2023_aus_brain"
#     input:
#         phenodata = config["phenodata"],
#         samples_rds = expand(os.path.join(config["downstream_dirs"]["normalization"],"{sample}.rds"),
#                sample = SAMPLES_DF.index)
#     output:
#         integrated = os.path.join(config["downstream_dirs"]["integration"],"merged.qs"),
#         features = os.path.join(config["downstream_dirs"]["integration"],"merged_features.qs"),
#         anchors = os.path.join(config["downstream_dirs"]["integration"],"merged_anchors.qs")
#     params:
#         samples = SAMPLES_DF.index
#     log:
#         "logs/integrate.log"
#     script:
#         "../scripts/integration.R"

# rule predict_merged_cell_types:
#     conda: "2023_aus_brain"
#     input: rules.integrate.output.integrated
#     output: os.path.join(config["downstream_dirs"]["garnett"],"merged.qs")
#     log:
#         "logs/predict_merged_cell_types.log"
#     script:
#         "../scripts/garnett_merged_predict.R"

# rule prepare_cellxgene:
#     conda: "2023_aus_brain"
#     input:
#         rules.predict_merged_cell_types.output
#     output:
#         os.path.join(config["cellxgene_dir"],"merged.h5ad")
#     log:
#         "logs/prepare_cellxgene.log"
#     script:
#         "../scripts/prepare_cellxgene.R"