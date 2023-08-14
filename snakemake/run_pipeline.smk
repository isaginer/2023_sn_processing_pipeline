from snakemake.utils import min_version
min_version("6.0")

include: "rules/load_samples.smk"
include: "rules/00_cellranger.smk"
include: "rules/01_decontx.smk"
include: "rules/02_doublets.smk"
include: "rules/03_qc_filtering.smk"
include: "rules/04_integration.smk"


all_input_list = [rules.cellranger_all.input,
                  rules.doublets_all.input,
                  rules.qc_all.input,
                  rules.integration_all.input]

if (not(config.get("skip_soup", False))):
    all_input_list.append([rules.decontx_all.input])


rule all:
    input:
        all_input_list
    default_target: True