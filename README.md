# Single nuclei RNA seq processing pipeline (2023 version)

Managed by [@aladyeva-wustl](https://github.com/aladyeva-wustl)

Repository with single nuclei RNAseq initial (QC + integration) data analysis pipeline
## Stages:
1. "Soup" removal
2. Doublets identification
3. QC filtering
4. Preprocessing

## Prerequisites:
### Create conda environment (use mamba for faster installation)
    CONDA_ENV=<conda_name>
    conda env create -n $CONDA_ENV --file snakemake/envs/pipeline_env.yml
    conda activate $CONDA_ENV
    Rscript snakemake/envs/install_packages.R

### Populate data (create symlinks)

    chmod +X scripts/populate_data.sh
    ./scripts/populate_data.sh

### Edit config.yaml
    vim snakemake/config/config.yaml

### Run pipeline
Dry-run
    snakemake -pr -s snakemake/run_pipeline.smk --configfile snakemake/config/config.yml \
        -c 4 --use-conda --conda-frontend mamba --rerun-incomplete --scheduler=greedy --dry-run

Real run
    snakemake -pr -s snakemake/run_pipeline.smk --configfile snakemake/config/config.yml \
        -c 4 --use-conda --conda-frontend mamba --rerun-incomplete --scheduler=greedy