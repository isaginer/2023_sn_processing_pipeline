# Single nuclei RNA seq processing pipeline (2023 version)

Managed by [@aladyeva-wustl](https://github.com/aladyeva-wustl)

Repository with single nuclei RNAseq initial (QC + integration) data analysis pipeline
## Stages:
1. "Soup" removal
2. Doublets identification
3. QC filtering
4. Preprocessing

## Prerequisites:
### Create conda environment
conda create -n <conda_name> --file envs/pipeline_env.yml
conda activate <conda_name>
Rscript scripts/installRPackages.R; touch {output.flag}

### Make 

### Populate data (create symlinks)