#!/usr/bin/bash

# Link AUS snRNA-seq data
incoming_dir="/home/data/WashU_Data/snRNAseq_Data/\
01.-Incoming_Data/2022_11_03_AUS_brain_bank/01.-MGIorigData/\
htcf.wustl.edu/files/5d9wJzXO/Karch_MGI3522_10X"

mkdir -p data/{fastq,cellranger}
ln -s ${incoming_dir}/FASTQ/*.gz data/fastq/
ln -s ${incoming_dir}/MGI* data/cellranger/