output_log <- file(snakemake@log[[2]], open="wt")
error_log <- file(snakemake@log[[1]], open="wt")
sink(output_log, type = "output")
sink(error_log, type = "message")

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(sceasy))
suppressPackageStartupMessages(library(reticulate))
suppressPackageStartupMessages(library(qs))
set.seed(123)

merged.integrated <- qread(snakemake@input[[1]])
sceasy::convertFormat(merged.integrated, from = "seurat", to = "anndata",
                       outFile = snakemake@output[[1]])

sink()
sink()