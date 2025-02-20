# make_plots_integration.R
output_log <- file(snakemake@log[[2]], open="wt")
error_log <- file(snakemake@log[[1]], open="wt")
sink(output_log, type = "output")
sink(error_log, type = "message")

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(qs))
set.seed(123)


integrated <- qread(snakemake@input[[1]])

#plot 1 UMAP with clusters
plt1 <- DimPlot(integrated)

#plot 2 UMAP with low quality cells
plt2 <- DimPlot(integrated, group.by = "Unknown_celltype")

#plot 3 UMAP with cell types predictions
plt3 <- DimPlot(integrated, group.by = "qc_sctype_prediction")

png(snakemake@output[[1]],
    height = 5, width = 6, units = "in", res = 300)
print(plt1)
dev.off()

png(snakemake@output[[2]],
    height = 5, width = 6, units = "in", res = 300)
print(plt2)
dev.off()

png(snakemake@output[[3]],
    height = 5, width = 6, units = "in", res = 300)
print(plt3)
dev.off()

sink()
sink()