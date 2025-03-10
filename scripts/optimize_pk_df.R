output_log <- file(snakemake@log[[2]], open="wt")
error_log <- file(snakemake@log[[1]], open="wt")
sink(output_log, type = "output")
sink(error_log, type = "message")

library(Seurat)
library(data.table)
library(dplyr)
library(DoubletFinder)
set.seed(123)

rename_columns <- function(pattern, columns_, new_name) {
    columns_[grepl(pattern, columns_)] <- new_name
    return(columns_)
}


min_cells <- snakemake@config[["min_cells"]]
dims_ <- snakemake@config[["doublet_dims"]]


counts <- Read10X(file.path(dirname(snakemake@input[[1]]),
                            snakemake@params[["mtx_location"]]))

if (length(names(counts)) > 1) {
  counts <- counts$`Gene Expression`
}
seu <- CreateSeuratObject(counts,
                          min.cells = min_cells) %>%
    SCTransform() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:dims_)

# Optimization of parameters
sweep_res_list <- paramSweep(seu, PCs = 1:dims_, sct = TRUE)
sweep_stats <- summarizeSweep(sweep_res_list, GT = FALSE)
bcmvn <- find.pK(sweep_stats)
saveRDS(bcmvn, file = snakemake@output[[1]])

sink()
