output_log <- file(snakemake@log[[2]], open="wt")
error_log <- file(snakemake@log[[1]], open="wt")
sink(output_log, type = "output")
sink(error_log, type = "message")

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(scDblFinder))
set.seed(123)

min_cells <- snakemake@config[["min_cells"]]

counts <- Read10X(file.path(dirname(snakemake@input[[1]]),
                            snakemake@params[["mtx_location"]]))

if (length(names(counts)) > 1) {
  counts <- counts$`Gene Expression`
}
seu <- CreateSeuratObject(counts,
                          min.cells = min_cells)
sce <- scDblFinder(GetAssayData(seu, slot = "counts"))
saveRDS(colData(sce), file = snakemake@output[[1]])

sink()
sink()