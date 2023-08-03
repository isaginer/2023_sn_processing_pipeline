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