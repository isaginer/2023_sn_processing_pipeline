output_log <- file(snakemake@log[[2]], open="wt")
error_log <- file(snakemake@log[[1]], open="wt")
sink(output_log, type = "output")
sink(error_log, type = "message")

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(celda))
suppressPackageStartupMessages(library(DropletUtils))
set.seed(123)

min_cells <- snakemake@config[["min_cells"]]
min_features <- snakemake@config[["min_features"]]
DECONTX_DIR <- snakemake@config["decontx_dir"]
SAMPLE <- snakemake@wildcards$sample

counts <- Read10X(snakemake@input[[1]])
counts.raw <- Read10X(snakemake@input[[2]])

if (length(names(counts)) > 1) {
  counts <- counts$`Gene Expression`
}

if (length(names(counts.raw)) > 1) {
  counts.raw <- counts.raw$`Gene Expression`
}

# Create a SingleCellExperiment object and run decontX
sce <- SingleCellExperiment(list(counts = counts))
sce.raw <- SingleCellExperiment(list(counts = counts.raw))
sce <- decontX(sce, background = sce.raw)
counts_new <- round(decontXcounts(sce))
saveRDS(CreateSeuratObject(counts_new),
        file.path(snakemake@output[[2]]))

nfeatures <- apply(counts_new, 2, sum)
nfeatures <- nfeatures[nfeatures > 0]
min.features <- round(quantile(nfeatures, seq(0, 1, 0.05))["5%"])
if (min.features < min_features) {
  min.features <- min_features
}

# Create a Seurat object from a SCE with decontX results
seu <- CreateSeuratObject(counts_new,
                        min.cells = min_cells,
                        min.features = min.features)

outs_path <- file.path(DECONTX_DIR, SAMPLE, "outs")
dir.create(outs_path, showWarnings = FALSE, recursive = TRUE)
write10xCounts(x = seu@assays$RNA@counts,
                overwrite = TRUE,
                version = "3",
                path = file.path(outs_path, "filtered_feature_bc_matrix"))

cat(NULL, file = snakemake@output[[1]])

sink()
sink()