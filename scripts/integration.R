suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(qs))
set.seed(123)
options(future.globals.maxSize = 150000 * 1024^2)

# TODO
# universal solution?
add_metadata <- function(x, metadata) {
    seu <- readRDS(x)
    pheno_name <- seu@project.name
    pheno_name <- gsub("-1$", "", gsub("^MGI3522_CKAE-", "", pheno_name))
    # That doesn't work
    # seu <- AddMetaData(seu, metadata[pheno_name, ])
    for (col_ in colnames(metadata)) {
        seu@meta.data[, col_] <- metadata[pheno_name, col_]
    }
    seu
}

metadata <- as.data.frame(read_excel(snakemake@input[["phenodata"]]))
rownames(metadata) <- lapply(metadata[, 1],
                            function(x) {
                                gsub("/", "-", gsub(" BA5$", "", x))
                            })
colnames(metadata) <- str_replace_all(colnames(metadata), " ", "_")

normalized_samples <- snakemake@input[["samples_rds"]]
names(normalized_samples) <- snakemake@params$samples

seu.list <- lapply(normalized_samples, function(x) {add_metadata(x, metadata)})

integration_nfeatures <- snakemake@config[["integration_nfeatures"]]
integration_ndims <- snakemake@config[["integration_ndims"]]
processing_ndims <- snakemake@config[["processing_ndims"]]
integration_reference <- snakemake@config[["integration_reference"]]
clustering_res <- snakemake@config[["clustering_res"]]

merged.features <- SelectIntegrationFeatures(object.list = seu.list, nfeatures = integration_nfeatures)
qsave(merged.features, file = snakemake@output[["features"]])

merged <- PrepSCTIntegration(object.list = seu.list, anchor.features = merged.features,
                            verbose = FALSE)

merged <- lapply(X = merged, FUN = function(x) {
  x <- RunPCA(x, features = merged.features)
})

merged.anchors <- FindIntegrationAnchors(object.list = merged, normalization.method = "SCT",
                                        reference = integration_reference, dims = 1:integration_ndims,
                                        anchor.features = merged.features, verbose = FALSE, reduction = 'rpca')
qsave(merged.anchors, file = snakemake@output[["anchors"]])

merged.integrated <- IntegrateData(anchorset = merged.anchors, normalization.method = "SCT",
                                  verbose = FALSE, dims = 1:integration_ndims)
merged.integrated <- RunPCA(merged.integrated, verbose = FALSE) %>%
        RunUMAP(dims = 1:processing_ndims, verbose = FALSE) %>%
        FindNeighbors(dims = 1:processing_ndims, verbose = FALSE) %>%
        FindClusters(verbose = FALSE, resolution = clustering_res)

qsave(merged.integrated, file = snakemake@output[["integrated"]])