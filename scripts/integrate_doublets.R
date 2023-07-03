suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(qs))
set.seed(123)
options(future.globals.maxSize = 150000 * 1024^2)

processed_samples <- snakemake@input[[1]]
names(processed_samples) <- snakemake@params$samples

seu.list <- lapply(processed_samples, function(x) {
    seu <- readRDS(x)
})

integration_nfeatures <- snakemake@config[["integration_nfeatures"]]
integration_ndims <- snakemake@config[["integration_ndims"]]
processing_ndims <- snakemake@config[["processing_ndims"]]
clustering_res <- snakemake@config[["clustering_res"]]

merged.features <- SelectIntegrationFeatures(object.list = seu.list, 
                                             nfeatures = integration_nfeatures)
qsave(merged.features, file = snakemake@output[["features"]])

merged <- PrepSCTIntegration(object.list = seu.list, anchor.features = merged.features,
                            verbose = FALSE)

merged <- lapply(X = merged, FUN = function(x) {
  x <- RunPCA(x, features = merged.features)
})

integration_reference <- snakemake@config[["integration_reference"]]
if (("integration_reference" %in% names(snakemake@params)) & (snakemake@config["doublets_integrated"])) {
    merged.anchors <- FindIntegrationAnchors(object.list = merged,
                                            normalization.method = "SCT",
                                            reference = integration_reference,
                                            dims = 1:integration_ndims,
                                            anchor.features = merged.features,
                                            verbose = FALSE, reduction = "rpca")
} else {
    merged.anchors <- FindIntegrationAnchors(object.list = merged,
                                             normalization.method = "SCT",
                                             dims = 1:integration_ndims, 
                                             anchor.features = merged.features,
                                             verbose = FALSE,
                                             reduction = "rpca")
}

qsave(merged.anchors, file = snakemake@output[["anchors"]])

merged.integrated <- IntegrateData(anchorset = merged.anchors,
                                   normalization.method = "SCT",
                                   verbose = FALSE, dims = 1:integration_ndims)
merged.integrated <- RunPCA(merged.integrated, verbose = FALSE) %>%
        RunUMAP(dims = 1:processing_ndims, verbose = FALSE) %>%
        FindNeighbors(dims = 1:processing_ndims, verbose = FALSE) %>%
        FindClusters(verbose = FALSE, resolution = clustering_res)

qsave(merged.integrated, file = snakemake@output[["integrated"]])