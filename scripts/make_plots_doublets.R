suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
set.seed(123)

seu <- readRDS(snakemake@input[[1]])
#plot 1 density doublet vs singlet scDblFinder Log UMI
plt1 <- ggplot(seu@meta.data,
                aes(x = log10(nCount_RNA), fill = scDblFinder.class)) +
geom_density(alpha = 0.7) + theme_minimal() +
guides(fill = guide_legend(title = "scDblFinder"))

#plot 2 density doublet vs singlet DoubletFinder Log UMI
plt2 <- ggplot(seu@meta.data, aes(x=log10(nCount_RNA), fill=tolower(doublet_prc))) +
geom_density(alpha = 0.7) + theme_minimal() + 
guides(fill = guide_legend(title = "DoubletFinder"))

#plot 3 density Unknown vs other Log UMI

plt3 <- ggplot(seu@meta.data, aes(x=log10(nCount_RNA), fill = Unknown_celltype)) +
geom_density(alpha = 0.7) + theme_minimal() +
guides(fill = guide_legend(title = "Garnett \nprediction"))

#plot 4 UMAP doublet overlap
plt4 <- DimPlot(seu, group.by = "merged_doublets",
                split.by = "merged_doublets", ncol = 2)

#plot 6 UMAP doublet QC pass
plt6 <- DimPlot(seu, group.by = "pass_doublets_QC", order = "FAIL")

#plot 7 UMAP clusters + UMAP garnett prediction
plt7_1 <- DimPlot(seu, group.by = "seurat_clusters", label = TRUE) + NoLegend()
plt7_2 <- DimPlot(seu, group.by = "sctype_prediction")
plt7 <- grid.arrange(plt7_1, plt7_2, ncol = 1)

#plot 8 Low quality cells with Unknown scType prediction
plt8 <- DimPlot(seu, group.by = "Unknown_celltype")

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

png(snakemake@output[[4]],
    height = 5, width = 6, units = "in", res = 300)
print(plt4)
dev.off()

png(snakemake@output[[5]],
    height = 5, width = 6, units = "in", res = 300)
print(plt6)
dev.off()

png(snakemake@output[[6]],
    height = 10, width = 6, units = "in", res = 300)
print(grid.arrange(plt7_1, plt7_2, ncol = 1))
dev.off()

png(snakemake@output[[7]],
    height = 5, width = 6, units = "in", res = 300)
print(plt8)
dev.off()