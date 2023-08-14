# make_plots_qc.R

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
set.seed(123)

nuclei_genes_hs <- c("MALAT1", "NEAT1", "FTX",
                  "FOXP1", "RBMS3", "ZBTB20", "LRMDA", "PBX1", "ITPR2", "AUTS2", "TTC28", "BNC2", "EXOC4", "RORA",
                  "PRKG1", "ARID1B", "PARD3B", "GPHN", "N4BP2L2", "PKHD1L1", "EXOC6B", "FBXL7", "MED13L",
                  "TBC1D5", "IMMP2L", "SYNE1", "RERE", "MBD5", "EXT1", "WWOX")
nuclei_genes_mm <- c("Malat1","Neat1","Ftx","Foxp1","Rbms3","Zbtb20","Lrmda",
                     "Pbx1","Itpr2","Auts2","Ttc28","Bnc2","Exoc4","Rora",
                     "Prkg1", "Arid1b", "Pard3b", "Gphn", "N4bp2l2","Pkhd1l1",
                     "Exoc6b", "Fbxl7", "Med13l", "Tbc1d5", "Immp2l", "Syne1",
                     "Rere", "Mbd5", "Ext1", "Wwox")
nuclei_genes <- c(nuclei_genes_hs,nuclei_genes_mm)

seu_sub <- readRDS(snakemake@input[["qc_step_1"]])
DefaultAssay(seu_sub) <- "RNA"

seu_sub <- PercentageFeatureSet(seu_sub, pattern = "^(MT|mt)-", col.name = "percent.mt")
nuclei_genes <- nuclei_genes[nuclei_genes %in% rownames(seu_sub)]
if (length(nuclei_genes)>0) {
  seu_sub <- PercentageFeatureSet(seu_sub, features = nuclei_genes, col.name = "percent.nuclei")  
} else {
  seu_sub$percent.nuclei <- NA
}

seu_sub$nCount_RNA_log <- log10(seu_sub$nCount_RNA)
seu_sub$nFeature_RNA_log <- log10(seu_sub$nFeature_RNA)

#plot 8 dot UMI vs nGenes
Idents(seu_sub) <- "orig.ident"
plt8 <- FeatureScatter(seu_sub, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()

#plot 9 dot MT vs Nuclei fraction (color by nFeatures)
plt9 <- ggplot(seu_sub@meta.data, aes(x=percent.nuclei, y=percent.mt, color=nFeature_RNA)) + 
geom_point(size=2) + theme_minimal()

#plot 9 dot MT vs nUMI (color by Nuclei fraction)
plt16 <- ggplot(seu_sub@meta.data, aes(x=nCount_RNA_log, y=percent.mt, color=percent.nuclei)) + 
geom_point(size=2) + theme_minimal()

#plot 10 UMAP Log UMI distribution
plt10 <- FeaturePlot(seu_sub, features = "nCount_RNA_log")

#plot 11 UMAP nGenes distribution
plt11 <- FeaturePlot(seu_sub, features = "nFeature_RNA")

#plot 12 UMAP UMI QC pass
plt12 <- DimPlot(seu_sub, group.by = "outlier_sum_log_6_lower", order = "TRUE")

#plot 13 UMAP Genes QC pass
plt13 <- DimPlot(seu_sub, group.by = "outlier_detected_log_2_lower", order = "TRUE")

#plot 14 UMAP clusters + UMAP garnett prediction
plt14_1 <- DimPlot(seu_sub, group.by = "seurat_clusters", label = T) + NoLegend()
plt14_2 <- DimPlot(seu_sub, group.by = "garnett_prediction") 
plt14 <- grid.arrange(plt14_1,plt14_2, ncol=1)

#plot 15 UMAP QC pass
plt15 <- DimPlot(seu_sub, group.by = "keep", order = "FALSE")

png(snakemake@output[[1]],
    height = 5, width = 6, units = "in", res = 300)
print(plt8)
dev.off()

png(snakemake@output[[2]],
    height = 5, width = 6, units = "in", res = 300)
print(plt9)
dev.off()

png(snakemake@output[[3]],
    height = 5, width = 6, units = "in", res = 300)
print(plt16)
dev.off()

png(snakemake@output[[4]],
    height = 5, width = 6, units = "in", res = 300)
print(plt10)
dev.off()

png(snakemake@output[[5]],
    height = 5, width = 6, units = "in", res = 300)
print(plt11)
dev.off()

png(snakemake@output[[6]],
    height = 5, width = 6, units = "in", res = 300)
print(plt12)
dev.off()

png(snakemake@output[[7]],
    height = 5, width = 6, units = "in", res = 300)
print(plt13)
dev.off()

png(snakemake@output[[8]],
    height = 10, width = 6, units = "in", res = 300)
print(grid.arrange(plt14_1,plt14_2, ncol=1))
dev.off()

png(snakemake@output[[9]],
    height = 5, width = 6, units = "in", res = 300)
print(plt15)
dev.off()
