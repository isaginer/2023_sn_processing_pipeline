
suppressPackageStartupMessages(library(Seurat))
set.seed(123)

seu_clean <- readRDS(snakemake@input[[1]])
## NORMALIZATION

seu_clean <- SCTransform(
  seu_clean,
  ncells = min(100000, ncol(seu_clean)),
  vars.to.regress = c("subsets_Mito_percent", "subsets_Ribo_percent"),
  verbose = TRUE,
  conserve.memory = TRUE
)

## CELL CYCLE SCORING

s.genes <- stringr::str_to_title(cc.genes$s.genes)
g2m.genes <- stringr::str_to_title(cc.genes$g2m.genes)


seu_clean <- CellCycleScoring(seu_clean, s.features = s.genes,
                              g2m.features = g2m.genes, set.ident = FALSE,
                              assay = "SCT")
saveRDS(seu_clean, file = snakemake@output[[1]])