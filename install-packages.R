#!/usr/bin/env Rscript

# devtools::install_github(lib='',repo = 'satijalab/seurat', ref = 'release/3.0', force=T)
devtools::install_github("satijalab/seurat", "seurat5", quiet = TRUE, force = T)

packages <- c("BASiCS", "celldex", "CellMixS", "clusterProfiler", "rrvgo",
              "DEsingle", "destiny", "DropletUtils", "enrichplot", "RUVSeq",
              "kBET", "Linnorm", "monocle3", "org.Hs.eg.db", "org.Mm.eg.db",
              "PCAtools", "rrvgo", "RUVSeq", "scater", "SCnorm", "singleCellTK",
              "SingleR", "slingshot","edgeR", "SingleCellExperiment", "scater",
              "singleCellTK", "slingshot", "scuttle", "scran", "SCnorm",
              "PCAtools", "org.Mm.eg.db", "org.Hs.eg.db" , "Linnorm", "enrichplot",
              "DropletUtils", "destiny","DEsingle", "clusterProfiler", "CellMixS",
              "celldex", "limma", "SummarizedExperiment", "sva", "DESeq2", "batchelor",
              "MAST")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

installed_packages <- installed.packages()

for (package in packages) {
  if (!package %in% installed_packages) {
    BiocManager::install(package)
  }
}

if(!require(limma)){
  BiocManager::install("limma")
  library(limma)
}
