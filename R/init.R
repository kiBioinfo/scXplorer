.onLoad <- function(libname, pkgname) {
  packages <- c("BASiCS", "celldex", "CellMixS", "clusterProfiler", "rrvgo",
                "DEsingle", "destiny", "DropletUtils", "enrichplot", "RUVSeq",
                 "Linnorm", "org.Hs.eg.db", "org.Mm.eg.db",
                "PCAtools", "rrvgo", "RUVSeq", "scater", "SCnorm",
                "SingleR", "slingshot","edgeR", "SingleCellExperiment", "scater",
                "slingshot", "scuttle", "scran", "SCnorm",
                "PCAtools", "org.Mm.eg.db", "org.Hs.eg.db" , "Linnorm", "enrichplot",
                "DropletUtils", "destiny","DEsingle", "clusterProfiler", "CellMixS",
                "celldex", "limma", "SummarizedExperiment", "sva", "DESeq2", "batchelor",
                "MAST", "HDF5Array", "ggrastr" , "terra", "S4Vectors", "lme4",
                'BiocGenerics', 'DelayedArray', 'DelayedMatrixStats', 'bluster')

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  installed_packages <- installed.packages()
  for (package in packages) {
    if (!package %in% installed_packages) {
      BiocManager::install(package)

    }
  }
}

#devtools::install_github('theislab/kBET')
