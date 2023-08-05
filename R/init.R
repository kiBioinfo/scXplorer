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
if (!requireNamespace("kBET", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::install_github("theislab/kBET")
}
#devtools::install_github('theislab/kBET')

## Function to parse arguments from yaml file
.parseConfig <- function(sctkConfig, arguments) {
  for (i in seq_along(arguments)) {
    arg <- arguments[i]
    assign(arg, sctkConfig[[arg]], envir = parent.frame())
  }
}
