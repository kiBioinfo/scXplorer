.onLoad <- function(libname, pkgname) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db", "limma"))
}

devtools::install_github('theislab/kBET')
