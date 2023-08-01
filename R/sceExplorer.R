#' @export
#' sceExplorer
sceExplorer <- function() {
  appDir <- system.file("shiny", package = "sceExplorer")

  shiny::runApp(appDir, display.mode = "normal", launch.browser = TRUE)
}

#' @export
#'
install_bioconductor_dependencies <- function() {
  packages <- c("BASiCS", "celldex", "CellMixS", "clusterProfiler",
                "DEsingle", "destiny", "DropletUtils", "enrichplot",
                "kBET", "Linnorm", "monocle3", "org.Hs.eg.db", "org.Mm.eg.db",
                "PCAtools", "rrvgo", "RUVSeq", "scater", "SCnorm", "singleCellTK", "SingleR", "slingshot","edgeR")

  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  installed_packages <- installed.packages()

  for (package in packages) {
    if (!package %in% installed_packages) {
      BiocManager::install(package)
    }
  }
}

