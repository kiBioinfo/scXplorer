#' @export
#' sceExplorer
sceExplorer <- function() {
  appDir <- system.file("shiny", package = "sceExplorer")
  shiny::runApp(appDir, display.mode = "normal", launch.browser = TRUE)
}

#' @export
#' install_bioconductor_dependencies
install_bioconductor_dependencies <- function() {
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
}


#' @export
#' install_cytotrace_package
install_cytotrace_package <- function() {

  #Check Installed packages
  installed_packages <- installed.packages()
  if(!'CytoTRACE' %in% installed_packages){
    if(!'reticulate' %in% installed_packages){
      install.packages("reticulate")
    }
    #Install python packages:
    reticulate::virtualenv_create(envname = "r-reticulate")
    reticulate::virtualenv_install(envname = "r-reticulate", packages = c("scanoramaCT", "numpy"))
    # Step 1: Download the file
    download_url <-  "https://cytotrace.stanford.edu/CytoTRACE_0.3.3.tar.gz"
    download_file <- "CytoTRACE_0.3.3.tar.gz"

    # Set the timeout value (in seconds)
    options(timeout=300)

    # Download the file with an increased timeout value
    download.file(download_url, download_file)

    # Step 2: Install the package
    install.packages(download_file, repos = NULL, type = "source")

    # Step 3: Clean up downloaded file
    file.remove(download_file)
  }
  else{
    #Print message if all packages are already installed
    message("All requirements are satisfied!")
  }
}
