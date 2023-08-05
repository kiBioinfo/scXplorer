# #!/usr/bin/env Rscript --vanilla
#
# ##Check to see if necessary packages are installed
# #CRAN packages
# cran.packages <- c("optparse", "yaml", "igraph", "Rtsne")
#
# cran.package.check <- lapply(cran.packages, FUN = function(x) {
#   if (!require(x, character.only = TRUE)) {
#     install.packages(x, dependencies = TRUE)
#   }
# })
#
# #Bioconductor packages
# bioc.packages <- c("singleCellTK",  "BiocParallel")
#
# bioc.package.check <- lapply(bioc.packages, FUN = function(x) {
#   if (!require(x, character.only = TRUE)) {
#     BiocManager::install(x)
#   }
# })
#
# devtools::install_github('cole-trapnell-lab/monocle3')
#
# ## Function to parse arguments from yaml file
# .parseConfig <- function(sctkConfig, arguments) {
#   for (i in seq_along(arguments)) {
#     arg <- arguments[i]
#     assign(arg, sctkConfig[[arg]], envir = parent.frame())
#   }
# }
