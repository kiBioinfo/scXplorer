#!/usr/bin/env Rscript

# devtools::install_github(lib='',repo = 'satijalab/seurat', ref = 'release/3.0', force=T)
devtools::install_github("satijalab/seurat", "seurat5", quiet = TRUE, force = T)
devtools::install_github(lib='',repo = 'cole-trapnell-lab/monocle3' , ref = 'release/1.3.1',  force = T)
