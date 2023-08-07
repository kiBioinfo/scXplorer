#!/usr/bin/env Rscript

#  devtools::install_github("satijalab/seurat", quiet = TRUE, force = T)
# #devtools::install_github(lib='',repo = 'cole-trapnell-lab/monocle3' , ref = 'release/1.3.1',  force = T, dependencies = T)
# devtools::install_github('cole-trapnell-lab/monocle3', force = T, dependencies = T,build = T)
# devtools::install_github('theislab/kBET', quiet = TRUE, force = T)

if(!require(shiny)){
  install.packages("shiny")
  library(shiny)
}
if(!require(purrr)){
  install.packages("purrr")
  library(purrr)
}
if(!require(shinyBS)){
  install.packages("shinyBS")
  library(shinyBS)
}
if(!require(shinycssloaders)){
  install.packages("shinycssloaders")
  library(shinycssloaders)
}
if(!require(DT)){
  install.packages("DT")
  library(DT)
}
if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}
if(!require(Matrix)){
  install.packages("Matrix")
  library(Matrix)
}
if(!require(tidyverse)){
  install.packages("tidyverse")
  library(tidyverse)
}



if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}



