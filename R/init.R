.onAttach <- function(libname, pkgname) {
  if(!requireNamespace("monocle3", quietly = TRUE)) {
    install.packages("remotes")
    remotes::install_github('cole-trapnell-lab/monocle3')
  }
}
