#' @export
#' sceExplorer
sceExplorer <- function() {
  appDir <- system.file("shiny", package = "sceExplorer")

  shiny::runApp(appDir, display.mode = "normal", launch.browser = TRUE)
}
