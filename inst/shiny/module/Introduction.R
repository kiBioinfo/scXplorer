

Intro_UI<- function(id) {
  ns <- NS(id)
  tagList(
    fluidPage(
      #includeMarkdown("fun/Readme.md"),
      htmlOutput(ns("my_md_output"))
      #uiOutput(ns("my_rmd_output"))
      # tags$head(tags$script(src = "www/mermaid.min.js")),
      # tags$div(ns("flowchart"))

    )
  )
}


Intro_Server <-function(id) {
  moduleServer(
    id,
    ## Below is the module function
    function(input, output, session) {
      # Render the RMD file
      output$my_md_output <- renderUI({
        mdFileContent <-readLines("fun/FlowChart.md") %>% paste(collapse = "\n")
        HTML(markdown::markdownToHTML(mdFileContent, fragment.only = TRUE))

      })

      # output$my_rmd_output <- renderUI({
      #   rmarkdown::render("fun/FlowChart.Rmd")
      #   includeHTML("fun/FlowChart.html")
      # })




    }
)
}
