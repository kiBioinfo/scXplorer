batchCorrct_UI<-  function(id) {
  ns <- NS(id)
  tagList(
    tabItem(tabName= "batch_correction",
            fluidPage(
              fluidRow(style = "height: 78vh; overflow-y: auto;",

              column(width = 10,
                    wellPanel(style = "border: 2px solidx blue;background-color: auto;",
                              h1(
                    tags$p("Batch effect correction on dimension reduction takes longer time")
                     ),
                    style = " background-color: #CEECF5; border: 3px solid #CEECF5;"
                    )
                     ),
              tags$br(),
              column(width = 2,

                     wellPanel(
                     #fileInput(ns("filter_data"),label = "Please Upload filterd .rds file", accept = ".rds",width = '100%'),
                     selectInput(ns("batch_correction_type"), label = "Batch correction on :", choices = c("Expression matrix" ="exprsn_mtrx", "Dimension reduction" = "dim_rdxn"), selected = "exprsn_mtrx"),
                     conditionalPanel(
                       condition = "input.batch_correction_type =='exprsn_mtrx' ",
                       selectInput(ns("corxn_method_exprsn"), label = "Choose batch effect correction method:", choices = c("BASiCS" = "bas",
                                                                                                                            "Limma"= "lima",
                                                                                                                            "ComBat" = "cbt",
                                                                                                                            "ComBat-seq" = "cbtseq",
                                                                                                                            "Mutual Nearest Neighbors (MNNs)" = "mnn",
                                                                                                                            "Scanorama" = "sca"), selected = "sca" ), ns = NS(id)

                     ),

                     conditionalPanel(
                       condition = "input.batch_correction_type =='dim_rdxn' ",
                       selectInput(ns("corxn_method_rdxn"), label = "Choose batch effect correction method:", choices = c("BEER" = "beer",
                                                                                                                            "Canonical Correlation Analysis (CCA)"= "cca",
                                                                                                                            "Batchelor_fastMNN" = "dmnn",
                                                                                                                            "Batchelor_harmony" = "hmy",
                                                                                                                            "LIGER" = "lig" ), selected = "dmnn" ), ns = NS(id)


                     ),#2nd conditionalPanel
                     selectInput(ns("expression_data"), label = "Choose count matrix", choices = c("Use Raw Counts" = "counts",
                                                                                                   "Use Normalized Counts"= "NMcounts"), selected = "NMcounts"),
                     actionBttn(ns("correct_batch"),label="EXEC",style = "jelly",color = "success",icon = icon("sliders")),
                     tags$br(),
                     tags$br(),
                     shinyjs::hidden(downloadBttn(ns('batch_corrct_data'), '.rds batch corrected data', color = "royal")),
                     style = " background-color: #CEECF5; border: 3px solid #CEECF5;"

                     )
              ) #input column


              ) #fluidPage
            )# main fluidPage
    ) # main TabItem
  )#tagList
}#UI_function


batchCorrct_Server <- function(id,normalization_data) {
  moduleServer(
    id,
    ## Below is the module function
    function(input, output, session) {

      vals=reactiveValues()
      #sce object
      scdata<-reactive({
        req(normalization_data())
        obj <- normalization_data()
        return(obj)

      })#scdata

      #batch correction
      batch_effect <- reactive({
        req(scdata())
        req(input$batch_correction_type)
        showNotification(paste("batch effect correction using", input$batch_correction_type ) , type = "message", duration = 3)
        if(input$batch_correction_type== "exprsn_mtrx"){
          data= BatchEffect_Matrix(scdata(), method = input$corxn_method_exprsn,used = 'VGcounts')
          data= scater::runPCA(data, exprs_values='BEVGcounts', altexp= 'BEVGcounts', name= "BEPCA")
          vals$data =data
          data
        }
        else if (input$batch_correction_type== "dim_rdxn") {
         data=BatchEffect_DimReduction(scdata(), method = input$corxn_method_rdxn, used = input$expression_data)
         vals$data =data
         data
        }

        return(data)
      })

      observeEvent(input$correct_batch,{
        req(batch_effect())
        shinyjs::show("batch_corrct_data")

        if(input$batch_correction_type=="dim_rdxn"){
          batch_effect()

          shinyjs::show(selector = "a[data-value=\"batch_evaluation\"]")
          shinyjs::show(selector = "a[data-value=\"non_linear_batch\"]")
          shinyjs::hide(selector = "a[data-value=\"variable_gene_batch\"]")
          shinyjs::hide(selector = "a[data-value=\"linear_batch\"]")
          shinyjs::runjs("$('a[data-value=\"batch_evaluation\"]').tab('show');")
        }
        else{
          batch_effect()
          shinyjs::show(selector = "a[data-value=\"batch_evaluation\"]")

          shinyjs::show(selector = "a[data-value=\"variable_gene_batch\"]")
          shinyjs::show(selector = "a[data-value=\"linear_batch\"]")
          shinyjs::hide(selector = "a[data-value=\"non_linear_batch\"]")
          shinyjs::runjs("$('a[data-value=\"batch_evaluation\"]').tab('show');")
        }
      })


      output$batch_corrct_data <- downloadHandler(
        filename = function() {
          paste0( "Batch_corrected",input$batch_correction_type, "-", Sys.Date(), ".rds")
        },
        content = function(file) {
          saveRDS(batch_effect(), file = file)
        }
      )


      return(batch_effect)

    }# I/O function
  )#ModuleServer
}#main function


### Server function to switch batch correction method between expression/reduction

rm_dim_data<-function(id,  batch_corrct_dmrd, batch_corrct_exprsn){
  moduleServer(
    id,
    function(input, output, session){
      ns<- session$ns

      data<- reactive(
        if(input$batch_correction_type== "exprsn_mtrx"){
          return(batch_corrct_exprsn())
        }
        else if(input$batch_correction_type=="dim_rdxn"){
          return(batch_corrct_dmrd())
        }
      )

      return(data)
    }
  )

}
