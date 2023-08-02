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
                       selectInput(ns("corxn_method_exprsn"),
                                   label = "Choose batch effect correction method:", choices = c("BASiCS" = "bas",
                                                                                                 "Limma"= "lima",
                                                                                                 "ComBat" = "cbt",
                                                                                                 "ComBat-seq" = "cbtseq",
                                                                                                 "Mutual Nearest Neighbors (MNNs)" = "mnn",
                                                                                                  "Scanorama" = "sca"), selected = "sca" ),
                       selectInput(ns("expression_data_exprsn"), label = "Choose count matrix",
                                   choices = c("Use Raw Counts" = "counts",
                                               "Use Normalized Counts"= "NMcounts",
                                               "Use Variable Genes Counts" = 'VGcounts'), selected = "VGcounts"),

                       ns = NS(id)

                     ),
                     conditionalPanel(
                       condition = "input.expression_data_exprsn =='counts' || input.expression_data_exprsn =='NMcounts' ",

                       numericInput(ns("ngene"), "Select No. Of Features :",2000, min=0, max=Inf),

                       ns = NS(id)

                     ),


                     conditionalPanel(
                       condition = "input.batch_correction_type =='dim_rdxn' ",
                       selectInput(ns("corxn_method_rdxn"),
                                   label = "Choose batch effect correction method:", choices = c("BEER" = "beer",
                                                                                                 "Canonical Correlation Analysis (CCA)"= "cca",
                                                                                                 "Batchelor_fastMNN" = "dmnn",
                                                                                                 "Batchelor_harmony" = "hmy",
                                                                                                  "LIGER" = "lig" ), selected = "dmnn" ),
                       selectInput(ns("expression_data"), label = "Choose count matrix",
                                   choices = c("Use Raw Counts" = "counts",
                                               "Use Normalized Counts"= "NMcounts"), selected = "NMcounts"),

                       ns = NS(id)


                     ),#2nd conditionalPanel

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
        errors <- c()
        req(input$batch_correction_type)
        tryCatch(
          {
            if(input$batch_correction_type== "exprsn_mtrx"){
              sce= BatchEffect_Matrix(scdata(), method = input$corxn_method_exprsn,used = input$expression_data_exprsn)
              if(input$expression_data_exprsn =='counts'){
                sce <- VariableGene_hvg(sce, used = 'BEcounts', ngene = input$ngene, batch = T)
                sce <- scater::runPCA(sce, exprs_values='BEVGcounts', altexp= 'BEVGcounts', name= "BEPCA")
              }
              else if (input$expression_data_exprsn=='NMcounts'){
                sce <- VariableGene_hvg(sce, used = 'BENMcounts', ngene = input$ngene, batch = T)
                sce <- scater::runPCA(sce, exprs_values='BENMVGcounts', altexp= 'BENMVGcounts', name= "BEPCA")
              }
              else if (input$expression_data_exprsn== 'VGcounts'){
                sce <- scater::runPCA(sce, exprs_values='BEVGcounts', altexp= 'BEVGcounts', name= "BEPCA")
              }

            }
            else if (input$batch_correction_type== "dim_rdxn") {
              sce <- BatchEffect_DimReduction(scdata(), method = input$corxn_method_rdxn, used = input$expression_data)
              # return(sce)
            }
            return(sce)
          },
          error = function(e) {
            errors <- c(errors, "Use another method")
            return(errors)
          }
        )
      })

      observeEvent(input$correct_batch,{
        withProgress(message = 'Batch Correction in progress', value = 0, {
          # for(i in 1:10) {
          #   incProgress(1/10)
          #   Sys.sleep(0.5)
          # }
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
