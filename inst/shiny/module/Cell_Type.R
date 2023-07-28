cell_type_analysis_UI<-function(id) {
  ns <- NS(id)
  tagList(
    fluidPage(
      useShinyjs(),
      fluidRow(style = "height: 78vh; overflow-y: auto;",
               column(width = 10,
                      wellPanel(
                      plotOutput(ns('cell_type'), height = "800px")%>% withSpinner(color="#0dc5c1",type = 6,size=0.9),
                      style = "height: auto; border: 3px solid #CEECF5; text-align: end;",
                      download_plot_UI(ns("Cell_type_p"))
                      )
                      ),


               column(width = 2,
                      wellPanel(
                        selectInput(ns('method'),
                                    label = "Chose cell type  annotation method:",
                                    choices = c("SingleR","sctype"), selected = "SingleR"),
                        conditionalPanel(
                          condition ="input.method=='SingleR'",
                          selectInput(ns("singleR"),
                                      label = "Choose the datasets :",
                                      choices = c("BlueprintEncodeData","DatabaseImmuneCellExpressionData",
                                                   "HumanPrimaryCellAtlasData","ImmGenData","MonacoImmuneData",
                                                   "MouseRNAseqData"), selected='BlueprintEncodeData'),  ns = NS(id)
                        ),
                        conditionalPanel(
                          condition ="input.method=='sctype'",
                            uiOutput(ns("scType")),ns = NS(id)
                        ),
                        uiOutput(ns("clust")),
                        actionBttn(ns("cell_bttn"),label="EXEC",style = "jelly",color = "success",icon = icon("sliders")),
                        tags$br(),
                        tags$br(),
                        shinyjs::hidden(actionBttn(ns("toAnalyze"), "Next Cell Development",style = "unite",color = "royal",
                                                   icon = icon("angles-right",class="fa-duotone fa-angles-right"))),
                        style = " background-color: #CEECF5; border: 3px solid #CEECF5;"
                        ) )



      )# fluidrow
    )#fluidPage
)#tagList
}#UI



cell_type_analysis_Server <-function(id,DGE) {
  moduleServer(
    id,
    ## Below is the module function
    function(input, output, session) {

      ns <- session$ns
      vals=reactiveValues()
      scdata <-reactive({
        req(DGE())
        obj <- DGE()
        return(obj)
      })
      output$clust<-renderUI({
        req(scdata())
        ns <- session$ns
        selectInput(ns("reduction"),
                    label = "select the dim reduction to plot",
                    choices = reducedDimNames(scdata()), selected = rev(reducedDimNames(scdata()))[1])
      })
      output$scType <- renderUI({
        sctype_data<-read_excel("scType/ScTypeDB_full.xlsx")
        selectizeInput(ns("SCtype"), label= "Choose the cell types:", choices=sctype_data$tissueType, selected="Immune system", multiple=T)
      })

      Cell_type<-reactive({
        req(scdata())
        req(input$cell_bttn)
        withProgress(message = 'Cell Type Annotation in progress...',
                     detail = 'This may take a while...', value = 0, {
          # for(i in 1:10) {
          #   incProgress(1/10)
          #   Sys.sleep(0.5)
          # }
        data= Cell_Type_Annotation(data= scdata(), dataset= input$singleR , TissueType= input$SCtype, method=input$method )
        data
        # setProgress(value = 1)
        })
      })

#
#       observeEvent(input$cell_bttn,{
#         withProgress(message = 'Calculation in progress', value = 0, {
#           for(i in 1:10) {
#             incProgress(1/10)
#             Sys.sleep(0.5)
#           }
#           req(Cell_type())
#
#         })
#
#       })
      output$cell_type <-renderPlot({
        req(Cell_type())
        p=plotReducedDim(Cell_type(),dimred = input$reduction, colour_by = "cellType", text_by = "cellType")
        shinyjs::show('toAnalyze')
        p
      }, res = 96)

      observeEvent(input$toAnalyze,
                   {
                     shinyjs::runjs("$('a[data-value=\"Cell_develope\"]').tab('show');")
                   })
      download_plot_Server("Cell_type_p", input_data = Cell_type , name ="Cell_type_plot")


      return(Cell_type)

    }#2nd Fun
  ) #module server
}#server
