Plot_Markers_UI<- function(id) {
  ns <- NS(id)
  tagList(
    fluidPage(
      useShinyjs(),
      fluidRow(style = "height: 78vh; overflow-y: auto;",
               column(width = 10,
                      wellPanel(
                        plotOutput(ns('marker_gene'), height = "800px")%>%
                          withSpinner(color="#0dc5c1",type = 6,size=0.9),
                        fluidRow(
                          column(width = 12, align = "right",
                          download_plot_UI(ns("marker_gene_p"))
                          )),
                          style = "height: auto; border: 3px solid #CEECF5; ;"

                      )
               ),
               column(width =2,
                       wellPanel(style = " background-color: #CEECF5; border: 3px solid #CEECF5;",

                                 selectizeInput( ns("select_method"),
                                              label = "Plot Markers By:",
                                              choices = c("Expression", "Reduction"), multiple = F,
                                              selected = NULL,
                                              options = list(placeholder = "Choose a value..."))
                                 ,
                                 selectizeInput(ns("features"), choices =NULL,
                                                label = "Select Markers of Interest:", multiple =TRUE),
                                 uiOutput(ns("genes")),
                                 conditionalPanel(
                                   condition  = "input.select_method == 'Expression' ",
                                   uiOutput(ns("colour"))
                                   ,ns = NS(id)
                                 ),

                                 conditionalPanel(
                                   condition = "input.select_method == 'Reduction' " ,
                                   uiOutput(ns("dimred")),

                                   ns = NS(id)
                                 ),
                                 # tags$br(),
                                 shinyjs::hidden(actionBttn(ns("toAnalyze"), "Next Cell Type",style = "unite",color = "royal",
                                                            icon = icon("angles-right",class="fa-duotone fa-angles-right")))
                       )
                       )
      )
    )
  )
}

Plot_Markers_Server <-function(id, DGE) {
  moduleServer(
    id,
    ## Below is the module function
    function(input, output, session) {

      ns <- session$ns

      scdata<-reactive({
        req(DGE())
        data <- DGE()
        return(data)
      })



     features_to_plot<-reactive({
       req(scdata())
       updateSelectizeInput( session, "features",

                             choices = rownames(scdata()), selected = head(rownames(scdata()),1), server  = TRUE)
     })


      output$colour <- renderUI({
        req(scdata())



        selectInput( ns("colour_by"),
                     label = "Colour By:",
                     choices = colnames(scdata()@colData), selected = rev(colnames(scdata()@colData))[1], selectize = F)
      })
      output$dimred <- renderUI({
        req(scdata())

        selectInput( ns("dimred_by"),
                     label = "Dim Reduction:",
                     choices = reducedDimNames(scdata()), selected = rev(reducedDimNames(scdata()))[1], selectize = F)
      })


        output$marker_gene <- renderPlot({
          req(scdata())

          req(features_to_plot())
          req(input$features)
          req(input$colour_by)
          req(input$dimred_by)

          p=plot_Markers(scdata(), method = input$select_method, colour_by = input$colour_by,
                           dimred = input$dimred_by, features = input$features,
                           by_exprs_values = rev(assayNames(scdata()))[1])
          shinyjs::show('toAnalyze')
          p
          }, res = 96)
        observeEvent(input$toAnalyze,
                     {
                       shinyjs::runjs("$('a[data-value=\"Cell_type\"]').tab('show');")
                     })

    }
  )
}
