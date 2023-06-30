cell_developement_analysis_UI<-function(id) {
  ns <- NS(id)
  tagList(
    fluidPage(
      useShinyjs(),
      fluidRow(style = "height: 78vh; overflow-y: auto;",
               column(width = 10,
                      wellPanel(
                        plotOutput(ns('cell_development'), height = "800px")%>% withSpinner(color="#0dc5c1",type = 6,size=0.9),
                        style = "height: auto; border: 3px solid #CEECF5;text-align: end;",
                        download_plot_UI(ns("cell_development_P"))
                        ),

                      tags$br(),
                      shinyjs::hidden(wellPanel(id=ns("curve"),
                        plotOutput(ns('cell_development_curve'), height = "800px")%>% withSpinner(color="#0dc5c1",type = 6,size=0.9),
                        style = "height: auto; border: 3px solid #CEECF5;text-align: end;",
                        download_plot_UI(ns("cell_development_curve_P"))
                        )

                       )),




               column(width = 2,
                      wellPanel(
                        selectInput(ns('method'), label = "Select method for cell development:",choices = c("Princurve_fit", "CytoTrace",
                                                                                                            "Slingshot", "Monocle"), selected = "Princurve_fit"),
                        conditionalPanel(
                          condition ="input.method=='Princurve_fit' || input.method=='Slingshot'",
                          uiOutput(ns("reduce_dim")),ns = NS(id)
                        ),
                        conditionalPanel(
                          condition ="input.method=='Monocle'",
                          uiOutput(ns("seedby")),ns = NS(id)
                        ),
                        uiOutput(ns("colour")),

                        style = " background-color: #CEECF5; border: 3px solid #CEECF5; "

                      ) )



      )# fluidrow
    )#fluidPage
  )#tagList
}#UI



cell_developement_analysis_Server <-function(id,cell_type) {
  moduleServer(
    id,
    ## Below is the module function
    function(input, output, session) {

      ns <- session$ns
      vals=reactiveValues()
      scdata<-reactive({
        cell_type()
      })

      output$colour <- renderUI({
        selectInput(ns("colour_by"),
                    label = "Colour by:",
                    choices = colnames(colData(scdata())), selected = rev(colnames(colData(scdata())))[1])
      })

      output$reduce_dim<-renderUI({
        ns <- session$ns
        selectInput(ns("reduction"),
                    label = "select the dim reduction to plot:",
                    choices = reducedDimNames(scdata()), selected = rev(reducedDimNames(scdata()))[1])
      })
      output$seedby <- renderUI({
        if(!("cellType" %in% colnames(colData(scdata())))){
          selectizeInput(ns("cell_type"), label= "Choose the cell types:", choices=unique((scdata()@colData$label)), selected=unique((scdata()@colData$label))[1], multiple=T)
        }
        else{
          selectizeInput(ns("cell_type"), label= "Choose the cell types:", choices=unique((scdata()@colData$cellType)), selected=unique((scdata()@colData$cellType))[1], multiple=T)
        }
      })

      Cell_develpment<-reactive({
        req(input$method)
        if(input$method=="Princurve_fit"){
          shinyjs::showElement("curve")
        }
        else{
          shinyjs::hideElement("curve")
        }
       p=CellDevelopment(scdata(),method=input$method, runWith = input$reduction, seedBy = input$cell_type)
       vals$p <- p
       p

      })





          #shinyjs::enable("cell_development_curve")

          output$cell_development <-renderPlot({

            plot_cell_development()

          }, res = 96)

          output$cell_development_curve <-renderPlot({
            plot_cell_development_curve()

          }, res = 96)




      #plot cell development

      plot_cell_development <- reactive({
        Plot_Development(CS.data = Cell_develpment() ,runWith = input$reduction, colorBy = input$colour_by )
      })

      plot_cell_development_curve <-reactive({

        Plot_Development_Princurve(CS.data =Cell_develpment(), runWith = input$reduction, colorBy = input$colour_by )
      })







      #Download plot
      download_plot_Server("cell_development_P", input_data = plot_cell_development , name ="Cell_development_plot")
      download_plot_Server("cell_development_curve_P", input_data = plot_cell_development_curve , name ="Cell_development_curve_plot")


    }#2nd Fun
  ) #module server
}#server
