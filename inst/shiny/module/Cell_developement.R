cell_developement_analysis_UI<-function(id) {
  ns <- NS(id)
  tagList(
    fluidPage(
      useShinyjs(),
      fluidRow(style = "height: 78vh; overflow-y: auto;",
               column(width = 10,
                      wellPanel(
                        plotOutput(ns('cell_development'), height = "800px")%>%
                          withSpinner(color="#0dc5c1",type = 6,size=0.9),
                        fluidRow(
                          column(width = 12, align = "right",
                        download_plot_UI(ns("cell_development_P"))
                          )),
                        style = "height: auto; border: 3px solid #CEECF5;"
                        ),

                      tags$br(),
                      shinyjs::hidden(wellPanel(id=ns("curve"),
                        plotOutput(ns('cell_development_curve'), height = "800px")%>%
                          withSpinner(color="#0dc5c1",type = 6,size=0.9),
                        fluidRow(
                          column(width = 12, align = "right",
                        download_plot_UI(ns("cell_development_curve_P"))
                          )),
                        style = "height: auto; border: 3px solid #CEECF5;",

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
                        downloadBttn(ns('final_result'), 'Download Result', color= 'royal'),

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
        req(cell_type())
        obj <- cell_type()
        return(obj)
      })

      output$colour <- renderUI({
        req(scdata())
        selectInput(ns("colour_by"),
                    label = "Colour by:",
                    choices = colnames(colData(scdata())), selected = rev(colnames(colData(scdata())))[1])
      })

      output$reduce_dim<-renderUI({
        req(scdata())
        ns <- session$ns
        selectInput(ns("reduction"),
                    label = "select the dim reduction to plot:",
                    choices = reducedDimNames(scdata()), selected = rev(reducedDimNames(scdata()))[1])
      })
      output$seedby <- renderUI({
        ns <- session$ns
        req(scdata())
        # if(!("cellType" %in% colnames(colData(scdata())))){
        #   selectizeInput(ns("seedBy"), label= "Choose the cell types:", choices=unique((scdata()@colData$label)), selected=unique((scdata()@colData$label))[1], multiple=T)
        # }
        # else{
        #   selectizeInput(ns("seedBy"), label= "Choose the cell types:", choices=unique((scdata()@colData$cellType)), selected=unique((scdata()@colData$cellType))[1], multiple=T)
        # }
        selectizeInput(ns("seedBy"), label= "Choose the cell types:",
                       choices=unique((scdata()@colData$cellType)), selected=unique((scdata()@colData$cellType))[1], multiple=T)
      })

      Cell_develpment<-reactive({
        req(scdata())
        # req(input$method)
        # req(input$reduction)
        # req(input$cell_type)
        if(input$method=="Princurve_fit"){
          shinyjs::showElement("curve")
        }
        else{
          shinyjs::hideElement("curve")
        }
        if(input$method=="Monocle"){
          req(input$seedBy)
          p=CellDevelopment(scdata(),method=input$method, seedBy = input$seedBy)
        }
        else{
          p=CellDevelopment(scdata(),method=input$method, runWith = input$reduction)
        # browser()
        }

       vals$p <- p
       p

      })

  #plot cell development

      plot_cell_development <- reactive({
        withProgress(message = 'Building single-cell trajectories...', value = 0, {
        req(Cell_develpment())
        Plot_Development(sce = Cell_develpment() ,runWith = input$reduction, colorBy = input$colour_by )
        })
        })

      plot_cell_development_curve <-reactive({
        req(Cell_develpment())

        Plot_Development_Princurve(sce =Cell_develpment(), runWith = input$reduction, colorBy = input$colour_by )
      })



      #shinyjs::enable("cell_development_curve")
      observeEvent(input$method,{
        if(input$method=='Princurve_fit'){
          output$cell_development <-renderPlot({
            req(plot_cell_development())

            plot_cell_development()

          }, res = 96)

          output$cell_development_curve <-renderPlot({

            plot_cell_development_curve()

          }, res = 96)
        }
        else{
          output$cell_development <-renderPlot({
            req(plot_cell_development())

            plot_cell_development()

          }, res = 96)
        }

        })


      #Download data
      output$final_result <- downloadHandler(
        filename = function() {
          paste0( "scRNA-Seq_analysis", "-", Sys.Date(), ".rds")
        },
        content = function(file) {
          saveRDS(Cell_develpment(), file = file)
        }
      )


      #Download plot
      download_plot_Server("cell_development_P", input_data = plot_cell_development , name ="Cell_development_plot")
      download_plot_Server("cell_development_curve_P", input_data = plot_cell_development_curve , name ="Cell_development_curve_plot")


    }#2nd Fun
  ) #module server
}#server
