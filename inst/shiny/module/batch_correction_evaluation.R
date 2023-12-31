batchCorrct_evalutn_UI<-  function(id) {
  ns <- NS(id)
  tagList(
    tabItem(tabName= "batch_correction",
            fluidPage(
              fluidRow(style = "height: 78vh; overflow-y: auto;",

                       column(width = 10,
                              wellPanel(
                                plotOutput(ns("batch_evaluation_plot"), height = "800px")%>%
                                  withSpinner(color="#0dc5c1",type = 6,size=0.9),
                                fluidRow(
                                  column(width = 12, align = "right",
                                         download_plot_UI(ns("batch_crrcn_evln_plot"))
                                  )),
                                style = "height: auto; border: 3px solid #CEECF5;",


                              )
                              ),

                       column(width = 2,

                              wellPanel(
                                selectInput(ns("batch_evaluation"), label = "Choose method to evaluate batch correction",
                                            choices = c("BatchEva_cmix" = "cmx",
                                                        "BatchEva_kbet" = "kbet",
                                                        "BatchEva_slt" = 'slt',
                                                        "BatchEva_lisi" = "lisi"), selected = "slt"),
                                style = " background-color: #CEECF5; border: 3px solid #CEECF5;"
                              )


                       )




              )#fluidRow
            )
    )
  )
}


batchCorrct_evalutn_Server <- function(id, batch_corrct) {
  moduleServer(
    id,
    function(input, output, session) {

      scdata<-reactive({
        req(batch_corrct())
        validate(
          need(class(batch_corrct()) == "SingleCellExperiment" && any(reducedDimNames(batch_corrct())=='BEPCA'),
               "Please use another approch for batch correction")

        )
        batch_corrct()
      })



      batch_evaluation<-reactive({
        req(scdata())
        req(input$batch_evaluation)

        withProgress(message = 'Batch evaluation in progress', value = 0, {
          # for(i in 1:10) {
          #   incProgress(1/10)
          #   Sys.sleep(0.5)
          # }
        d= BatchEva(sce = scdata(), method = input$batch_evaluation,used= "BEPCA")

        return(d)
        })
      })

      output$batch_evaluation_plot <-renderPlot({
        batch_evaluation()
      })

      download_plot_Server("batch_crrcn_evln_plot", input_data = batch_evaluation ,
                           name =paste0("batch_correction_evaluation_by_", input$batch_evaluation))

      }
  )
}
