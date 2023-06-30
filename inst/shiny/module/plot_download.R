download_plot_UI<-function(id) {
  ns <- NS(id)
 
    boxDropdown(
      numericInput(ns("w"), label = "Figure Width", value = 14),
      numericInput(ns("h"), label = "Figure Height", value = 12),
      numericInput(ns("ppi"), label = "Figure Resolution", value = 350),
      
      downloadBttn(
        outputId = ns("pdf"),
        label="PDF figure",
        style = "jelly",
        color = "success",
        size='sm',
        block=TRUE
      ),
      downloadBttn(
        outputId = ns("png"),
        label="PNG figure",
        style = "jelly",
        color = "success",
        size='sm',
        block=TRUE
      ),
      downloadBttn(
        outputId = ns("jpeg"),
        label="JPEG figure",
        style = "jelly",
        color = "success",
        size='sm',
        block=TRUE
      ),
      downloadBttn(
        outputId = ns("tiff"),
        label="TIFF figure",
        style = "jelly",
        color = "success",
        size='sm',
        block=TRUE
      )
      , icon=icon("download"))
    

}

download_plot_Server <- function(id, input_data, name) {
  moduleServer(
    id,
    ## Below is the module function
    function(input, output, session) {
      ns <- session$ns
      
      output$pdf <- downloadHandler(
        filename=paste0(name, ".pdf"),
        content = function(file){
          pdf(file,width=input$w,height=input$h)
          
          print(input_data())
          dev.off()
        }
      )
      
      output$tiff <- downloadHandler(
        filename=paste0(name, ".tiff"),
        content = function(file){
          tiff(file,width=input$w,height=input$h,units="in",res=input$ppi)
          
          print(input_data())
          dev.off()
        }
      )
      output$png <- downloadHandler(
        filename=paste0(name,".png"),
        content = function(file){
          png(file,width=input$w,height=input$h,units="in",res=input$ppi)
          print(input_data())
          dev.off()
        }
      )
      output$jpeg <- downloadHandler(
        filename=paste0(name,".jpeg"),
        content = function(file){
          jpeg(file,width=input$w,height=input$h,units="in",res=input$ppi)
          print(input_data())
          dev.off()
        })

      
      
    }
  )
}