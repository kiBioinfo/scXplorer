# variable gene selection:
batch_variable_gene_UI<-  function(id) {
  ns <- NS(id)
  tagList(
    tabItem(tabName= "variable_gene_batch",
            fluidPage(
              fluidRow(style = "height: 78vh; overflow-y: auto;",
                       column(width = 10,
                              box(title = "Most variable features:",solidHeader=TRUE,status= 'primary',width = 12,
                                  column(width = 12, plotOutput(ns('QC_MeanExp_vs_SD'), height = "600px", width = "100%") %>% withSpinner(color="#0dc5c1",type = 6,size=0.9)),
                                  dropdownMenu = boxDropdown(
                                    boxDropdownItem(
                                      downloadButton(ns("plot_download"),style="color: #fff; background-color: purple; border-color: purple;", class = "primary", label = "Download Plot")
                                    ), icon=icon("download"))),
                              
                              
                              box(title ="Stats:",solidHeader=TRUE,status= 'primary',width = 12,
                                  style='overflow-x: scroll;height:400px;overflow-y: scroll;',
                                  
                                  column(width = 12, DT::dataTableOutput(ns('QC_MeanExp_stats'),width = "100%", height = "auto") %>% withSpinner(color="#0dc5c1",type = 6,size=0.9 )),
                                  dropdownMenu = boxDropdown(
                                    boxDropdownItem(
                                      downloadButton(ns("stats_download"),style="color: #fff; background-color: purple; border-color: purple;", class = "primary", label = "Download stats")
                                    ), icon=icon("download"))  )),#Output column
                       column(width = 2,
                              
                              box(title = "Initial Parameters", solidHeader = T, status = "primary", width = 12, collapsible = T,
                                   numericInput("range", "Select No. Of Features :",2000, min=0, max=Inf),
                                  selectInput( ns("INPUT_feature_selection_Method"),
                                               label = "Select variable feature selection Method :",
                                               choices = c("ModelGeneVar"  = 'mvg'), selected='mvg')
                              ),
                              actionBttn(ns("submit"),label="Apply",style = "jelly",color = "success",icon = icon("sliders")),
                              tags$br(),
                              tags$br(),
                              
                              downloadBttn(ns('normalized_data'), '.rds normalized batch corrected data')
                       )#input column
                       



)#1st fluidrow
)#fluidpage
)#tabItem
)#taglist
}#ui


batch_variable_gene_Server <- function(id,batch_corrct) {
  moduleServer(
    id,
    ## Below is the module function
    function(input, output, session) {
      
      
      vals=reactiveValues()
      scdata<-reactive({
        batch_corrct()
      })
      
      batch_VG<- reactive({
        data= VariableGene_hvg_plot(scdata(), used = 'counts',ngene = input$range , batch = TRUE)

        vals$data= data
        data
      })
      
      observeEvent(input$submit,{
        
        output$QC_MeanExp_vs_SD <-renderPlot({
          
          batch_VG()[[1]]
          
        }, res = 96)
        
        output$QC_MeanExp_stats <- DT::renderDataTable({
          
          
          datatable(as.data.frame( batch_VG()[[2]]))
          
        })
        
        
        output$normalized_data <- downloadHandler(
          filename = function() {
            paste0(input$INPUT_Nrmlz_Method, "_filtered_normalized_batch_crrected_VG", "-", Sys.Date(), ".rds")
          },
          content = function(file) {
            saveRDS(batch_VG()[[3]], file = file)
          }
        )
      }) #observe Actionbutton
      
      return(list(scdata=reactive({batch_VG()[[3]]})))  
    }#function session
  )}#server