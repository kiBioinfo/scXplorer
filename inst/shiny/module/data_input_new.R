library(shinythemes)
source("helper.R")
datainputUI<-  function(id) {
  ns <- NS(id)
  tagList(
    tabItem(tabName= "dataInput",
            fluidPage(#theme = shinytheme("darkly"),
              
              fluidRow(
                column(width = 2,
                       box(title="Upload Raw Data:",solidHeader=TRUE, status='primary',width = 12,
                           
                           column(width = 12,
                                  radioButtons(ns("file_type"), NULL,choices = c(
                                    
                                    "10X Data (1 .mtx and 2 .tsv)"="10X",
                                    "HDF5 Data (*.h5)"="10X_f5",
                                    "Raw Count Matrix"="count_Mtrx",
                                    "R Object"="rds_obj"), selected = "rds_obj"),
                                  
                                  conditionalPanel(condition  = "input.file_type == '10X'",
                                                   tags$b("Select 10x files folder (2 .tsv and 1 .mtx files):"),
                                                   tags$p(),
                                                   
                                                   shinyDirButton(ns("normal_10x"), "Select 10x file folder", "Please select a folder",buttonType = "success",icon=icon("folder",lib = "font-awesome",class="fa-solid fa-folder"),style = "bold" ),
                                                   ns = NS(id)
                                                   
                                  ),
                                  
                                  
                                  conditionalPanel(
                                    condition = "input.file_type=='10X_f5'",
                                    tags$b("Please select a .f5 data file:"),
                                    tags$p(),
                                    
                                    fileInput(ns("File_f5_Obj"), "Select file", accept=".h5", multiple = FALSE ),
                                    
                                    ns = NS(id)
                                  ),
                                  
                                  conditionalPanel(
                                    condition = "input.file_type=='count_Mtrx'",
                                    tags$b("Please select a file containg raw count data (.txt/.csv/.gz):"),
                                    tags$p(),
                                    
                                    fileInput(ns("Cmtrx"),"Please select file",accept=c(".gz", ".txt", ".csv"), multiple = FALSE),
                                    
                                    ns = NS(id)
                                  ),
                                  
                                  conditionalPanel(
                                    condition = "input.file_type=='rds_obj'",
                                    tags$b("Please select a .rds data file:"),
                                    tags$p(),
                                    
                                    fileInput(ns("rdsObj"), "Please Select file", accept=".rds", multiple = FALSE ),
                                    ns = NS(id)
                                    
                                  ),
                                  tags$br(),
                                  #uiOutput(ns("process_raw_data")),
                                  shinyjs::hidden(actionBttn(ns("sce_object_init"),label="Process raw data", style = "jelly",color = "success",icon = icon("arrows-spin", class="fa-solid fa-arrows-spin"),class = "btn-primary")),
                                  tags$br(),
                                  tags$br(),
                                  shinyjs::hidden(materialSwitch(inputId = ns("plot_data_switch"), label = strong("Plot raw data"), status = "danger",value = FALSE )),
                                  shinyjs::hidden(uiOutput(ns("colour_plot"))),
                                  
                                  
                                  tags$br(),
                                  shinyjs::hidden(actionBttn(ns("toAnalyze"), "Continue to QC Tab",class = "button button-3d button-block button-pill button-caution",
                                                             style = "unite",color = "royal", icon = icon("angles-right",class="fa-duotone fa-angles-right")))
                                  
                           )#Column
                           
                       )#box
                ), #column
                
                column(width = 10,
                       box(title = " Basic Stats",width=12, solidHeader = TRUE, status='primary',collapsible = T,
                           style = "height: 78vh; overflow-y: auto;",
                           tags$style(HTML(".shiny-output-error-validation {color: #fa0052;}")),
                           column(width=12,
                                  
                                  withSpinner(DT::dataTableOutput(ns("basic_stats")), color="#0dc5c1",type = 6,size=0.9),
                                  #verbatimTextOutput(ns("object")),
                                  shinyjs::hidden(wellPanel(id=ns("plot"),
                                                            withSpinner(plotOutput(ns("raw_data_plot"), height = "900px", width = "100%"), color="#0dc5c1",type = 6,size=0.9)
                                                            ,style = "height= auto;width: 100%; border: 3px solid #CEECF5; text-align: end;",
                                                            download_plot_UI(ns("plt_raw_data"))
                                  ))
                                  
                                  
                                  
                                  
                                  
                           )
                       ))
                
                
              )
              
            )
    )
  )
}


datainputServer <- function(id) {
  moduleServer(
    id,
    ## Below is the module function
    function(input, output, session, fileRoot = NULL) {
      ns <- session$ns
      
      volumes <- (c(Home = fs::path_home(), getVolumes()()))
      
      shinyDirChoose(input, "normal_10x", roots = volumes, session = session, restrictions = system.file(package = "base"))
      
      
      #Load Raw Data
      # Create a reactive object for each uploaded file
      rv <- reactiveValues(data = NULL)
      observe({
        req(input$normal_10x)
        
        Path_10x_normal <- as.character(parseDirPath(volumes, input$normal_10x))
        validate(
          need(Path_10x_normal != "", "Please select a data set")
          
        )
        data <- Seurat::Read10X(Path_10x_normal)
        rv$data<- SingleCellExperiment::SingleCellExperiment(assay = list(counts = data))
      })
      observe({
        req(input$rdsObj)
        validate(
          need(input$rdsObj != "", "Please select a data set")
          
        )
        rv$data<-load_sce_obj(input$rdsObj$datapath)
      })
      observe({
        req(input$Cmtrx)
        validate(
          need(input$Cmtrx != "", "Please select a data set")
          
        )
        
        rv$data<- load_sce_obj(input$Cmtrx$datapath)
      })
      observe({
        req(input$File_f5_Obj)
        validate(
          need(input$File_f5_Obj != "", "Please select a data set")
          
        )
        data<-Seurat::Read10X_h5(input$File_f5_Obj$datapath)
        rv$data<- SingleCellExperiment::SingleCellExperiment(assay = list(counts = data))
        
      })
      

      # Clear the uploaded file(s) when the user selects a different radio button option
      observeEvent(input$file_type, {
        shinyjs::hide("sce_object_init")
        shinyjs::hide("plot_data_switch")
        shinyjs::hide("toAnalyze")
        shinyjs::hide("colour_plot")
        shinyjs::hide("plot")
        if (input$file_type == "10X") {
          rv$data=NULL
          if (!is.null(input$File_f5_Obj)) {
            shinyjs::reset("File_f5_Obj")
          }
          if (!is.null(input$Cmtrx)) {
            shinyjs::reset("Cmtrx")
          }
          if(!is.null(input$rdsObj)){
            shinyjs::reset("rdsObj")
          }
        } else if (input$file_type == "10X_f5") {
          rv$data=NULL
          if (!is.null(input$normal_10x)) {
            shinyjs::reset("normal_10x")
          }
          if (!is.null(input$Cmtrx)) {
            shinyjs::reset("Cmtrx")
          }
          if(!is.null(input$rdsObj)){
            shinyjs::reset("rdsObj")
          }
        } else if (input$file_type == "count_Mtrx") {
          rv$data=NULL
          if (!is.null(input$normal_10x)) {
            shinyjs::reset("normal_10x")
          }
          if (!is.null(input$File_f5_Obj)) {
            shinyjs::reset("File_f5_Obj")
          }
          if(!is.null(input$rdsObj)){
            shinyjs::reset("rdsObj")
          }
        }else if (input$file_type == "rds_obj") {
          rv$data=NULL
          if (!is.null(input$normal_10x)) {
            shinyjs::reset("normal_10x")
          }
          if (!is.null(input$File_f5_Obj)) {
            shinyjs::reset("File_f5_Obj")
          }
          if (!is.null(input$Cmtrx)) {
            shinyjs::reset("Cmtrx")
          }
          
        }
      })
      observeEvent(rv$data,{
        if(class(rv$data)=='SingleCellExperiment'){
          shinyjs::show("sce_object_init")
          
        }
      })
      
      
      #instruction to load metadata
      dataModal <- function(failed = FALSE) {
        #Example dataset
        
        modalDialog(
          span(h5("Please upload the meta data file"),
               "The first column should be the CELLs included in the experiment"),
          tags$hr(),
          materialSwitch(inputId = ns("modal_metadata_switch"), label = "Show example annotation file", status = "danger",value = FALSE ),
          DT::dataTableOutput(ns("modal_metadata")),
          downloadButton(ns('download_anno_data'), 'Download'),
          #renderTable(tail(d,20),rownames=FALSE),
          easyClose=TRUE,
          footer = tagList(
            fileInput(ns("metad"), 'Choose a ".txt" Meta data file',accept = ".txt", multiple = F),
            tags$br(),
            Switch.shinyInput(
              inputId = ns("skip_metadata"),
              value = FALSE,
              label = "Skip adding  Metadata"
            ),
            modalButton("Close")
          )
        )
      }
      
      #Stats of sce object
      scdata<-reactive({
        
        req(input$sce_object_init)
        
        validate(
          need(inherits(rv$data,"SingleCellExperiment"), "Please select a data set")
        )
        
        #Loading metadata file
        #manipulate metadata and count matrix to keep matched CELLs in sce object
        # req(input$skip_metadata)
        if(input$skip_metadata==TRUE){
          metadata=NULL
        }
        else{
          file<- input$metad
          ext <- tools::file_ext(file$datapath)
          
          req(file)
          validate(need(ext == "txt", "Please upload a '.txt' metadata file"))
          metadata=file$datapath
        }
        
        
        removeModal()
        
        scdata=add_QC_matrix(rv$data, metadata= metadata )
        
        
        shinyjs::show("plot_data_switch")
        shinyjs::show("toAnalyze")
        scdata
      })
      
      
      
      #Process Raw data and add annotation
      observeEvent(input$sce_object_init,{
        req(input$sce_object_init)
        
        updateMaterialSwitch(session, 'plot_data_switch', value = FALSE)
        
        
        if(inherits(rv$data,"SingleCellExperiment")){
          
          showModal(dataModal())
        }
        
      })
      
      
      output$basic_stats<-DT::renderDataTable({
        validate(
          need(inherits(scdata(),"SingleCellExperiment"), "Please add Metadata")
        )
        QC_Stats(scdata())
        
      })
      
      
      #Modal annotation data
      
      modal_data <- reactive({
        if(input$modal_metadata_switch==TRUE){
          d=metadata_modal(current_data=NULL)
        }
        else{
          d=metadata_modal(current_data=rv$data)
        }
        return(d)
      })
      
      
      output$modal_metadata<- DT::renderDataTable(tail(modal_data(),20),rownames=FALSE)
      
      output$download_anno_data<-downloadHandler(
        filename = "annotation.txt",
        content = function(file) {
          write.table(modal_data(), file = file, row.names = FALSE, sep = "\t")
        })
      
      
      # Choose Color for plot
      output$colour_plot<- renderUI({
        req(scdata())
        tagList(
          h6("Colour By:"),
          Select.shinyInput(inputId=ns("colour_by"),
                            
                            items = colnames(scdata()@colData) , selected = rev(colnames(scdata()@colData))[1],
                            noResults= "No Result:")
        )
      })
      
      observe({
        if(input$plot_data_switch==TRUE & class(rv$data)=='SingleCellExperiment' )
        {
          shinyjs::show("plot")
          shinyjs::show("raw_data_plot")
          shinyjs::show("colour_plot")
          
          shinyjs::hide("basic_stats")
          
        }
        else{
          shinyjs::hide("raw_data_plot")
          shinyjs::hide("colour_plot")
          shinyjs::hide("plot")
          shinyjs::show("basic_stats")
        }
      })
      
      
      
      #plot Raw data
      raw_data_plt<- reactive({
        
        p=Plot_raw_data(sce= scdata(), colour_by= input$colour_by$text)
        p
      })
      
      output$raw_data_plot <-renderPlot({
        raw_data_plt()
        
      }, res = 90)
      
      download_plot_Server("plt_raw_data", input_data = raw_data_plt , name ="Raw_data_plot")
      
      
      #Disable batch correction tab
      shinyjs::hide(selector = "a[data-value=\"batch_correction\"]")
      shinyjs::hide(selector = "a[data-value=\"variable_gene_batch\"]")
      shinyjs::hide(selector = "a[data-value=\"linear_batch\"]")
      shinyjs::hide(selector = "a[data-value=\"non_linear_batch\"]")
      shinyjs::hide(selector = "a[data-value=\"batch_evaluation\"]")
      
      observeEvent(input$toAnalyze,
                   {
                     shinyjs::runjs("$('a[data-value=\"PreQC\"]').tab('show');")
                   })
      return(scdata)
      
    })#moduleServer
}#dataInputServer