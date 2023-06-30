library(shinythemes)
source("helper.R")
datainputUI<-  function(id) {
  ns <- NS(id)
  tagList(
    tabItem(tabName= "dataInput",
      fluidPage(#theme = shinytheme("darkly"),
      fluidRow(style = "height: 78vh; overflow-y: auto;", 
        column(width = 4,
    box(title="Upload Raw Data:",solidHeader=TRUE, status='primary',
        column(width = 12,
        radioButtons(ns("data"), NULL,choices = c(
          
                                              "10X Data (1 .mtx and 2 .tsv)"="10X",
                                              "HDF5 Data (*.h5)"="10X_f5",
                                              "Raw Count Matrix"="count_Mtrx",
                                              "R Object"="seurat_obj"), selected = "10X"),
        conditionalPanel(condition  = "input.data == '10X'",
                         tags$b("Select 10x files folder:"),
                         tags$p(),
                         verbatimTextOutput(ns("input_10x_normal"), placeholder = TRUE),
                         shinyDirButton(ns("normal_10x"), "Select 10x file folder", "Please select a folder",buttonType = "success",icon=icon("folder",lib = "font-awesome",class="fa-solid fa-folder"),style = "bold" ),
                         ns = NS(id) 
                         
        ),
        
        
        conditionalPanel(
          condition = "input.data=='10X_f5'",
          tags$b("Please select a .f5 data file:"),
          tags$p(),
          verbatimTextOutput(ns("input_f5"), placeholder = TRUE),
          shinyFilesButton(ns("f5"), "Select file",
                           "Please select file",buttonType = "success",icon=my_icon("file-upload",lib = "font-awesome", class="fa-solid fa-file-upload"),style = "bold", multiple = FALSE),
          ns = NS(id)
        ),
        
        conditionalPanel(
          condition = "input.data=='count_Mtrx'",
          tags$b("Please select a file containg raw count data:"),
          tags$p(),
          verbatimTextOutput(ns("input_Cmtrx"), placeholder = TRUE),
          shinyFilesButton(ns("Cmtrx"), "Select file",
                           "Please select file",buttonType = "success",icon=my_icon("file-upload",lib = "font-awesome", class="fa-solid fa-file-upload"),style = "bold" , multiple = FALSE),
          ns = NS(id)
        ),
    
    conditionalPanel(
      condition = "input.data=='seurat_obj'",
      tags$b("Please select a .rds data file:"),
      tags$p(),
      verbatimTextOutput(ns("input_seuObj"), placeholder = TRUE),
      shinyFilesButton(ns("seuObj"), "Select file",
                       "Please select file", buttonType = "success",icon=my_icon("file-upload",lib = "font-awesome", class="fa-solid fa-file-upload"),style = "bold", multiple = FALSE),
      ns = NS(id)
    
        ))),
    
    box(title = "2. Initialize SingleCellExperiment object", height = 350, solidHeader = TRUE, status='primary',
        collapsible = T, id = "object_init",
        
        textInput(ns("project_name"), value = "Project1", label = "Project Name"),
        
        fileInput(ns("metad"), 'Choose Meta data file', multiple = F),
        withBusyIndicatorUI(icon_name = "sinit",
                            actionBttn(ns("seurat_init"),label="Process raw data", style = "jelly",color = "success",icon = icon("arrows-spin", class="fa-solid fa-arrows-spin"),class = "btn-primary")
        )
        
    )),
    column(width = 8,
    box(title = "3. Basic Stats",width=8, solidHeader = TRUE, status='primary',collapsible = T,
      column(width=12,
             
             withSpinner(tableOutput(ns("basic_stats")), color="#0dc5c1",type = 6,size=0.9),
             
       
             tags$br(),
             
             
               
             
      )
    )),
    column(width = 2,
           actionBttn(ns("toAnalyze"), "Continue to QC Tab",class = "button button-3d button-block button-pill button-caution",style = "unite",color = "royal", icon = icon("angles-right",class="fa-duotone fa-angles-right")))
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
      shinyFileChoose(input, "f5", roots = volumes, session = session)
      shinyFileChoose(input, "Cmtrx", roots = volumes, session = session)
      shinyFileChoose(input, "seuObj", roots = volumes, session = session)
      
      output$input_10x_normal <- renderPrint({
        parseDirPath(volumes, input$normal_10x)
      })
      output$input_Cmtrx <- renderPrint({
        as.character(parseFilePaths(volumes, input$Cmtrx)[4])
      })
      output$input_f5 <- renderPrint({
        as.character(parseFilePaths(volumes, input$f5)[4])
      })
      output$input_seuObj <- renderPrint({
        as.character(parseFilePaths(volumes, input$seuObj)[4])
      })
      
      
      observeEvent(input$seurat_init, {
        
        withProgress(message = "Initializing Seurat Object ...",{
          shinyjs::show("load_sinit")
          
          if(input$data == "seurat_obj"){
            File_R_Obj <- isolate(as.character(parseFilePaths(volumes, input$seuObj)[4]))
            rawdata <- readRDS(File_R_Obj)
            if(!class(rawdata)=="Seurat"){
              
              rawdata=as.Seurat(rawdata)
            }
            
            #rownames(rawdata) <- unlist(map(gsub(".*_","",rownames(rawdata)),1))
          }
          else if(input$data == "10X"){
            Path_10x_normal <- isolate(as.character(parseDirPath(volumes, input$normal_10x)))
            rawdata <- Read10X(data.dir = Path_10x_normal,gene.column = 1)
            rownames(rawdata) <- unlist(map(gsub(".*_","",rownames(rawdata)),1))
            if(is.null(input$metad)){
            rawdata <- CreateSeuratObject(counts = rawdata, project = input$project_name, min.cells = 3, min.features = 200)
            }
            else
            {
              rawdata <- CreateSeuratObject(counts = rawdata, meta.data=input$metad,project = input$project_name, min.cells = 3, min.features = 200)
            }
          }
          else if(input$data =="count_Mtrx"){
            File_c_mtrx <- isolate(as.character(parseFilePaths(volumes, input$Cmtrx)[4]))
            rawdata <- as.data.frame(vroom::vroom(File_c_mtrx))
            rownames(rawdata) <- rawdata[,1]
            rownames(rawdata) <- unlist(map(gsub(".*_","",rownames(rawdata)),1))
            if(is.null(input$metad)){
              rawdata <- CreateSeuratObject(counts = rawdata, project = input$project_name, min.cells = 3, min.features = 200)
            }
            else
            {
              rawdata <- CreateSeuratObject(counts = rawdata, meta.data=input$metad,project = input$project_name, min.cells = 3, min.features = 200)
            }
          }
          else if(input$data == "10X_f5"){
            File_f5_Obj <- isolate(as.character(parseFilePaths(volumes, input$f5)[4]))
            rawdata <- Read10X_h5(File_f5_Obj)
            rownames(rawdata) <- unlist(map(gsub(".*_","",rownames(rawdata)),1))
            if(is.null(input$metad)){
              rawdata <- CreateSeuratObject(counts = rawdata, project = input$project_name, min.cells = 3, min.features = 200)
            }
            else
            {
              rawdata <- CreateSeuratObject(counts = rawdata, meta.data=input$metad,project = input$project_name, min.cells = 3, min.features = 200)
            }
          }
          else{
            rawdata <- Read10X(data.dir = "useful/normal_analysis/",gene.column = 1)
            rownames(rawdata) <- unlist(map(gsub(".*_","",rownames(rawdata)),1))
          }
          
          setProgress(value = 0.5)
          
         
          
          
          rawdata[["percent.mt"]] <- PercentageFeatureSet(object = rawdata, pattern = "(?i)^MT-") #(?i) for case insensitive
          CS.data <<- as.SingleCellExperiment(rawdata)
          names(colData(CS.data))[c(1:3)]<-c("orig.ident"  , "nCount_RNA" ,  "nFeature_RNA")
          #data to use in other modules
          CS.data <<- CS.data
          rm(rawdata)
          setProgress(value = 0.7)
          
          median_gene <- median(CS.data@colData$nFeature_RNA)
          min_gene <- min(CS.data@colData$nFeature_RNA)
          max_gene <- max(CS.data@colData$nFeature_RNA)
          
          min_mito <- round(min(CS.data@colData$percent.mt), 2)
          max_mito <- round(max(CS.data@colData$percent.mt), 2)
          median_mito <- round(median(CS.data@colData$percent.mt), 2)
          
          min_UMI <- min(CS.data@colData$nCount_RNA)
          max_UMI <- max(CS.data@colData$nCount_RNA)
          median_UMI <- median(CS.data@colData$nCount_RNA)
          
          nGene_summary <- data.frame(c(min_gene, median_gene, max_gene),
                                      c(min_UMI, median_UMI, max_UMI),
                                      c(min_mito,median_mito,max_mito),
                                      row.names = c("Min", "Median", "Max"))
          nGene_summary <- t(nGene_summary)
          row.names(nGene_summary) <- c("nGene", "nUMI","nMito")
          output$basic_stats <-renderTable(spacing = "xs", align = "c", rownames = TRUE,
                                           striped = TRUE, hover = TRUE, bordered = TRUE, digits = 0,
                                           { nGene_summary }
          )
          
          setProgress(value = 1)
          
          shinyjs::hide("load_sinit")
          shinyjs::show("check_sinit")
        })#Wsprogress
      
      })#observeE_seurat
     
      #Disable batch correction tab
      #shinyjs::hide(selector = "a[data-value=\"batch_correction\"]")
      shinyjs::hide(selector = "a[data-value=\"variable_gene_batch\"]")
      shinyjs::hide(selector = "a[data-value=\"linear_batch\"]")
      shinyjs::hide(selector = "a[data-value=\"non_linear_batch\"]")
      shinyjs::hide(selector = "a[data-value=\"batch_evaluation\"]")
     
   
      observeEvent(input$toAnalyze,
                   {
                     shinyjs::runjs("$('a[data-value=\"PreQC\"]').tab('show');")
                   })
       
    })#moduleServer
}#dataInputServer
