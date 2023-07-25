normalizationUI<-  function(id) {
  ns <- NS(id)
  tagList(
    tabItem(tabName= "normalization",
            fluidPage(



                 fluidRow( style = "height: 78vh; overflow-y: auto;",

                   column(width = 10,
                box(title = "Most variable features:",solidHeader=TRUE,status= 'primary',width = 12,
                    tags$style(HTML(".table.dataTable tbody td.active, .table.dataTable tbody tr.active td {
    background-color: #A020F0;color: white;}")),
                    column(width = 12, plotOutput(ns('QC_MeanExp_vs_SD'), height = "600px", width = "100%") %>% withSpinner(color="#0dc5c1",type = 6,size=0.9)),
                    dropdownMenu =download_plot_UI(ns("normalization_plot"))),


                box(title ="Stats:",solidHeader=TRUE,status= 'primary',width = 12,
                    style='overflow-x: scroll;height:400px;overflow-y: scroll;',

                    column(width = 12, DT::dataTableOutput(ns('QC_MeanExp_stats'),width = "100%", height = "auto") %>% withSpinner(color="#0dc5c1",type = 6,size=0.9 )),
                    dropdownMenu = boxDropdown(

                      downloadBttn(ns("stats_download"),style="jelly", label = "Download stats", size = 'sm', color = 'success', block = TRUE)
                      , icon=icon("download"))  )),

                        column(width = 2,

                       box(title = "Initial Parameters", solidHeader = T, status = "primary", width = 12, collapsible = T,
                           #fileInput(ns("filter_data"),label = "Please Upload filterd .rds file", accept = ".rds",width = '100%'),

                           column(width=12, selectInput( ns("INPUT_Nrmlz_Method"),
                                                         label = "Select Method for Normalization :",
                                                         choices = c(
                                                           "Log normalization" ='LogNormalize',
                                                           "SCTransform"  ='sct',
                                                           "DESeq2::vst"  ='vst',
                                                           "DESeq2::lib_size" = 'lib',
                                                          "SCnorm" = 'scn',
                                                          "Linnorm" = 'lim',
                                                           'ruvseq'), selected='LogNormalize')),
                           column(width=12, numericInput(ns("range"), "Select No. Of Features :",2000, min=0, max=Inf)),
                           column(width=12,
                                  conditionalPanel(
                                    condition ="input.INPUT_Nrmlz_Method=='LogNormalize' || input.INPUT_Nrmlz_Method=='sct'|| input.INPUT_Nrmlz_Method=='vst' || input.INPUT_Nrmlz_Method=='lib' || input.INPUT_Nrmlz_Method=='scn' || input.INPUT_Nrmlz_Method=='lim' || input.INPUT_Nrmlz_Method=='ruvseq'",

                                  selectInput( ns("INPUT_feature_selection_Method"),
                                                         label = "Select variable feature selection Method :",
                                                         choices = c("modelGeneVar"  = 'mvg'), selected='mvg'),  ns = NS(id) )#Conditional panel1
                                  #
                                  # conditionalPanel(
                                  #   condition ="input.INPUT_Nrmlz_Method=='seu'",
                                  #   selectInput( ns("seurat_feature_selection_Method"),
                                  #                label = "Select variable feature selection Method :",
                                  #                choices = c('vst', 'dispersion',"mean.var.plot"), selected='vst'),  ns = NS(id)
                                  # )#conditionalPanel2



                                  ),


                           actionBttn(ns("QC_submit_sidebarPanel"),label="EXEC",style = "jelly",color = "success",icon = icon("sliders")),
                           tags$br(),
                           tags$br(),
                           shinyjs::hidden(wellPanel(id=ns("download_data"),
                           downloadBttn(ns('normalized_data'), 'Normalized data', color= 'royal'),
                           tags$br(),
                           tags$br(),
                           actionBttn(ns("toAnalyze"), "Next Dim reduction",
                                      style = "unite",color = "royal", icon = icon("angles-right",class="fa-duotone fa-angles-right"))))
                          ))#


                       #mainbox
              )#1st fluidrow
            )#fluidpage
    )#tabItem
  )#taglist
}#ui


normalizationServer <- function(id,filtered_data) {
  moduleServer(
    id,
    ## Below is the module function
    function(input, output, session) {
      ns <- session$ns

      vals=reactiveValues()

      # scdata<-reactive({
      #   errors <- c()
      #   req(filtered_data())
      #
      #   # try to read in file
      #   tryCatch(
      #     {
      #       obj <- filtered_data()
      #     },
      #     error = function(e) {
      #       errors <- c(errors, "Invalid rds file.")
      #       return(errors)
      #     }
      #   )
      #
      #   # Validate obj is a SingleCellExperiment object
      #   if (!inherits(obj, "SingleCellExperiment")){
      #     errors <- c(errors, "File is not a SingleCellExperiment object")
      #     return(errors)
      #   }
      #
      #   return(obj)
      # })

      scdata<-reactive({
        req(filtered_data())
        obj <- filtered_data()
        return(obj)
      })
      #Normalization
      normalize <- reactive({
        req(scdata())
        showNotification(paste("Normalization using", input$INPUT_Nrmlz_Method ) , type = "message", duration = 3)
        req(input$INPUT_Nrmlz_Method)

         data= Normalize_Matrix(scdata(),method =input$INPUT_Nrmlz_Method)

        #data@colData = subset(data@colData, select = -c(nCount_RNA,nFeature_RNA, Mito_gene_percent, Hemoglobin_gene_percent, Ribosomal_gene_percent, sizeFactor))
        data
      })#normalization reactive
      #call variable feature plot function from visualizatin.R, return three inputs(1:plot,2:variance and men,3:singlecell object)
      plotlog<-reactive({
        req(normalize())
        showNotification(paste("Variable gene Selection started using", input$INPUT_feature_selection_Method), type= "message", duration = 3 )

        p=VariableGene_hvg_plot(normalize(),used = "NMcounts", ngene = input$range)
        vals$p =p


        # else if(input$INPUT_Nrmlz_Method=='seu')
        # {
        #   object=runSeuratFindHVG(normalize(),'seuratNormData',method = input$seurat_feature_selection_Method)
        #   object=setTopHVG(object, method = input$seurat_feature_selection_Method, hvgNumber = input$range)
        #   hvgs <- getTopHVG(object, hvgNumber = input$range)
        #   p = plotTopHVG(object,method = input$seurat_feature_selection_Method,hvgNumber=input$range,labelsCount = 20)
        #   vals$p
        #   variable_genes_data=as.data.frame(rowData(object))
        #   return(list(p,variable_genes_data,object,hvgs))
        # }

        shinyjs::show('download_data')

        return(p)

      })#plot log


      observeEvent(input$QC_submit_sidebarPanel,{


        output$QC_MeanExp_vs_SD <-renderPlot({


          plotlog()[[1]]

        }, res = 96, bg=useAutoColor())

        output$QC_MeanExp_stats <- DT::renderDataTable({
          req(plotlog())
          datatable(as.data.frame( plotlog()[[2]]))

        })



       # download_plot_Server(ns("normalization_plot"), input=)

      }) #observe Actionbutton

      #Download data and plot
      output$normalized_data <- downloadHandler(
        filename = function() {
          paste0(input$INPUT_Nrmlz_Method, "_filtered_normalized", "-", Sys.Date(), ".rds")
        },
        content = function(file) {
          saveRDS(plotlog()[[3]], file = file)
        }
      )

      output$stats_download<-downloadHandler(
        filename = "Highly_Variable_genes.txt",
        content = function(file) {
          write.table(plotlog()[[2]], file = file, row.names = FALSE, sep = "\t")
        })

      # output$pdf <- downloadHandler(
      #   filename="PreQC.pdf",
      #   content = function(file){
      #     pdf(file,width=input$w,height=input$h)
      #
      #     print(plotlog()[[1]])
      #     dev.off()
      #   }
      # )
      normalizationPlot <- reactive({
        plotlog()[[1]]
      })

      download_plot_Server("normalization_plot", input_data = normalizationPlot , name ="normalization_plot")

      scd<-  reactive({
        data=plotlog()[[3]]
        if(any(c("nCount_RNA","nFeature_RNA", "Mito_gene_percent", "Hemoglobin_gene_percent", "Ribosomal_gene_percent") %in% colnames(colData(data)))){
          data@colData = subset(data@colData, select = -c(nCount_RNA,nFeature_RNA, Mito_gene_percent, Hemoglobin_gene_percent, Ribosomal_gene_percent))
        }

        return(data)
        })

      observeEvent(input$toAnalyze,
                   {
                     shinyjs::runjs("$('a[data-value=\"dim_reduction\"]').tab('show');")
                   })
      return(scd)

    }#function session
  )}#server
