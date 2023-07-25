preQCUI<-  function(id) {
  ns <- NS(id)
  tagList(

            fluidPage(

              fluidRow(
                style = "height: 78vh; overflow-y: auto;",
                box(width = 10,status = 'primary',
                    fluidRow(
                      box(title="Total Gene Plot",solidHeader=TRUE,status='primary',closable=T,
                          column(width = 12,


                                 plotOutput(ns("beforeQC_totalGene"),width="100%") %>% withSpinner(color="#0dc5c1",type = 5,size=0.5),

                                 maximizable = T)),

                      box(title = "Total Count",solidHeader=TRUE,status='primary', closable=T,
                          column(width = 12,
                                 plotOutput(ns('beforeQC_totalCount'), width = '100%') %>% withSpinner(color="#0dc5c1",type = 5,size=0.5),

                                 maximizable = T

                          )),

                      box(title="%MT",solidHeader=TRUE,status='primary',closable=T,
                          column(width = 12,
                                 plotOutput(ns('beforeQC_totalMT'),   width = '100%')%>% withSpinner(color="#0dc5c1",type = 5,size=0.5),
                                 maximizable = T)),

                      box(title="%Ribo",solidHeader=TRUE,status='primary',closable=T,
                          column(width = 12,
                                 plotOutput(ns('QC_RiboGene'),   width = '100%')%>% withSpinner(color="#0dc5c1",type = 5,size=0.5), maximizable = T)),



                      box(title="TotalCount Vs TotalGene",solidHeader=TRUE,status='primary',closable=T,
                          column(width = 12,
                                 plotOutput(ns('beforeQC_totalCount_vs_totalGene'),   width = '100%')%>% withSpinner(color="#0dc5c1",type = 5,size=0.5),maximizable = T)),

                      box(title="Summary",solidHeader=TRUE,status='primary',closable=T,
                          column(width = 12,
                                 tableOutput(ns('sce_obj')),   maximizable = T))

                      )),


                column(width=2,
                           box(title = "QC Filter Parameters", solidHeader = T, status = "primary", width = 12, collapsible = T,id = "uploadbox",
                               tags$style(HTML(".irs--shiny .irs-bar  {background:  #A020F0; border:  #A020F0}")),

                               tags$style(HTML(".irs--shiny .irs-single {background:  #A020F0}")),
                           #sliderInput(ns('QC_TotalGene_Down'), label = 'Min-genes', min=0, max=5000,value = 100),
                           uiOutput(ns("min_gene_per_cell")),

                           #sliderInput(ns('QC_TotalGene_Up'),   label = 'Max-gene', min=500,max=15000, value = 2500),
                           uiOutput(ns("max_gene_per_cell")),

                           #sliderInput(ns('QC_TotalCount_Down'), label = 'Min-count',min=5, max=10000, value = 100),
                           uiOutput(ns("min_count")),

                           #sliderInput(ns('QC_TotalCount_Up'),   label = 'Max-count', min=500,max=20000000,value = 10000000),
                           uiOutput(ns("max_count")),

                           #sliderInput(ns('QC_MT'), label = '%MT cut-off', min=0, max=100,value = 5)  ,
                           uiOutput(ns("mt_perc")),

                           #sliderInput(ns('QC_ribo'), label = '%Ribo cut-off', min=0, max=100,value = 5),
                           uiOutput(ns("ribo_perc")),
                           uiOutput(ns("beforeQC_selectFeatureVector")),

                            actionBttn(ns("beforeQC_submit_sidebarPanel"),label="EXEC",style = "jelly",color = "success",icon = icon("sliders"))),

                           box(width = 12, status = "primary", collapsible = TRUE, solidHeader = TRUE, id=ns("download_box"),title = "Download",
                           numericInput(ns("w"), label = "Figure Width", value = 14),
                           numericInput(ns("h"), label = "Figure Height", value = 12),
                           numericInput(ns("ppi"), label = "Figure Resolution", value = 350),
                           dropdownButton(
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
                             ),
                             circle=FALSE,
                             label="Download Figure",
                             status="success"
                           ),
                           tags$br(),
                           downloadBttn(ns('filterd_data'), 'QC filtered data', color = "royal"),
                           tags$br(),
                           tags$br(),
                           actionBttn(ns("toAnalyze"), "Next Normalize",
                                      style = "unite",color = "royal",
                                      icon = icon("angles-right",class="fa-duotone fa-angles-right"))
                           )
                       )#parameter column
              )
            )
    )
}

preQCServer <- function(id,raw_data) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns


      scdata<-reactive({
        # errors <- c()
         req(raw_data())
         obj <- raw_data()
        #
        # # try to read in file
        # tryCatch(
        #   {
        #     obj <- raw_data()
        #   },
        #   error = function(e) {
        #     errors <- c(errors, "Invalid rds file.")
        #     return(errors)
        #   }
        # )
        #
        # # Validate obj is a SingleCellExperiment object
        # if (!inherits(obj, "SingleCellExperiment")){
        #   errors <- c(errors, "File is not a SingleCellExperiment object")
        #   return(errors)
        # }


        return(obj)
      })

      output$min_gene_per_cell <- renderUI({
        ns <- session$ns
        req(scdata())
        sliderInput(ns('QC_TotalGene_Down'), label = "Min-genes", min= min(scdata()@colData$nFeature_RNA) , max=max(scdata()@colData$nFeature_RNA)/2,value = min(scdata()@colData$nFeature_RNA))
      })

      output$max_gene_per_cell <- renderUI({
        ns <- session$ns
        req(scdata())
        sliderInput(ns('QC_TotalGene_Up'),   label = 'Max-genes', min= min(scdata()@colData$nFeature_RNA) , max=max(scdata()@colData$nFeature_RNA),value = max(scdata()@colData$nFeature_RNA))
      })

      output$min_count <- renderUI({
        ns <- session$ns
        req(scdata())
        sliderInput(ns('QC_TotalCount_Down'), label = 'Min-count',min= min(scdata()@colData$nCount_RNA) , max=max(scdata()@colData$nCount_RNA)/2,value = min(scdata()@colData$nCount_RNA))
      })

      output$max_count <- renderUI({
        ns <- session$ns
        req(scdata())
        sliderInput(ns('QC_TotalCount_Up'), label = 'Max-count',min= min(scdata()@colData$nCount_RNA) , max=max(scdata()@colData$nCount_RNA),value = max(scdata()@colData$nCount_RNA))

      })

      output$mt_perc <- renderUI({
        ns <- session$ns
        req(scdata())
        sliderInput(ns('QC_MT'), label = '%MT cut-off',min= round(min(scdata()@colData$Mito_gene_percent)) , max=round(max(scdata()@colData$Mito_gene_percent))+5,value = round(min(scdata()@colData$Mito_gene_percent))+1)
      })
      output$ribo_perc<- renderUI({
        ns <- session$ns
        req(scdata())
        sliderInput(ns('QC_ribo'), label = '%Ribo cut-off',min= round(min(scdata()@colData$Ribosomal_gene_percent)) , max=round(max(scdata()@colData$Ribosomal_gene_percent))+5,value = round(min(scdata()@colData$Ribosomal_gene_percent))+1)
      })






      vals=reactiveValues()
      scdata_filt<-reactive({
        req(scdata())

        req(input$QC_TotalGene_Down)
        req(input$QC_TotalGene_Up)
        req(input$QC_TotalCount_Down)
        req(input$QC_TotalCount_Up)
        req(input$QC_MT)
        req(input$QC_ribo)


              data=filter_raw_data(sce = scdata(), min_genes_per_cell = input$QC_TotalGene_Down, max_genes_per_cell = input$QC_TotalGene_Up, min_count_cell = input$QC_TotalCount_Down,
                             max_count_cell =input$QC_TotalCount_Up, mt_percnt = input$QC_MT, ribsoml_percnt = input$QC_ribo)

              vals$data
              return(data)

      })
      output$beforeQC_selectFeatureVector = renderUI({
        ns <- session$ns
        req(scdata_filt())
        selectInput(ns("colour_by"),
                    label = "Colour By:",
                    choices = c('All', colnames(scdata_filt()@colData))[1:4],selected='All')
      })


      plot_gene<- reactive({
        p=scater::plotColData(scdata_filt(), y = "nFeature_RNA", x = input$colour_by, colour_by = input$colour_by) + theme(legend.position = "none")
        p+labs(x = input$colour_by , title = "nFeature_RNA", y=NULL) + theme(plot.title = element_text(hjust = 0.5, size=16) , axis.title = element_text(size = 14), axis.text = element_text(size = 12))
      })
      plot_count<-reactive({

        p=scater::plotColData(scdata_filt(), y = "nCount_RNA", x = input$colour_by, colour_by = input$colour_by) +  theme(legend.position = "none")
        p+labs(x = input$colour_by , title = "nCount_RNA", y=NULL) + theme(plot.title = element_text(hjust = 0.5, size=16) , axis.title = element_text(size = 14), axis.text = element_text(size = 12))
      })
      plot_MT<-reactive({

          p=scater::plotColData(scdata_filt(), y="Mito_gene_percent", x=input$colour_by, color_by =input$colour_by )  + theme(legend.position = "none")
          p+labs(x = input$colour_by , title = "Mito_gene_percent", y=NULL) + theme(plot.title = element_text(hjust = 0.5, size=16) , axis.title = element_text(size = 14), axis.text = element_text(size = 12))

      })
      plot_Ribo<- reactive({
        p=scater::plotColData(scdata_filt(), y="Ribosomal_gene_percent",  x=input$colour_by, colour_by=input$colour_by) + theme(legend.position = "none")
        p+labs(x = input$colour_by , title = "Ribosomal_gene_percent", y=NULL) + theme(plot.title = element_text(hjust = 0.5, size=16) , axis.title = element_text(size = 14), axis.text = element_text(size = 12))
      })
      plot_count_vs_gene<-reactive({
        p= scater::plotColData(scdata_filt(), x = "nCount_RNA", y = "nFeature_RNA", colour_by = input$colour_by) +  theme(legend.position = "none")
        p + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12))
      })

      shinyjs::hide("download_box")

      observeEvent(input$beforeQC_submit_sidebarPanel,{
        req(scdata_filt())
        if(class(scdata_filt())!='SingleCellExperiment'){
          shinyjs::hide("download_box")
        }

        #withProgress(message = "Creating plot  ...",{

         # waiter_show(html = spin_ball())
        withProgress(message = 'Calculation in progress', value = 0, {
          for(i in 1:10) {
            incProgress(1/10)
            Sys.sleep(0.5)
          }

          output$beforeQC_totalGene <- renderPlot({
            plot_gene()
          })

          output$beforeQC_totalCount<-renderPlot({
            plot_count()
          })

          output$beforeQC_totalMT<-renderPlot({
            plot_MT()

          })

          output$QC_RiboGene<-renderPlot({
            plot_Ribo()
          })

          output$beforeQC_totalCount_vs_totalGene<-renderPlot({
            plot_count_vs_gene()
          })



        output$sce_obj <- renderTable({
          req(scdata_filt())
          table= c(ncol(scdata_filt()), nrow(scdata_filt())) %>% as.data.frame()
          colnames(table)  <- c("Summary")
          row.names(table) <- c("Total Cells", "Total Genes")

          table

        }, rownames = T, bordered=TRUE, align='c')

        output$pdf <- downloadHandler(
          filename="PreQC.pdf",
          content = function(file){
            pdf(file,width=input$w,height=input$h)

            print(ggpubr::ggarrange(plotlist = list(plot_gene(),plot_count(),plot_MT(),plot_Ribo(),plot_count_vs_gene() )))
            dev.off()
          }
        )
        shinyjs::show("download_box")
        output$tiff <- downloadHandler(
          filename="PreQC.tiff",
          content = function(file){
            tiff(file,width=input$w,height=input$h,units="in",res=input$ppi)

            print(ggpubr::ggarrange(plotlist = list(plot_gene(),plot_count(),plot_MT(),plot_Ribo(),plot_count_vs_gene())))
            dev.off()
          }
        )
        output$png <- downloadHandler(
          filename="PreQC.png",
          content = function(file){
            png(file,width=input$w,height=input$h,units="in",res=input$ppi)
            print(ggpubr::ggarrange(plotlist = list(plot_gene(),plot_count(),plot_MT(),plot_Ribo(),plot_count_vs_gene())))
            dev.off()
          }
        )
        output$jpeg <- downloadHandler(
          filename="PreQC.jpeg",
          content = function(file){
            jpeg(file,width=input$w,height=input$h,units="in",res=input$ppi)
            print(ggpubr::ggarrange(plotlist = list(plot_gene(),plot_count(),plot_MT(),plot_Ribo(),plot_count_vs_gene())))
            dev.off()
          })
        output$filterd_data <- downloadHandler(
          filename = function() {
            paste0( "QC_filtered_data", "-", Sys.Date(), ".rds")
          },
          content = function(file) {
            saveRDS(scdata_filt(), file = file)
          }
        )


        setProgress(value = 1)


        })
      })#End observe plot

      observeEvent(input$toAnalyze,
                   {
                     shinyjs::runjs("$('a[data-value=\"normalization\"]').tab('show');")
                   })


      return(scdata_filt)
       })#main

}

