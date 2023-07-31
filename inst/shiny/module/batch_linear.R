batch_dim_reduction_UI<-function(id) {
  ns <- NS(id)
  tagList(
    fluidPage(
      fluidRow(style = "height: 78vh; overflow-y: auto;",
               tabBox(width = 12,
                      tabPanel("Linear Dimensional Reduction", value="linear",
                               tabsetPanel(id = "subTab1",
                                     tabPanel("Dim Loadings",value=1,
                                                fluidRow(
                                                  tags$br(),
                                                  tags$br(),
                                                  column(width = 10,
                                                         wellPanel( plotOutput(ns('dotp'),height = "1000px") %>%
                                                                     withSpinner(color="#0dc5c1",type = 6,size=0.9),
                                                                    #Align download button to the bottom left corner
                                                                    fluidRow(
                                                                      column(width = 12, align = "right",
                                                                             download_plot_UI(ns("dim_loadings"))
                                                                      )),
                                                                    style = "height: auto; width: 100%; border: 3px solid #CEECF5;",

                                                                    )),

                                                  column(width=2,
                                                         wellPanel( numericInput(ns("npc"), "Select No. Of PC :",6, min=1, max=50),
                                                                    style = " background-color: #CEECF5; border: 3px solid #CEECF5;"))

                                                )),

                                       tabPanel("PCA Plot",value=2,
                                                fluidRow(
                                                  tags$br(),
                                                  tags$br(),

                                                  column(width = 10,
                                                         wellPanel(
                                                           plotOutput(ns('pcaP'), height=900) %>%
                                                             withSpinner(color="#0dc5c1",type = 6,size=0.9),
                                                           #Align download button to the bottom left corner
                                                           fluidRow(
                                                             column(width = 12, align = "right",
                                                                    download_plot_UI(ns("PCA"))
                                                             )),

                                                            style = "height= auto;width: 100%; border: 3px solid #CEECF5;",

                                                         )),

                                                  column(width=2,
                                                         wellPanel(numericInput(ns("xpc"), "PC (X-Axis)",1, min=1, max=50),numericInput(ns("ypc"), "PC (Y-Axis)",2, min=1, max=50),
                                                                   uiOutput(ns("colPC")),
                                                                   style = " background-color: #CEECF5; border: 3px solid #CEECF5;"
                                                         ))

                                                )),

                                       tabPanel("Heatmap",value=3,
                                                fluidRow(
                                                  tags$br(),
                                                  tags$br(),
                                                  column(width = 10,
                                                         wellPanel(plotOutput(ns('heat'),height = "750px") %>%
                                                                     withSpinner(color="#0dc5c1",type = 6,size=0.9),
                                                                   #Align download button to the bottom left corner
                                                                   fluidRow(
                                                                     column(width = 12, align = "right",
                                                                            download_plot_UI(ns("Heatmap_P"))
                                                                     )),
                                                                   style = "height: auto; width: 100%; border: 3px solid #CEECF5;",

                                                         )),
                                                  column(width=2,
                                                         wellPanel(
                                                           numericInput(ns("dms"), "Select No. Of DIMs to show:",1, min=1, max=50),
                                                           numericInput(ns("ngenes"), "Select No. Of Features to show:",20, min=1, max=Inf),
                                                           style = " background-color: #CEECF5; border: 3px solid #CEECF5;"
                                                         )
                                                  ))),

                                       tabPanel("Pairwise PCA Plot",value=4,
                                                fluidRow(
                                                  tags$br(),
                                                  tags$br(),
                                                  column(width = 10,
                                                         wellPanel(plotOutput(ns('Pairwise_PCA'),height = "600px") %>%
                                                                     withSpinner(color="#0dc5c1",type = 6,size=0.9),
                                                                   #Align download button to the bottom left corner
                                                                   fluidRow(
                                                                     column(width = 12, align = "right",
                                                                            download_plot_UI(ns("Pair_PCA"))
                                                                     )),
                                                                   style = "height: auto; width: 100%; border: 3px solid #CEECF5;",

                                                         )),

                                                  column(width=2,
                                                         wellPanel(numericInput(ns("n_pc"), "Select No. Of PC :",3, min=1, max=50),uiOutput(ns("col_PC")),
                                                                   style = " background-color: #CEECF5; border: 3px solid #CEECF5;"
                                                         ))

                                                )),

                                       tabPanel("Elbow Plot",value=5,
                                                fluidRow(
                                                  tags$br(),
                                                  tags$br(),

                                                  column(width = 10,
                                                         wellPanel(plotOutput( ns('elbow'), height = "600px") %>%
                                                                     withSpinner(color="#0dc5c1",type = 6,size=0.9),
                                                                   #Align download button to the bottom left corner
                                                                   fluidRow(
                                                                     column(width = 12, align = "right",
                                                                            download_plot_UI(ns("elbow_plot"))
                                                                            )),
                                                                   style = "height: auto; width: 100%; border: 3px solid #CEECF5;",

                                                         )),

                                                  column(width=2,
                                                         wellPanel(numericInput(ns("elb"), "Select No. Of PC :",20, min=1, max=50),
                                                                   style = " background-color: #CEECF5; border: 3px solid #CEECF5;"
                                                         ))
                                                )),

                                       selected = 1
                      )),

                      tabPanel("Non Linear Dimensional Reduction", value="non linear",

                               fluidRow(

                                 column(width = 10,
                                        wellPanel(plotOutput(ns('non_linear_reduction'), width = "100%", height = "600px")%>%
                                                     withSpinner(color="#0dc5c1",type = 6,size=0.9),
                                                  #Align download button to the bottom left corner
                                                  fluidRow(
                                                    column(width = 12, align = "right",
                                                           download_plot_UI(ns("Dim_reduction"))
                                                    )),
                                                   style = "height: auto; width: 100%; border: 3px solid #CEECF5;",

                                        )),



                                 #batch correctionoption


                                 column(width=2,
                                        wellPanel(selectInput(ns("dimM"), "Select Dim reduction Method :",choices = c("tSNE", "Umap"),selected="Umap"),
                                                  conditionalPanel(
                                                    condition = "input.dimM=='tSNE'",
                                                    numericInput(ns("np"), "Select No. Of PC :",20, min=1, max=50),
                                                    numericInput(ns("perp"), "Select Perplexity Value :",30, min=1, max=Inf)
                                                    ,  ns = NS(id) ),
                                                  conditionalPanel(
                                                    condition = "input.dimM=='Umap'",
                                                    numericInput(ns("Unp"), "Select No. Of PC :",20, min=1, max=50),  ns = NS(id) ),

                                                  conditionalPanel(
                                                    condition = "input.dimM=='DFmap'",
                                                    numericInput(ns("Dnp"), "Select No. Of PC :",20, min=1, max=50),  ns = NS(id) ),


                                                  selectInput(ns("C_M"),"Select The Cluster Find Method",
                                                              choices = c("K-Means","DBscan","Mclust","Hclust",
                                                                          "Graph", "Hclust_scran", "Kmean_scran"),selected = "Graph"),

                                                  conditionalPanel(
                                                    condition="input.C_M=='DBscan'",  ns = NS(id) ),
                                                  conditionalPanel(
                                                    condition = "input.C_M=='K-Means' || input.C_M=='Mclust' || input.C_M=='Hclust'",
                                                    numericInput(ns("K_cln"),"Choose number of clusters:", 5,min=1, max=Inf),  ns = NS(id)),

                                                  conditionalPanel(
                                                    condition="input.C_M=='Kmean_scran'",
                                                    numericInput(ns("K"),"Choose number of K:", 30,min=1, max=Inf),ns = NS(id) ),

                                                  conditionalPanel(
                                                    condition="input.C_M=='Graph'",
                                                    selectInput(ns("algo"),label = "Choose the Algorithm:",
                                                                choices = c("walktrap", "louvain", "infomap", 'fast_greedy', "label_prop", 'leading_eigen'), selected = "fast_greedy"),ns = NS(id) ),
                                                  uiOutput(ns('color_non_linear')),
                                                  shinyjs::hidden(actionBttn(ns("go"),label="EXEC",style = "jelly",color = "success",icon = icon("sliders"))),
                                                  tags$br(),
                                                  tags$br(),
                                                  shinyjs::hidden(actionBttn(ns("DGE_analsis_tab"), label="Next DE Analysis",
                                                             block = TRUE,style = "unite",color = "royal", icon = icon("angles-right",class="fa-duotone fa-angles-right"))),
                                                  style = " background-color: #CEECF5; border: 3px solid #CEECF5;"
                                        )
                                            )








                               )# nonlinear fluidRow
                      )#nonlinear tabPanel
               )#tabBox
      ))
  )}#UI end

batch_dim_reduction_Server <- function(id,batch_corrct) {
  moduleServer(
    id,
    ## Below is the module function
    function(input, output, session) {
      ns <- session$ns
      vals=reactiveValues()
      scdata<-reactive({
        req(batch_corrct())
        return(batch_corrct())

      })

      #colour_input
      output$colPC<-renderUI({
        ns <- session$ns
        req(scdata())
        selectInput(ns("colour_by"),
                    label = "Colour By:",
                    choices =  colnames(scdata()@colData)[!colnames(scdata()@colData) %in% names],
                    selected =  rev(names(scdata()@colData))[1])
      })
      output$col_PC<-renderUI({
        ns <- session$ns
        req(scdata())
        selectInput(ns("colourP"),
                    label = "Colour By:",
                    choices =  colnames(scdata()@colData)[!colnames(scdata()@colData) %in% names],
                    selected =  rev(names(scdata()@colData))[1])
      })

      #loadings plot
      lodings_p<- reactive({

        validate(
          need(!is.na(input$npc) , "It should  be positive number")
        )
        req(scdata())
        p<-dim_lodingsViz(scdata(),ndims = input$npc)
        p
      })

      #PCA Plot
      pca_plot<- reactive({
        req(scdata())
        req(input$colour_by)
        validate(
          need(!is.na(input$xpc) , "It should  be positive number")
        )
        validate(
          need(!is.na(input$ypc), "It should  be positive number")
        )
        p=singleCellTK::plotPCA(scdata(),pcX=paste0("PC",input$xpc),
                                pcY = paste0("PC",input$ypc), reducedDimName = "BEPCA",colorBy = input$colour_by) + theme_cowplot()
        p
      })

      #Heatmap
      heatmap_plot<- reactive({
        req(scdata())
        validate(
          need(!is.na(input$ngenes), "It should  be positive number")
        )
        validate(
          need(!is.na(input$dms), "It should  be positive number")
        )
        p=dim_heatmap(scdata(), ndims=input$dms,nfeatures=input$ngenes)
        p
      })
      #Parwise PCA
      pairwise_pca_plot <- reactive({
        req(scdata())
        req(input$colourP)
        validate(
          need(!is.na(input$n_pc), "It should at least 2 dimension")

        )
        p=scater::plotReducedDim(scdata(),dimred="BEPCA", ncomponents=input$n_pc,colour_by=input$colourP)
        p
      })
      #Elbow
      elbow_plt <- reactive({
        req(scdata())
        validate(
          need(!is.na(input$elb), "It should  be positive number")

        )
        p= ElbowPlot(scdata(),ndims = input$elb)
        p
      })

      output$dotp<-renderPlot({
        lodings_p()

      }, res = 96)#dotp

      output$pcaP<-renderPlot({
        pca_plot()
      },res = 96)

      output$heat<-renderPlot({
        heatmap_plot()
      },res = 96)

      output$Pairwise_PCA<- renderPlot({
        pairwise_pca_plot()
      },res = 96)

      output$elbow<- renderPlot({
        elbow_plt()
      },res = 96)



      #Download plots (called download_plot_Server() from plot_download.R )

      download_plot_Server("dim_loadings", input_data = lodings_p , name ="Batch_Dim_loadings")
      download_plot_Server("PCA", input_data = pca_plot , name ="PCA_Plot")
      download_plot_Server("Heatmap_P", input_data = heatmap_plot , name ="Batch_Heatmap")
      download_plot_Server("Pair_PCA", input_data = pairwise_pca_plot , name ="Batch_pairwise_pca_plot")
      download_plot_Server("elbow_plot", input_data = elbow_plt , name ="Elbow_plot")
      download_plot_Server("Dim_reduction", input_data = non_linear_plot , name =paste0("Batch_Dim_reduction_",input$dimM,"_",input$C_M ))



      #Non Linear

      cs_data<- reactive({
        req(scdata())
        withProgress(message = 'Dimension Reduction in progress......',
                     detail = 'This may take a while...', value = 0, {
        if(input$dimM=="tSNE")
        {
          CS.data = scater::runTSNE(scdata(),dimred= "BEPCA", n_dimred=input$np, perplexity = input$perp, name= "tSNE")
          # return(CS.data)
        }
        else
        {
          CS.data = scater::runUMAP(scdata(),  dimred= "BEPCA", n_dimred=input$np, name= "Umap")
          # return(CS.data)
        }
        shinyjs::show("go")
        return(CS.data)
        })
      })

      # Modal to plot optimum number of clusters in the dataset
      dataModal <- function(failed = FALSE) {

        modalDialog(
          span(h5("Optimal number of clusters")),
          tags$hr(),
          uiOutput(ns("data_for_clusterSearch")),
          numericInput(ns("n_PCs_for_clusterSearch"), "Select No. Of PC :",20, min=1, max=50),
          plotOutput(ns("plotClusterSearch"))%>% withSpinner(color="#0dc5c1",type = 6,size=0.9),
          #renderTable(tail(d,20),rownames=FALSE),
          easyClose=TRUE,
          footer = tagList(
            modalButton("Close")
          )
        )
      }
      #Select data as input
      output$data_for_clusterSearch <- renderUI({
        req(cs_data())
        shiny::selectInput(ns("runwith"),
                           label = "Run with:",
                           choices = c("VGcounts", SingleCellExperiment::reducedDimNames(cs_data())),
                           selected = rev(reducedDimNames(cs_data()))[1], selectize = F)
      })
      #Print optimum clusters in dataset
      output$plotClusterSearch <- renderPlot({
        req(cs_data())
        req(input$runwith)
        req(input$n_PCs_for_clusterSearch)
        ClusterNumSearch(cs_data(),runWith = input$runwith , PCNum = input$n_PCs_for_clusterSearch)
      })


      #Show modal for optimum cluster identification
      observeEvent(input$C_M,{
        req(input$C_M)
        if(input$C_M=='K-Means' || input$C_M=='Mclust' || input$C_M=='Hclust'){
          shiny::showModal(dataModal())
        }

      })
      # Cluster find in dataset
      cluster_find_data<- reactive({
        req(cs_data())
        withProgress(message = 'Finding clusters in progress......',
                     detail = 'This may take a while...', value = 0, {

                       CS.data = ClusterFind(cs_data(),method = input$C_M,runWith = input$dimM,
                                             k = input$K, ClusterNum = input$K_cln,PCNum= input$np,cluster.fun = input$algo)
                       if(any(c("sizeFactor") %in% colnames(colData(CS.data)))){
                         CS.data@colData <- subset(CS.data@colData, select = -c(All, sizeFactor))
                       }

                       # shinyjs::show("go")
                       return(CS.data)
                     })
      })
      output$color_non_linear<-renderUI({
        ns <- session$ns
        req(cluster_find_data())
        selectInput(ns("colourby"),
                    label = "Colour By:",
                    choices =  colnames(cluster_find_data()@colData)[!colnames(cluster_find_data()@colData) %in% names],  selected = rev(names(cluster_find_data()@colData))[1])
      })

      non_linear_plot<-reactive({
        req(cluster_find_data())
        obj= cluster_find_data()
        obj@colData$cluster <- as.factor(obj@colData$cluster)
        obj@colData$cluster=factor(obj@colData$cluster, levels = levels( obj@colData$cluster) %>% str_sort(numeric = TRUE))

        p=plotReducedDim(obj, dimred =  input$dimM,
                         colour_by = input$colourby,text_by= input$colourby) + ggplot2::ggtitle(label = input$dimM )  + ggplot2::theme(plot.title = element_text(hjust = 0.5))

        p

      })

      observeEvent(input$go,{


        output$non_linear_reduction<-renderPlot({
          non_linear_plot()
        },res = 96)

        shinyjs::show("DGE_analsis_tab")
      })
      #show/hide DGE_Analysis tab

      observeEvent(input$DGE_analsis_tab,
                   {

                     showModal(modalDialog(
                       title = "Somewhat important message",
                       "By default data without batch correction has been selected
                       Please click on the use batch correction button!",
                       easyClose = TRUE,
                       footer = tagList(
                         modalButton("Ok")
                        )
                     ))
                     shinyjs::show(selector = "a[data-value=\"DGE_Analysis\"]")
                     shinyjs::runjs("$('a[data-value=\"DGE_Analysis\"]').tab('show');")


                     })



      # isolate( return(list(scdata=reactive({cluster_find_data()}))))

return(cluster_find_data)
    })}
