dim_reduction_UI<-function(id) {
  ns <- NS(id)
  tagList(
    fluidPage(
      fluidRow(style = "height: 78vh; overflow-y: auto;",
               shinyFeedback::useShinyFeedback(),
      tabBox(width = 12,
             tabPanel("Linear Dimensional Reduction", value="linear",
                      tabsetPanel(id = "subTab1",

              tabPanel("Dim Loadings",value=1,
                         fluidRow(
                          tags$br(),
                          tags$br(),
                          column(width = 10,
                                 wellPanel(
                                   plotOutput(ns('dotp'), height=1000, width = '100%')%>%
                                     withSpinner(color="#0dc5c1",type = 6,size=0.9),
                                   #Align download button to the bottom left corner
                                   fluidRow(
                                     column(width = 12, align = "right",
                                           download_plot_UI(ns("dim_loadings"))
                                   )),
                                  style = "height= auto;width: 100%; border: 3px solid #CEECF5;"
                                 )),

                          column(width=2,
                            wellPanel( numericInput(ns("npc"), "Select No. Of PC :",1, min=1, max=50),
                                      style = " background-color: #CEECF5; border: 3px solid #CEECF5;"))

                               )),

                tabPanel("PCA Plot",value=2,
                         fluidRow(
                           tags$br(),
                           tags$br(),

                           column(width = 10,
                                  wellPanel(
                                    plotOutput(ns('pcaP'), height="1000px", width = '100%')%>%
                                      withSpinner(color="#0dc5c1",type = 6,size=0.9),
                                    #Align download button to the bottom left corner
                                    fluidRow(
                                      column(width = 12, align = "right",
                                             download_plot_UI(ns("PCA"))
                                      )),
                                        # style = "height= auto;width: 100%; border: 3px solid #CEECF5; text-align: end;",
                                    style = "height= auto;width: 100%; border: 3px solid #CEECF5;"

                                  )),

                           column(width=2,
                                  wellPanel(numericInput(ns("xpc"), "PC (X-Axis)",1, min=1, max=50),
                                            numericInput(ns("ypc"), "PC (Y-Axis)",2, min=1, max=50),
                                            uiOutput(ns("colPC")),
                                             style = " background-color: #CEECF5; border: 3px solid #CEECF5;"
                                      ))

                               )),

                tabPanel("Heatmap",value=3,
                         fluidRow(
                           tags$br(),
                           tags$br(),
                         column(width = 10,
                                wellPanel(
                                  plotOutput(ns('heat'),height = "1000px", width = '100%')%>%
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
                                wellPanel(
                                  plotOutput(ns('Pairwise_PCA'),height = "1000px", width = '100%') %>%
                                          withSpinner(color="#0dc5c1",type = 6,size=0.9),
                                  fluidRow(
                                    column(width = 12, align = "right",
                                           download_plot_UI(ns("Pair_PCA"))
                                    )),
                                    style = "height: auto; width: 100%; border: 3px solid #CEECF5;",

                                  )),

                             column(width=2,
                                    wellPanel(numericInput(ns("n_pc"), "Select No. Of PC :",3, min=1, max=50),
                                              uiOutput(ns("col_PC")),
                                              style = " background-color: #CEECF5; border: 3px solid #CEECF5;"
                                              ))

                         )),

                tabPanel("Elbow Plot",value=5,
                         fluidRow(
                         tags$br(),
                         tags$br(),

                         column(width = 10,
                                wellPanel(plotOutput( ns('elbow'), height = "1000px", width = '100%') %>%
                                            withSpinner(color="#0dc5c1",type = 6,size=0.9),
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
                                   wellPanel( plotOutput(ns('non_linear_reduction'), width = "100%", height = "1000px") %>%
                                              withSpinner(color="#0dc5c1",type = 6,size=0.9),
                                              fluidRow(
                                                column(width = 12, align = "right",
                                                       download_plot_UI(ns("Dim_reduction"))
                                                )),
                                              style = "height: auto; width: 100%; border: 3px solid #CEECF5;",

                                              )),



                        #batch correctionoption


                            column(width=2,
                                  wellPanel(selectInput(ns("dimM"), "Select Dim reduction Method :",choices = c("tSNE", "Umap"),
                                                        selected="tSNE"),
                                            numericInput(ns("np"), "Select No. Of PC :",20, min=1, max=50),
                                   conditionalPanel(
                                     condition = "input.dimM=='tSNE'",
                                     numericInput(ns("perp"), "Select Perplexity Value :",30, min=1, max=Inf)
                                     ,  ns = NS(id) ),

                                     selectInput(ns("C_M"),"Select The Cluster Find Method",
                                               choices = c("K-Means","DBscan","Mclust","Hclust","Graph", "Hclust_scran", "Kmean_scran"),
                                               selected = "Graph"),

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
                                                 choices = c("walktrap", "louvain", "infomap",
                                                             'fast_greedy', "label_prop", 'leading_eigen'),
                                                 selected = "fast_greedy"),ns = NS(id) ),
                                   uiOutput(ns('color_non_linear')),
                                   shinyjs::hidden(actionBttn(ns("go"),label="EXEC",style = "jelly",color = "success",icon = icon("sliders"))),
                                   style = " background-color: #CEECF5; border: 3px solid #CEECF5;"
                                   ),
                                  tags$br(),
                                  shinyjs::hidden(wellPanel(id=ns("Skip_box"),
                                                            actionBttn(ns("batch_corcn_tab"), label="Continue to Batch Correction ",block = TRUE,
                                                                                                             style = "unite",color = "royal", icon = icon("angles-right",class="fa-duotone fa-angles-right")),
                                                            tags$br(),
                                                            actionBttn(ns("skip_corcn_tab"),label= "Skip Batch Correction ", block = TRUE,
                                                                                                               style = "unite",color = "royal", icon = icon("angles-right",class="fa-duotone fa-angles-right")),
                                                            tags$br(),
                                                            downloadBttn(ns('dim_reduction_data'), '.rds dim rdxn data', block = TRUE,style = "unite",color = "royal" ),
                                                            style = " background-color: #CEECF5; border: 3px solid #CEECF5;"
                                  ))


                                   )

                            )# nonlinear fluidRow
                      )#nonlinear tabPanel
      )#tabBox
      ))
    )}#UI end

dim_reduction_Server <- function(id,normalization_data) {
  moduleServer(
    id,
    ## Below is the module function
    function(input, output, session) {
      ns <- session$ns
      vals=reactiveValues()
      scdata<-reactive({
        req(normalization_data())
        obj <- normalization_data()
        return(obj)
      })



    #PCA Diemsional reduction
    sc_data<- reactive({
      withProgress(message = 'Dimension Reduction in progress......',
                   detail = 'This may take a while...', value = 0, {
      req(scdata())
      scater::runPCA(scdata(), exprs_values=  "VGcounts", altexp= "VGcounts")
     })
    })
      #colour_input
     output$colPC<-renderUI({
       req(scdata())
        ns <- session$ns
        selectInput(ns("colour_by"),
                    label = "Colour By:",
                    choices = colnames(scdata()@colData)[!colnames(scdata()@colData) %in% names],
                    selected = colnames(scdata()@colData)[1])
      })
      output$col_PC<-renderUI({
        req(scdata())
        ns <- session$ns
        selectInput(ns("colourP"),
                    label = "Colour By:",
                    choices = colnames(scdata()@colData)[!colnames(scdata()@colData) %in% names],
                    selected =  rev(names(scdata()@colData))[1])
      })

      #loadings plot
      lodings_p<- reactive({
        # req(input$npc)
        # If the length of the input is 0
        # (i.e. nothing is selected),we show
        # a feedback to the user in the form of a text
        # # If the length > 0, we remove the feedback.
        # if (length(input$npc) == 0){
        #   showFeedbackWarning(
        #     inputId = ns("npc"),
        #     text = "Select at least one PC"
        #   )
        # } else {
        #   shinyFeedback::hideFeedback("npc")
        # }
        # req() allows to stop further code execution
        # if the condition is not a truthy.
        # Hence if input$npc is NULL, the computation
        # will be stopped here.
        validate(
          need(!is.na(input$npc) , "It should  be positive number")
        )
        req(sc_data())
        p<-dim_lodingsViz(sc_data(),ndims = input$npc)
        p
      })

      #PCA Plot
      pca_plot<- reactive({
      req(sc_data())
          validate(
          need(!is.na(input$xpc) , "It should  be positive number")
          )
          validate(
            need(!is.na(input$ypc), "It should  be positive number")
          )
        p=singleCellTK::plotPCA(sc_data(),pcX=paste0("PC",input$xpc),
                                pcY = paste0("PC",input$ypc), reducedDimName = "PCA",colorBy = input$colour_by) +
                                theme_cowplot()
        p
      })

      #Heatmap
      heatmap_plot<- reactive({
        req(sc_data())
        validate(
          need(!is.na(input$ngenes), "It should  be positive number")
        )
        validate(
          need(!is.na(input$dms), "It should  be positive number")
        )
        p=dim_heatmap(sc_data(), ndims=input$dms,nfeatures=input$ngenes)
        p
      })
      #Parwise PCA
      pairwise_pca_plot <- reactive({
        req(sc_data())
        validate(
          need(!is.na(input$n_pc), "It should at least 2 dimension")

        )
        p=scater::plotReducedDim(sc_data(),dimred="PCA", ncomponents=input$n_pc,colour_by=input$colourP)
        p
      })
      #Elbow
      elbow_plt <- reactive({
        req(sc_data())
        validate(
          need(!is.na(input$elb), "It should  be positive number")

        )
        p= ElbowPlot(sc_data(),ndims = input$elb)
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

      download_plot_Server("dim_loadings", input_data = lodings_p , name ="Dim_loadings")
      download_plot_Server("PCA", input_data = pca_plot , name ="PCA_Plot")
      download_plot_Server("Heatmap_P", input_data = heatmap_plot , name ="Heatmap")
      download_plot_Server("Pair_PCA", input_data = pairwise_pca_plot , name ="pairwise_pca_plot")
      download_plot_Server("elbow_plot", input_data = elbow_plt , name ="Elbow_plot")
      download_plot_Server("Dim_reduction", input_data = non_linear_plot , name =paste0("Dim_reduction_",input$dimM,"_",input$C_M ))


      #Non Linear

      # Umap/tSNE dimensional reduction
      cs_data<- reactive({
        req(sc_data())
        # withProgress(message = 'Dimension Reduction in progress......',
        #              detail = 'This may take a while...', value = 0, {

        validate(
          need(!is.na(input$np), "It should  be positive number")

        )
        if(input$dimM=="tSNE")
        {
          validate(
            need(!is.na(input$perp), "It should  be positive number")

          )

          CS.data = DimReduction_tsne(sc_data(),  PCNum = input$np, perplexity = input$perp)
          # return(CS.data)
        }
        else if((input$dimM=="Umap"))
        {

          CS.data = DimReduction_umap(sc_data(),  PCNum = input$np)
          # return(CS.data)
        }
        # shinyjs::show("go")
        return(CS.data)
        #})
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


      # find Clusters in dataset
      cluster_find_data<- reactive({
        req(sc_data())
        withProgress(message = 'Finding clusters in progress......',
                     detail = 'This may take a while...', value = 0, {
         validate(
           need(!is.na(input$K_cln) , "It should  be positive number")

         )


        CS.data = ClusterFind(cs_data(),method = input$C_M,runWith = input$dimM,
                              k = input$K, ClusterNum = input$K_cln,PCNum= input$np,cluster.fun = input$algo)
        if(any(c("sizeFactor") %in% colnames(colData(CS.data)))){
          CS.data@colData <- subset(CS.data@colData, select = -c(All, sizeFactor))
        }

          shinyjs::show("go")
          return(CS.data)
        })
      })
      #colour by column
      output$color_non_linear<-renderUI({
        req(cluster_find_data())
        ns <- session$ns
        selectInput(ns("colourby"),
                    label = "Colour By:",
                    choices = colnames(cluster_find_data()@colData)[!colnames(cluster_find_data()@colData) %in% names],
                    selected = rev(names(cluster_find_data()@colData))[1])
      })
      #Plot reduce dim after clustering
      non_linear_plot<-reactive({
        req(cluster_find_data())
        validate(
          need(inherits(cluster_find_data(),"SingleCellExperiment"), "Please select a data set")
        )
        if(class(cluster_find_data())!= 'SingleCellExperiment'){
          shinyjs::hide("go")

        }
        obj= cluster_find_data()
        obj@colData$cluster <- as.factor(obj@colData$cluster)
        obj@colData$cluster=factor(obj@colData$cluster, levels = levels( obj@colData$cluster) %>% str_sort(numeric = TRUE))


          p=scater::plotReducedDim(obj, dimred =  input$dimM,
                         colour_by = input$colourby,text_by= input$colourby) +
                         ggplot2::ggtitle(label = input$dimM ) + ggplot2::theme(plot.title = element_text(hjust = 0.5))

        p

      })

      observeEvent(input$go,{
        req(input$C_M)
        req(non_linear_plot())
        if(any(typeof(non_linear_plot())=='list')){
          shinyjs::show("Skip_box")
        }
        else{
          shinyjs::hide("Skip_box")
        }

        output$non_linear_reduction<-renderPlot({
          non_linear_plot()
        },res = 96)

      })


      #Download dataset with clustering information
      output$dim_reduction_data <- downloadHandler(
        filename = function() {
          paste0(input$dimM,"_",input$C_M,".rds")
        },
        content = function(file) {
          saveRDS(cluster_find_data(), file = file)
        }
      )



      #show/hide batch_correction tab

      observeEvent(input$batch_corcn_tab,
                   {
                     shinyjs::show(selector = "a[data-value=\"batch_correction\"]")
                     shinyjs::runjs("$('a[data-value=\"batch_correction\"]').tab('show');")
                   })

      observeEvent(input$skip_corcn_tab,
                   {
                     shinyjs::hide(selector = "a[data-value=\"batch_correction\"]")
                     shinyjs::runjs("$('a[data-value=\"DGE_Analysis\"]').tab('show');")

                   })

     # return(list(scdata=reactive({c(cluster_find_data())})))
    return(cluster_find_data)
      })}
