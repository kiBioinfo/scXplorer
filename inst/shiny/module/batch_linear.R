batch_dim_reduction_UI<-function(id) {
  ns <- NS(id)
  tagList(
    fluidPage(
      fluidRow(style = "height: 78vh; overflow-y: auto;",
               tabBox(width = 12,
                      tabPanel("Linear Dimensional Reduction", value="linear", tabsetPanel(id = "subTab1",
                                                                                           tabPanel("Dim Loadings",value=1,
                                                                                                    fluidRow(
                                                                                                      tags$br(),
                                                                                                      tags$br(),
                                                                                                      column(width = 10,
                                                                                                             wellPanel( plotOutput(ns('dotp'),height = "1000px")%>% withSpinner(color="#0dc5c1",type = 6,size=0.9)
                                                                                                                        , style = "height: auto; width: 100%; border: 3px solid #CEECF5; text-align: end;",

                                                                                                                        download_plot_UI(ns("dim_loadings")))),

                                                                                                      column(width=2,
                                                                                                             wellPanel( numericInput(ns("npc"), "Select No. Of PC :",6, min=1, max=50),
                                                                                                                        style = " background-color: #CEECF5; border: 3px solid #CEECF5;"))

                                                                                                    )),

                                                                                           tabPanel("PCA Plot",value=2,
                                                                                                    fluidRow(
                                                                                                      tags$br(),
                                                                                                      tags$br(),

                                                                                                      column(width = 10,
                                                                                                             wellPanel(plotOutput(ns('pcaP'), height=900),
                                                                                                                       style = "height= auto;width: 100%; border: 3px solid #CEECF5;text-align: end;",
                                                                                                                       download_plot_UI(ns("PCA"))
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
                                                                                                             wellPanel(plotOutput(ns('heat'),height = "750px"),
                                                                                                                       style = "height: auto; width: 100%; border: 3px solid #CEECF5; text-align: end;",
                                                                                                                       download_plot_UI(ns("Heatmap_P"))
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
                                                                                                             wellPanel(plotOutput(ns('Pairwise_PCA'),height = "600px") %>% withSpinner(color="#0dc5c1",type = 6,size=0.9),
                                                                                                                       style = "height: auto; width: 100%; border: 3px solid #CEECF5; text-align: end;",
                                                                                                                       download_plot_UI(ns("Pair_PCA"))
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
                                                                                                             wellPanel(plotOutput( ns('elbow'), height = "600px") %>% withSpinner(color="#0dc5c1",type = 6,size=0.9),
                                                                                                                       style = "height: auto; width: 100%; border: 3px solid #CEECF5; text-align: end;",
                                                                                                                       download_plot_UI(ns("elbow_plot"))
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
                                        wellPanel( plotOutput(ns('non_linear_reduction'), width = "100%", height = "600px")%>% withSpinner(color="#0dc5c1",type = 6,size=0.9),
                                                   style = "height: auto; width: 100%; border: 3px solid #CEECF5; text-align: end;",
                                                   download_plot_UI(ns("Dim_reduction"))
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


                                                  selectInput(ns("C_M"),"Select The Cluster Find Method",choices = c("K-Means","DBscan","Mclust","Hclust","Graph", "Hclust_cran", "Kmean_scran"),selected = "Graph"),

                                                  conditionalPanel(
                                                    condition="input.C_M=='DBscan'",  ns = NS(id) ),
                                                  conditionalPanel(
                                                    condition = "input.C_M=='K-Means' || input.C_M=='Mclust' || input.C_M=='Hclust'",
                                                    numericInput(ns("K_cln"),"Choose number of clusters:", 5,min=1, max=Inf),  ns = NS(id) ),

                                                  conditionalPanel(
                                                    condition="input.C_M=='Kmean_scran'",
                                                    numericInput(ns("K"),"Choose number of K:", 30,min=1, max=Inf),ns = NS(id) ),

                                                  conditionalPanel(
                                                    condition="input.C_M=='Graph'",
                                                    selectInput(ns("algo"),label = "Choose the Algorithm:",
                                                                choices = c("walktrap", "louvain", "infomap", 'fast_greedy', "label_prop", 'leading_eigen'), selected = "fast_greedy"),ns = NS(id) ),
                                                  uiOutput(ns('color_non_linear')),
                                                  actionBttn(ns("go"),"submit", block = TRUE, style = "unite",color = "royal"),
                                                  tags$br(),
                                                  shinyjs::hidden(actionBttn(ns("DGE_analsis_tab"), label="Continue to DE Analysis",
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
        return(batch_corrct())

      })





      # #PCA
      # scdata<- reactive({
      #   used= "BEVGcounts"
      #   scater::runPCA(scdata(), exprs_values=used, altexp= used, name= "BEPCA")
      # })
      #colour_input
      output$colPC<-renderUI({
        ns <- session$ns
        selectInput(ns("colour_by"),
                    label = "Colour By:",
                    choices =  colnames(scdata()@colData)[!colnames(scdata()@colData) %in% names],selected =  rev(names(scdata()@colData))[1])
      })
      output$col_PC<-renderUI({
        ns <- session$ns
        selectInput(ns("colourP"),
                    label = "Colour By:",
                    choices =  colnames(scdata()@colData)[!colnames(scdata()@colData) %in% names],selected =  rev(names(scdata()@colData))[1])
      })

      #loadings plot
      lodings_p<- reactive({
        # If the length of the input is 0
        # (i.e. nothing is selected),we show
        # a feedback to the user in the form of a text
        # If the length > 0, we remove the feedback.
        if (length(input$npc) == 0){
          showFeedbackWarning(
            inputId = "npc",
            text = "Select at least one PC"
          )
        } else {
          shinyFeedback::hideFeedback("npc")
        }
        # req() allows to stop further code execution
        # if the condition is not a truthy.
        # Hence if input$npc is NULL, the computation
        # will be stopped here.
        req(input$npc)
        p<-dim_lodingsViz(scdata(),ndims = input$npc)
        p
      })

      #PCA Plot
      pca_plot<- reactive({

        p=singleCellTK::plotPCA(scdata(),pcX=paste0("PC",input$xpc),pcY = paste0("PC",input$ypc), reducedDimName = "BEPCA",colorBy = input$colour_by) + theme_cowplot()
        p
      })

      #Heatmap
      heatmap_plot<- reactive({
        p=dim_heatmap(scdata(), ndims=input$dms,nfeatures=input$ngenes)
        p
      })
      #Parwise PCA
      pairwise_pca_plot <- reactive({
        p=scater::plotReducedDim(scdata(),dimred="BEPCA", ncomponents=input$n_pc,colour_by=input$colourP)
        p
      })
      #Elbow
      elbow_plt <- reactive({
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
        if(input$dimM=="tSNE")
        {
          CS.data = scater::runTSNE(scdata(),dimred= "BEPCA", n_dimred=input$np, perplexity = input$perp, name= "tSNE")
          return(CS.data)
        }
        else
        {
          CS.data = scater::runUMAP(scdata(),  dimred= "BEPCA", n_dimred=input$np, name= "Umap")
          return(CS.data)
        }
      })

      cluster_find_data<- reactive({
        CS.data = ClusterFind(cs_data(),method = input$C_M,runWith = input$dimM,k = input$K, ClusterNum = input$K_cln,PCNum= input$np,cluster.fun=input$algo)
        if(any(c("All", "sizeFactor") %in% colnames(colData(CS.data)))){
          CS.data@colData <- subset(CS.data@colData, select = -c(All, sizeFactor))
        }

        CS.data

      })# Cluster find data
      output$color_non_linear<-renderUI({
        ns <- session$ns
        selectInput(ns("colourby"),
                    label = "Colour By:",
                    choices =  colnames(cluster_find_data()@colData)[!colnames(cluster_find_data()@colData) %in% names],  selected = rev(names(cluster_find_data()@colData))[1])
      })

      non_linear_plot<-reactive({
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


        output$stat <- renderPrint({

          cluster_find_data()

        })
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
