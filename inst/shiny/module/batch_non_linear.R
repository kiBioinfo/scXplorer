batch_non_linear_UI<-  function(id) {
  ns <- NS(id)
  tagList(
    tabItem(tabName= "normalization",
            fluidPage(
              fluidRow(

                column(width = 10,
                       wellPanel( plotOutput(ns('non_linear_reduction'), width = "100%", height = "600px")%>% withSpinner(color="#0dc5c1",type = 6,size=0.9),
                                  style = "height: auto; width: 100%; border: 3px solid #CEECF5;"
                       )),



                #batch correctionoption


                column(width=2,
                       wellPanel(selectInput(ns("dimM"), "Select Dim reduction Method :",choices = c("tSNE", "Umap"),selected="tSNE"),
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
                                 actionBttn(ns("go"),label="EXEC",style = "jelly",color = "success",icon = icon("sliders")),
                                 tags$br(),
                                 tags$br(),
                                 shinyjs::hidden(actionBttn(ns("DGE_analsis_tab"), label="Continue to DE Analysis",
                                                            block = TRUE,style = "unite",color = "royal", icon = icon("angles-right",class="fa-duotone fa-angles-right"))),
                                 style = " background-color: #CEECF5; border: 3px solid #CEECF5;"
                       )
                )

            )#1st fluidrow
    )#fluidpage
  )#tabItem
  )#taglist
}#ui



batch_non_linear_Server <- function(id,batch_corrct) {
  moduleServer(
    id,
    ## Below is the module function
    function(input, output, session) {


      vals=reactiveValues()
      scdata<-reactive({
        req(batch_corrct())
        batch_corrct()
      })

      cs_data<- reactive({
        req(scdata())
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
        req(cs_data())
        CS.data = ClusterFind(cs_data(),method = input$C_M,runWith = input$dimM,k = input$K, ClusterNum = input$K_cln,PCNum= input$np, cluster.fun=input$algo)
        if(any(c("All", "sizeFactor") %in% colnames(colData(CS.data)))){
          CS.data@colData <- subset(CS.data@colData, select = -c(All, sizeFactor))
        }
        CS.data

      })# Cluster find data
      output$color_non_linear<-renderUI({
        ns <- session$ns
        req(cluster_find_data())
        selectInput(ns("colourby"),
                    label = "Colour By:",
                    choices = colnames(cluster_find_data()@colData)[!colnames(cluster_find_data()@colData) %in% names],  selected = rev(names(cluster_find_data()@colData))[1])
      })

      non_linear_plot<-reactive({
        req(cluster_find_data())
        obj= cluster_find_data()
        obj@colData$cluster <- as.factor(obj@colData$cluster)
        obj@colData$cluster=factor(obj@colData$cluster, levels = levels( obj@colData$cluster) %>% str_sort(numeric = TRUE))

        p=plotReducedDim(obj, dimred =  input$dimM,
                         colour_by = input$colourby,text_by= input$colourby) + ggplot2::ggtitle(label = input$dimM )

        p

      })

      observeEvent(input$go,{


        output$non_linear_reduction<-renderPlot({
          non_linear_plot()
        },res = 96)
        shinyjs::show("DGE_analsis_tab")
      })

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



      return(cluster_find_data)


    }#function session
  )}#server
