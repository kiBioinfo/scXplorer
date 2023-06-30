de_analysis_UI<-function(id) {
  ns <- NS(id)
  tagList(
    fluidPage(
      useShinyjs(),
      fluidRow(style = "height: 78vh; overflow-y: auto;",
               tabBox(width = 12,
                       tabsetPanel(id = "subTab1",

                           tabPanel("Differential Gene Expression Analysis",value=1,
                                    tags$br(),
                                    fluidRow(

                                      column(width = 10,
                                            wellPanel(
                                             DT::dataTableOutput(ns('DEG_stats'),width = "100%", height = "auto") %>% withSpinner(color="#0dc5c1",type = 6,size=0.9 ),
                                             style = "height: auto; border: 3px solid #CEECF5;text-align: end;",
                                             boxDropdown(downloadBttn(
                                               outputId = ns("DEG_table"),
                                               label="DEG Table",
                                               style = "jelly",
                                               color = "success",
                                               size='sm',
                                               block=TRUE
                                             )
                                             , icon=icon("download"))
                                             )


                                            ),

                                      column(2,
                                             wellPanel(
                                               actionBttn(ns("Go_to_batch_correction_tab"), strong("Re-do Batch Correction"),style = "material-flat", color = "royal"), # To control previous files
                                               tags$br(),
                                               tags$br(),
                                               materialSwitch(ns("use_batch_corrected_data"), "Use batch corrected data:", value = FALSE, status = "warning"), #Control data input


                                                    selectInput(ns("DE_Analysis"), choices = c("pairwise",
                                                                                               "One Vs Rest" ="All"), label = "Select DEG Analysis", selected = "pairwise"),


                                                    uiOutput(ns("comparision")),

                                               conditionalPanel(
                                                 condition = "input.DE_Analysis=='pairwise' " ,
                                                 uiOutput(ns("group_1")),
                                                 uiOutput(ns("group_2")),  ns = NS(id)
                                               ),





                                                    selectInput(ns("method_DEGA"), choices = c("limma","edger","deseq","","mast","monocle","desingle"), label = "Choose method to identify DEG:", selected = "limma"),


                                               conditionalPanel(
                                                 condition = "input.DE_Analysis=='All' " ,

                                                 uiOutput(ns("deg_list")),  ns = NS(id)
                                               ),


                                                    numericInput(ns("fold_change"), label = "Choose the log2 fold change value:", 1.2, min = 0.1, max = 20),
                                                    numericInput(ns("adj_p_value"), label = "Choose the adjuste P value:  ", 0.05, min = 0.00000000000000001, max = 0.049),

                                               style = " background-color: #CEECF5; border: 3px solid #CEECF5; "
                                                    #actionBttn(ns("plot_option"), label = "Volcano Plot", status = "success"),
                                             # column 3
                                             )#welpanel
                                      )# main column

                                    )
                                      ),

                           tabPanel("Volcano Plot",value=2,
                                    column(width = 12,

                                           plotOutput(ns('pairwise_plot'), height = "800px")%>% withSpinner(color="#0dc5c1",type = 6,size=0.9),
                                           style = "height: auto; width: 100%; border: 3px solid #CEECF5; text-align: end;",
                                           download_plot_UI(ns("Volcano_Plt"))

                                    )),

                           tabPanel("Heatmap",value=3,
                                    fluidRow(
                                    column(width = 10,

                                           plotOutput(ns('heatmap'), height = "800px")%>% withSpinner(color="#0dc5c1",type = 6,size=0.9),
                                           style = "height: auto; width: 100%; border: 3px solid #CEECF5; text-align: end;",
                                           download_plot_UI(ns("deg_heatmap"))

                                           ),

                                       column(width = 2  ,
                                              uiOutput(ns("group_by")),
                                              uiOutput(ns("top_genes")),
                                              uiOutput(ns("used"))
                                    ))),


                           selected = 1
                      )
               )#tab Box

      ) #main fluidRow

    )# main fluidPage
  )#tagList
}#UI

de_analysis_Server <-function(id,dim_reduction_data , dim_reduction_batch_correction) {
  moduleServer(
    id,
    ## Below is the module function
    function(input, output, session) {
      ns <- session$ns
      vals=reactiveValues()

      #switch to batch correction tab
      observeEvent(input$Go_to_batch_correction_tab,
                   {
                     showModal(modalDialog(title = "Redo the batch correction or use batch corrected data",

                                           easyClose = TRUE,     footer = tagList(
                                             materialSwitch(ns("start_batch_correction"), label = "Redo Batch Correction:", value = FALSE, status="warning"),
                                             modalButton("Cancel"))

                     ))

                   })

      observeEvent(input$start_batch_correction, {
        if(input$start_batch_correction==TRUE)
        {
          removeModal()
          shinyjs::show(selector = "a[data-value=\"batch_correction\"]")
          shinyjs::runjs("$('a[data-value=\"batch_correction\"]').tab('show');")

        }

        else{

          shinyjs::hide(selector = "a[data-value=\"batch_correction\"]")
          shinyjs::runjs("$('a[data-value=\"batch_correction\"]').tab('hide');")
        }
      })




     # Get the dim reduction(batch corrected/not corrected) data
      scdata<-reactive({
        req(dim_reduction_data())
      if(input$use_batch_corrected_data){

        data = dim_reduction_batch_correction()
        }
      else{
          data= dim_reduction_data()
      }
        if(!'batch' %in% colnames(data@colData)){
          shinyjs::hide("use_batch_corrected_data")
          shinyjs::hide('Go_to_batch_correction_tab')
        }
        else{
          shinyjs::show("use_batch_corrected_data")
          shinyjs::show('Go_to_batch_correction_tab')
        }
        return(data)
      }
      )



      output$comparision<-renderUI({
        ns <- session$ns
        req(scdata())
        validate(
          need(inherits(scdata(),"SingleCellExperiment"), "Please select a data set")
        )
        selectInput(ns("comparision_vector"),
                    label = "Comparision By:",
                    choices = c(colnames(scdata()@colData)), selected = rev(names(scdata()@colData))[1])
      })
      # comparision vector
      group<-reactive({
        req(scdata())
        req(input$comparision_vector)
        g=colData(scdata())[,input$comparision_vector]
        g=unique(sort(g))
        g
      })

       #Pairwise comparision vector
     output$group_1<-renderUI({
       req(group())
       selectInput( ns("group1"),
                    label = "Select group 1 :",
                    choices = group(), selected = group()[1], selectize = FALSE)
     })

     output$group_2<-renderUI({
       req(group())
       selectInput( ns("group2"),
                    label = "Select group 2 :",
                    choices = group(),  selected = rev(group())[1] , selectize = FALSE)
     })






      output$group_by <- renderUI({
        selectInput(ns("group_By"),
                    label = "Group By:",
                    choices = c(colnames(scdata()@colData)), selected = rev(names(scdata()@colData))[1], multiple = T)
      })
      output$top_genes <- renderUI({
        numericInput(ns("Top_genes"), label = "# Top DEG to plot", 20, min = 1, max = Inf)
      })

      output$used <- renderUI({
        materialSwitch(
          inputId = ns("used_data"),
          label = "Normalized count",
          value = FALSE,
          status = "warning"
        )
      })

      # pairwise DEG analysis
      DEG_P<-reactive({
        req(scdata())
        validate(
          need(inherits(scdata(),"SingleCellExperiment"), "Please select a data set")
        )
        if(input$DE_Analysis== "pairwise"){
          deg<-DEG_Analyzer(scdata(),cmp = input$comparision_vector, method = input$method_DEGA, group1 = input$group1, group2 = input$group2)
        }
       else{
         deg<-DEG_Analyzer_c2c(scdata(),cmp = input$comparision_vector, method = input$method_DEGA)
       }


        return(deg)
      })

      #One Vs Rest
      #Select DEG table


      output$deg_list <-renderUI({
        if(input$DE_Analysis== "pairwise"){
          return()
        }
        else {
          selectInput(ns("deg_table"),
                      label = "Select DEG Table:",
                      choices = c(names(DEG_P())), selected = names(DEG_P())[1])
        }

      })

      DEG_select_table<- reactive({
        if(input$DE_Analysis== "pairwise"){
          deg<- DEG_P()
        }
        else {
          deg<- DEG_P()[[input$deg_table]] %>% select(-c(DEgene))
        }

        shinyalert(
          html = TRUE,
          type = 'success',
          text = tagList(
            div(
            tableOutput(ns('stat')),
            style="text-align: -webkit-center; font-size: 27px; font-family: Bradley Hand;"
          ))
        )

        return(deg)
      })

      output$DEG_table<-downloadHandler(
        filename = "Differentially_expressed_genes.txt",
        content = function(file) {
          write.table( DEG_select_table(), file = file, row.names = TRUE, sep = "\t")
        })

      ###Plot DEG



      volcano_p<- reactive({
        if(input$DE_Analysis== "pairwise"){
          p= volcano_plot(DEG_P(), fc =input$fold_change, pval = input$adj_p_value )
        }
        else{
          p= volcano_plot(DEG_select_table(), fc =input$fold_change, pval = input$adj_p_value )
        }

        vals$p=p
        p
      })

      heatmap_P<- reactive({

        #DEG Table
        if(input$DE_Analysis== "pairwise"){
          deg= DEG_P()
        }
        else{
          deg= DEG_select_table()
        }

        #Count matrix
        if(input$used_data){
          for(k in assayNames(scdata())){
            if(k=="BEcounts"){
              count= "BEcounts"
            }
            else if (k == "BENMcounts") {
              count= "BENMcounts"
            }
            else{
              count = "NMcounts"
            }
          }
          p=Plot_HeatmapDEG(scdata(),  used = count, degTable = deg, degTop = input$Top_genes, savePlot = FALSE, by=input$group_By, group1 = input$group1, group2 = input$group2)
        }
        else{
          p=Plot_HeatmapDEG(scdata(),  used = "counts", degTable = deg, degTop = input$Top_genes, savePlot = FALSE, by=input$group_By, group1 = input$group1, group2 = input$group2)
        }
        vals$p=p
        p
      })


      #plot
        output$pairwise_plot<-renderPlot({
          volcano_p()
        })
        output$heatmap<-renderPlot({
          heatmap_P()

        })


        #Download plot
        download_plot_Server("Volcano_Plt", input_data = volcano_p , name ="DEG_Volcano_plot")
        download_plot_Server("deg_heatmap", input_data = heatmap_P , name ="DEG_heatmap_plot")

      output$DEG_stats <- DT::renderDataTable({


          datatable(as.data.frame( DEG_select_table()))

      })




       output$stat <- renderTable({

         if(input$DE_Analysis== "pairwise"){
           p= DEG_table(DEG_P(), fc =input$fold_change, pval = input$adj_p_value )
         }
         else{
           p= DEG_table(DEG_select_table(), fc =input$fold_change, pval = input$adj_p_value )
         }

         p


       }, rownames = TRUE, bordered=TRUE, align='c')
       return(scdata)


    }# inputFun
  )#modl Srvr
}#Server
