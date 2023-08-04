#require(shinydashboard)


packages=c('shiny', 'devtools','tidyverse', 'rrvgo','gprofiler2', 'org.Mm.eg.db', 'enrichplot','clusterProfiler','readxl','openxlsx','SingleR','celldex','HGNChelper','DEsingle','bluster','scran','scater',
           'singleCellTK','ggfortify','cowplot','ggrepel','bslib','scuttle','sass','patchwork','Seurat','purrr','DropletUtils','spsComps','SingleCellExperiment','waiter','ggplotify','grid','viridis',
           'pheatmap','shinyBS','shinyjs','DT','shinyFiles','shinycssloaders','shinyWidgets','thematic','bs4Dash','fresh', 'shiny.blueprint','CytoTRACE')
# Use lapply to load each package
suppressPackageStartupMessages(lapply(packages, require, character.only = TRUE))

#xfun::pkg_attach(all_packages)
# for NIDCR shiny server
# reticulate::use_python("/home/shiny/miniconda3/bin/python", required= TRUE)

options(shiny.autoreload = TRUE)
#
options(shiny.maxRequestSize=10000*1024^2)
#Import Functions
#source("fun/cs.R",encoding = "utf-8")
#source("fun/visualization.R",encoding = "utf-8")

source("fun/Preprocessing_Fun.R",encoding = "utf-8")
source("fun/Normalization_Fun.R",encoding = "utf-8")
source("fun/Dimensional_Reduction_Fun.R",encoding = "utf-8")
source("fun/Clustering_Fun.R",encoding = "utf-8")
source("fun/Batch_correction_Fun.R",encoding = "utf-8")
source("fun/Differential_Gene_Expression_Analysis_Fun.R",encoding = "utf-8")
source("fun/Cell_type_annotation_Fun.R",encoding = "utf-8")
source("fun/Cell_development_Fun.R",encoding = "utf-8")
source("fun/Other_fun.R",encoding = "utf-8")


#import modules

source("module/sidebar.R",encoding = "utf-8")
source("module/data_input_new.R",encoding = "utf-8")
source("module/Pre_QC.R",encoding = "utf-8")
source("module/Normalization.R",encoding = "utf-8")
source("module/PCA_Dimension_Reduction.R",encoding = "utf-8")
source("module/marker_identification_DEG.R",encoding = "utf-8")
source("module/Cell_Type.R",encoding = "utf-8")
source("module/Batch_correction.R",encoding = "utf-8")
source("module/batch_non_linear.R",encoding = "utf-8")
source("module/variable_gene_batch.R",encoding = "utf-8")
source("module/batch_linear.R",encoding = "utf-8")
source("module/batch_correction_evaluation.R",encoding = "utf-8")
source("module/Cell_developement.R",encoding = "utf-8")
source("module/plot_download.R", encoding = "utf-8")
source("module/plot_markers.R", encoding = "utf-8")
source("module/Introduction.R")



ui <- bs4DashPage( preloader = list(html = tagList( spin_google(), "Please Wait Loading ..."), color = 	"#A020F0"),


                   fullscreen = TRUE,
                   scrollToTop = TRUE,


                header =dashboardHeader(
                  titleWidth = "thik",
                title = dashboardBrand(

                  title = strong("Single cell analysis tool"),
                  color = "purple",
                  opacity=1

                ),

                border = TRUE,
                sidebarIcon = icon("bars"),
                controlbarIcon = icon("th"),
                fixed = FALSE

),
sidebar=bs4DashSidebar(
  width = 350,
  skin = "dark",
  status = "purple",
  elevation = 6,
  sidebarUserPanel(
    image = "hex-ScTat.png",
    name = "Welcome to sceExplorer"
  ),
  sidebarUI("sidebar")
   ),


  body=bs4DashBody(
	useShinyjs(),

	# use_theme(create_theme(
	#   theme = "default",
	#   bs_vars_wells(
	#     bg = "#CEECF5",
	#     border = "#3f2d54"
	#   )
	# )),
  #include styling
	tags$head( includeCSS("www/dark_mode.css")),
	#tags$style("@import url(https://use.fontawesome.com/releases/v6.3.0/css/all.css);"),
	#tags$style(includeCSS("~/Desktop/Tapan/fontawesome/fontawesome-free-6.4.0-web/css/all.css")),
	tags$style(includeCSS("www/fontawesome/v6.4.0/css/all.css")),


	bs4TabItems(
	  bs4TabItem(tabName="intro",Intro_UI("intro")),
	  bs4TabItem(tabName="dataInput",datainputUI("dataInput")),
	  bs4TabItem(tabName="PreQC",preQCUI("preQC")),
	  bs4TabItem(tabName="normalization",normalizationUI("normalization")),
	  bs4TabItem(tabName="batch_correction",batchCorrct_UI("batchcorrection")),
	  bs4TabItem(tabName="batch_evaluation",batchCorrct_evalutn_UI("batchCorrct_evalutn")),
	  #bs4TabItem(tabName="variable_gene_batch",batch_variable_gene_UI("batch_VG")),
	  bs4TabItem(tabName="linear_batch", batch_dim_reduction_UI("linear_batch")),
	  bs4TabItem(tabName="non_linear_batch",batch_non_linear_UI("batch_nonLinear")),
	  bs4TabItem(tabName="dim_reduction",dim_reduction_UI("dimreduction")),
	  bs4TabItem(tabName="DGE_Analysis",de_analysis_UI("de_analysis")),
	  bs4TabItem(tabName="plot_markers",Plot_Markers_UI("marker_plot")),
	  bs4TabItem(tabName="Cell_type",cell_type_analysis_UI("cell_type")),
	  bs4TabItem(tabName="Cell_develope",cell_developement_analysis_UI("Cell_develope"))


      )

    ),

footer=dashboardFooter(
    wellPanel(
    HTML(
      '
       <p align="center" width="4">Bioinformatics Facility, Karolinska Institute, Ming Wai Lau Centre for Reparative Medicine
, Hong Kong</p>

       <p align="center" width="4">Using SingleCellExperiment object </p>'
  )
  )
  )

)

server <- function(input, output, session) {

  observeEvent(input$dark_mode, {
    toast(
      title = if (input$dark_mode) "Dark theme on!" else "Light theme on",
      options = list(position = "topRight", class = "bg-warning", autohide = TRUE)
    )
  })
  useAutoColor()

  waiter_hide()
  Intro_Server('intro')
  sidebarServer("sidebar")
  raw_data <-datainputServer("dataInput")
  filtered_data<-preQCServer("preQC", raw_data)

  normalization_data<-normalizationServer("normalization",filtered_data)
  dim_reduction_data<-dim_reduction_Server("dimreduction",normalization_data)

  batch_corrct<-batchCorrct_Server("batchcorrection",normalization_data)
  batch_corrct_dmrd<-batch_non_linear_Server("batch_nonLinear",batch_corrct)
  #batch_crrct_VG<-batch_variable_gene_Server("batch_VG",batch_corrct )
  batch_corrct_exprsn<-batch_dim_reduction_Server("linear_batch", batch_corrct)
  batchCorrct_evalutn_Server("batchCorrct_evalutn", batch_corrct)

  dim_reduction_batch_correction<-rm_dim_data("batchcorrection", batch_corrct_dmrd, batch_corrct_exprsn)
 # Dim_reduction_data <- dim_rdxn_data("dim_rdxn_data")


  DGE<-de_analysis_Server("de_analysis", dim_reduction_data , dim_reduction_batch_correction )
  Plot_Markers_Server("marker_plot", DGE)
  cell_type<-cell_type_analysis_Server("cell_type",DGE)

  cell_developement_analysis_Server("Cell_develope", cell_type)

}

shinyApp(ui, server)
