sidebarUI <- function(id) {
  ns <- NS(id)
  bs4SidebarMenu(
  id="tabselected",
  bs4SidebarMenuItem("Data Input", tabName = "dataInput", icon = icon("database",lib = "font-awesome",class="fa-solid fa-database")),
  bs4SidebarMenuItem("QC", tabName = "PreQC", icon = icon("filter",lib = "font-awesome",class="fa-solid fa-filter")),
 
  bs4SidebarMenuItem("Normalization", tabName = "normalization", icon = icon("scale-balanced",lib = "font-awesome",class="fa-duotone fa-scale-balanced")),
  bs4SidebarMenuItem("Dimension Reduction", tabName = "dim_reduction", icon = icon("minimize",lib = "font-awesome",class="fa-duotone fa-minimize")),
  bs4SidebarMenuItem("Batch Correction", tabName = "batch_correction", icon = icon("check-double",lib = "font-awesome",class="fa-duotone fa-check-double")),
  bs4SidebarMenuItem("Batch correction evaluation",tabName = "batch_evaluation", icon = icon("check-double",lib = "font-awesome",class="fa-solid fa-check-double")),
  #bs4SidebarMenuItem("Variable Features Batch",tabName = "variable_gene_batch", icon = icon("check-double",lib = "font-awesome",class="fa-duotone fa-check-double")),
  bs4SidebarMenuItem("Dimension Reduction Batch",tabName = "linear_batch", icon = icon("check-double",lib = "font-awesome",class="fa-duotone fa-check-double")),
  bs4SidebarMenuItem("Non-Linear Dim Reduction Batch",tabName = "non_linear_batch", icon = icon("check-double",lib = "font-awesome",class="fa-duotone fa-check-double")),
  bs4SidebarMenuItem("DEG Analysis", tabName = "DGE_Analysis", icon = icon("retweet",lib = "font-awesome",class="fa-duotone fa-retweet")),
  bs4SidebarMenuItem("Plot Markers", tabName = "plot_markers", icon = icon("retweet",lib = "font-awesome",class="fa-duotone fa-retweet")),
  bs4SidebarMenuItem("Cell Type", tabName = "Cell_type", icon = icon("circle-nodes",lib = "font-awesome",class="fa-duotone fa-circle-nodes")),
  bs4SidebarMenuItem("Cell Developement", tabName = "Cell_develope", icon = icon("circle-nodes",lib = "font-awesome",class="fa-duotone fa-circle-nodes")) 
                     
                     
                   
  
  
)
}

sidebarServer <- function(input, output, session) {
}
