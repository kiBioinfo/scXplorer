
################################################################################
##################### METHODS FOR TRAJECTORY ANALYSIS ##########################
################################################################################


### Find pseudo trajectory ----

#' A function to Find pseudo trajectory
#'
#' @param sce a SingleCellExperiment object with cell development trajectory analysis info
#' @param runWith run with dimension reduction Umap/tSNE/PCA
#' @param batch batch T/F
#' @param method a method to Find pseudo trajectory
#' @return  a plot shows cell development in Pseudo Trajectory
#'
#' Plot_Development(sce, runWith = 'Umap', colorBy = 'cellType' )
#'

CellDevelopment <- function(sce,  method='Princurve_fit', runWith = 'PCA', dimSe = c(1,2),
                            used = 'counts', batch = F , clusterSe= 'cluster', seedBy){

  sce <- switch(EXPR = method,
                "Princurve_fit" = CellDev_Princurve_fit(sce, runWith = runWith, dimSe = dimSe),
                "CytoTrace" = CellDev_CytoTrace(sce, used = used, batch = batch),
                "Slingshot" = CellDev_Slingshot(sce, clusterSe= clusterSe, runWith = runWith),
                "Monocle" = CellDev_Monocle(sce, seedBy))
  return(sce)
}



CellDev_Princurve_fit = function(sce, runWith = 'PCA', dimSe = c(1,2)) {

  fit1 = princurve::principal_curve(SingleCellExperiment::reducedDim(sce, runWith)[,dimSe], plot = F, smoother = "smooth_spline")
  cell.order = fit1$lambda
  SummarizedExperiment::colData(sce)$cellOd = cell.order
  SingleCellExperiment::reducedDim(sce, 'PrinFit') = fit1$s
  return(sce)
}

CellDev_CytoTrace = function(sce, used = 'counts', batch = F) {
  edata = SummarizedExperiment::assay(sce, used)
  #suppressPackageStartupMessages(library(CytoTRACE, quietly = T))
  if(!batch) {
    results2 = CytoTRACE::CytoTRACE(mat = log2(edata + 1))
    SummarizedExperiment::colData(sce)$cellOd = results2$CytoTRACE
  } else {
    batch = SingleCellExperiment::colData(sce)[,'batch']
    results2 = CytoTRACE::CytoTRACE(mat = log2(edata + 1), batch = batch)
    SummarizedExperiment::colData(sce)$cellOd = results2$CytoTRACE
  }
  return(sce)
}

CellDev_Slingshot = function(sce, clusterSe= 'label', runWith = 'PCA') {
  #clusterSe = 'ClusterGraph'
  sim = slingshot::slingshot(sce, clusterLabels = clusterSe, reducedDim = runWith)
  SummarizedExperiment::colData(sce)$cellOd = sim@colData$slingPseudotime_1
  #plot(SingleCellExperiment::reducedDim(sce, runWith))
  #lines(slingshot::SlingshotDataSet(sim), lwd=2, col='black')
  return(sce)

}

#seedBy = Cell Type or cluster:(Single or multiple)
#E.g. cluster: 1:5 or 6
CellDev_Monocle = function(sce, seedBy) {

  cds = monocle3::new_cell_data_set(expression_data = SummarizedExperiment::assay(sce, 'counts'),
                                    cell_metadata = SingleCellExperiment::colData(sce),
                                    gene_metadata = data.frame(gene_short_name = rownames(SummarizedExperiment::assay(sce, 'counts')),
                                                               row.names = rownames(SummarizedExperiment::assay(sce, 'counts'))))
  col.name = colnames(SingleCellExperiment::colData(sce))[apply(SingleCellExperiment::colData(sce), 2, function(col) any(col == seedBy))]
  cds = monocle3::preprocess_cds(cds, num_dim = 50)
  if('batch' %in% colnames(SingleCellExperiment::colData(sce))) { cds = monocle3::align_cds(cds, alignment_group = "batch") }
  cds = monocle3::reduce_dimension(cds)
  cds = monocle3::cluster_cells(cds)
  cds = monocle3::learn_graph(cds)
  cell_ids = which(SingleCellExperiment::colData(cds)[, col.name] == seedBy)
  closest_vertex = cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex = as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes = igraph::V(monocle3::principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  cds = monocle3::order_cells(cds, root_pr_nodes=root_pr_nodes)
  SummarizedExperiment::colData(sce)$cellOd = monocle3::pseudotime(cds)[colnames(SummarizedExperiment::assay(sce, 'counts'))]
  return(sce)
}


################################################################################
########################## VISUALIZATION #######################################
################################################################################


#' A function to plot cell cycle development IN Pseudo time
#'
#' @param sce a SingleCellExperiment object with cell development trajectory analysis info
#' @param runWith run with dimension reduction Umap/tSNE
#' @param colorBy name of the column to use for color
#'
#' @return  a plot shows cell development in Pseudo Trajectory
#'
#' Plot_Development(sce, runWith = 'Umap', colorBy = 'cellType' )
#'
Plot_Development = function(sce, runWith = 'Umap', dimSe = c(1,2), colorBy = 'default', rev = F, devWindowNum = 5) {

  if (colorBy == 'default') {
    colorBy = colnames(colData(sce))[1]
  }

  cell.order = sce$cellOd
  if(rev) {
    cell.order = max(cell.order) - cell.order
  }

  df.plot = data.frame(reducedDim(sce, runWith)[,dimSe],
                       colorBy = as.factor(colData(sce)[,colorBy]),
                       color = as.character(colData(sce)[,colorBy]),
                       cellOd = cell.order,
                       stringsAsFactors = F)

  df.plot.od = df.plot[order(df.plot$cellOd, decreasing = F),]
  df.plot.od.tmp = data.frame()

  #### method 2
  x.ave = .myAverageWindow(VectorX = df.plot.od[,1], windowNum = round(nrow(df.plot.od)/devWindowNum))
  y.ave = .myAverageWindow(VectorX = df.plot.od[,2], windowNum = round(nrow(df.plot.od)/devWindowNum))
  for(i in c(1:nrow(df.plot.od)-1)) {
    #i = 1
    df.plot.od.tmp[i,1] = 50*(x.ave[i+1] - df.plot.od[i, 1])/(nrow(df.plot.od) - 1) + df.plot.od[i, 1]
    df.plot.od.tmp[i,2] = 50*(y.ave[i+1] - df.plot.od[i, 2])/(nrow(df.plot.od) - 1) + df.plot.od[i, 2]
  }
  df.plot.od.tmp[nrow(df.plot.od),] = df.plot.od[nrow(df.plot.od),c(1,2)]
  #### method 2
  df.plot.arrow = data.frame(x = df.plot.od[,1],
                             xend = df.plot.od.tmp[,1],
                             y = df.plot.od[,2],
                             yend = df.plot.od.tmp[,2])
  library(ggplot2, quietly = T)
  g.arrow = ggplot() +
    geom_point(data = df.plot.od, mapping = aes(x = df.plot.od[,1], y = df.plot.od[,2], color = colorBy)) +
    geom_segment(data = df.plot.arrow, aes(x = x,
                                           xend = xend,
                                           y = y,
                                           yend = yend),
                 arrow = arrow(length = unit(0.1,"cm"))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black")) +
    theme(legend.position="bottom",
          legend.title = element_blank()) +
    xlab(colnames(df.plot.od)[1]) +
    ylab(colnames(df.plot.od)[2])

  library(colorspace, quietly = T)
  library(RColorBrewer, quietly = T)

  for(i in 1:nrow(df.plot)) {
    df.plot$color[i] = qualitative_hcl(length(unique(df.plot$colorBy)))[as.vector(as.integer(df.plot$colorBy))[i]]
  }

  library(ggplot2, quietly = T)
  g.density.x = ggplot(data = df.plot, mapping = aes(x = cellOd, y = after_stat(density), color = colorBy)) +
    geom_density() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    theme(legend.position = "none") +
    xlab('Cell order') +
    ylab('Density')
  #g.density.x

  figures_without_legend = cowplot::plot_grid(
    g.arrow + theme(legend.position = "none"),
    g.density.x,

    align = "hv",
    nrow = 2,
    #rel_widths = c(0.5, 1),
    rel_heights = c(1, 0.5)
  )

  legends = cowplot::get_legend(g.arrow)

  g.all = cowplot::plot_grid(figures_without_legend,
                             legends,
                             rel_heights = c(10,2), ncol = 1)
  return(g.all)
}

Plot_Development_Princurve = function(sce, runWith = 'Umap', dimSe = c(1,2), colorBy = 'default', rev = F) {
  if (colorBy == 'default') {
    colorBy = colnames(colData(sce))[1]
  }

  cell.order = sce$cellOd
  if(rev) {
    cell.order = max(cell.order) - cell.order
  }

  suppressPackageStartupMessages(library(ggplot2, quietly = T))
  df.plot = data.frame(reducedDim(sce, runWith)[,dimSe],
                       fit_x = reducedDim(sce, 'PrinFit')[,1],
                       fit_y = reducedDim(sce, 'PrinFit')[,2],
                       colorBy = as.factor(colData(sce)[,colorBy]),
                       cellOd = cell.order,
                       stringsAsFactors = F)
  df.plot.od = df.plot[order(df.plot$cellOd, decreasing = T),]
  g = ggplot() +
    geom_point(data = df.plot.od, mapping = aes(x = df.plot.od[,1], y = df.plot.od[,2], color = colorBy)) +
    geom_path(data = df.plot.od, mapping = aes(x = fit_x, y = fit_y)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    theme(legend.position="bottom",
          legend.title = element_blank()) +
    xlab(colnames(df.plot.od)[1]) +
    ylab(colnames(df.plot.od)[2])

  g.density.x = ggplot(data = df.plot, mapping = aes(x = cellOd, y = ..density.., color = colorBy)) +
    geom_density() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    theme(legend.position = "none") +
    xlab('Cell order') +
    ylab('Density')

  figures_without_legend = cowplot::plot_grid(
    g + theme(legend.position = "none"),
    g.density.x,

    align = "hv",
    nrow = 2,
    rel_heights = c(1, 0.5)
  )

  legends = cowplot::get_legend(g)

  g.all = cowplot::plot_grid(figures_without_legend,
                             legends,
                             rel_heights = c(10,2), ncol = 1)
  return(g.all)
}
