################################################################################
######################### DIMENSION REDUCTION METHODS ##########################
################################################################################

#' A function to perform PCA dimension reduction
#'
#' @param sce a SingleCellExperiment object
#' @param runWith PCA dimension reduction Variable gene counts/ batch corrected variable gene counts/Normalized counts/ normalized
#'
#' @return  a SummarizedExperiment object with PCA dimension reduction
#' @examples
#'
#' DimReduction_pca(sce)
#'
DimReduction_pca = function(CS.data, runWith = 'VGcounts') {

  matrix.PJ = switch (EXPR = runWith,
                      'VGcounts'   = log2(SummarizedExperiment::assay(SingleCellExperiment::altExp(CS.data), 'VGcounts')+1),
                      'BEVGcounts' = SummarizedExperiment::assay(SingleCellExperiment::altExp(CS.data, 'BEVGcounts'), 'BEVGcounts'),
                      'NMcounts'   = SummarizedExperiment::assay(CS.data, 'NMcounts')[rownames(SummarizedExperiment::assay(SingleCellExperiment::altExp(CS.data), 'VGcounts')),],
                      'BEcounts'   = SummarizedExperiment::assay(CS.data, 'BEcounts')[rownames(SummarizedExperiment::assay(SingleCellExperiment::altExp(CS.data), 'VGcounts')),])
  emb.pca = prcomp(t(matrix.PJ))
  SingleCellExperiment::reducedDim(CS.data, 'PCA') = emb.pca$x[,c(1:50)]
  return(CS.data)
}

#' A function to perform ICA dimension reduction
#'
#' @param sce a SingleCellExperiment object
#' @param runWith ICA dimension reduction Variable gene counts/ batch corrected variable gene counts/Normalized counts/ normalized
#'
#' @return  a SummarizedExperiment object with ICA dimension reduction
#' @examples
#'
#' DimReduction_ica(sce)
#'
DimReduction_ica = function(CS.data, runWith = 'VGcounts') {

  matrix.PJ = switch (EXPR = runWith,
                      'VGcounts' = log2(SummarizedExperiment::assay(SingleCellExperiment::altExp(CS.data), 'VGcounts')+1),
                      'BEVGcounts' = SummarizedExperiment::assay(SingleCellExperiment::altExp(CS.data, 'BEVGcounts'), 'BEVGcounts'),
                      'NMcounts' = SummarizedExperiment::assay(CS.data, 'NMcounts')[rownames(SummarizedExperiment::assay(SingleCellExperiment::altExp(CS.data), 'VGcounts')),],
                      'BEcounts' = SummarizedExperiment::assay(CS.data, 'BEcounts')[rownames(SummarizedExperiment::assay(SingleCellExperiment::altExp(CS.data), 'VGcounts')),])

  emb.ica = fastICA::fastICA(X = matrix.PJ, n.comp = 50)#nrow(matrix.PJ))

  SingleCellExperiment::reducedDim(CS.data, 'PCA') = emb.ica$A

  return(CS.data)
}

#' A function to perform tSNE dimension reduction
#'
#' @param sce a SingleCellExperiment object
#' @param runWith use dimension reduction data PCA/ICA
#' @param perplexity the perplexity parameter determines how many nearest
#' neighbors are considered for each point when computing the pairwise similarities between points
#' @param PCNum number of dimension to consider
#' @return  a SummarizedExperiment object with tSNE dimension reduction
#' @examples
#'
#' DimReduction_tsne(sce)
#'
DimReduction_tsne = function(sce, runWith = 'PCA', PCNum = 20, perplexity = 30, parallelRun = TRUE) {

  matrix.PJ = switch (EXPR = runWith,
                      'VGcounts' = log2(assay(altExp(sce),'VGcounts') + 1),
                      'PCA' = t(reducedDim(sce, runWith)[,c(1:PCNum)]),
                      'BEPCA' = t(reducedDim(sce, runWith)[,c(1:PCNum)]))

  library(Rtsne, quietly = T)
  if (parallelRun) {
    set.seed(15555)
    emb.tsne = Rtsne(t(as.matrix(matrix.PJ)),
                     is_distance = FALSE,
                     perplexity = perplexity,
                     num_threads = parallel::detectCores(),
                     pca = F, check_duplicates = F,
                     verbose = FALSE)$Y
    rownames(emb.tsne) = colnames(matrix.PJ)
    colnames(emb.tsne) = c('tSNE_1', 'tSNE_2')
  } else {
    set.seed(15555)
    emb.tsne = Rtsne(t(as.matrix(matrix.PJ)),
                     is_distance = FALSE,
                     perplexity = perplexity,
                     num_threads = 1,
                     pca = F, check_duplicates = F,
                     verbose = FALSE)$Y
    rownames(emb.tsne) = colnames(matrix.PJ)
    colnames(emb.tsne) = c('tSNE_1', 'tSNE_2')
  }

  reducedDim(sce, 'tSNE') = emb.tsne

  return(sce)
}

#' A function to perform Umap dimension reduction
#'
#' @param sce a SingleCellExperiment object
#' @param runWith use dimension reduction data PCA/ICA
#'
#' @param PCNum number of dimension to consider
#' @return  a SummarizedExperiment object with Umap dimension reduction
#' @examples
#'
#' DimReduction_umap(sce)
#'
DimReduction_umap = function(sce, runWith = 'PCA', PCNum = 20) {

  matrix.PJ = switch (EXPR = runWith,
                      'VGcounts' = log2(assay(altExp(sce),'VGcounts') + 1),
                      'PCA' = t(reducedDim(sce, runWith)[,c(1:PCNum)]),
                      'BEPCA' = t(reducedDim(sce, runWith)[,c(1:PCNum)]))

  library(umap, quietly = T)
  set.seed(15555)
  emb.umap = umap(t(as.matrix(matrix.PJ)))

  emb.umap = emb.umap$layout
  colnames(emb.umap) = paste0('Umap_', c(1:ncol(emb.umap)))
  reducedDim(sce, 'Umap') = emb.umap

  return(sce)
}


#' A function to perform dfmap dimension reduction
#'
#' @param sce a SingleCellExperiment object
#' @param runWith use dimension reduction data PCA/ICA
#'
#' @param PCNum number of dimension to consider
#' @return  a SummarizedExperiment object with dfmap dimension reduction
#' @examples
#'
#' DimReduction_dfmap(sce)
#'

DimReduction_dfmap = function(sce, runWith = 'PCA', PCNum = 20) {

  matrix.PJ = switch (EXPR = runWith,
                      'VGcounts' = log2(assay(altExp(sce),'VGcounts') + 1),
                      'PCA' = t(reducedDim(sce, runWith)[,c(1:PCNum)]),
                      'BEPCA' = t(reducedDim(sce, runWith)[,c(1:PCNum)]))

  library(destiny, quietly = T)
  set.seed(15555)
  emb.dm = DiffusionMap(t(matrix.PJ))
  emb.dm = as.data.frame(emb.dm@eigenvectors)
  row.names(emb.dm) = colnames(matrix.PJ)

  reducedDim(sce, 'DFmap') = emb.dm

  return(sce)
}

# A function to perform PCA using scater
PCA_run<-function(sce,topGenes)
{
  data=scater::runPCA(sce,subset_row=topGenes)
  return(data)
}


################################################################################
########################## VISUALIZATION #######################################
################################################################################


#' A function to plot dimension reduction
#'
#' @param sce a SingleCellExperiment object
#' @param runWith dimension reduction PCA/tSNE/Umap
#' @param colorBy select a column in colData
#' @param showDensity T/F to show density on plot
#' @return  a ggplot on dimension reduction
#' @examples
#'
#' Plot_DimReduction(sce)
#'
Plot_DimReduction = function(sce, runWith = 'tSNE', colorBy = 'default', dim.se = c(1:2), showDensity = TRUE) {
  if (colorBy == 'default') {
    colorBy = colnames(colData(sce))[1]
  }

  #sce = DimReduction_pca()

  #plot
  df.plot = data.frame(reducedDim(sce, runWith)[,dim.se],
                       colorBy = as.factor(colData(sce)[,colorBy]),
                       color = as.character(colData(sce)[,colorBy]),

                       stringsAsFactors = F)

  library(colorspace, quietly = T)
  library(RColorBrewer, quietly = T)

  for(i in 1:nrow(df.plot)) {
    df.plot$color[i] = qualitative_hcl(length(unique(df.plot$colorBy)))[as.vector(as.integer(df.plot$colorBy))[i]]
  }

  library(ggplot2, quietly = T)
  g = ggplot(data = df.plot, mapping = aes(x = df.plot[,1], y = df.plot[,2], color = colorBy)) +
    geom_point() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    theme(legend.position="bottom",
          legend.title = element_blank()) +
    xlab(colnames(df.plot)[1]) +
    ylab(colnames(df.plot)[2])

  g.density.x = ggplot(data = df.plot, mapping = aes(x = df.plot[,1], y = ..density.., color = colorBy)) +
    geom_density() +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black")) +
    theme(#axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()) +
    theme(legend.position = "none") +
    xlab(colnames(df.plot)[1]) +
    ylab('Density') #+

  g.density.y = ggplot(data = df.plot, mapping = aes(x = ..density.., y = df.plot[,2], color = colorBy)) +
    geom_density() +
    scale_x_reverse() +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black")) +
    theme(#axis.title.y=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()) +
    theme(legend.position = "none") +
    ylab(colnames(df.plot)[2]) +
    xlab('Density')#

  g.legend = ggplot(data = df.plot, mapping = aes(x = 0, y = 0, color = colorBy)) +
    geom_point() +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black")) +
    theme(legend.position="bottom",
          legend.title = element_blank())


  if(!showDensity) {
    g.return = g
  } else {
    figures_without_legend = cowplot::plot_grid(
      g.density.y,
      g + theme(legend.position = "none"),
      NULL,
      g.density.x,

      align = "hv",
      nrow = 2, rel_widths = c(0.5, 1), rel_heights = c(1, 0.5))

    legends = cowplot::get_legend(g.legend)

    g.all = cowplot::plot_grid(figures_without_legend,
                               legends,
                               rel_heights = c(10,2), nrow = 2)

    g.return = g.all
  }

  return(g.return)
}

#' A function Visualize Dimensional Reduction genes
#'
#'
#'
#' @param sce a SingleCellExperiment object as input
#' @param ndims number of PCs to visualize
#' @inputs sce a SingleCellExperiment object as input
#' @return a ggplot
#' @examples
#'
#' dim_lodingsViz(sce, ndims = 2)
#'
#Dimension loadings plot (PCA)
dim_lodingsViz<-function(sce,ndims, ngenes){

  #condition to select ngenes to plot
  if(ndims<=6) ngenes = 20
  if(ndims >6 & ndims <=12) ngenes = 10
  if(ndims >12 & ndims <=24) ngenes =5
  if(ndims >24) ngenes = 2

  plots <- lapply(
    1:ndims,
    FUN=function(i) {
      #check and select PCA name
      for(j in reducedDimNames(sce)){
        if(j=="BEPCA"){
          rdim_name= "BEPCA"
        }
        else if (j == "PCA") {
          rdim_name= "PCA"
        }
      }
      loading<-reducedDim(sce,rdim_name)
      loading<-attr(loading,"rotation")
      loading<-as.data.frame(loading)


      data.plot=head(loading[order(abs(loading[,i]),decreasing=TRUE),],ngenes)[i]
      data.plot$feature <- factor(x = rownames(x = data.plot), levels = rownames(x = data.plot))
      data.plot$feature<-with(data.plot,reorder(data.plot[,2],data.plot[,1]))
      p=ggplot(
        data = data.plot,
        mapping = aes_string(x = colnames(x = data.plot)[1], y = 'feature')) +

        geom_point(col = "blue") +
        labs(y = NULL) + theme_cowplot()
      return(p)
    })
  #condition to devide plot space
  if(ndims==1) ncol=1
  if(ndims == 2) ncol = 2
  if(ndims >= 3) ncol = 3
  plots <- wrap_plots(plots, ncol = ncol)
  return(plots)
}

#' A function Visualize Dimensional reduction heatmap
#'
#' @param sce a SingleCellExperiment object as input
#' @param ndims number of PCs to visualize
#' @param nfeatures number of genes to plot
#' @inputs sce a SingleCellExperiment object as input
#' @return a a ggplot
#' @examples
#'
#' dim_heatmap(sce, ndims = 2, nfeatures =10)
#'
#Dim Heatmap
dim_heatmap=function(sce,ndims,nfeatures){
  plots <- lapply(
    1:ndims,
    FUN=function(i) {
      for(j in SingleCellExperiment::reducedDimNames(sce)){
        if(j=="BEPCA"){
          rdim_name= "BEPCA"
        }
        else if (j == "PCA") {
          rdim_name= "PCA"
        }
      }

      loading<-reducedDim(sce,rdim_name)
      loading<-attr(loading,"rotation")
      loading<-as.data.frame(loading)
      data.plot=head(loading[order(abs(loading[,i]),decreasing=TRUE),],nfeatures)[i]
      for(k in assayNames(sce)){
        if(k=="BEcounts"){
          count= "BEcounts"
        }
        else if (k == "BENMcounts") {
          count= "BENMcounts"
        }
        else if(k == 'NMcounts'){
          count = "NMcounts"
        }
        else if(k == 'logcounts'){
          count = "logcounts"
        }
        else{
          count = 'counts'
        }
      }
      edata<-assay(sce, count)
      ss=edata[row.names(edata) %in% rownames(data.plot),]
      df_num_scale = scale(ss)
      df_num_scale[is.na(df_num_scale)] = 0.1
      p=ggplotify::as.ggplot(pheatmap::pheatmap(df_num_scale[,1:100],main = paste0("PC",i),treeheight_row=0,treeheight_col=0,show_colnames=F))
      p
    }
  )
  #condition to devide plot space
  if(ndims==1) ncol=1
  if(ndims == 2) ncol = 2
  if(ndims >= 3) ncol = 3
  plots <- patchwork::wrap_plots(plots, ncol = ncol)
  return(plots)
}


#' A function Visualize Dimensions to pick the relevent Dimensions
#'
#' @param sce a SingleCellExperiment object as input
#' @param ndims number of PCs to visualize
#' @param reduction PCA
#' @inputs sce a SingleCellExperiment object as input
#' @return a a ggplot
#' @examples
#'
#' ElbowPlot(sce, ndims = 2, reduction = 'PCA')
#'
ElbowPlot <- function(sce, ndims = 20, reduction = 'PCA') {
  for(j in SingleCellExperiment::reducedDimNames(sce)){
    if(j=="BEPCA"){
      rdim_name= "BEPCA"
    }
    else if (j == "PCA") {
      rdim_name= "PCA"
    }
  }

  data.use<-reducedDim(sce,rdim_name)
  data.use<-attr(data.use,"percentVar")

  #data.use=reducedDim(sce, 'PCA')
  if (length(x = data.use) == 0) {
    stop(paste("No standard deviation info stored for", reduction))
  }
  if (ndims > length(x = data.use)) {
    warning("The object only has information for ", length(x = data.use), " reductions")
    ndims <- length(x = data.use)
  }
  stdev <- 'Standard Deviation'
  plot <- ggplot(data = data.frame(dims = 1:ndims, stdev = data.use[1:ndims])) +
    geom_point(mapping = aes_string(x = 'dims', y = 'stdev')) +
    labs(
      x="Principal Component #", y="Percent Variance"
    ) +
    theme_cowplot()
  return(plot)
}


#' A function to search the optimum dimension to select for the further analysis
#'
#' @param sce a SingleCellExperiment object
#' @param runWith count matrix variable genes count/ raw count matrix
#'
#'
#' @return  a ggplot
#' @examples
#'
#' PCNumSearch(sce)
#'
PCNumSearch = function(sce, runWith = 'VGcounts') {
  if (runWith == 'VGcounts') {
    matrix.PJ = assay(altExp(sce), runWith)
    matrix.PJ = log2(matrix.PJ + 1)
  }
  else{
    matrix.PJ = assay(sce, 'counts')
  }


  suppressPackageStartupMessages(library(PCAtools, quietly = T))
  p = PCAtools::pca(matrix.PJ)
  #PCAtools::screeplot(p, axisLabSize = 18, titleLabSize = 22)
  #PCAtools::biplot(p)
  #PCAtools::plotloadings(p, labSize = 3)
  elbow = PCAtools::findElbowPoint(p$variance)

  horn = PCAtools::parallelPCA(matrix.PJ)

  ymax = max(elbow, horn$n)
  df = data.frame(pc = factor(PCAtools::getComponents(p, 1: (ymax+5)), levels = c(PCAtools::getComponents(p, 1: (ymax+5)))),
                  ev = p$variance[PCAtools::getComponents(p, 1: (ymax+5))])

  g = ggplot2::ggplot(df) +
    ggplot2::geom_col(aes(x = pc, y = ev), color = "dodgerblue", fill = "white") +
    ggplot2::geom_line(aes(x = pc, y = ev), color = "red", group = 1, alpha = 0.5) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          #panel.border = element_rect(colour = "black", size=1),
          panel.background = element_blank()) +
    ggplot2::xlab('PCs') +
    ggplot2::ylab('Explained variation') +
    ggplot2::theme(legend.position="bottom",
                   legend.title = element_blank()) +
    ggplot2::theme(axis.text.y = element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.text.x = element_text(angle = 90)) +

    ggplot2::geom_vline(xintercept = elbow, linetype = "longdash",
                        colour = "black") +
    ggplot2::geom_label(x = elbow + 1, y = max(df$ev)/2,
                        label = paste0('Elbow'), vjust = 0, size = 4) +

    ggplot2::geom_vline(xintercept = horn$n, linetype = "longdash",
                        colour = "black") +
    ggplot2::geom_label(x = horn$n + 1, y = (max(df$ev) + 4)/2,
                        label = paste0('Horn'), vjust = 0, size = 4)
  #g
  return(g)
}




