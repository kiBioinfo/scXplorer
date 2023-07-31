################################################################################
######################## CLUSTERING METHODS ####################################
################################################################################


#' A function to find clusters in the dataset
#'
#' @param sce a SingleCellExperiment object
#' @param runWith run with VGcounts/PCA/ICA/tSNE/Umap
#' @param method a method to predict number of clusters in the dataset
#' @param ClusterNum user defind number of clusters
#' @param PCNum number of PCs to consider only in PCA/ICA condition
#'
#' @return  a SummarizedExperiment object with clusters
#' @examples
#'
#' ClusterFind(sce)
#'

ClusterFind = function(sce, method = 'Hclust', ClusterNum = 3, runWith = 'PCA', PCNum = 10, k=100,cluster.fun='walktrap') {

  sce = switch(EXPR = method,
               'Hclust'	 = ClusterFind_hclust(   sce = sce,  ClusterNum = ClusterNum, runWith = runWith, PCNum = PCNum),
               'K-Means'     = ClusterFind_kmean(  sce = sce,  ClusterNum = ClusterNum, runWith = runWith, PCNum = PCNum),
               'Mclust'	 = ClusterFind_mclust(    sce = sce,  ClusterNum = ClusterNum, runWith = runWith, PCNum = PCNum),
               'Graph'	 = ClusterFind_graph( sce = sce,   runWith = runWith, PCNum = PCNum, cluster.fun= cluster.fun),
               'DBscan'      = ClusterFind_dbscan(sce=sce, runWith = runWith, PCNum = PCNum),
               'Hclust_scran'	= ClusterFind_hclust_scran(sce=sce, runWith = runWith),
               'Kmean_scran' = ClusterFind_kmean_scran(sce = sce,  k = k, runWith = runWith)
  )


  return(sce)
}

# Clustering with k-mean
ClusterFind_kmean = function(sce, ClusterNum, runWith = 'PCA', PCNum = 10){

  matrix.run = switch(EXPR       = runWith,
                      'VGcounts' = assay(altExp(sce),'VGcounts'),
                      'PCA'      = reducedDim(x = sce, type = 'PCA')[,1:PCNum],
                      'ICA'      = reducedDim(x = sce, type = 'ICA')[,1:PCNum],
                      'tSNE'     = reducedDim(x = sce, type = 'tSNE')[,1:2],
                      'Umap'     = reducedDim(x = sce, type = 'Umap')[,1:2])

  set.seed(0)
  res.km = kmeans(x = matrix.run, centers = ClusterNum)
  res.km.cluster = res.km$cluster
  #plot(c(1:nrow(matrix.run)), match(rownames(matrix.run), names(res.km.cluster)))

  colData(sce)$cluster = as.character(res.km.cluster)
  #colLabels(sce)<-as.factor(res.km.cluster)
  return(sce)
}

# Clustering with dbscan
ClusterFind_dbscan = function(sce, runWith = 'PCA', PCNum = 10){

  matrix.run = switch(EXPR       = runWith,
                      'VGcounts' = assay(altExp(sce),'VGcounts'),
                      'PCA'      = reducedDim(x = sce, type = 'PCA')[,1:PCNum],
                      'ICA'      = reducedDim(x = sce, type = 'ICA')[,1:PCNum],
                      'tSNE'     = reducedDim(x = sce, type = 'tSNE')[,1:2],
                      'Umap'     = reducedDim(x = sce, type = 'Umap')[,1:2])

  #library(dbscan, quietly = T)
  set.seed(0)
  res.db = dbscan::dbscan(x = matrix.run, eps = 0.4, minPts = 5)
  res.db.cluster = res.db$cluster

  colData(sce)$cluster = as.character(res.db.cluster)
  #colLabels(sce)<- as.factor(res.db.cluster)
  return(sce)
}

# Clustering with mclust

ClusterFind_mclust = function(sce, ClusterNum = 3, runWith = 'PCA', PCNum = 10){

  matrix.run = switch(EXPR       = runWith,
                      'VGcounts' = assay(altExp(sce),'VGcounts'),
                      'PCA'      = reducedDim(x = sce, type = 'PCA')[,1:PCNum],
                      'ICA'      = reducedDim(x = sce, type = 'ICA')[,1:PCNum],
                      'tSNE'     = reducedDim(x = sce, type = 'tSNE')[,1:2],
                      'Umap'     = reducedDim(x = sce, type = 'Umap')[,1:2])

  suppressPackageStartupMessages(library(mclust, quietly = T))
  res.mc = mclust::Mclust(data = matrix.run, G =  c(1:ClusterNum))
  res.mc.cluster = res.mc$classification
  #plot(c(1:nrow(matrix.run)), match(rownames(matrix.run), names(res.mc.cluster)))

  colData(sce)$cluster = as.character(res.mc.cluster)
  #colLabels(sce) <- as.factor(res.mc.cluster)
  return(sce)
}

# Clustering with hierarchical

ClusterFind_hclust = function(sce, ClusterNum = 3, runWith = 'PCA', PCNum = 10) {

  matrix.run = switch(EXPR       = runWith,
                      'VGcounts' = assay(altExp(sce),'VGcounts'),
                      'PCA'      = reducedDim(x = sce, type = 'PCA')[,1:PCNum],
                      'ICA'      = reducedDim(x = sce, type = 'ICA')[,1:PCNum],
                      'tSNE'     = reducedDim(x = sce, type = 'tSNE')[,1:2],
                      'Umap'     = reducedDim(x = sce, type = 'Umap')[,1:2])

  suppressPackageStartupMessages(library(dplyr, quietly = T))
  res.hc = matrix.run %>%
    scale() %>%                    # Scale the data
    dist(method = "euclidean") %>% # Compute dissimilarity matrix
    hclust(method = "ward.D2")     # Compute hierachical clustering

  res.hc.cluster = cutree(tree = res.hc, k = ClusterNum)
  #plot(c(1:nrow(matrix.run)), match(rownames(matrix.run), names(res.hc.cluster)))

  colData(sce)$cluster = as.character(res.hc.cluster)
  #colLabels(sce)<- as.factor(res.hc.cluster)
  return(sce)
}

# Graph basecd clustering
# cluster.fun = one of walktrap/louvain/infomap/fast_greedy/label_prop/leading_eigen
ClusterFind_graph = function(sce, runWith = 'PCA', PCNum = 10,cluster.fun) {

  suppressPackageStartupMessages(library(scran, quietly = T))
  suppressPackageStartupMessages(library(bluster, quietly = T))
  #  matrix.run = switch(EXPR       = runWith,
  #                      'VGcounts' = assay(altExp(sce),'VGcounts'),
  #                      'PCA'      = reducedDim(x = sce, type = 'PCA')[,1:PCNum],
  #                      'ICA'      = reducedDim(x = sce, type = 'ICA')[,1:PCNum],
  #                      'BEPCA'    = reducedDim(x = sce, type = 'BEPCA')[,1:PCNum],
  #                      'tSNE'     = reducedDim(x = sce, type = 'tSNE')[,1:2],
  #                      'Umap'     = reducedDim(x = sce, type = 'Umap')[,1:2])
  #
  #  com.graph = MUDAN::getComMembership(matrix.scerun,
  #                                       k=k,
  #                                       method=igraph::cluster_infomap,
  #                                       verbose=FALSE)
  nn.clusters2 <- clusterCells(x=sce, use.dimred= runWith,
                               BLUSPARAM=NNGraphParam( cluster.fun=cluster.fun))


  colData(sce)$cluster = as.character(nn.clusters2)
  #colLabels(sce)<-as.factor(nn.clusters2)
  return(sce)
}

# hierarchical clustering from scran package

ClusterFind_hclust_scran<- function(sce, runWith)
{
  set.seed(1111)
  khclust.info <- clusterCells(x=sce, use.dimred=runWith,
                               BLUSPARAM=TwoStepParam(
                                 first=KmeansParam(centers=1000),
                                 second=HclustParam(method="ward.D2", cut.dynamic=TRUE,
                                                    cut.param=list(deepSplit=3)) # for higher resolution.
                               ),
                               full=TRUE
  )
  colData(sce)$cluster = as.character(khclust.info$clusters)
  #colLabels(sce)<- as.factor(khclust.info$clusters)
  return(sce)

}
# Setting the seed due to the randomness of k-means.
# kmean clustering from scran package
ClusterFind_kmean_scran<- function(sce, runWith,k=30){
  set.seed(0101010)
  kgraph.clusters <- clusterCells(x=sce, use.dimred= runWith,
                                  BLUSPARAM=TwoStepParam(
                                    first=KmeansParam(centers=1000),
                                    second=NNGraphParam(k=k)
                                  )
  )
  colData(sce)$cluster = as.character(kgraph.clusters)
  #colLabels(sce)<- as.factor(kgraph.clusters)
  return(sce)
}




################################################################################
########################## VISUALIZATION #######################################
################################################################################

#' A function to find optimum number of clusters in the dataset
#'
#' @param sce a SingleCellExperiment object
#' @param runWith run with VGcounts/PCA/ICA/tSNE/Umap
#'
#' @param PCNum number of PCs to consider only in PCA/ICA condition
#'
#' @return  a ggplot represents number of clusters
#' @examples
#'
#' ClusterNumSearch(sce)
#'

ClusterNumSearch = function(sce, runWith = 'VGcounts', PCNum = 10) {

  matrix.run = switch(EXPR       = runWith,
                      'VGcounts' = assay(altExp(sce),'VGcounts'),
                      'PCA'      = reducedDim(x = sce, type = 'PCA')[,1:PCNum],
                      'BEPCA'    = reducedDim(x = sce, type = 'PCA')[,1:PCNum],
                      'ICA'      = reducedDim(x = sce, type = 'ICA')[,1:PCNum],
                      'tSNE'     = reducedDim(x = sce, type = 'tSNE')[,1:2],
                      'Umap'     = reducedDim(x = sce, type = 'Umap')[,1:2])

  #library(factoextra, quietly = T)
  #Partitioning methods, such as k-means clustering require the users to specify the number of clusters to be generated.
  res.cluster.num = factoextra::fviz_nbclust(x = matrix.run, FUNcluster = kmeans, method = 'silhouette', k.max = 15)

  plot(res.cluster.num)
}
