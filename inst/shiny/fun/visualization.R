#' Plot highly variable genes
#'
#' This function takes an object of SingleCellExperiment class
#'
#'
#' @param used can be the raw count matrix or the normalized count matrix
#' @param ngene number of genes of interest
#' @param batch true or false
#' @return a list containg a ggplot, a highly variable gene information table and the SingleCellExperiment object
#' @examples
#' sce <- VariableGene_hvg_plot(CS.data, used = 'counts', ngene = 1000, batch = F)
#' the plot can be visualize sce[[1]]
#' the table can be shown sce[[2]]
#' the SingleCellExperiment object can be accessed sce[[3]]
#'
#'
 VariableGene_hvg_plot = function(CS.data, used = 'counts', ngene = 1000, batch = F) {
  # use assay of interest to find top HVG

    edata = assay(CS.data, used) %>% as.data.frame()
    #edata = log2(edata+1)

  if(batch) {
    batch = CS.data@colData$batch

    var = as.data.frame(scran::modelGeneVar(edata, block = batch))
    var=var %>% drop_na()
    top.hvgs2 <- scran::getTopHVGs(var, n=ngene)
    var$Status<-is.element(rownames(var),top.hvgs2)
    top20<-var %>% filter(Status=="TRUE") %>% arrange(FDR) %>% arrange(desc(bio)) %>% head(20)
    selected_genes<-var %>% filter(Status=="TRUE") %>% arrange(FDR)
    colnames(selected_genes)[c(1,4)]<-c("Average_expression","Variance")
    edata.od = edata[top.hvgs2,]


  } else {
    var = as.data.frame(scran::modelGeneVar(edata))
    #var=var %>% drop_na()
    top.hvgs2 <- scran::getTopHVGs(var, n=ngene)
    var$Status<-is.element(rownames(var),top.hvgs2)
    top20<-var %>% filter(Status=="TRUE") %>% arrange(FDR) %>% arrange(desc(bio)) %>% head(20)
    selected_genes<-var %>% filter(Status=="TRUE") %>% arrange(FDR)
    colnames(selected_genes)[c(1,4)]<-c("Average_expression","Variance")
    edata.od = edata[top.hvgs2,]

  }
  altExp(CS.data, 'VGcounts') = SingleCellExperiment(assay = list(VGcounts = edata.od))
  if(any(c("nCount_RNA","nFeature_RNA", "Mito_gene_percent", "Hemoglobin_gene_percent", "Ribosomal_gene_percent") == colnames(colData(CS.data)))){
    CS.data@colData = subset(CS.data@colData, select = -c(nCount_RNA,nFeature_RNA, Mito_gene_percent, Hemoglobin_gene_percent, Ribosomal_gene_percent))
  }
  # top hvg count matrix saved in altExp
  # if(used == 'counts') {
  # altExp(CS.data, 'VGcounts') = SingleCellExperiment(assay = list(VGcounts = edata.od))
  # }
  # if(used == 'BENMcounts') {
  #   SingleCellExperiment::altExp(CS.data, 'BENMVGcounts') = SingleCellExperiment::SingleCellExperiment(assays = list(BENMVGcounts = edata.od))
  # }
  # if(used == 'NMcounts') {
  #   altExp(CS.data, 'NMVGcounts') = SingleCellExperiment(assay = list(NMVGcounts = edata.od))
  # }
  # if(used == 'BEcounts') {
  #   altExp(CS.data, 'BEVGcounts') = SingleCellExperiment(assay = list(BEVGcounts = edata.od))
  # }
  #

  p <-ggplot(var, aes(x = mean, y = bio)) +
  geom_point(colour = ifelse(var$Status=="TRUE","red","black"), size = 1.5, alpha = 1.5) +
    ggrepel::geom_text_repel(data = top20, mapping = aes(label = rownames(top20),x = mean,  y = bio), box.padding=unit(1, "lines"),
                             point.padding=unit(0.5, "lines"),
                             segment.colour = "purple",segment.size = 0.5,segment.alpha = 0.5,max.overlaps = Inf) +
    geom_point(data = top20, mapping = aes(label = rownames(top20)), color = "purple") + cowplot::theme_cowplot()+
    labs(x="Average Expression",y="Standardized Variance")

  return(list(p,selected_genes,CS.data))
}

#Dimension loadings plot (PCA)
dim_lodingsViz<-function(CS.data,ndims, ngenes=2){

  #condition to select ngenes to plot
  if(ndims<=6) ngenes = 20
  if(ndims >6 & ndims <=12) ngenes = 10
  if(ndims >12 & ndims <=24) ngenes =5
  if(ndims >24) ngenes = 2

  plots <- lapply(
    1:ndims,
    FUN=function(i) {
      for(j in reducedDimNames(CS.data)){
        if(j=="BEPCA"){
          rdim_name= "BEPCA"
        }
        else if (j == "PCA") {
          rdim_name= "PCA"
        }
      }
      loading<-reducedDim(CS.data,rdim_name)
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

#Dim Heatmap
dim_heatmap=function(CS.data,ndims,nfeatures){
  plots <- lapply(
    1:ndims,
    FUN=function(i) {
      for(j in SingleCellExperiment::reducedDimNames(CS.data)){
        if(j=="BEPCA"){
          rdim_name= "BEPCA"
        }
        else if (j == "PCA") {
          rdim_name= "PCA"
        }
      }

      loading<-reducedDim(CS.data,rdim_name)
      loading<-attr(loading,"rotation")
      loading<-as.data.frame(loading)
      data.plot=head(loading[order(abs(loading[,i]),decreasing=TRUE),],nfeatures)[i]
      for(k in assayNames(CS.data)){
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
      edata<-assay(CS.data, count)
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

ElbowPlot <- function(CS.data, ndims = 20, reduction = 'PCA') {
  for(j in SingleCellExperiment::reducedDimNames(CS.data)){
    if(j=="BEPCA"){
      rdim_name= "BEPCA"
    }
    else if (j == "PCA") {
      rdim_name= "PCA"
    }
  }

  data.use<-reducedDim(CS.data,rdim_name)
  data.use<-attr(data.use,"percentVar")

  #data.use=reducedDim(CS.data, 'PCA')
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

PCA_run<-function(CS.data,topGenes)
{
  data=scater::runPCA(CS.data,subset_row=topGenes)
  return(data)
}

volcano_plot<- function(CS.data, fc=2, pval= 0.005){
  library(ggpubr)
  library(ggrepel)
  d=CS.data
  names(d)[1]=c("log2FC")
  d$label= rownames(d)


  d=d[!is.na(d[,2]),]
  d$log10=-log10(d[,2])  # col 2 fdr
  d$class="none"
  if(nrow(d[d[,2] <= pval & d[,1] >= fc,])!=0) d[d[,2] <= pval & d[,1] >= fc,]$class="UP"
  if(nrow(d[d[,2] <= pval & d[,1] <= -fc,])!=0) d[d[,2] <= pval & d[,1] <= -fc,]$class="DOWN"
  UP <- d[grepl("UP", d$class),]
  up=head(UP[order(UP$log2FC, decreasing = T), ],10)
  DOWN <- d[grepl("DOWN", d$class),]
  down=head(DOWN[order(DOWN$log2FC), ],10)
  top=rbind(up,down)
  up_num=nrow(d[d$class == "UP",])
  down_num=nrow(d[d$class == "DOWN",])
  d$class <- as.factor(d$class) # col 3 class
  xval=ceiling(max(abs(d[,1])))
  colors <- c("UP"="#FC4E07", "none"="#E7B800", "DOWN"="#00AFBB")
  p=ggplot(data=d, aes(x=log2FC,y=log10,color=class,label=label)) +
    geom_point(data = d[d$class=="UP",],aes(y=log10,color="UP"))+
    geom_point(data = d[d$class=="none",],aes(y=log10,color="none"))+
    geom_point(data = d[d$class=="DOWN",],aes(y=log10,color="DOWN"))+
    scale_color_manual(values = colors)+
    ylab(paste0("-log10 ",names(d)[2]))+xlab("Log2FoldChange")+
    theme_pubr(base_size = 12,border=TRUE)+geom_hline(yintercept=-log10(pval), linetype="dashed",
                                                      color = "black", linewidth=0.5)+geom_vline(xintercept=c(-fc,fc), linetype="dashed",
                                                                                            color = "black", linewidth=0.5)+geom_label_repel(data=top,

                                                                                              fontface="bold",
                                                                                              color="purple",
                                                                                              box.padding=unit(1, "lines"),
                                                                                              point.padding=unit(0.5, "lines"),
                                                                                              segment.colour = "purple",segment.size = 0.5,segment.alpha = 0.5,max.overlaps = Inf)+
    annotate(geom = 'text', label = paste0('UP_Number: ', up_num), x = Inf, y = Inf, hjust = 1.1, vjust = 1.5)+
    annotate(geom = 'text', label = paste0('DOWN_Number: ', down_num), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5)+labs(color = "class")
return(p)
}


#DEG Table
DEG_table<- function(CS.data, fc=2, pval= 0.005){
  library(ggpubr)
  library(ggrepel)
  d=CS.data
  names(d)[1]=c("log2FC")
  d$label= rownames(d)


  d=d[!is.na(d[,2]),]
   # col 2 fdr
  d$class="none"
  if(nrow(d[d[,2] <= pval & d[,1] >= fc,])!=0) d[d[,2] <= pval & d[,1] >= fc,]$class="UP"
  if(nrow(d[d[,2] <= pval & d[,1] <= -fc,])!=0) d[d[,2] <= pval & d[,1] <= -fc,]$class="DOWN"
  UP <- d[grepl("UP", d$class),]
  DOWN <- d[grepl("DOWN", d$class),]
  NONE <- d[grepl("none", d$class),]
  up_num=nrow(d[d$class == "UP",])
  down_num=nrow(d[d$class == "DOWN",])
  nont_reg<- nrow(d[d$class == "none",])
  table<- cbind(up_num,down_num, nont_reg ) %>% as.data.frame()
  colnames(table) <- c("UP", "Down", "None")
  rownames(table) <- "DEG"
  return(table)
 }


Cell_Type_Annotation<- function(data=CS.data,method='SingleR', dataset="BlueprintEncodeData", TissueType=TissueType,  clusters='', used='counts')
{
  clusters=data@colData$cluster
  # CS.data= SingleCellExperiment object
  # dataset = for singleR method: celldex::
  # TissueType = input for sctype (Lung , liver ...)
  # method= SingleR/sctype
  # cluster= clusters column
  # reduction = PCA/Umap/tSNE

  if(method=="SingleR")
  {
    if(dataset=="BlueprintEncodeData")
    {
      ref=celldex::BlueprintEncodeData()
    }
    else if(dataset== "DatabaseImmuneCellExpressionData"){
      ref= celldex::DatabaseImmuneCellExpressionData()
    }

    else if(dataset == "HumanPrimaryCellAtlasData"){
      ref = celldex::HumanPrimaryCellAtlasData()
    }
    else if(dataset == "ImmGenData"){
      ref= celldex::ImmGenData()
    }
    else if(dataset =="MonacoImmuneData"){
      ref = celldex::MonacoImmuneData()
    }
    else if (dataset== "MouseRNAseqData") {
      ref = celldex::MouseRNAseqData()
    }


    pred <- SingleR::SingleR(test=data, ref=ref, labels= ref$label.main, assay.type.test = used)
    tab <- table(Assigned=pred$pruned.labels, Cluster=clusters)
    #p0=pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101))

    pred1 <- SingleR::SingleR(test=data, ref=ref, labels= ref$label.main, clusters=clusters, assay.type.test = used)
    pred1$Var1=rownames(pred1)
    data@colData$cellType=""
    data@colData$cellType<-pred1$pruned.labels[match(clusters, pred1$Var1)]

   # p= gridExtra::grid.arrange(
   #    plotReducedDim(data, dimred = reduction, colour_by="cellType", text_by="cellType"),
   #    plotReducedDim(data, dimred = reduction, colour_by="cluster", text_by = "cluster")
   #  )

  }
  else if (method== "sctype") {

    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
    # load cell type annotation function
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

    db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";


    tissue = TissueType
    gs_list = gene_sets_prepare(db_, tissue)

    scRNAseqData= as.matrix(assay(data, used))

    es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE,
                          gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
    cL_resutls = do.call("rbind", lapply(unique(data@colData$cluster), function(cl){
      es.max.cl = sort(rowSums(es.max[ ,rownames(data@colData[data@colData$cluster==cl, ])]), decreasing = !0)
      head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(data@colData$cluster==cl)), 10)
    }))

    sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

    sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
    print(sctype_scores[,1:3])
    data@colData$cellType = ""


    for(j in unique(sctype_scores$cluster)){
      cl_type = sctype_scores[sctype_scores$cluster==j,];
      data@colData$cellType[data@colData$cluster == j] = as.character(cl_type$type[1])
    }

   # p= plotReducedDim(data,dimred = reduction, colour_by = "cellType", text_by = "cellType")

    }

  return(data)
}


#Plot_raw data
Plot_raw_data<-function(sce, colour_by="All"){

  cowplot::plot_grid(plotColData(sce, y = "nFeature_RNA", x = colour_by, colour_by = colour_by),
                     scater::plotColData(sce, y = "nCount_RNA", x = colour_by, colour_by = colour_by),
                      plotColData(sce, y = "Mito_gene_percent",x = colour_by, colour_by = colour_by),
                      plotColData(sce, y = "Hemoglobin_gene_percent",x = colour_by, colour_by = colour_by),
                      plotColData(sce, y = "Ribosomal_gene_percent",x = colour_by, colour_by = colour_by),
                      plotColData(sce, x = "nCount_RNA", y = "nFeature_RNA", colour_by = colour_by),
                      plotColData(sce, x = "nFeature_RNA", y = "Mito_gene_percent", colour_by = colour_by),
                      plotColData(sce, x = "nFeature_RNA", y = "Ribosomal_gene_percent", colour_by = colour_by),
                      plotColData(sce, x = "nFeature_RNA", y = "Hemoglobin_gene_percent", colour_by = colour_by)
                     ,ncol = 3)

}
#p + scale_color_manual(values=c("#999999"))
#Plot QC after filteration

#
# ggplot(edata, aes(y=nFeature_RNA, x=Type, fill=Type)) + geom_violin( alpha=0.9) + cowplot::theme_cowplot() + geom_point(position = position_jitter(seed = 1, width = 0.2), size= 1, shape=19)
# p2=edata %>% ggplot(aes(y=nCount_RNA, x=featureVector, fill=featureVector)) +
#   geom_violin( alpha=0.9, position = position_dodge(width = 0.9)) + cowplot::theme_cowplot() +
#   geom_point(position = position_jitter(seed = 1, width = 0.2), size= 1, shape=19) +
#   theme(legend.position = "none") +
#   labs(x = "PBMC", title = "nCount_RNA")+
#   theme(axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5))
#
#
# edata_filt %>% filter(nFeature_RNA> 200 & nFeature_RNA<15000 & nCount_RNA> 1000 & nCount_RNA<1000000 & Ribosomal_gene_percent>5) %>%
#   ggplot(aes(y=nCount_RNA, x=featureVector, fill=featureVector)) +
#   geom_violin( alpha=0.9, show.legend = FALSE) + cowplot::theme_cowplot() +
#   geom_point(position = position_jitter(seed = 1, width = 0.2), size= 1, shape=19) +
#
#   labs(x = "PBMC", title = "nCount_RNA", y=NULL)+
#   theme(axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5))
#
#Plot Filterd data
plot_filterd_data<- function(scdata,  y=NULL, colour_by=NULL, x=NULL){
  edata=scdata%>% colData() %>% data.frame()
  if(is.null(x)){
    if(colour_by=="All"){

      p= edata %>%  ggplot(aes(y=edata[,y], x=All, fill=All)) +
        geom_violin( alpha=0.3, show.legend = FALSE) + cowplot::theme_cowplot() +
        # geom_point(position = position_jitter(seed = 1, width = 0.2), size= 1, shape=19, show.legend = FALSE) +

        labs(x = "All", title =y, y=NULL)+
        theme(axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5))
    }

    else {

      p=edata %>% ggplot(aes(y=edata[,y], x=edata[,colour_by], fill=edata[,colour_by])) +
        geom_violin( alpha=0.3, show.legend = FALSE) + cowplot::theme_cowplot() +
        #geom_point(position = position_jitter(seed = 1, width = 0.2), size= 1, shape=19, show.legend = FALSE) +

        labs(x =colour_by , title = y, y=NULL)+
        theme(plot.title = element_text(hjust = 0.5))
    }
  }


  else{
    if(colour_by=="All"){

      p= edata %>%  ggplot(aes(y=edata[,y], x=edata[,x], col=All)) +
        geom_point( alpha=0.3, show.legend = TRUE, size= 3) + cowplot::theme_cowplot() +
        # geom_point(position = position_jitter(seed = 1, width = 0.2), size= 1, shape=19, show.legend = FALSE) +

        labs(x = x, title =NULL, y=y)+
        theme(axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5))
    }

    else {

      p=edata %>% ggplot(aes(y=edata[,y], x=edata[,x], col=edata[,colour_by])) +
        geom_point( alpha=0.3, show.legend = TRUE, size =3) + cowplot::theme_cowplot() +
        #geom_point(position = position_jitter(seed = 1, width = 0.2), size= 1, shape=19, show.legend = FALSE) +

        labs(x = x, title =NULL, y=y, color = colour_by)+
        theme(plot.title = element_text(hjust = 0.5))

    }
  }
  return(p)
}


### Plot Marker genes
plot_Markers <- function(sce, method, colour_by, dimred = "PCA", features, by_exprs_values, x, text_by= 'cluster' ){
  plot = switch(EXPR = method,
               Expression = scater::plotExpression(object = sce, exprs_values= by_exprs_values, features = features, colour_by = colour_by, x= colour_by ),
               Reduction = plot_Markers_reduction(sce, features = features, by_exprs_values = by_exprs_values, dimred = dimred, text_by =text_by)
  )
  return(plot)
}

plot_Markers_reduction <- function(sce, dimred, colour_by, by_exprs_values, features, text_by ){
  plotlist <- list()
  for (i in features) {
    plotlist[[i]] <- scater::plotReducedDim(sce, dimred =dimred , colour_by = i, by_exprs_values = by_exprs_values, text_by=text_by) +
      scale_fill_gradientn(colours = colorRampPalette(c("grey90", "orange3", "firebrick",
                                                        "firebrick", "red", "red"))(10)) + ggtitle(label = i) + theme(plot.title = element_text(size = 20))
  }
  if(length(plotlist)==1){
    plotlist
  }
  else if(length(plotlist)==2){
    cowplot::plot_grid(ncol = 2, plotlist = plotlist)
  }
  else{
    cowplot::plot_grid(ncol = 3, plotlist = plotlist)
  }

}




















































