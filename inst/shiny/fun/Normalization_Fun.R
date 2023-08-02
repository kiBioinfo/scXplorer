################################################################################
######################### NORMALIZATION METHODS ################################
################################################################################

# Normalize the raw count

#' A function to normalize data
#'
#' @param sce a SingleCellExperiment object
#' @param method method for narmalization
#' @param used count matrix
#' @param batch batch True/false
#' @return  a SingleCellExperiment object with normalized data
#' @examples
#'
#' Normalize_Matrix(sce)
#'

Normalize_Matrix = function(sce, method = 'lim', used = 'counts', batch = F) {

  sce = switch(EXPR = method,
               'sct'    = Normalize_scVST(   sce = sce, used = used),
               'vst'    = Normalize_dqVST(   sce = sce, used = used),
               'lib'    = Normalize_libsize( sce = sce, used = used),
               'scn'    = Normalize_scNorm(  sce = sce, used = used),
               'lim'    = Normalize_linnorm( sce = sce, used = used),
               'LogNormalize' = Normalize_LogNormalize(sce = sce, used = used),
               'ruvseq' = Normalize_ruvseq(sce = sce, used = used)
  )
  #sce = Normalize_scVST(sce)
  #sce = Normalize_dqVST(sce)
  #sce = Normalize_libsize(sce)
  #sce = Normalize_scNorm(sce)
  #sce = Normalize_linnorm(sce)

  return(sce)
}

.Normalize_Matrix = function(sce, method = 'lim', used = 'counts') {

  sce = switch(EXPR = method,
               'sct'    = .Normalize_scVST(   sce = sce, used = used),
               'vst'    = .Normalize_dqVST(   sce = sce, used = used),
               'lib'    = .Normalize_libsize( sce = sce, used = used),
               'scn'    = .Normalize_scNorm(  sce = sce, used = used),
               'lim'    = .Normalize_linnorm( sce = sce, used = used))

  return(sce)
}

# Log normalization
Normalize_LogNormalize= function(sce, used = 'counts') {
  sce <- scran::computeSumFactors(sce)
  sce <- scater::logNormCounts(sce)
  return(sce)
}
#Variance stabilizing transformation using sctransform
Normalize_scVST = function(sce, used = 'counts') {
  #View(assay(sce, used))

  library(SingleCellExperiment, quietly = T)
  edata = assay(sce, used)
  #batch = sce@colData$batch
  #pheno = colData(sce)
  #
  library(sctransform, quietly = T)
  edata_vst = sctransform::vst(as.matrix(edata))$y

  edata_vst = as.data.frame(edata_vst)[rownames(edata),]
  row.names(edata_vst) = rownames(edata)
  edata_vst[is.na(edata_vst)] = 0

  assay(sce, 'NMcounts') = as.matrix(edata_vst)
  return(sce)
}

.Normalize_scVST = function(sce, used = 'counts') {
  #View(assay(sce, used))

  library(SingleCellExperiment, quietly = T)
  edata = assay(sce, used)
  batch = sce@colData$batch
  #pheno = colData(sce)

  ###
  unique.batch = unique(batch)

  list.counts = list()
  for(i in c(1 : length(unique.batch))) {
    list.counts[[i]] = edata[, batch == unique.batch[i]]
  }

  suppressPackageStartupMessages(library(sctransform, quietly = T))
  list.counts.norm = lapply(list.counts, function(x) sctransform::vst(umi = x)$y)

  count.matrix = as.matrix(list.counts.norm[[1]])
  for (i in 2:length(list.counts.norm)) {
    count.matrix = merge(count.matrix, as.matrix(list.counts.norm[[i]]), by = "row.names", all=T)
    rownames(count.matrix) = count.matrix$Row.names
    count.matrix = count.matrix[ , !(names(count.matrix) %in% "Row.names")]
  }
  count.matrix[is.na(count.matrix)] = 0

  count.matrix = count.matrix[,colnames(edata)]
  ###

  edata_vst = count.matrix

  edata_vst = as.data.frame(edata_vst)[rownames(edata),]
  row.names(edata_vst) = rownames(edata)
  edata_vst[is.na(edata_vst)] = 0

  assay(sce, 'NMcounts') = as.matrix(edata_vst)
  return(sce)
}

#Variance stabilizing transformation using DESeq2
Normalize_dqVST = function(sce, used = 'counts') {
  #View(assay(sce, used))

  library(SingleCellExperiment, quietly = T)
  edata = assay(sce, used)
  #batch = sce@colData$batch
  #pheno = colData(sce)

  library(DESeq2, quietly = T)
  edata_vst = DESeq2::vst(as.matrix(edata))

  assay(sce, 'NMcounts') = as.matrix(edata_vst)
  return(sce)
}

.Normalize_dqVST = function(sce, used = 'counts') {
  #View(assay(sce, used))

  library(SingleCellExperiment, quietly = T)
  edata = assay(sce, used)
  batch = sce@colData$batch
  #pheno = colData(sce)

  ###
  unique.batch = unique(batch)

  list.counts = list()
  for(i in c(1 : length(unique.batch))) {
    list.counts[[i]] = edata[, batch == unique.batch[i]]
  }

  suppressPackageStartupMessages(library(DESeq2, quietly = T))
  list.counts.norm = lapply(list.counts, function(x) DESeq2::varianceStabilizingTransformation(x))

  count.matrix = as.matrix(list.counts.norm[[1]])
  for (i in 2:length(list.counts.norm)) {
    count.matrix = merge(count.matrix, as.matrix(list.counts.norm[[i]]), by = "row.names", all=T)
    rownames(count.matrix) = count.matrix$Row.names
    count.matrix = count.matrix[ , !(names(count.matrix) %in% "Row.names")]
  }
  count.matrix[is.na(count.matrix)] = 0

  count.matrix = count.matrix[,colnames(edata)]
  ###

  edata_vst = count.matrix

  edata_vst = as.data.frame(edata_vst)[rownames(edata),]
  row.names(edata_vst) = rownames(edata)

  assay(sce, 'NMcounts') = as.matrix(edata_vst)
  return(sce)
}


Normalize_libsize = function(sce, used = 'counts') {
  #View(assay(sce, used))
  edata = assay(sce, used)
  #batch = sce@colData$batch
  #pheno = colData(sce)

  library(DESeq2, quietly = T)
  lib.size = estimateSizeFactorsForMatrix(edata)
  edata_ed = t(t(edata)/lib.size)

  assay(sce, 'NMcounts') = as.matrix(edata_ed)
  return(sce)
}

.Normalize_libsize = function(sce, used = 'counts') {
  #View(assay(sce, used))
  edata = assay(sce, used)
  batch = sce@colData$batch
  #pheno = colData(sce)

  ###
  unique.batch = unique(batch)

  list.counts = list()
  for(i in c(1 : length(unique.batch))) {
    list.counts[[i]] = edata[, batch == unique.batch[i]]
  }

  suppressPackageStartupMessages(library(DESeq2, quietly = T))
  list.counts.norm = lapply(list.counts, function(x) t(t(x)/estimateSizeFactorsForMatrix(x)))

  count.matrix = as.matrix(list.counts.norm[[1]])
  for (i in 2:length(list.counts.norm)) {
    count.matrix = merge(count.matrix, as.matrix(list.counts.norm[[i]]), by = "row.names", all=T)
    rownames(count.matrix) = count.matrix$Row.names
    count.matrix = count.matrix[ , !(names(count.matrix) %in% "Row.names")]
  }
  count.matrix[is.na(count.matrix)] = 0

  count.matrix = count.matrix[,colnames(edata)]
  ###
  edata_ed = count.matrix

  edata_ed = as.data.frame(edata_ed)[rownames(edata),]
  row.names(edata_ed) = rownames(edata)

  assay(sce, 'NMcounts') = as.matrix(edata_ed)
  return(sce)
}


Normalize_scNorm = function(sce, used = 'counts') {
  edata = assay(sce, used)
  batch = sce@colData$batch
  #pheno = colData(sce)

  if(is.null(batch)) {
    library(SCnorm, quietly = T)
    #UMI count
    set.seed(0)
    DataNorm = SCnorm(Data = edata, Conditions = rep('x', ncol(edata)),
                      PrintProgressPlots = F,
                      FilterCellNum = 10,
                      PropToUse = .1,
                      NCores = parallel::detectCores(),
                      Thresh = .1,
                      ditherCounts = TRUE)
  } else {
    Conditions = batch

    library(SCnorm, quietly = T)
    #UMI count
    set.seed(0)
    DataNorm = SCnorm(Data = edata, Conditions= Conditions,
                      PrintProgressPlots = F,
                      FilterCellNum = 10,
                      PropToUse = .1,
                      NCores = parallel::detectCores(),
                      Thresh = .1,
                      ditherCounts = TRUE)
  }

  edata_scnorm = assay(DataNorm,'normcounts')

  assay(sce, 'NMcounts') = as.matrix(edata_scnorm)
  return(sce)
}

.Normalize_scNorm = function(sce, used = 'counts') {
  edata = assay(sce, used)
  batch = sce@colData$batch
  #pheno = colData(sce)

  suppressPackageStartupMessages(library(SCnorm, quietly = T))

  ###
  unique.batch = unique(batch)

  list.counts = list()
  for(i in c(1 : length(unique.batch))) {
    list.counts[[i]] = edata[, batch == unique.batch[i]]
  }

  set.seed(0)
  list.counts.norm = lapply(list.counts, function(x) SCnorm::SCnorm(Data = x, Conditions = rep('x', ncol(x)),
                                                                    PrintProgressPlots = F,
                                                                    FilterCellNum = 10,
                                                                    PropToUse = .1,
                                                                    NCores = parallel::detectCores(),
                                                                    Thresh = .1,
                                                                    ditherCounts = TRUE))

  count.matrix = assay(list.counts.norm[[1]],'normcounts')
  for (i in 2:length(list.counts.norm)) {
    count.matrix = merge(count.matrix, assay(list.counts.norm[[i]], 'normcounts'), by = "row.names", all=T)
    rownames(count.matrix) = count.matrix$Row.names
    count.matrix = count.matrix[ , !(names(count.matrix) %in% "Row.names")]
  }
  count.matrix[is.na(count.matrix)] = 0

  count.matrix = count.matrix[,colnames(edata)]
  ###

  edata_scnorm = count.matrix
  edata_scnorm = as.data.frame(edata_scnorm)[rownames(edata),]
  row.names(edata_scnorm) = rownames(edata)

  assay(sce, 'NMcounts') = as.matrix(edata_scnorm)
  return(sce)
}

#Normalization using Linnorm algorithm.
Normalize_linnorm = function(sce, used = 'counts') {
  #View(assay(sce, used))
  edata = assay(sce, used)
  #batch = sce@colData$batch
  #pheno = colData(sce)

  library(Linnorm, quietly = T)
  edata_lin = Linnorm::Linnorm.Norm(datamatrix = edata)


  assay(sce, 'NMcounts') = as.matrix(edata_lin)
  return(sce)
}

.Normalize_linnorm = function(sce, used = 'counts') {
  #View(assay(sce, used))
  edata = assay(sce, used)
  batch = sce@colData$batch
  #pheno = colData(sce)

  ###
  unique.batch = unique(batch)

  list.counts = list()
  for(i in c(1 : length(unique.batch))) {
    list.counts[[i]] = edata[, batch == unique.batch[i]]
  }

  suppressPackageStartupMessages(library(Linnorm, quietly = T))
  list.counts.norm = lapply(list.counts, function(x) Linnorm::Linnorm.Norm(datamatrix = x))

  count.matrix = as.matrix(list.counts.norm[[1]])
  for (i in 2:length(list.counts.norm)) {
    count.matrix = merge(count.matrix, as.matrix(list.counts.norm[[i]]), by = "row.names", all=T)
    rownames(count.matrix) = count.matrix$Row.names
    count.matrix = count.matrix[ , !(names(count.matrix) %in% "Row.names")]
  }
  count.matrix[is.na(count.matrix)] = 0

  count.matrix = count.matrix[,colnames(edata)]
  ###

  edata_lin = count.matrix
  edata_lin = as.data.frame(edata_lin)[rownames(edata),]
  row.names(edata_lin) = rownames(edata)

  assay(sce, 'NMcounts') = as.matrix(edata_lin)
  return(sce)
}

# Normalization using ruvseq package
Normalize_ruvseq = function(sce, used = 'counts', spkins) {
  edata = assay(sce, used)
  #batch = sce@colData$batch
  spkins=rownames(sce)[grep("(?i)^ERCC", rownames(sce))]

  spkins = spkins[spkins %in% rownames(edata)]

  library(RUVSeq, quietly = T)
  edata_ruv = RUVSeq::RUVg(x = as.matrix(edata), cIdx = spkins, k = 2)

  edata_ruv = edata_ruv$normalizedCounts

  assay(sce, 'NMcounts') = edata_ruv

  return(sce)
}


################################################################################
#################### VARIABLE GENE SELECTION METHODS ###########################
################################################################################


#================#
# Most variable genes
#================#


#' A function to find most variable genes in the dataset
#'
#' @param sce a SingleCellExperiment object
#' @param method method to predict most variable genes
#' @param ngene number of genes to consider
#' @param used count matrix
#' @param batch batch T/F
#' @return  a SummarizedExperiment object with vgcounts
#' @examples
#'
#' SelectVariableGene(sce)
#'

SelectVariableGene = function(sce, method = 'mvg', used = 'counts', ngene = 1000, batch = F) {

  if(!batch) {sce = switch(EXPR = method,

                           'seu'   = VariableGene_seu(   sce = sce, used = used, ngene = ngene),
                           'mvg'   = VariableGene_mvg(   sce = sce, used = used, ngene = ngene),
                           'hvg'   = VariableGene_hvg(   sce = sce, used = used, ngene = ngene))
  } else {
    #print("xxxxx")
    sce = VariableGene_hvg(   sce = sce, used = used, ngene = ngene, batch = T)
  }

  return(sce)
}



# Variable Gene selection using Seurat

VariableGene_seu = function(sce, used = 'counts', ngene = 1000) {
  library(Seurat, quietly = T)
  edata = assay(sce, used)

  ctrl = CreateSeuratObject(counts = edata, min.cells = 0, min.features = 0)
  ctrl = NormalizeData(ctrl)
  ctrl = FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = ngene)

  edata.od = edata[match(VariableFeatures(ctrl), rownames(ctrl@assays$RNA@counts)),]

  altExp(sce, 'VGcounts') = SingleCellExperiment(assay = list(VGcounts = edata[rownames(edata.od),]))

  return(sce)
}

# Variable Gene selection using Deseq2 mvg

VariableGene_mvg = function(sce, used = 'counts', ngene = 1000) {
  edata = assay(sce, used)

  library(DESeq2, quietly = T)
  lib.size = estimateSizeFactorsForMatrix(edata)
  ed = t(t(edata)/lib.size)

  cv2Cutoff = 0.3
  geneNumForPlot = 500

  means = rowMeans(ed)
  vars = apply(ed,1,var)
  cv2 = vars/(means^2)

  library(statmod, quietly = T)
  minMeanForFit = unname( quantile( means[ which( cv2 >cv2Cutoff) ], .95 ) )
  useForFit = means >= minMeanForFit # & spikeins
  fit = glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ),cv2[useForFit] )
  a0 = unname( fit$coefficients["a0"] )
  a1 = unname( fit$coefficients["a1tilde"])
  fit$coefficients
  xg = exp(seq( min(log(means[means>0])), max(log(means)), length.out=1000 ))
  vfit = a1/xg + a0
  df = ncol(ed) - 1
  afit = a1/means+a0
  varFitRatio = vars/(afit*means^2)
  varorder = order(varFitRatio,decreasing=T)
  oed = ed[varorder,]
  pdf(file = '~/mvg.pdf', width = 5, height = 5, useDingbats = FALSE)
  par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(means),log(cv2)); lines( log(xg), log(vfit), col="black", lwd=3 ); lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black"); lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black");
  points(log(means[varorder[1:geneNumForPlot]]),log(cv2[varorder[1:geneNumForPlot]]),col="red", pch = 16)
  dev.off()

  edata.od = edata[varorder,]

  altExp(sce, 'VGcounts') = SingleCellExperiment(assay = list(VGcounts = edata.od[1:ngene,]))
  return(sce)
}

# Variable Gene selection using scran
VariableGene_hvg = function(sce, used = 'counts', ngene = 1000, batch = F) {



    edata = assay(sce, 'counts') %>% as.data.frame()
    edata = log2(edata+1)



  if(batch) {
    batch = sce@colData$batch

    dec = scran::modelGeneVar(edata, block = batch)
    top.hvgs2 <- scran::getTopHVGs(dec, n=ngene)

    if(used == 'counts'){
      edata.od = edata[top.hvgs2,]
      altExp(sce, 'VGcounts') = SingleCellExperiment(assay = list(VGcounts = edata.od))
    }
    if(used == 'BEcounts'){
      edata = assay(sce,'BEcounts')
      edata.od = edata[top.hvgs2,]
      altExp(sce, 'BEVGcounts') = SingleCellExperiment(assay = list(BEVGcounts = edata.od))
    }
    if(used == 'BENMcounts'){
      edata = assay(sce,'BENMcounts')
      edata.od = edata[top.hvgs2,]
      altExp(sce, 'BENMVGcounts') = SingleCellExperiment(assay = list(BENMVGcounts = edata.od))
    }

  } else {
    dec = scran::modelGeneVar(edata)
    top.hvgs2 <- scran::getTopHVGs(dec, n=ngene)

    edata.od = edata[top.hvgs2,]

    altExp(sce, 'VGcounts') = SingleCellExperiment(assay = list(VGcounts = edata.od))
  }

  return(sce)
}

################################################################################
########################## VISUALIZATION #######################################
################################################################################

#' A function to plot highly variable genes
#'
#' This function takes an object of SingleCellExperiment class
#'
#'
#' @param used can be the raw count matrix or the normalized count matrix
#' @param ngene number of genes of interest
#' @param batch true or false
#' @return a list containg a ggplot, a highly variable gene information table and the SingleCellExperiment object
#' @examples
#' sce <- VariableGene_hvg_plot(sce, used = 'counts', ngene = 1000, batch = F)
#' the plot can be visualize sce[[1]]
#' the table can be shown sce[[2]]
#' the SingleCellExperiment object can be accessed sce[[3]]
#'
#'
VariableGene_hvg_plot = function(sce, used = 'counts', ngene = 1000, batch = F, convert_to_log = F) {
  # use assay of interest to find top HVG
  if(used == 'counts'){
    edata = assay(sce, used) %>% as.data.frame()
    edata = log2(edata+1)
  }
  else if(used == 'logcounts'){
    edata = assay(sce, used) %>% as.data.frame()
  }
  else if(used == 'NMcounts'){

    if(convert_to_log){
      edata = assay(sce, used) %>% as.data.frame()
      edata = log2(edata+1)
    }
    else{
      edata = assay(sce, used) %>% as.data.frame()
    }
  }
  else if(used == 'BEcounts'){
    edata = assay(sce, used) %>% as.data.frame()
  }
  else if(used == 'BENMcounts'){
    edata = assay(sce, used) %>% as.data.frame()
  }

  if(batch) {
    batch = sce@colData$batch

    var = as.data.frame(scran::modelGeneVar(edata, block = batch))
    var=var %>% drop_na()
    top.hvgs2 <- scran::getTopHVGs(var, n=ngene)
    var$Status<-is.element(rownames(var),top.hvgs2)
    top20<-var %>% filter(Status=="TRUE") %>% arrange(FDR) %>% arrange(desc(bio)) %>% head(20)
    selected_genes<-var %>% filter(Status=="TRUE") %>% arrange(FDR)
    colnames(selected_genes)[c(1,4)]<-c("Average_expression","Variance")
    #Select count matrix for variable genes
    if(used == 'counts'){
      edata = assay(sce, used) %>% as.data.frame()
      edata = log2(edata+1)
    }
    else{
      edata = assay(sce, used) %>% as.data.frame()
    }
    edata.od = edata[top.hvgs2,]


  } else {
    var = as.data.frame(scran::modelGeneVar(edata))
    #var=var %>% drop_na()
    top.hvgs2 <- scran::getTopHVGs(var, n=ngene)
    var$Status<-is.element(rownames(var),top.hvgs2)
    top20<-var %>% filter(Status=="TRUE") %>% arrange(FDR) %>% arrange(desc(bio)) %>% head(20)
    selected_genes<-var %>% filter(Status=="TRUE") %>% arrange(FDR)
    colnames(selected_genes)[c(1,4)]<-c("Average_expression","Variance")
    #Select count matrix for variable genes
    if(used == 'counts'){
      edata = assay(sce, used) %>% as.data.frame()
      edata = log2(edata+1)
    }
    else{
      edata = assay(sce, used) %>% as.data.frame()
    }

    edata.od = edata[top.hvgs2,]

  }

  # altExp(sce, 'VGcounts') = SingleCellExperiment(assay = list(VGcounts = edata.od))
  if(any(c("nCount_RNA","nFeature_RNA", "Mito_gene_percent", "Hemoglobin_gene_percent", "Ribosomal_gene_percent") == colnames(colData(sce)))){
    sce@colData = subset(sce@colData, select = -c(nCount_RNA,nFeature_RNA, Mito_gene_percent, Hemoglobin_gene_percent, Ribosomal_gene_percent))
  }
  # top hvg count matrix saved in altExp
  if(used == 'counts') {
    SingleCellExperiment::altExp(sce, 'VGcounts') = SingleCellExperiment::SingleCellExperiment(assay = list(VGcounts = edata.od))
  }
  if(used == 'NMcounts') {
    SingleCellExperiment::altExp(sce, 'VGcounts') = SingleCellExperiment::SingleCellExperiment(assay = list(VGcounts = edata.od))
  }
  if(used == 'logcounts') {
    SingleCellExperiment::altExp(sce, 'VGcounts') = SingleCellExperiment::SingleCellExperiment(assay = list(VGcounts = edata.od))
  }
  if(used == 'BENMcounts') {
    SingleCellExperiment::altExp(sce, 'BEVGcounts') = SingleCellExperiment::SingleCellExperiment(assays = list(BEVGcounts = edata.od))
  }
  if(used == 'BEcounts') {
    SingleCellExperiment::altExp(sce, 'BEVGcounts') = SingleCellExperiment::SingleCellExperiment(assay = list(BEVGcounts = edata.od))
  }


  p <-ggplot(var, aes(x = mean, y = bio)) +
    geom_point(colour = ifelse(var$Status=="TRUE","red","black"), size = 1.5, alpha = 1.5) +
    ggrepel::geom_text_repel(data = top20, mapping = aes(label = rownames(top20),x = mean,  y = bio), box.padding=unit(1, "lines"),
                             point.padding=unit(0.5, "lines"),
                             segment.colour = "purple",segment.size = 0.5,segment.alpha = 0.5,max.overlaps = Inf) +
    geom_point(data = top20, mapping = aes(label = rownames(top20)), color = "purple") + cowplot::theme_cowplot()+
    labs(x="Average Expression",y="Standardized Variance")

  return(list(p,selected_genes,sce))
}
