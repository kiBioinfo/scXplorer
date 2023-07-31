
################################################################################
##################### METHODS FOR BATCH CORRECTION #############################
################################################################################


#================#
# Batch effect correction -- on expression matrix
#================#
#' A function for batch correction on expression matrix
#'
#' @param sce a SingleCellExperiment object
#' @param method method for batch correction
#' @param used  count matrix raw/Normalized
#'
#' @return  a batch corrected SingleCellExperiment object
#' @examples
#'
#' BatchEffect_Matrix(sce)
#'
BatchEffect_Matrix = function(sce, method = 'cbtseq', used = 'counts') {

  sce = switch(EXPR = method,
               'bas'      = RmBatch_denoise(   sce = sce, used = used),
               'lima'     = RmBatch_limma(     sce = sce, used = used),
               'cbt'      = RmBatch_combat(    sce = sce, used = used),
               'cbtseq'   = RmBatch_combatseq( sce = sce, used = used),
               'mnn'      = RmBatch_mnn(       sce = sce, used = used),
               'sca'      = RmBatch_scano(     sce = sce, used = used))

  return(sce)
}


RmBatch_denoise = function(sce, used = 'counts') {
  cat(crayon::red('===Running denoise!===\n'))
  #View(SummarizedExperiment::assay(sce, used))
  #edata = SummarizedExperiment::assay(sce, used)
  edata = switch(EXPR = used,
                 'counts' = SummarizedExperiment::assay(sce, used),
                 'NMcounts' = SummarizedExperiment::assay(sce, used),
                 'VGcounts' = SummarizedExperiment::assay(SingleCellExperiment::altExp(sce, used), used))

  batch = SingleCellExperiment::colData(sce)[,'batch']
  #library(BASiCS, quietly = T)
  Data = SingleCellExperiment::SingleCellExperiment(assays = list(counts = edata),
                                                    colData = S4Vectors::DataFrame(BatchInfo = batch))
  Chain = BASiCS::BASiCS_MCMC(Data = Data,
                              #N = 120000, Thin = 20, Burn = 10000,
                              N = 1000, Thin = 10, Burn = 500,
                              PrintProgress = FALSE, Regression = TRUE, WithSpikes = F,
                              Threads = 2)#parallel::detectCores())
  edata_Denoised = BASiCS::BASiCS_DenoisedCounts(Data = Data, Chain = Chain)

  if(used == 'counts') {
    SummarizedExperiment::assay(sce, 'BEcounts') = edata_Denoised
  }
  if(used == 'VGcounts') {
    SingleCellExperiment::altExp(sce, 'BEVGcounts') = SingleCellExperiment::SingleCellExperiment(assays = list(BEVGcounts = edata_Denoised))
  }

  return(sce)
}


RmBatch_limma = function(sce, used = 'counts') {
  cat(crayon::red('===Running limma!===\n'))
  edata = switch(EXPR = used,
                 'counts' = log2(SummarizedExperiment::assay(sce, used)+1),
                 'NMcounts' = SummarizedExperiment::assay(sce, used),
                 'VGcounts' = log2(SummarizedExperiment::assay(SingleCellExperiment::altExp(sce), used)+1))
  edata_rbe = limma::removeBatchEffect(x = edata, batch = SingleCellExperiment::colData(sce)[,'batch'])
  edata_rbe[edata_rbe < 0] = 0
  #SummarizedExperiment::assay(sce, 'BEcounts') = edata_rbe
  if(used == 'counts') {
    SummarizedExperiment::assay(sce, 'BEcounts') = edata_rbe
  }
  if(used == 'VGcounts') {
    SingleCellExperiment::altExp(sce, 'BEVGcounts') = SingleCellExperiment::SingleCellExperiment(assays = list(BEVGcounts = edata_rbe))
  }
  if(used== 'NMcounts' ){
    SummarizedExperiment::assay(sce, 'BENMcounts') = edata_rbe
  }
  return(sce)

}


RmBatch_combat = function(sce, used = 'counts') {
  cat(crayon::red('===Running combat!===\n'))
  edata = switch(EXPR = used,
                 'counts' = log2(SummarizedExperiment::assay(sce, used)+1),
                 'NMcounts' = SummarizedExperiment::assay(sce, used),
                 'VGcounts' = log2(SummarizedExperiment::assay(SingleCellExperiment::altExp(sce), used)+1))
  batch = SingleCellExperiment::colData(sce)[,'batch']

  var.data = apply(edata, 1, var)
  if(sum(var.data == 0) > 0) {
    rm.edata = edata[which(var.data == 0),,drop = F]

    edata = edata[-which(var.data == 0 ),]
    modcombat = stats::model.matrix(~1, data = SummarizedExperiment::colData(sce))
    edata_combat = sva::ComBat(dat = as.matrix(edata), batch = batch, mod = modcombat, par.prior = TRUE)

    edata_combat = rbind(edata_combat, rm.edata)
    edata_combat = edata_combat[rownames(SummarizedExperiment::assay(sce)),]
  } else {
    modcombat = stats::model.matrix(~1, data = SummarizedExperiment::colData(sce))
    edata_combat = sva::ComBat(dat = as.matrix(edata), batch = batch, mod = modcombat, par.prior=TRUE)
  }

  if(used == 'counts') {
    SummarizedExperiment::assay(sce, 'BEcounts') = edata_combat
  }
  if(used == 'VGcounts') {
    SingleCellExperiment::altExp(sce, 'BEVGcounts') = SingleCellExperiment::SingleCellExperiment(assays = list(BEVGcounts = edata_combat))
  }
  if(used== 'NMcounts' ){
    SummarizedExperiment::assay(sce, 'BENMcounts') = edata_combat
  }
  return(sce)
}

RmBatch_combatseq = function(sce, used = 'counts') {
  cat(crayon::red('===Running combatseq!===\n'))
  edata = switch(EXPR = used,
                 'counts' = SummarizedExperiment::assay(sce, used),
                 'NMcounts' = SummarizedExperiment::assay(sce, used),
                 'VGcounts' = SummarizedExperiment::assay(SingleCellExperiment::altExp(sce,used), used))
  batch = SingleCellExperiment::colData(sce)[,'batch']
  edata_combatseq = sva::ComBat_seq(counts = edata, batch = batch)
  if(used == 'counts') {
    SummarizedExperiment::assay(sce, 'BEcounts') = log2(edata_combatseq + 1)
  }
  if(used == 'VGcounts') {
    SingleCellExperiment::altExp(sce, 'BEVGcounts') = SingleCellExperiment::SingleCellExperiment(assays = list(BEVGcounts = log2(edata_combatseq + 1)))
  }
  if(used== 'NMcounts' ){
    SummarizedExperiment::assay(sce, 'BENMcounts') = edata_combatseq
  }
  return(sce)
}

#RmBatch_mnn = function(sce, used = 'counts') {
#  edata = SummarizedExperiment::assay(sce, used)
#  batch = SingleCellExperiment::colData(sce)[,'batch']
#
#  combined = SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(edata)))
#  combined = batchelor::multiBatchNorm(combined, batch = batch)
#  combined$batch = batch
#
#  f.out2 = batchelor::fastMNN(combined, batch = combined$batch)
#
#  edata_mnn = as.matrix(SummarizedExperiment::assay(f.out2))
#
#  SummarizedExperiment::assay(sce, 'BEcounts') = edata_mnn
#  return(sce)
#}

RmBatch_scano = function(sce, used = 'counts') {
  cat(crayon::red('===Running scanoramaCT!===\n'))
  edata = switch(EXPR = used,
                 'counts' = log2(SummarizedExperiment::assay(sce, used)+1),
                 'NMcounts' = SummarizedExperiment::assay(sce, used),
                 'VGcounts' = log2(SummarizedExperiment::assay(SingleCellExperiment::altExp(sce, used), used)+1))
  batch = SingleCellExperiment::colData(sce)[,'batch']

  list.matrix = list()
  for(i in c(1: length(unique(batch)))) {
    edata.tmp = as.data.frame(edata[,unique(batch)[i] == batch])
    list.matrix[[i]] = edata.tmp
  }

  processedMats = list.matrix
  cell_list = unlist(lapply(processedMats, colnames))
  genes_list = lapply(processedMats, rownames)

  processedMats = lapply(processedMats, t)
  #countsOrig = lapply(processedMats, function(x) colSums(x > 0))#countsNorm = lapply(countsOrig, function(x) { (x - min(x))/(max(x) - min(x)) })#countsNorm = lapply(countsNorm, function(x) cbind(x, x))
  scanoramaCT = reticulate::import("scanoramaCT")
  np = reticulate::import("numpy")
  #merge = scanoramaCT$merge_datasets(processedMats, genes_list)#process = scanoramaCT$process_data(merge[[1]], merge[[2]], #hvg = 0L, dimred = 100L)#gcalign = scanoramaCT$find_alignments_gc(process[[1]], countsNorm)
  # Batch correction
  integrated.corrected.data = scanoramaCT$correct(processedMats,
                                                  genes_list,
                                                  return_dimred = TRUE,
                                                  return_dense = T,verbose = F)
  #countsNormCorrected = unlist(lapply(gcalign, function(x) x[,1]))
  matCorrected = t(do.call(rbind, integrated.corrected.data[[2]]))
  rownames(matCorrected) = integrated.corrected.data[[3]]
  m = do.call(rbind, integrated.corrected.data[[1]])
  colnames(matCorrected) = cell_list

  matCorrected = matCorrected[rownames(edata), colnames(edata)]
  if(used == 'counts') {
    SummarizedExperiment::assay(sce, 'BEcounts') = matCorrected
  }
  if(used == 'VGcounts') {
    SingleCellExperiment::altExp(sce, 'BEVGcounts') = SingleCellExperiment::SingleCellExperiment(assays = list(BEVGcounts = matCorrected))
  }
  if(used== 'NMcounts' ){
    SummarizedExperiment::assay(sce, 'BENMcounts') = matCorrected
  }

  return(sce)
}


#================#
# Batch effect correction -- on dimension reduction
#================#

#' A function for batch correction on dimension reduction
#'
#' @param sce a SingleCellExperiment object
#' @param method method for batch correction onon dimension reduction
#' @param used  count matrix raw/Normalized
#'
#' @return  a batch corrected SingleCellExperiment object
#' @examples
#'
#' BatchEffect_DimReduction(sce)
#'
BatchEffect_DimReduction = function(sce, method = 'beer', used = 'counts') {
  sce = switch(EXPR = method,
               'beer'  = RmBatch_beer( sce = sce, used = used),
               'cca'   = RmBatch_cca(  sce = sce, used = used),
               'dmnn'  = RmBatch_dmnn( sce = sce, used = used),
               'hmy'   = RmBatch_hmy(  sce = sce, used = used),
               'lig'   = RmBatch_lig(  sce = sce, used = used))

  return(sce)

}


RmBatch_beer = function(sce, used = 'counts', showBy = 'PCA') {

  #BeerPath = readline("Path of BEER source file:")
  BeerPath = 'fun/BEER.R'
  #BeerPath = gsub(pattern = '\"', replacement = '', x = BeerPath)
  #BeerPath = gsub(pattern = '\'', replacement = '', x = BeerPath)

  source(BeerPath)
  #source("~/T1/BE/BEER/BEER-0.1.7/BEER.R")

  edata = assay(sce, used)
  batch = sce@colData$batch

  mybeer=BEER(DATA = edata,
              BATCH = batch,
              GNUM=30,
              PCNUM=50,
              ROUND=1,
              GN=2000,
              SEED=1,
              COMBAT=TRUE,
              RMG=NULL)

  pbmc = mybeer$seurat
  PCUSE = mybeer$select

  pcs_beer = pbmc@reductions$pca@cell.embeddings[,PCUSE]

  #if(showBy == 'PCA') {
  reducedDim(sce, 'BEPCA') = pcs_beer
  #}
  return(sce)
}


RmBatch_cca = function(sce, used = 'counts', showby = 'PCA') {
  edata = assay(sce, used)
  batch = sce@colData$batch

  batch.factor = unique(batch)

  suppressPackageStartupMessages(library(Seurat, quietly = T))
  i=1
  list.data = list(edata[,which(batch == batch.factor[i])])
  for(i in c(2:length(batch.factor))) {
    list.data = append(list.data, list(edata[,which(batch == batch.factor[i])]))
  }
  names(list.data) = batch.factor

  list.data.obj = mapply(function(x, y) CreateSeuratObject(counts = x, project = y), list.data, names(list.data))

  immune.anchors = FindIntegrationAnchors(object.list = list.data.obj, dims = 1:20, k.filter = round(min(table(batch))/2) )
  immune.combined = IntegrateData(anchorset = immune.anchors, dims = 1:20)
  DefaultAssay(immune.combined) = "integrated"

  immune.combined = ScaleData(immune.combined, verbose = FALSE)
  immune.combined = RunPCA(immune.combined, npcs = 30, verbose = FALSE)

  reducedDim(sce, 'BEPCA') = as.matrix(immune.combined@reductions$pca@cell.embeddings)

  return(sce)
}


RmBatch_dmnn = function(sce, used = 'counts', showby = 'PCA') {
  edata = assay(sce, used)
  batch = sce@colData$batch

  #suppressPackageStartupMessages(library(SingleCellExperiment, quietly = T))
  combined = SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(edata)))

  #suppressPackageStartupMessages(library(batchelor, quietly = T))
  combined = batchelor::multiBatchNorm(combined, batch = batch)

  #suppressPackageStartupMessages(library(scran, quietly = T))
  chosen.hvgs = scran::modelGeneVar(combined)
  chosen.hvgs <- chosen.hvgs$bio > 0

  #suppressPackageStartupMessages(library(scater, quietly = T))
  combined <- scater::runPCA(combined, subset_row=chosen.hvgs)
  combined$batch = batch

  set.seed(0)
  f.out2 <- batchelor::fastMNN(combined, batch=combined$batch, subset.row=chosen.hvgs)

  reducedDim(sce, 'BEPCA') = reducedDim(f.out2, "corrected")
  #plot(reducedDim(f.out2, "corrected")[,c(1,2)], col = rainbow(2)[as.factor(groupVector)])
  return(sce)
}


RmBatch_hmy = function(sce, used = 'counts', showby = 'PCA') {
  edata = assay(sce, used)
  batch = sce@colData$batch

  #suppressPackageStartupMessages(library(SingleCellExperiment, quietly = T))
  combined = SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(edata)))

  #suppressPackageStartupMessages(library(batchelor, quietly = T))
  combined = batchelor::multiBatchNorm(combined, batch = batch)

  #suppressPackageStartupMessages(library(scran, quietly = T))
  chosen.hvgs = scran::modelGeneVar(combined)
  chosen.hvgs <- chosen.hvgs$bio > 0

  #suppressPackageStartupMessages(library(scater, quietly = T))
  combined <- scater::runPCA(combined, subset_row=chosen.hvgs)

  #plot(reducedDim(combined, 'PCA')[,c(1,2)], col = as.factor(meta.table$Patients))
  #plot(reducedDim(sce, 'PCA')[,c(1,2)], col = as.factor(meta.table$Patients))

  harmony_embeddings = harmony::HarmonyMatrix(data_mat = reducedDim(combined, 'PCA')[,c(1:5)],
                                              meta_data = data.frame(batch = batch), vars_use = "batch",
                                              do_pca=FALSE)

  #plot(harmony_embeddings[,c(1,2)], col = as.factor(meta.table$Patients))

  #harmony_embeddings2 <- harmony::HarmonyMatrix(data_mat = reducedDim(sce, 'PCA')[,c(1:5)],
  #                                             meta_data = data.frame(batch = batch), vars_use = "batch",
  #                                             do_pca=FALSE)
  #plot(harmony_embeddings2[,c(1,2)], col = as.factor(meta.table$Patients))

  reducedDim(sce, 'BEPCA') = harmony_embeddings
  return(sce)
}


RmBatch_lig = function(sce, used = 'counts', showby = 'PCA') {
  edata = assay(sce, used)
  batch = sce@colData$batch

  batch.factor = unique(batch)

  #suppressPackageStartupMessages(library(liger, quietly = T))
  i=1
  list.data = list(edata[,which(batch == batch.factor[i])])
  for(i in c(2:length(batch.factor))) {
    list.data = append(list.data, list(edata[,which(batch == batch.factor[i])]))
  }
  names(list.data) = batch.factor

  ifnb_liger = liger::createLiger(list.data)
  ifnb_liger = liger::normalize(ifnb_liger)
  ifnb_liger = liger::selectGenes(ifnb_liger)
  ifnb_liger = liger::scaleNotCenter(ifnb_liger)

  ifnb_liger <- liger::optimizeALS(ifnb_liger, k = 5)
  ifnb_liger <- liger::quantile_norm(ifnb_liger)

  data.use <- ifnb_liger@H.norm
  #plot(data.use[,c(1,2)], col = as.factor(meta.table$Patients))

  data.use.pca = prcomp(data.use)$x
  #plot(data.use.pca[,c(1,2)], col = as.factor(meta.table$Patients))

  reducedDim(sce, 'BEPCA') = data.use.pca
  return(sce)
}


################################################################################
########################## VISUALIZATION #######################################
################################################################################

#================#
# Batch effect correction -- Evaluation
#================#
#' A function for batch correction Evaluation
#'
#' @param sce a SingleCellExperiment object
#' @param method method for batch correction evaluation
#' @param used batch corrected count matrix/ batch corrected PCA
#' @param PCNum number of PCs to use
#' @return  a ggplot shows batch correction evaluation score
#' @examples
#'
#' BatchEva(sce)
#'
BatchEva = function(sce, method = 'cmx', used = NULL, PCNum = 3) {

  # Condition to select count matrix
  if(!is.null(used)){
    used=used
  }
  else if ("BENMcounts" %in% assayNames(sce)) {
    used= "BENMcounts"
  }
  else{
    used= "BEPCA"
  }
  evaPlot = switch(EXPR = method,
                   'cmx' = BatchEva_cmix( sce = sce, used = used, PCNum = PCNum),
                   'kbet' = BatchEva_kbet( sce = sce, used = used, PCNum = PCNum),
                   'slt'  = BatchEva_slt(  sce = sce, used = used, PCNum = PCNum))

  return(evaPlot)
}



BatchEva_kbet = function(sce, used = 'BEPCA', PCNum = 3) {

  if(used == 'BEPCA') {
    emb.bepca = reducedDim(sce, used)
    if(ncol(emb.bepca) < PCNum) {
      PCNum = ncol(emb.bepca)
    }
  }

  if(used == 'PCA') {
    emb.bepca = reducedDim(sce, used)
    if(ncol(emb.bepca) < PCNum) {
      PCNum = ncol(emb.bepca)
    }
  }

  edata = switch (EXPR = used,
                  'counts' = assay(sce, 'counts'),
                  'VGcounts' = log2(assay(altExp(sce),'VGcounts') + 1),
                  'PCA' = t(reducedDim(sce, used)[,c(1:PCNum)]),
                  'BEPCA' = t(reducedDim(sce, used)[,c(1:PCNum)]),
                  'BEcounts' = assay(sce, 'BEcounts'),
                  'BENMcounts' =assay(sce, 'BENMcounts'))

  batch = sce@colData$batch
  pheno = colData(sce)

  batch.estimate = kBET::kBET(df = edata, batch = batch, plot = F)

  plot.data = data.frame(class=rep(c('observed', 'expected'),
                                   each=length(batch.estimate$stats$kBET.observed)),
                         data =  c(batch.estimate$stats$kBET.observed,
                                   batch.estimate$stats$kBET.expected))
  g = ggplot(plot.data, mapping = aes(x = class, y = data)) + geom_boxplot() +
    labs(x='Test', y='Rejection rate',title='kBET test results') +
    theme_bw() +
    scale_y_continuous(limits=c(0,1))
  return(g)
}

.BatchEva_kbet = function(sce, used = 'BEPCA', PCNum = 3) {

  if(used == 'BEPCA') {
    emb.bepca = reducedDim(sce, used)
    if(ncol(emb.bepca) < PCNum) {
      PCNum = ncol(emb.bepca)
    }
  }

  if(used == 'PCA') {
    emb.bepca = reducedDim(sce, used)
    if(ncol(emb.bepca) < PCNum) {
      PCNum = ncol(emb.bepca)
    }
  }

  edata = switch (EXPR = used,
                  'counts' = assay(sce, 'counts'),
                  'VGcounts' = log2(assay(altExp(sce),'VGcounts') + 1),
                  'PCA' = t(reducedDim(sce, used)[,c(1:PCNum)]),
                  'BEPCA' = t(reducedDim(sce, used)[,c(1:PCNum)]),
                  'BEcounts' = assay(sce, 'BEcounts'),
                  'BENMcounts' =assay(sce, 'BENMcounts'))

  batch = sce@colData$batch
  pheno = colData(sce)

  batch.estimate = kBET::kBET(df = edata, batch = batch, plot = F)

  plot.data = data.frame(class=rep(c('observed', 'expected'),
                                   each=length(batch.estimate$stats$kBET.observed)),
                         data =  c(batch.estimate$stats$kBET.observed,
                                   batch.estimate$stats$kBET.expected))
  return(plot.data)
}


BatchEva_cmix = function(sce, used = 'PCA', PCNum = 3) {

  if(used == 'BEPCA') {
    emb.bepca = reducedDim(sce, used)
    if(ncol(emb.bepca) < PCNum) {
      PCNum = ncol(emb.bepca)
    }
  }

  if(used == 'PCA') {
    emb.bepca = reducedDim(sce, used)
    if(ncol(emb.bepca) < PCNum) {
      PCNum = ncol(emb.bepca)
    }
  }
  # use batch effect normalized count for batch effect evaluation
  if(used == 'BENMcounts') {
    used = "BEPCA"
    sce <- scater::runPCA(sce, exprs_values="BENMcounts", name ="BEPCA")
    emb.bepca = reducedDim(sce, used)
    if(ncol(emb.bepca) < PCNum) {
      PCNum = ncol(emb.bepca)
    }
  }

  edata = switch (EXPR = used,
                  #'counts' = assay(sce, 'counts'),
                  #'VGcounts' = log2(assay(altExp(sce),'VGcounts') + 1),
                  'PCA' = reducedDim(sce, used)[,c(1:PCNum)],
                  'BEPCA' = reducedDim(sce, used)[,c(1:PCNum)])

  batch = sce@colData$batch
  pheno = colData(sce)

  #suppressPackageStartupMessages(library(CellMixS, quietly = T))
  sce = SingleCellExperiment::SingleCellExperiment(assay = list(counts = assay(sce, 'counts'),
                                                                logcounts = log2(assay(sce, 'counts') + 1)),
                                                   colData = data.frame(batch = batch))

  reducedDim(sce, 'PCA') = edata

  sce = CellMixS::cms(sce, k = 30, group = "batch",
                      dim_red = "PCA",
                      res_name = "unaligned",
                      n_dim = 3, cell_min = 4)

  suppressPackageStartupMessages(library(dplyr, quietly = T))
  #cms_mat = sce$cms.unaligned %>% bind_cols() %>% magrittr::set_colnames(used)
  cms_mat = data.frame(sce$cms.unaligned)
  colnames(cms_mat) = used


  #plot.g = CellMixS::visHist(cms_mat, n_col = 1)
  #CellMixS::visHist(sce$cms.unaligned)
  plot.g = ggplot2::ggplot(data = cms_mat, mapping = aes(x = cms_mat[,1])) +
    ggplot2::geom_histogram(aes(y = ..density..), color = 'deepskyblue', fill = 'deepskyblue', bins=10) +
    ggplot2::geom_density(alpha = .2, color = 'red') +
    ggplot2::theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.line = element_line(colour = "black"),
                   panel.border = element_blank(),
                   #panel.border = element_rect(colour = "black", size=1),
                   panel.background = element_blank()) +
    ggplot2::xlab('cms') +
    ggplot2::ylab('Cell counts') +
    ggplot2::theme(legend.position="bottom",
                   legend.title = element_blank()) +
    ggplot2::theme(axis.text.y = element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.text.x = element_text(angle = 90)) +

    ggplot2::geom_vline(xintercept = mean(cms_mat[,1]), linetype = "longdash",
                        colour = "black") +
    ggplot2::geom_label(aes(y = 2, x = mean(cms_mat[,1]),
                            label = paste0(format(mean(cms_mat[,1]), digits = 2))), vjust = 0, size = 4)
  return(plot.g)
}

.BatchEva_cmix = function(sce, used = 'PCA', PCNum = 3) {

  if(used == 'BEPCA') {
    emb.bepca = reducedDim(sce, used)
    if(ncol(emb.bepca) < PCNum) {
      PCNum = ncol(emb.bepca)
    }
  }

  if(used == 'PCA') {
    emb.bepca = reducedDim(sce, used)
    if(ncol(emb.bepca) < PCNum) {
      PCNum = ncol(emb.bepca)
    }
  }

  edata = switch (EXPR = used,
                  #'counts' = assay(sce, 'counts'),
                  #'VGcounts' = log2(assay(altExp(sce),'VGcounts') + 1),
                  'PCA' = reducedDim(sce, used)[,c(1:PCNum)],
                  'BEPCA' = reducedDim(sce, used)[,c(1:PCNum)])

  batch = sce@colData$batch
  pheno = colData(sce)

  #suppressPackageStartupMessages(library(CellMixS, quietly = T))
  sce = SingleCellExperiment::SingleCellExperiment(assay = list(counts = assay(sce, 'counts'),
                                                                logcounts = log2(assay(sce, 'counts') + 1)),
                                                   colData = data.frame(batch = batch))

  reducedDim(sce, 'PCA') = edata

  sce = CellMixS::cms(sce, k = 30, group = "batch",
                      dim_red = "PCA",
                      res_name = "unaligned",
                      n_dim = 3, cell_min = 4)

  suppressPackageStartupMessages(library(dplyr, quietly = T))
  #cms_mat = sce$cms.unaligned %>% bind_cols() %>% magrittr::set_colnames(used)
  cms_mat = data.frame(sce$cms.unaligned)
  colnames(cms_mat) = 'cms'

  return(cms_mat)
}


BatchEva_slt = function(sce, used = 'BEPCA', PCNum = 3) {
  edata = switch (EXPR = used,
                  'counts' = assay(sce, 'counts'),
                  'VGcounts' = log2(assay(altExp(sce),'VGcounts') + 1),
                  'PCA' = t(reducedDim(sce, used)[,c(1:PCNum)]),
                  'BEPCA' = t(reducedDim(sce, used)[,c(1:PCNum)]),
                  'BEcounts' = assay(sce, 'BEcounts'),
                  'BENMcounts' =assay(sce, 'BENMcounts'))

  batch = sce@colData$batch
  pheno = colData(sce)

  #stats::dist(edata)
  ss = cluster::silhouette(x = c(1:length(unique(batch)))[as.factor(batch)],
                           dist = dist(t(edata)))
  df = data.frame(y = ss[,3],
                  batch = factor(batch, levels = unique(batch)),
                  row.names = rownames(t(edata)))

  #boxp = boxplot(df$y,
  #               col = c('lightblue'),
  #               notch = TRUE,
  #               main = "Silhouette scores")

  df = df[order(df$y, decreasing = T),]
  df = df[order(df$batch, decreasing = T),]
  df$x = c(1:nrow(df))

  sltg = ggplot2::ggplot() +
    ggplot2::geom_col(data = df, mapping = aes(x = x, y = y, fill = batch, color = batch)) +
    #ggplot2::scale_color_manual(values = rainbow(length(unique(batch)))) +
    #ggplot2::scale_fill_manual(values = rainbow(length(unique(batch)))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          #panel.border = element_rect(colour = "black", size=1),
          panel.background = element_blank()) +
    ggplot2::coord_flip() +
    ggplot2::xlab('') +
    ggplot2::ylab('Silhouette scores') +
    ggplot2::theme(legend.position="bottom",
                   legend.title = element_blank()) +
    ggplot2::theme(axis.text.y = element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.text.x = element_text(angle = 90)) +

    ggplot2::geom_hline(yintercept = mean(df$y), linetype = "longdash",
                        colour = "black") +
    ggplot2::geom_label(aes(x = nrow(df)/2, y = mean(df$y),
                            label = paste0(format(mean(df$y), digits = 2))), vjust = 0, size = 4)
  #g
  return(sltg)
}

.BatchEva_slt = function(sce, used = 'BEPCA', PCNum = 3) {
  edata = switch (EXPR = used,
                  'PCA' = t(reducedDim(sce, used)[,c(1:PCNum)]),
                  'BEPCA' = t(reducedDim(sce, used)[,c(1:PCNum)]))

  batch = sce@colData$batch
  pheno = colData(sce)

  ss = cluster::silhouette(x = c(1:length(unique(batch)))[as.factor(batch)],
                           dist = dist(t(edata)))
  df = data.frame(y = ss[,3],
                  batch = factor(batch, levels = unique(batch)),
                  row.names = rownames(t(edata)))

  df = df[order(df$y, decreasing = T),]
  df = df[order(df$batch, decreasing = T),]
  df$x = c(1:nrow(df))

  return(df)
}
