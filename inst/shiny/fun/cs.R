#================#
#Pre-analysis
#================#
Plot_TotalCount = function(CS.data = NULL, used = 'counts', by = NULL, outIndex = getwd(), width = 4, height = 4.2, savePlot = T) {

  edata = assay(CS.data, used)
  batch = CS.data@colData$batch
  pheno = colData(CS.data)

  matrix = edata
  totalReadsForCell = apply(matrix, 2, sum)

  if (is.null(by)) {
    by = colnames(pheno)[1]
    metaVector = pheno[,by]
  } else {
    if(by == 'batch') {
      metaVector = batch
    } else {
      metaVector = pheno[,by]
    }
  }

  library(ggplot2, quietly = T)
  df = data.frame(TotalCountsForCell = totalReadsForCell, Type = metaVector)
  g = ggplot(df, mapping = aes(x = TotalCountsForCell, color = Type)) +
    geom_density() +
    scale_color_manual(values = rainbow(length(unique(metaVector)))) + theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          #panel.border = element_rect(colour = "black", size=1),
          panel.background = element_blank()) +
    xlab('Total read counts') +
    ylab('Density') +
    theme(legend.position="bottom",
          legend.title = element_blank()) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y=element_blank())
  #g

  if(savePlot) {
    if (!(file.exists(paste0(outIndex, '/Figures_ExpFeature/')))){
      dir.create(paste0(outIndex, '/Figures_ExpFeature/'))
    }
    ggsave(filename = paste0(paste0(outIndex, '/Figures_ExpFeature/'), "Plot_TotalCount_With_", used, "_On_", by, ".pdf"), plot = g, width = width, height = height, useDingbats=FALSE)
  }
  return(g)
}


Plot_TotalGene = function(CS.data = NULL, used = 'counts', by = NULL, cutOff = 1, outIndex = getwd(), width = 4, height = 4.2, savePlot = T) {

  edata = assay(CS.data, used)
  batch = CS.data@colData$batch
  pheno = colData(CS.data)

  matrix = edata

  if (is.null(by)) {
    by = colnames(pheno)[1]
    metaVector = pheno[,by]
  } else {
    if(by == 'batch') {
      metaVector = batch
    } else {
      metaVector = pheno[,by]
    }
  }


  NumberOfExpressedGeneForCell = NULL
  for (i in 1:ncol(matrix)) {
    NumberOfExpressedGeneForCell[i] = sum(matrix[,i] >= cutOff)
  }

  library(ggplot2, quietly = T)
  df = data.frame(NumberOfExpressedGeneForCell = NumberOfExpressedGeneForCell, Type = metaVector)
  g = ggplot(df, mapping = aes(x = NumberOfExpressedGeneForCell, color = Type)) +
    geom_density() +
    scale_color_manual(values = rainbow(length(unique(metaVector))))+ theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          #panel.border = element_rect(colour = "black", size=1),
          panel.background = element_blank()) +
    xlab(paste0('Number of expressed genes','\n',"(Count cut-off = ", cutOff,")")) +
    ylab('Density') +
    theme(legend.position="bottom",
          legend.title = element_blank()) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y=element_blank()) +
    ##
    geom_vline(data = df, mapping = aes(xintercept = min(NumberOfExpressedGeneForCell)), colour="darkgray", linetype = 'dashed') +
    geom_text(data = df, mapping = aes(x=min(NumberOfExpressedGeneForCell), label=min(NumberOfExpressedGeneForCell), y= 0), colour="black", angle=90, vjust = 1.2)
    ##
  #g
  if(savePlot) {
    if (!(file.exists(paste0(outIndex, '/Figures_ExpFeature/')))){
      dir.create(paste0(outIndex, '/Figures_ExpFeature/'))
    }
    ggsave(filename = paste0(paste0(outIndex, '/Figures_ExpFeature/'), "Plot_TotalGene_With_", used, "_On_", by, ".pdf"), plot = g, width = width, height = height, useDingbats=FALSE)
  }

  return(g)
}


Plot_MeanCount = function(CS.data = NULL, used = 'counts', by = NULL, outIndex = getwd(), width = 4, height = 4.2, savePlot = T) {

  edata = assay(CS.data, used)
  batch = CS.data@colData$batch
  pheno = colData(CS.data)

  matrix = edata

  if (is.null(by)) {
    by = colnames(pheno)[1]
    metaVector = pheno[,by]
  } else {
    if(by == 'batch') {
      metaVector = batch
    } else {
      metaVector = pheno[,by]
    }
  }


  totalReadsForCell = apply(matrix, 2, mean)

  library(ggplot2, quietly = T)
  df = data.frame(TotalCountsForCell = totalReadsForCell, Type = metaVector)
  g = ggplot(df, mapping = aes(x = TotalCountsForCell, color = Type)) +
    geom_density() +
    scale_color_manual(values = rainbow(length(unique(metaVector)))) + theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          #panel.border = element_rect(colour = "black", size=1),
          panel.background = element_blank()) +
    xlab('Mean read counts') +
    ylab('Density') +
    theme(legend.position="bottom",
          legend.title = element_blank()) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y=element_blank())
  #g
  if(savePlot) {
    if (!(file.exists(paste0(outIndex, '/Figures_ExpFeature/')))){
      dir.create(paste0(outIndex, '/Figures_ExpFeature/'))
    }
    ggsave(filename = paste0(paste0(outIndex, '/Figures_ExpFeature/'), "Plot_MeanCount_With_", used, "_On_", by, ".pdf"), plot = g, width = width, height = height, useDingbats=FALSE)
  }

  return(g)
}


#================#
# QC
#================#
Filter_Matrix = function(CS.data = NULL, matrix = NULL, totalCount = 200, totalGene = 200, geneCapt = 2, MT = 0.1, cutOff = 1) {
  if(!is.null(CS.data)) {
    print("Filtering by S4!")
    edata = assay(CS.data, 'counts')
    batch = CS.data@colData$batch
    pheno = colData(CS.data)

    matrix = as.matrix(edata)
    totalCount.num = apply(matrix, 2, sum)
    totalGene.num  = apply(matrix >= cutOff, 2, sum)
    geneCapt.num   = apply(matrix >= 1, 1, sum)

    edata = edata[geneCapt.num >= geneCapt, ((totalCount.num >= totalCount) & (totalGene.num >= totalGene))]

    #batch = batch[((totalCount.num >= totalCount) & (totalGene.num >= totalGene))]

    pheno = pheno[((totalCount.num >= totalCount) & (totalGene.num >= totalGene)),]

    if(length(grep('^MT-', rownames(edata))) > 0) {
      mt.exp = edata[grep('^MT-', rownames(edata)),]
      mt.exp.pt = apply(mt.exp, 2, sum) / apply(edata, 2, sum)
      edata = edata[,mt.exp.pt > MT]
    }
    if(length(grep('^mt-', rownames(edata))) > 0) {
      mt.exp = edata[grep('^MT-', rownames(edata)),]
      mt.exp.pt = apply(mt.exp, 2, sum) / apply(edata, 2, sum)
      edata = edata[,mt.exp.pt > MT]
    }

    CS.data = SingleCellExperiment(assay = list(counts = edata, logcounts = log2(edata + 1)),
                                   colData = pheno,
                                   metadata = list(study = Study.Name))

    return(CS.data)

  } else {

    matrix = as.matrix(matrix)
    totalCount.num = apply(matrix, 2, sum)
    totalGene.num  = apply(matrix >= cutOff, 2, sum)
    geneCapt.num   = apply(matrix >= 1, 1, sum)

    matrix = matrix[geneCapt.num >= geneCapt, ((totalCount.num >= totalCount) & (totalGene.num >= totalGene))]

    return(matrix)
  }


}


#================#
# Normalization
#================#
Normalize_Matrix = function(CS.data, method = 'lim', used = 'counts', batch = F) {

  CS.data = switch(EXPR = method,
                   'sct'    = Normalize_scVST(   CS.data = CS.data, used = used),
                   'vst'    = Normalize_dqVST(   CS.data = CS.data, used = used),
                   'lib'    = Normalize_libsize( CS.data = CS.data, used = used),
                   'scn'    = Normalize_scNorm(  CS.data = CS.data, used = used),
                   'lim'    = Normalize_linnorm( CS.data = CS.data, used = used),
                   'LogNormalize' = Normalize_LogNormalize(CS.data = CS.data, used = used))
  #CS.data = Normalize_scVST(CS.data)
  #CS.data = Normalize_dqVST(CS.data)
  #CS.data = Normalize_libsize(CS.data)
  #CS.data = Normalize_scNorm(CS.data)
  #CS.data = Normalize_linnorm(CS.data)

  return(CS.data)
}

.Normalize_Matrix = function(CS.data, method = 'lim', used = 'counts') {

  CS.data = switch(EXPR = method,
                   'sct'    = .Normalize_scVST(   CS.data = CS.data, used = used),
                   'vst'    = .Normalize_dqVST(   CS.data = CS.data, used = used),
                   'lib'    = .Normalize_libsize( CS.data = CS.data, used = used),
                   'scn'    = .Normalize_scNorm(  CS.data = CS.data, used = used),
                   'lim'    = .Normalize_linnorm( CS.data = CS.data, used = used))

  return(CS.data)
}

Normalize_LogNormalize= function(CS.data, used = 'counts') {
  CS.data <- scran::computeSumFactors(CS.data)
  CS.data <- scater::logNormCounts(CS.data,name="NMcounts")
  return(CS.data)
}
Normalize_scVST = function(CS.data, used = 'counts') {
  #View(assay(CS.data, used))

  library(SingleCellExperiment, quietly = T)
  edata = assay(CS.data, used)
  #batch = CS.data@colData$batch
  #pheno = colData(CS.data)

  library(sctransform, quietly = T)
  edata_vst = sctransform::vst(as.matrix(edata))$y

  edata_vst = as.data.frame(edata_vst)[rownames(edata),]
  row.names(edata_vst) = rownames(edata)
  edata_vst[is.na(edata_vst)] = 0

  assay(CS.data, 'NMcounts') = as.matrix(edata_vst)
  return(CS.data)
}

.Normalize_scVST = function(CS.data, used = 'counts') {
  #View(assay(CS.data, used))

  library(SingleCellExperiment, quietly = T)
  edata = assay(CS.data, used)
  batch = CS.data@colData$batch
  #pheno = colData(CS.data)

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

  assay(CS.data, 'NMcounts') = as.matrix(edata_vst)
  return(CS.data)
}


Normalize_dqVST = function(CS.data, used = 'counts') {
  #View(assay(CS.data, used))

  library(SingleCellExperiment, quietly = T)
  edata = assay(CS.data, used)
  #batch = CS.data@colData$batch
  #pheno = colData(CS.data)

  library(DESeq2, quietly = T)
  edata_vst = DESeq2::vst(as.matrix(edata))

  assay(CS.data, 'NMcounts') = as.matrix(edata_vst)
  return(CS.data)
}

.Normalize_dqVST = function(CS.data, used = 'counts') {
  #View(assay(CS.data, used))

  library(SingleCellExperiment, quietly = T)
  edata = assay(CS.data, used)
  batch = CS.data@colData$batch
  #pheno = colData(CS.data)

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

  assay(CS.data, 'NMcounts') = as.matrix(edata_vst)
  return(CS.data)
}


Normalize_libsize = function(CS.data, used = 'counts') {
  #View(assay(CS.data, used))
  edata = assay(CS.data, used)
  #batch = CS.data@colData$batch
  #pheno = colData(CS.data)

  library(DESeq2, quietly = T)
  lib.size = estimateSizeFactorsForMatrix(edata)
  edata_ed = t(t(edata)/lib.size)

  assay(CS.data, 'NMcounts') = as.matrix(edata_ed)
  return(CS.data)
}

.Normalize_libsize = function(CS.data, used = 'counts') {
  #View(assay(CS.data, used))
  edata = assay(CS.data, used)
  batch = CS.data@colData$batch
  #pheno = colData(CS.data)

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

  assay(CS.data, 'NMcounts') = as.matrix(edata_ed)
  return(CS.data)
}


Normalize_scNorm = function(CS.data, used = 'counts') {
  edata = assay(CS.data, used)
  batch = CS.data@colData$batch
  #pheno = colData(CS.data)

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

  assay(CS.data, 'NMcounts') = as.matrix(edata_scnorm)
  return(CS.data)
}

.Normalize_scNorm = function(CS.data, used = 'counts') {
  edata = assay(CS.data, used)
  batch = CS.data@colData$batch
  #pheno = colData(CS.data)

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

  assay(CS.data, 'NMcounts') = as.matrix(edata_scnorm)
  return(CS.data)
}


Normalize_linnorm = function(CS.data, used = 'counts') {
  #View(assay(CS.data, used))
  edata = assay(CS.data, used)
  #batch = CS.data@colData$batch
  #pheno = colData(CS.data)

  library(Linnorm, quietly = T)
  edata_lin = Linnorm::Linnorm.Norm(datamatrix = edata)


  assay(CS.data, 'NMcounts') = as.matrix(edata_lin)
  return(CS.data)
}

.Normalize_linnorm = function(CS.data, used = 'counts') {
  #View(assay(CS.data, used))
  edata = assay(CS.data, used)
  batch = CS.data@colData$batch
  #pheno = colData(CS.data)

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

  assay(CS.data, 'NMcounts') = as.matrix(edata_lin)
  return(CS.data)
}


Normalize_ruvseq = function(CS.data, used = 'counts', spkins) {
  edata = assay(CS.data, used)
  #batch = CS.data@colData$batch

  spkins = spkins[spkins %in% rownames(edata)]

  library(RUVSeq, quietly = T)
  edata_ruv = RUVSeq::RUVg(x = as.matrix(edata), cIdx = spkins, k = 2)

  edata_ruv = edata_ruv$normalizedCounts

  assay(CS.data, 'NMcounts') = edata_ruv

  return(CS.data)
}


#================#
# Batch effect correction -- on expression matrix
#================#

BatchEva = function(CS.data, method = 'cmx', used = NULL, PCNum = 3) {

  # Condition to select "used"
  if(!is.null(used)){
    used=used
  }
  else if ("BENMcounts" %in% assayNames(CS.data)) {
    used= "BENMcounts"
  }
  else{
    used= "BEPCA"
  }



  evaPlot = switch(EXPR = method,
               'cmx' = BatchEva_cmix( CS.data = CS.data, used = used, PCNum = PCNum),
               'kbet' = BatchEva_kbet( CS.data = CS.data, used = used, PCNum = PCNum),
               'slt'  = BatchEva_slt(  CS.data = CS.data, used = used, PCNum = PCNum))

  return(evaPlot)
}


BatchEva_kbet = function(CS.data, used = 'BEPCA', PCNum = 3) {

  if(used == 'BEPCA') {
    emb.bepca = reducedDim(CS.data, used)
    if(ncol(emb.bepca) < PCNum) {
      PCNum = ncol(emb.bepca)
    }
  }

  if(used == 'PCA') {
    emb.bepca = reducedDim(CS.data, used)
    if(ncol(emb.bepca) < PCNum) {
      PCNum = ncol(emb.bepca)
    }
  }

  edata = switch (EXPR = used,
                  'counts' = assay(CS.data, 'counts'),
                  'VGcounts' = log2(assay(altExp(CS.data),'VGcounts') + 1),
                  'PCA' = t(reducedDim(CS.data, used)[,c(1:PCNum)]),
                  'BEPCA' = t(reducedDim(CS.data, used)[,c(1:PCNum)]),
                  'BEcounts' = assay(CS.data, 'BEcounts'),
                  'BENMcounts' =assay(CS.data, 'BENMcounts'))

  batch = CS.data@colData$batch
  pheno = colData(CS.data)

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

.BatchEva_kbet = function(CS.data, used = 'BEPCA', PCNum = 3) {

  if(used == 'BEPCA') {
    emb.bepca = reducedDim(CS.data, used)
    if(ncol(emb.bepca) < PCNum) {
      PCNum = ncol(emb.bepca)
    }
  }

  if(used == 'PCA') {
    emb.bepca = reducedDim(CS.data, used)
    if(ncol(emb.bepca) < PCNum) {
      PCNum = ncol(emb.bepca)
    }
  }

  edata = switch (EXPR = used,
                  'counts' = assay(CS.data, 'counts'),
                  'VGcounts' = log2(assay(altExp(CS.data),'VGcounts') + 1),
                  'PCA' = t(reducedDim(CS.data, used)[,c(1:PCNum)]),
                  'BEPCA' = t(reducedDim(CS.data, used)[,c(1:PCNum)]),
                  'BEcounts' = assay(CS.data, 'BEcounts'),
                  'BENMcounts' =assay(CS.data, 'BENMcounts'))

  batch = CS.data@colData$batch
  pheno = colData(CS.data)

  batch.estimate = kBET::kBET(df = edata, batch = batch, plot = F)

  plot.data = data.frame(class=rep(c('observed', 'expected'),
                                   each=length(batch.estimate$stats$kBET.observed)),
                         data =  c(batch.estimate$stats$kBET.observed,
                                   batch.estimate$stats$kBET.expected))
  return(plot.data)
}


BatchEva_cmix = function(CS.data, used = 'PCA', PCNum = 3) {

  if(used == 'BEPCA') {
    emb.bepca = reducedDim(CS.data, used)
    if(ncol(emb.bepca) < PCNum) {
      PCNum = ncol(emb.bepca)
    }
  }

  if(used == 'PCA') {
    emb.bepca = reducedDim(CS.data, used)
    if(ncol(emb.bepca) < PCNum) {
      PCNum = ncol(emb.bepca)
    }
  }
  # use batch effect normalized count for batch effect evaluation
  if(used == 'BENMcounts') {
    used = "BEPCA"
    CS.data <- scater::runPCA(CS.data, exprs_values="BENMcounts", name ="BEPCA")
    emb.bepca = reducedDim(CS.data, used)
    if(ncol(emb.bepca) < PCNum) {
      PCNum = ncol(emb.bepca)
    }
  }

  edata = switch (EXPR = used,
                  #'counts' = assay(CS.data, 'counts'),
                  #'VGcounts' = log2(assay(altExp(CS.data),'VGcounts') + 1),
                  'PCA' = reducedDim(CS.data, used)[,c(1:PCNum)],
                  'BEPCA' = reducedDim(CS.data, used)[,c(1:PCNum)])

  batch = CS.data@colData$batch
  pheno = colData(CS.data)

  #suppressPackageStartupMessages(library(CellMixS, quietly = T))
  sce = SingleCellExperiment::SingleCellExperiment(assay = list(counts = assay(CS.data, 'counts'),
                                                                logcounts = log2(assay(CS.data, 'counts') + 1)),
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

.BatchEva_cmix = function(CS.data, used = 'PCA', PCNum = 3) {

  if(used == 'BEPCA') {
    emb.bepca = reducedDim(CS.data, used)
    if(ncol(emb.bepca) < PCNum) {
      PCNum = ncol(emb.bepca)
    }
  }

  if(used == 'PCA') {
    emb.bepca = reducedDim(CS.data, used)
    if(ncol(emb.bepca) < PCNum) {
      PCNum = ncol(emb.bepca)
    }
  }

  edata = switch (EXPR = used,
                  #'counts' = assay(CS.data, 'counts'),
                  #'VGcounts' = log2(assay(altExp(CS.data),'VGcounts') + 1),
                  'PCA' = reducedDim(CS.data, used)[,c(1:PCNum)],
                  'BEPCA' = reducedDim(CS.data, used)[,c(1:PCNum)])

  batch = CS.data@colData$batch
  pheno = colData(CS.data)

  #suppressPackageStartupMessages(library(CellMixS, quietly = T))
  sce = SingleCellExperiment::SingleCellExperiment(assay = list(counts = assay(CS.data, 'counts'),
                                                                logcounts = log2(assay(CS.data, 'counts') + 1)),
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


BatchEva_slt = function(CS.data, used = 'BEPCA', PCNum = 3) {
  edata = switch (EXPR = used,
                  'counts' = assay(CS.data, 'counts'),
                  'VGcounts' = log2(assay(altExp(CS.data),'VGcounts') + 1),
                  'PCA' = t(reducedDim(CS.data, used)[,c(1:PCNum)]),
                  'BEPCA' = t(reducedDim(CS.data, used)[,c(1:PCNum)]),
                  'BEcounts' = assay(CS.data, 'BEcounts'),
                  'BENMcounts' =assay(CS.data, 'BENMcounts'))

  batch = CS.data@colData$batch
  pheno = colData(CS.data)

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

.BatchEva_slt = function(CS.data, used = 'BEPCA', PCNum = 3) {
  edata = switch (EXPR = used,
                  'PCA' = t(reducedDim(CS.data, used)[,c(1:PCNum)]),
                  'BEPCA' = t(reducedDim(CS.data, used)[,c(1:PCNum)]))

  batch = CS.data@colData$batch
  pheno = colData(CS.data)

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




#================#
# Batch effect correction -- on expression matrix
#================#

BatchEffect_Matrix = function(CS.data, method = 'cbtseq', used = 'counts') {

  CS.data = switch(EXPR = method,
                   'bas'      = RmBatch_denoise(   CS.data = CS.data, used = used),
                   'lima'     = RmBatch_limma(     CS.data = CS.data, used = used),
                   'cbt'      = RmBatch_combat(    CS.data = CS.data, used = used),
                   'cbtseq'   = RmBatch_combatseq( CS.data = CS.data, used = used),
                   'mnn'      = RmBatch_mnn(       CS.data = CS.data, used = used),
                   'sca'      = RmBatch_scano(     CS.data = CS.data, used = used))

  return(CS.data)
}


RmBatch_denoise = function(CS.data, used = 'counts') {
  cat(crayon::red('===Running denoise!===\n'))
  #View(SummarizedExperiment::assay(CS.data, used))
  #edata = SummarizedExperiment::assay(CS.data, used)
  edata = switch(EXPR = used,
                 'counts' = SummarizedExperiment::assay(CS.data, used),
                 'NMcounts' = SummarizedExperiment::assay(CS.data, used),
                 'VGcounts' = SummarizedExperiment::assay(SingleCellExperiment::altExp(CS.data, used), used))

  batch = SingleCellExperiment::colData(CS.data)[,'batch']
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
    SummarizedExperiment::assay(CS.data, 'BEcounts') = edata_Denoised
  }
  if(used == 'VGcounts') {
    SingleCellExperiment::altExp(CS.data, 'BEVGcounts') = SingleCellExperiment::SingleCellExperiment(assays = list(BEVGcounts = edata_Denoised))
  }

  return(CS.data)
}


RmBatch_limma = function(CS.data, used = 'counts') {
  cat(crayon::red('===Running limma!===\n'))
  edata = switch(EXPR = used,
                 'counts' = log2(SummarizedExperiment::assay(CS.data, used)+1),
                 'NMcounts' = SummarizedExperiment::assay(CS.data, used),
                 'VGcounts' = log2(SummarizedExperiment::assay(SingleCellExperiment::altExp(CS.data), used)+1))
  edata_rbe = limma::removeBatchEffect(x = edata, batch = SingleCellExperiment::colData(CS.data)[,'batch'])
  edata_rbe[edata_rbe < 0] = 0
  #SummarizedExperiment::assay(CS.data, 'BEcounts') = edata_rbe
  if(used == 'counts') {
    SummarizedExperiment::assay(CS.data, 'BEcounts') = edata_rbe
  }
  if(used == 'VGcounts') {
    SingleCellExperiment::altExp(CS.data, 'BEVGcounts') = SingleCellExperiment::SingleCellExperiment(assays = list(BEVGcounts = edata_rbe))
  }
  if(used== 'NMcounts' ){
    SummarizedExperiment::assay(CS.data, 'BENMcounts') = edata_rbe
  }
  return(CS.data)

}


RmBatch_combat = function(CS.data, used = 'counts') {
  cat(crayon::red('===Running combat!===\n'))
  edata = switch(EXPR = used,
                 'counts' = log2(SummarizedExperiment::assay(CS.data, used)+1),
                 'NMcounts' = SummarizedExperiment::assay(CS.data, used),
                 'VGcounts' = log2(SummarizedExperiment::assay(SingleCellExperiment::altExp(CS.data), used)+1))
  batch = SingleCellExperiment::colData(CS.data)[,'batch']

  var.data = apply(edata, 1, var)
  if(sum(var.data == 0) > 0) {
    rm.edata = edata[which(var.data == 0),,drop = F]

    edata = edata[-which(var.data == 0 ),]
    modcombat = stats::model.matrix(~1, data = SummarizedExperiment::colData(CS.data))
    edata_combat = sva::ComBat(dat = as.matrix(edata), batch = batch, mod = modcombat, par.prior = TRUE)

    edata_combat = rbind(edata_combat, rm.edata)
    edata_combat = edata_combat[rownames(SummarizedExperiment::assay(CS.data)),]
  } else {
    modcombat = stats::model.matrix(~1, data = SummarizedExperiment::colData(CS.data))
    edata_combat = sva::ComBat(dat = as.matrix(edata), batch = batch, mod = modcombat, par.prior=TRUE)
  }

  if(used == 'counts') {
    SummarizedExperiment::assay(CS.data, 'BEcounts') = edata_combat
  }
  if(used == 'VGcounts') {
    SingleCellExperiment::altExp(CS.data, 'BEVGcounts') = SingleCellExperiment::SingleCellExperiment(assays = list(BEVGcounts = edata_combat))
  }
  if(used== 'NMcounts' ){
    SummarizedExperiment::assay(CS.data, 'BENMcounts') = edata_combat
  }
  return(CS.data)
}

RmBatch_combatseq = function(CS.data, used = 'counts') {
  cat(crayon::red('===Running combatseq!===\n'))
  edata = switch(EXPR = used,
                 'counts' = SummarizedExperiment::assay(CS.data, used),
                 'NMcounts' = SummarizedExperiment::assay(CS.data, used),
                 'VGcounts' = SummarizedExperiment::assay(SingleCellExperiment::altExp(CS.data,used), used))
  batch = SingleCellExperiment::colData(CS.data)[,'batch']
  edata_combatseq = sva::ComBat_seq(counts = edata, batch = batch)
  if(used == 'counts') {
    SummarizedExperiment::assay(CS.data, 'BEcounts') = log2(edata_combatseq + 1)
  }
  if(used == 'VGcounts') {
    SingleCellExperiment::altExp(CS.data, 'BEVGcounts') = SingleCellExperiment::SingleCellExperiment(assays = list(BEVGcounts = log2(edata_combatseq + 1)))
  }
  if(used== 'NMcounts' ){
    SummarizedExperiment::assay(CS.data, 'BENMcounts') = edata_combatseq
  }
  return(CS.data)
}

#RmBatch_mnn = function(CS.data, used = 'counts') {
#  edata = SummarizedExperiment::assay(CS.data, used)
#  batch = SingleCellExperiment::colData(CS.data)[,'batch']
#
#  combined = SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(edata)))
#  combined = batchelor::multiBatchNorm(combined, batch = batch)
#  combined$batch = batch
#
#  f.out2 = batchelor::fastMNN(combined, batch = combined$batch)
#
#  edata_mnn = as.matrix(SummarizedExperiment::assay(f.out2))
#
#  SummarizedExperiment::assay(CS.data, 'BEcounts') = edata_mnn
#  return(CS.data)
#}

RmBatch_scano = function(CS.data, used = 'counts') {
  cat(crayon::red('===Running scanoramaCT!===\n'))
  edata = switch(EXPR = used,
                 'counts' = log2(SummarizedExperiment::assay(CS.data, used)+1),
                 'NMcounts' = SummarizedExperiment::assay(CS.data, used),
                 'VGcounts' = log2(SummarizedExperiment::assay(SingleCellExperiment::altExp(CS.data, used), used)+1))
  batch = SingleCellExperiment::colData(CS.data)[,'batch']

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
  integrated.corrected.data = scanoramaCT$correct(processedMats,
                                                  genes_list,
                                                  return_dimred = TRUE,
                                                  return_dense = T)
  #countsNormCorrected = unlist(lapply(gcalign, function(x) x[,1]))
  matCorrected = t(do.call(rbind, integrated.corrected.data[[2]]))
  rownames(matCorrected) = integrated.corrected.data[[3]]
  m = do.call(rbind, integrated.corrected.data[[1]])
  colnames(matCorrected) = cell_list

  matCorrected = matCorrected[rownames(edata), colnames(edata)]
  if(used == 'counts') {
    SummarizedExperiment::assay(CS.data, 'BEcounts') = matCorrected
  }
  if(used == 'VGcounts') {
    SingleCellExperiment::altExp(CS.data, 'BEVGcounts') = SingleCellExperiment::SingleCellExperiment(assays = list(BEVGcounts = matCorrected))
  }
  if(used== 'NMcounts' ){
    SummarizedExperiment::assay(CS.data, 'BENMcounts') = matCorrected
  }

  return(CS.data)
}


#================#
# Batch effect correction -- on dimension reduction
#================#
BatchEffect_DimReduction = function(CS.data, method = 'beer', used = 'counts') {
  CS.data = switch(EXPR = method,
                   'beer'  = RmBatch_beer( CS.data = CS.data, used = used),
                   'cca'   = RmBatch_cca(  CS.data = CS.data, used = used),
                   'dmnn'  = RmBatch_dmnn( CS.data = CS.data, used = used),
                   'hmy'   = RmBatch_hmy(  CS.data = CS.data, used = used),
                   'lig'   = RmBatch_lig(  CS.data = CS.data, used = used))

  return(CS.data)

}


RmBatch_beer = function(CS.data, used = 'counts', showBy = 'PCA') {

  #BeerPath = readline("Path of BEER source file:")
  BeerPath = 'fun/BEER.R'
  #BeerPath = gsub(pattern = '\"', replacement = '', x = BeerPath)
  #BeerPath = gsub(pattern = '\'', replacement = '', x = BeerPath)

  source(BeerPath)
  #source("~/T1/BE/BEER/BEER-0.1.7/BEER.R")

  edata = assay(CS.data, used)
  batch = CS.data@colData$batch

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
  reducedDim(CS.data, 'BEPCA') = pcs_beer
  #}
  return(CS.data)
}


RmBatch_cca = function(CS.data, used = 'counts', showby = 'PCA') {
  edata = assay(CS.data, used)
  batch = CS.data@colData$batch

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

  reducedDim(CS.data, 'BEPCA') = as.matrix(immune.combined@reductions$pca@cell.embeddings)

  return(CS.data)
}


RmBatch_dmnn = function(CS.data, used = 'counts', showby = 'PCA') {
  edata = assay(CS.data, used)
  batch = CS.data@colData$batch

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

  reducedDim(CS.data, 'BEPCA') = reducedDim(f.out2, "corrected")
  #plot(reducedDim(f.out2, "corrected")[,c(1,2)], col = rainbow(2)[as.factor(groupVector)])
  return(CS.data)
}


RmBatch_hmy = function(CS.data, used = 'counts', showby = 'PCA') {
  edata = assay(CS.data, used)
  batch = CS.data@colData$batch

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
  #plot(reducedDim(CS.data, 'PCA')[,c(1,2)], col = as.factor(meta.table$Patients))

  harmony_embeddings = harmony::HarmonyMatrix(data_mat = reducedDim(combined, 'PCA')[,c(1:5)],
                                              meta_data = data.frame(batch = batch), vars_use = "batch",
                                              do_pca=FALSE)

  #plot(harmony_embeddings[,c(1,2)], col = as.factor(meta.table$Patients))

  #harmony_embeddings2 <- harmony::HarmonyMatrix(data_mat = reducedDim(CS.data, 'PCA')[,c(1:5)],
  #                                             meta_data = data.frame(batch = batch), vars_use = "batch",
  #                                             do_pca=FALSE)
  #plot(harmony_embeddings2[,c(1,2)], col = as.factor(meta.table$Patients))

  reducedDim(CS.data, 'BEPCA') = harmony_embeddings
  return(CS.data)
}


RmBatch_lig = function(CS.data, used = 'counts', showby = 'PCA') {
  edata = assay(CS.data, used)
  batch = CS.data@colData$batch

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

  reducedDim(CS.data, 'BEPCA') = data.use.pca
  return(CS.data)
}


#================#
# Combinations of normalization and batch correction
#================#



#================#
# Dimension reduction
#================#
PCNumSearch = function(CS.data, runWith = 'VGcounts') {
  if (runWith == 'VGcounts') {
    matrix.PJ = assay(altExp(CS.data), runWith)
    matrix.PJ = log2(matrix.PJ + 1)
  }

  if (runWith == 'PCA') {
    matrix.PJ = t(reducedDim(CS.data, runWith))
  }

  if (runWith == 'BEPCA') {
    matrix.PJ = t(reducedDim(CS.data, runWith))
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


#================#
# Dimension reduction
#================#
Plot_DimReduction = function(CS.data, runWith = 'tSNE', colorBy = 'default', dim.se = c(1:2), showDensity = TRUE) {
  if (colorBy == 'default') {
    colorBy = colnames(colData(CS.data))[1]
  }

  #CS.data = DimReduction_pca()

  #plot
  df.plot = data.frame(reducedDim(CS.data, runWith)[,dim.se],
                       colorBy = as.factor(colData(CS.data)[,colorBy]),
                       color = as.character(colData(CS.data)[,colorBy]),

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


DimReduction_pca = function(CS.data, runWith = 'VGcounts') {
  matrix.PJ = assay(altExp(CS.data), 'VGcounts')
  matrix.PJ = log2(matrix.PJ+1)

  if(runWith == 'NMcounts') {
    matrix.PJ = assay(CS.data, 'NMcounts')[rownames(assay(altExp(CS.data), 'VGcounts')),]
  }
  if(runWith == 'BEcounts') {
    matrix.PJ = assay(CS.data, 'BEcounts')[rownames(assay(altExp(CS.data), 'VGcounts')),]
  }

  #dim(matrix.PJ)
  emb.pca = prcomp(t(matrix.PJ))

  reducedDim(CS.data, 'PCA') = emb.pca$x[,c(1:50)]

  return(CS.data)
}


DimReduction_ica = function(CS.data, runWith = 'VGcounts') {
  matrix.PJ = assay(altExp(CS.data), 'VGcounts')
  matrix.PJ = log2(matrix.PJ+1)

  if('NMcounts' %in% names(CS.data@assays)) {
    matrix.PJ = assay(CS.data, 'NMcounts')[rownames(assay(altExp(CS.data), 'VGcounts')),]
  }
  if('BEcounts' %in% names(CS.data@assays)) {
    matrix.PJ = assay(CS.data, 'BEcounts')[rownames(assay(altExp(CS.data), 'VGcounts')),]
  }

  library(fastICA, quietly = T)
  emb.ica = fastICA(X = matrix.PJ, n.comp = 50)#nrow(matrix.PJ))

  reducedDim(CS.data, 'PCA') = emb.ica$A

  return(CS.data)
}


DimReduction_tsne = function(CS.data, runWith = 'PCA', PCNum = 20, perplexity = 30, parallelRun = TRUE) {

  matrix.PJ = switch (EXPR = runWith,
                      'VGcounts' = log2(assay(altExp(CS.data),'VGcounts') + 1),
                      'PCA' = t(reducedDim(CS.data, runWith)[,c(1:PCNum)]),
                      'BEPCA' = t(reducedDim(CS.data, runWith)[,c(1:PCNum)]))

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

  reducedDim(CS.data, 'tSNE') = emb.tsne

  return(CS.data)
}


DimReduction_umap = function(CS.data, runWith = 'PCA', PCNum = 20) {

  matrix.PJ = switch (EXPR = runWith,
                      'VGcounts' = log2(assay(altExp(CS.data),'VGcounts') + 1),
                      'PCA' = t(reducedDim(CS.data, runWith)[,c(1:PCNum)]),
                      'BEPCA' = t(reducedDim(CS.data, runWith)[,c(1:PCNum)]))

  library(umap, quietly = T)
  set.seed(15555)
  emb.umap = umap(t(as.matrix(matrix.PJ)))

  emb.umap = emb.umap$layout
  colnames(emb.umap) = paste0('Umap_', c(1:ncol(emb.umap)))
  reducedDim(CS.data, 'Umap') = emb.umap

  return(CS.data)
}


DimReduction_dfmap = function(CS.data, runWith = 'PCA', PCNum = 20) {

  matrix.PJ = switch (EXPR = runWith,
                      'VGcounts' = log2(assay(altExp(CS.data),'VGcounts') + 1),
                      'PCA' = t(reducedDim(CS.data, runWith)[,c(1:PCNum)]),
                      'BEPCA' = t(reducedDim(CS.data, runWith)[,c(1:PCNum)]))

  library(destiny, quietly = T)
  set.seed(15555)
  emb.dm = DiffusionMap(t(matrix.PJ))
  emb.dm = as.data.frame(emb.dm@eigenvectors)
  row.names(emb.dm) = colnames(matrix.PJ)

  reducedDim(CS.data, 'DFmap') = emb.dm

  return(CS.data)
}


DimReduction_ddtree = function(CS.data, runWith = 'PCA', PCNum = 20) {
  matrix.PJ = switch (EXPR = runWith,
                      'VGcounts' = log2(assay(altExp(CS.data),'VGcounts') + 1),
                      'PCA' = t(reducedDim(CS.data, runWith)[,c(1:PCNum)]),
                      'BEPCA' = t(reducedDim(CS.data, runWith)[,c(1:PCNum)]))

  library(DDRTree)
  DDRTree_res = DDRTree::DDRTree(X = matrix.PJ, dimensions = 2)
  Z <- DDRTree_res$Z #obatain matrix
  Y <- DDRTree_res$Y
  stree <- DDRTree_res$stree
  library(ggplot2)
  d <- data.frame(t(Z))
  row.names(d) = colnames(x)
  out_coordinate <- paste(outPath,".txt",sep = "")
  write.table(file = out_coordinate,x = d,sep = "\t")
  if (metaNum == 1) {
    g <- ggplot(data = d, mapping = aes(x=X1, y=X2, color=as.factor(metaData[,1])))
    g <- g + geom_point()
    g <- g + theme_bw()
    g <- g + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.line = element_line(colour = "black"))
  } else {
    g <- ggplot(data = d, mapping = aes(x=X1, y=X2, color=as.factor(metaData[,1]), shape = as.factor(metaData[,2])))
    g <- g + geom_point()
    g <- g + theme_bw()
    g <- g + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.line = element_line(colour = "black"))
  }
  #plot(Y[1, ], Y[2, ], col = 'blue', pch = 17) #center of the Z
  #title(main="DDRTree smooth principal curves", col.main="red", font.main=4)
  ggsave(plot = g,height = height, width = width, dpi = 600, filename = outPath, useDingbats=FALSE)

}


#================#
# Most variable genes
#================#
SelectVariableGene = function(CS.data, method = 'mvg', used = 'counts', ngene = 1000, batch = F) {

  if(!batch) {CS.data = switch(EXPR = method,

                               'seu'   = VariableGene_seu(   CS.data = CS.data, used = used, ngene = ngene),
                               'mvg'   = VariableGene_mvg(   CS.data = CS.data, used = used, ngene = ngene),
                               'hvg'   = VariableGene_hvg(   CS.data = CS.data, used = used, ngene = ngene))
  } else {
    #print("xxxxx")
    CS.data = VariableGene_hvg(   CS.data = CS.data, used = used, ngene = ngene, batch = T)
  }

  return(CS.data)
}





VariableGene_seu = function(CS.data, used = 'counts', ngene = 1000) {
  library(Seurat, quietly = T)
  edata = assay(CS.data, used)

  ctrl = CreateSeuratObject(counts = edata, min.cells = 0, min.features = 0)
  ctrl = NormalizeData(ctrl)
  ctrl = FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = ngene)

  edata.od = edata[match(VariableFeatures(ctrl), rownames(ctrl@assays$RNA@counts)),]

  altExp(CS.data, 'VGcounts') = SingleCellExperiment(assay = list(VGcounts = edata[rownames(edata.od),]))

  return(CS.data)
}


VariableGene_mvg = function(CS.data, used = 'counts', ngene = 1000) {
  edata = assay(CS.data, used)

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

  altExp(CS.data, 'VGcounts') = SingleCellExperiment(assay = list(VGcounts = edata.od[1:ngene,]))
  return(CS.data)
}


VariableGene_hvg = function(CS.data, used = 'counts', ngene = 1000, batch = F) {


    edata = assay(CS.data, 'counts')
    edata = log2(edata+1)



  if(batch) {
    batch = CS.data@colData$batch

    dec = scran::modelGeneVar(edata, block = batch)
    top.hvgs2 <- scran::getTopHVGs(dec, n=ngene)

    edata.od = edata[top.hvgs2,]

    altExp(CS.data, 'VGcounts') = SingleCellExperiment(assay = list(VGcounts = edata.od))
  } else {
    dec = scran::modelGeneVar(edata)
    top.hvgs2 <- scran::getTopHVGs(dec, n=ngene)

    edata.od = edata[top.hvgs2,]

    altExp(CS.data, 'VGcounts') = SingleCellExperiment(assay = list(VGcounts = edata.od))
  }

  return(CS.data)
}


#================#
# Clustering
#================#
ClusterNumSearch = function(CS.data, runWith = 'VGcounts', PCNum = 10) {

  matrix.run = switch(EXPR       = runWith,
                      'VGcounts' = assay(altExp(CS.data),'VGcounts'),
                      'PCA'      = reducedDim(x = CS.data, type = 'PCA')[,1:PCNum],
                      'ICA'      = reducedDim(x = CS.data, type = 'ICA')[,1:PCNum],
                      'tSNE'     = reducedDim(x = CS.data, type = 'tSNE')[,1:2],
                      'Umap'     = reducedDim(x = CS.data, type = 'Umap')[,1:2])

  #library(factoextra, quietly = T)
  #Partitioning methods, such as k-means clustering require the users to specify the number of clusters to be generated.
  res.cluster.num = factoextra::fviz_nbclust(x = matrix.run, FUNcluster = kmeans, method = 'silhouette', k.max = 15)

  plot(res.cluster.num)
}


ClusterFind_kmean = function(CS.data, ClusterNum, runWith = 'PCA', PCNum = 10){

  matrix.run = switch(EXPR       = runWith,
                      'VGcounts' = assay(altExp(CS.data),'VGcounts'),
                      'PCA'      = reducedDim(x = CS.data, type = 'PCA')[,1:PCNum],
                      'ICA'      = reducedDim(x = CS.data, type = 'ICA')[,1:PCNum],
                      'tSNE'     = reducedDim(x = CS.data, type = 'tSNE')[,1:2],
                      'Umap'     = reducedDim(x = CS.data, type = 'Umap')[,1:2])

  set.seed(0)
  res.km = kmeans(x = matrix.run, centers = ClusterNum)
  res.km.cluster = res.km$cluster
  #plot(c(1:nrow(matrix.run)), match(rownames(matrix.run), names(res.km.cluster)))

  colData(CS.data)$cluster = as.character(res.km.cluster)
  #colLabels(CS.data)<-as.factor(res.km.cluster)
  return(CS.data)
}


ClusterFind_dbscan = function(CS.data, runWith = 'PCA', PCNum = 10){

  matrix.run = switch(EXPR       = runWith,
                      'VGcounts' = assay(altExp(CS.data),'VGcounts'),
                      'PCA'      = reducedDim(x = CS.data, type = 'PCA')[,1:PCNum],
                      'ICA'      = reducedDim(x = CS.data, type = 'ICA')[,1:PCNum],
                      'tSNE'     = reducedDim(x = CS.data, type = 'tSNE')[,1:2],
                      'Umap'     = reducedDim(x = CS.data, type = 'Umap')[,1:2])

  #library(dbscan, quietly = T)
  set.seed(0)
  res.db = dbscan::dbscan(x = matrix.run, eps = 0.4, minPts = 5)
  res.db.cluster = res.db$cluster

  colData(CS.data)$cluster = as.character(res.db.cluster)
  #colLabels(CS.data)<- as.factor(res.db.cluster)
  return(CS.data)
}


ClusterFind_mclust = function(CS.data, ClusterNum = 3, runWith = 'PCA', PCNum = 10){

  matrix.run = switch(EXPR       = runWith,
                      'VGcounts' = assay(altExp(CS.data),'VGcounts'),
                      'PCA'      = reducedDim(x = CS.data, type = 'PCA')[,1:PCNum],
                      'ICA'      = reducedDim(x = CS.data, type = 'ICA')[,1:PCNum],
                      'tSNE'     = reducedDim(x = CS.data, type = 'tSNE')[,1:2],
                      'Umap'     = reducedDim(x = CS.data, type = 'Umap')[,1:2])

  suppressPackageStartupMessages(library(mclust, quietly = T))
  res.mc = mclust::Mclust(data = matrix.run, G =  c(1:ClusterNum))
  res.mc.cluster = res.mc$classification
  #plot(c(1:nrow(matrix.run)), match(rownames(matrix.run), names(res.mc.cluster)))

  colData(CS.data)$cluster = as.character(res.mc.cluster)
  #colLabels(CS.data) <- as.factor(res.mc.cluster)
  return(CS.data)
}


ClusterFind_hclust = function(CS.data, ClusterNum = 3, runWith = 'PCA', PCNum = 10) {

  matrix.run = switch(EXPR       = runWith,
                      'VGcounts' = assay(altExp(CS.data),'VGcounts'),
                      'PCA'      = reducedDim(x = CS.data, type = 'PCA')[,1:PCNum],
                      'ICA'      = reducedDim(x = CS.data, type = 'ICA')[,1:PCNum],
                      'tSNE'     = reducedDim(x = CS.data, type = 'tSNE')[,1:2],
                      'Umap'     = reducedDim(x = CS.data, type = 'Umap')[,1:2])

  suppressPackageStartupMessages(library(dplyr, quietly = T))
  res.hc = matrix.run %>%
    scale() %>%                    # Scale the data
    dist(method = "euclidean") %>% # Compute dissimilarity matrix
    hclust(method = "ward.D2")     # Compute hierachical clustering

  res.hc.cluster = cutree(tree = res.hc, k = ClusterNum)
  #plot(c(1:nrow(matrix.run)), match(rownames(matrix.run), names(res.hc.cluster)))

  colData(CS.data)$cluster = as.character(res.hc.cluster)
  #colLabels(CS.data)<- as.factor(res.hc.cluster)
  return(CS.data)
}


ClusterFind_graph = function(CS.data, runWith = 'PCA', PCNum = 10,cluster.fun) {
  library(bluster)
library(scran)
#  matrix.run = switch(EXPR       = runWith,
#                      'VGcounts' = assay(altExp(CS.data),'VGcounts'),
#                      'PCA'      = reducedDim(x = CS.data, type = 'PCA')[,1:PCNum],
#                      'ICA'      = reducedDim(x = CS.data, type = 'ICA')[,1:PCNum],
#                      'BEPCA'    = reducedDim(x = CS.data, type = 'BEPCA')[,1:PCNum],
#                      'tSNE'     = reducedDim(x = CS.data, type = 'tSNE')[,1:2],
#                      'Umap'     = reducedDim(x = CS.data, type = 'Umap')[,1:2])
#
#  com.graph = MUDAN::getComMembership(matrix.CS.datarun,
#                                       k=k,
#                                       method=igraph::cluster_infomap,
#                                       verbose=FALSE)
nn.clusters2 <- clusterCells(x=CS.data, use.dimred= runWith,
    BLUSPARAM=NNGraphParam( cluster.fun=cluster.fun))


  colData(CS.data)$cluster = as.character(nn.clusters2)
  #colLabels(CS.data)<-as.factor(nn.clusters2)
  return(CS.data)
}

ClusterFind_hclust_scran<- function(CS.data, runWith)
{
set.seed(1111)
khclust.info <- clusterCells(x=CS.data, use.dimred=runWith,
    BLUSPARAM=TwoStepParam(
        first=KmeansParam(centers=1000),
        second=HclustParam(method="ward.D2", cut.dynamic=TRUE,
            cut.param=list(deepSplit=3)) # for higher resolution.
    ),
    full=TRUE
)
colData(CS.data)$cluster = as.character(khclust.info$clusters)
#colLabels(CS.data)<- as.factor(khclust.info$clusters)
 return(CS.data)

}
# Setting the seed due to the randomness of k-means.
ClusterFind_kmean_scran<- function(CS.data, runWith,k=30){
set.seed(0101010)
kgraph.clusters <- clusterCells(x=CS.data, use.dimred= runWith,
    BLUSPARAM=TwoStepParam(
        first=KmeansParam(centers=1000),
        second=NNGraphParam(k=k)
    )
)
colData(CS.data)$cluster = as.character(kgraph.clusters)
#colLabels(CS.data)<- as.factor(kgraph.clusters)
return(CS.data)
}


ClusterFind = function(CS.data, method = 'Hclust', ClusterNum = 3, runWith = 'PCA', PCNum = 10, k=100,cluster.fun='walktrap') {

  CS.data = switch(EXPR = method,
                   'Hclust'	 = ClusterFind_hclust(   CS.data = CS.data,  ClusterNum = ClusterNum, runWith = runWith, PCNum = PCNum),
                   'K-Means'     = ClusterFind_kmean(  CS.data = CS.data,  ClusterNum = ClusterNum, runWith = runWith, PCNum = PCNum),
                   'Mclust'	 = ClusterFind_mclust(    CS.data = CS.data,  ClusterNum = ClusterNum, runWith = runWith, PCNum = PCNum),
                   'Graph'	 = ClusterFind_graph( CS.data = CS.data,   runWith = runWith, PCNum = PCNum, cluster.fun= cluster.fun),
                   'DBscan'      = ClusterFind_dbscan(CS.data=CS.data, runWith = runWith, PCNum = PCNum),
		   'Hclust_cran'	= ClusterFind_hclust_scran(CS.data=CS.data, runWith = runWith),
		   'Kmean_scran' = ClusterFind_kmean_scran(CS.data = CS.data,  k = k, runWith = runWith)
)


  return(CS.data)
}

#================#
# Differentially expressed analysis
#================#
Plot_HeatmapDEG = function(CS.data, used = 'counts', by = NULL, group1 = NULL, group2 = NULL, degTable, degTop = 10, expScale = F, expLog = T, outIndex = getwd(), savePlot = T) {
  #by: will take the first one as major to present the data
  edata = assay(CS.data, used)
  batch = CS.data@colData$batch
  pheno = colData(CS.data)

  if(is.null(batch)) {
    annotationForCol = data.frame(pheno)
    colnames(annotationForCol) = colnames(pheno)
  } else {
    annotationForCol = data.frame(pheno,
                                  batch)
    colnames(annotationForCol) = c(colnames(pheno), 'batch')
  }

  if (!is.null(by)) {
    annotationForCol = annotationForCol[,by, drop = F]
  }

  annotationForCol = annotationForCol[sample(nrow(annotationForCol)),,drop = F]

  if(!is.null(group1)) {
    if(!is.null(group2)) {
      od = c(grep(paste0("^", group1, '$'), annotationForCol[,1]),
             grep(paste0("^", group2, '$'), annotationForCol[,1]))

      annotationForCol = annotationForCol[od,, drop = F]
    } else {
      od.1 = grep(paste0("^", group1, '$'), annotationForCol[,1])
      od.2 = c(1:nrow(annotationForCol))[-grep(paste0("^", group1, '$'), annotationForCol[,1])]
      od = c(od.1, od.2)

      annotationForCol = annotationForCol[od,, drop = F]
    }
  } else {
    od = order(annotationForCol[,1], decreasing = F)
    annotationForCol = annotationForCol[od,, drop = F]
  }

  #if(is.null(degTop)) {
    degTable = .mySelectTopTailDEGs(degTable, FCVector = degTable$log2FC, TopTailNum = degTop)
  #}


  if(expLog) {
    #print('Expression matrix is log-transformed!')
    matrix = log2(edata[rownames(degTable), rownames(annotationForCol)] + 1)
  } else {
    matrix = edata[rownames(degTable), rownames(annotationForCol)]
  }

  if(expScale) {
    #print('Expression matrix is scaled!')
    x = matrix
    xz = t(apply(x, 1, function(x) { (x-mean(x))/sd(x) }))
    #apply(xz, 1, mean); apply(xz, 1, sd); max(xz); min(xz)
    xz[xz <= -3] = -3 ;xz[xz >=  3] =  3 #apply(xz, 1, mean); apply(xz, 1, sd)
    matrix = xz
  }


  bk = unique(c(seq(min(matrix), max(matrix), length=50))) #to set 1 as white
  colors = colorRampPalette(c("blue","white", "red"))(length(bk)-1)
# colors=colorRampPalette(c("#3182bd","#272222","#f2f609"))(1000)
  if(savePlot) {
    if (!(file.exists(paste0(outIndex, '/Figures_DEGHeatmap/')))){
      dir.create(paste0(outIndex, '/Figures_DEGHeatmap/'))
    }
    suppressPackageStartupMessages(library(pheatmap, quietly = T))
    ip = pheatmap::pheatmap(mat = matrix,
                           cluster_rows=F, #cutree_cols = 7,
                           cluster_cols=F,
                           clustering_method = "ward.D2",
                           show_rownames = T, fontsize_row = 8,
                           show_colnames = F,#fontsize_col = 3,
                           annotation_col = annotationForCol,
                           #annotation_row = annotationForRow,
                           #annotation_colors = anno_col,
                           color = colors, breaks = bk,
                           border_color = 'black',
                           nan_col = "gray",
                           #cellwidth = 8,
                           #cellheight = 8,
                           height = 10,
                           width = 10,
                           filename = paste0(paste0(outIndex, '/Figures_DEGHeatmap/'), "DEGHeatmap_With_", used, "_On_", paste(by, collapse="_"), ".pdf"))
  } else {
    p = pheatmap::pheatmap(mat = matrix,
                           cluster_rows=F, #cutree_cols = 7,
                           cluster_cols=F,
                           clustering_method = "ward.D2",
                           show_rownames = T, fontsize_row = 8,
    	                       show_colnames = F,#fontsize_col = 3,
                           annotation_col = annotationForCol,
                           #annotation_row = annotationForRow,
                           #annotation_colors = anno_col,
                           color = colors, breaks = bk,
                           border_color = NA,
                           nan_col = "gray")
  }
  #p
  #p = plot(p[[4]])
  return(p)

}

DEG_Analyzer_All = function(CS.data, used = 'counts', cmp, useBatch = FALSE, method = 'limma', groupVector = NULL, expt = 0.2, parallel_TF = FALSE, ncores = 1) {
  #library(SingleCellExperiment)
  #cmp = any column name of interest to find deg
  edata = assay(CS.data, used)
  #batch = CS.data@colData$batch
  pheno = colData(CS.data)

  if(useBatch) {
    batch = CS.data@colData$batch
  } else {
    batch = NULL
  }

  degList = NULL

  group1s = unique(pheno[,cmp][order(pheno[,cmp], decreasing = F)])

  for(i in c(1:length(group1s))) {
    #i = 2
    group1 = group1s[i]

    group2 = 'XXotherXX'
    groupVector = pheno[,cmp]
    groupVector[groupVector != group1] = group2
    groupVector = factor(groupVector, levels = c(group1, group2))

    res = switch(EXPR = method,
                 'limma'    = DEG_limma(   edata = edata, groupVector = groupVector,                group1 = group1, group2 = group2, expt = expt),
                 'edger'    = DEG_edger(   edata = edata, groupVector = groupVector, batch = batch, group1 = group1, group2 = group2, expt = expt),
                 'deseq'    = DEG_deseq(   edata = edata, groupVector = groupVector, batch = batch, group1 = group1, group2 = group2, expt = expt, parallel_TF = parallel_TF, ncores = ncores),
                 'mast'     = DEG_mast(    edata = edata, groupVector = groupVector, batch = batch, group1 = group1, group2 = group2, expt = expt, parallel_TF = parallel_TF, ncores = ncores),
                 'monocle'  = DEG_monocle( edata = edata, groupVector = groupVector, batch = batch, group1 = group1, group2 = group2, expt = expt, parallel_TF = parallel_TF, ncores = ncores),
                 'desingle' = DEG_desingle(edata = edata, groupVector = groupVector,                group1 = group1, group2 = group2, expt = expt, parallel_TF = parallel_TF, ncores = ncores)
    )

    degList[[i]] = data.frame(DEgene = rownames(res),
                              res)
  }

  names(degList) = group1s[c(1:length(group1s))]

  return(degList)
}

### For all factors ----
DEG_Analyzer_c2c = function(CS.data, used = 'counts', cmp, useBatch = FALSE, method = 'limma', groupVector = NULL, expt = 0.2, parallel_TF = FALSE, ncores = 1) {
  #library(SingleCellExperiment)
  edata = SummarizedExperiment::assay(CS.data, used)
  #batch = SingleCellExperiment::colData(CS.data)[,'batch']
  pheno = SummarizedExperiment::colData(CS.data)

  if(useBatch) {
    batch = SingleCellExperiment::colData(CS.data)[,'batch']
  } else {
    batch = NULL
  }

  degList = NULL

  group1s = unique(pheno[,cmp][order(pheno[,cmp], decreasing = F)])

  for(i in c(1:length(group1s))) {
    #i = 2
    group1 = group1s[i]

    group2 = 'XXotherXX'
    groupVector = pheno[,cmp]
    groupVector[groupVector != group1] = group2
    groupVector = factor(groupVector, levels = c(group1, group2))

    res = switch(EXPR = method,
                 'limma'    = DEG_limma(   edata = edata, groupVector = groupVector,                group1 = group1, group2 = group2, expt = expt),
                 'edger'    = DEG_edger(   edata = edata, groupVector = groupVector, batch = batch, group1 = group1, group2 = group2, expt = expt),
                 'deseq'    = DEG_deseq(   edata = edata, groupVector = groupVector, batch = batch, group1 = group1, group2 = group2, expt = expt, parallel_TF = parallel_TF, ncores = ncores),
                 'mast'     = DEG_mast(    edata = edata, groupVector = groupVector, batch = batch, group1 = group1, group2 = group2, expt = expt, parallel_TF = parallel_TF, ncores = ncores),
                 'monocle'  = DEG_monocle( edata = edata, groupVector = groupVector, batch = batch, group1 = group1, group2 = group2, expt = expt, parallel_TF = parallel_TF, ncores = ncores),
                 'desingle' = DEG_desingle(edata = edata, groupVector = groupVector,                group1 = group1, group2 = group2, expt = expt, parallel_TF = parallel_TF, ncores = ncores)
    )

    degList[[i]] = data.frame(DEgene = rownames(res),
                              res)
  }

  names(degList) = group1s[c(1:length(group1s))]

  return(degList)
}

DEG_Analyzer = function(CS.data, used = 'counts', cmp, useBatch = FALSE, method = 'limma', groupVector = NULL, group1, group2 = NULL, expt = 0.2, parallel_TF = FALSE, ncores = 1) {
  edata = assay(CS.data, used)
  #batch = CS.data@colData$batch
  pheno = colData(CS.data)

  if(useBatch) {
    batch = CS.data@colData$batch
  } else {
    batch = NULL
  }

  if (is.null(groupVector)) {
    if (is.null(group2)) {
      group2 = 'XXotherXX'
      groupVector = pheno[,cmp]
      groupVector[groupVector != group1] = group2
      groupVector = factor(groupVector, levels = c(group1, group2))
    } else {
      groupVector = factor(pheno[,cmp], levels = c(group1, group2))
      edata = edata[,!is.na(groupVector)]

      if(useBatch) {
        batch = CS.data@colData$batch[!is.na(groupVector)]
      }

      groupVector = groupVector[!is.na(groupVector)]
    }
  } else {
    if (is.null(group2)) {
      group2 = 'XXotherXX'
      groupVector[groupVector != group1] = group2
      groupVector = factor(groupVector, levels = c(group1, group2))
    } else {
      groupVector = factor(groupVector, levels = c(group1, group2))
      edata = edata[,!is.na(groupVector)]

      if(useBatch) {
        batch = CS.data@colData$batch[!is.na(groupVector)]
      } else {
        batch = NULL
      }

      groupVector = groupVector[!is.na(groupVector)]
    }
  }


  res = switch(EXPR = method,
               'limma'    = DEG_limma(   edata = edata, groupVector = groupVector,                group1 = group1, group2 = group2, expt = expt),
               'edger'    = DEG_edger(   edata = edata, groupVector = groupVector, batch = batch, group1 = group1, group2 = group2, expt = expt),
               'deseq'    = DEG_deseq(   edata = edata, groupVector = groupVector, batch = batch, group1 = group1, group2 = group2, expt = expt, parallel_TF = parallel_TF, ncores = ncores),
               'mast'     = DEG_mast(    edata = edata, groupVector = groupVector, batch = batch, group1 = group1, group2 = group2, expt = expt, parallel_TF = parallel_TF, ncores = ncores),
               'monocle'  = DEG_monocle( edata = edata, groupVector = groupVector, batch = batch, group1 = group1, group2 = group2, expt = expt, parallel_TF = parallel_TF, ncores = ncores),
               'desingle' = DEG_desingle(edata = edata, groupVector = groupVector,                group1 = group1, group2 = group2, expt = expt, parallel_TF = parallel_TF, ncores = ncores)
               )

  return(res)
}


#DEG_limma = function(CS.data, used = 'counts', cmp, groupVector = NULL, group1, group2 = NULL, expt = NULL) {
DEG_limma = function(edata, groupVector = NULL, group1, group2 = NULL, expt = NULL) {
  #edata = assay(CS.data, used)
  ##batch = CS.data@colData$batch
  #pheno = colData(CS.data)

  #if (is.null(groupVector)) {
  #  if (is.null(group2)) {
  #    group2 = 'XXotherXX'
  #    groupVector = pheno[,cmp]
  #    groupVector[groupVector != group1] = group2
  #    groupVector = factor(groupVector, levels = c(group1, group2))
  #  } else {
  #    groupVector = factor(pheno[,cmp], levels = c(group1, group2))
  #    edata = edata[,!is.na(groupVector)]
  #    groupVector = groupVector[!is.na(groupVector)]
  #  }
  #} else {
  #  if (is.null(group2)) {
  #    group2 = 'XXotherXX'
  #    groupVector[groupVector != group1] = group2
  #    groupVector = factor(groupVector, levels = c(group1, group2))
  #  } else {
  #    groupVector = factor(groupVector, levels = c(group1, group2))
  #    edata = edata[,!is.na(groupVector)]
  #    groupVector = groupVector[!is.na(groupVector)]
  #  }
  #}

  #
  groupVector = factor(groupVector, levels = c(group1, group2))
  #


  if(is.null(expt)) {
    edata = edata[apply(edata > 0, 1, sum) > 2,]
  } else {
    edata = edata[apply(edata > 0, 1, sum) > (expt * ncol(edata)),]
  }

  suppressPackageStartupMessages(library(limma, quietly = T))
  dge = edgeR::DGEList(counts=edata)
  design = model.matrix(~groupVector)

  dge = edgeR::calcNormFactors(dge)

  #v = voom(dge, design, plot=TRUE)
  #fit = lmFit(v, design)
  #fit = eBayes(fit)

  logCPM = edgeR::cpm(dge, log=TRUE, prior.count=1)
  fit = limma::lmFit(logCPM, design)
  fit = limma::eBayes(fit, trend=TRUE)

  top.table = limma::topTable(fit, coef=ncol(design), sort.by = 'P', n = Inf)

  res = top.table[,c(1,5)]; colnames(res) = c('log2FC', 'adjP');
  res$log2FC = res$log2FC * (-1); res = res[order(res$log2FC, decreasing = T),]
  res = res[res$adjP < 0.05,]
  #myShowDEGBoxplot(DEGlist = res, DEGlistFC = res$log2FC, ReadCount = edata, GroupVector = groupVector, TopTailNum = 10, outIndex = '~/')

  edata.group1 = edata[rownames(res),groupVector == group1] #dim(edata.group1)
  ratioG1 = apply(edata.group1 > 0, 1, sum) / ncol(edata.group1)
  edata.group2 = edata[rownames(res),groupVector == group2] #dim(edata.group2)
  ratioG2 = apply(edata.group2 > 0, 1, sum) / ncol(edata.group2)

  res$PTG1 = ratioG1
  res$PTG2 = ratioG2

  return(res)
}


DEG_edger = function(edata, groupVector = NULL, batch = NULL, group1, group2 = NULL, expt = NULL) {
  #
  groupVector = factor(groupVector, levels = c(group1, group2))
  #

  if(is.null(expt)) {
    edata = edata[apply(edata > 0, 1, sum) > 2,]
  } else {
    edata = edata[apply(edata > 0, 1, sum) > (expt * ncol(edata)),]
  }

  suppressPackageStartupMessages(library(edgeR, quietly = T))

  y = DGEList(counts = edata,
              group = groupVector)

  y = calcNormFactors(y)

  if(is.null(batch)) {
    y.design = model.matrix(~groupVector)
  } else {
    y.design = model.matrix(~batch + groupVector)
  }

  rownames(y.design) = colnames(y)

  y = estimateDisp(y, design = y.design)

  y = estimateGLMCommonDisp(y, y.design)
  y = estimateGLMTrendedDisp(y, y.design)
  y = estimateGLMTagwiseDisp(y, y.design)

  fitY = glmFit(y, y.design)
  qlfY = glmLRT(fitY)

  reDEG = topTags(object = qlfY, n = nrow(fitY), adjust.method = "BH", sort.by = "PValue")
  #print(reDEG$comparison)

  res = reDEG$table
  res = res[res$FDR < 0.05,]
  res$logFC = res$logFC * (-1); res = res[order(res$logFC, decreasing = T),]
  res = res[,c(1,5)]
  colnames(res) = c('log2FC','padj')
  #myShowDEGBoxplot(DEGlist = res, DEGlistFC = res$log2FC, ReadCount = edata, GroupVector = groupVector, TopTailNum = 10, outIndex = '~/')

  edata.group1 = edata[rownames(res),groupVector == group1] #dim(edata.group1)
  ratioG1 = apply(edata.group1 > 0, 1, sum) / ncol(edata.group1)
  edata.group2 = edata[rownames(res),groupVector == group2] #dim(edata.group2)
  ratioG2 = apply(edata.group2 > 0, 1, sum) / ncol(edata.group2)

  res$PTG1 = ratioG1
  res$PTG2 = ratioG2

  return(res)
}


DEG_deseq = function(edata, groupVector = NULL, batch = NULL, group1, group2 = NULL, expt = NULL, parallel_TF = FALSE, ncores = 1) {
  suppressPackageStartupMessages(library(DESeq2, quietly = T))

  if(is.null(expt)) {
    edata = edata[apply(edata > 0, 1, sum) > 2,]
  } else {
    edata = edata[apply(edata > 0, 1, sum) > (expt * ncol(edata)),]
  }

  cts = edata

  if(is.null(batch)) {
    coldata = data.frame(group = factor(groupVector, levels = c(group1, group2)))
    colnames(coldata) = c("Group")

    dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts,
                                          colData = coldata,
                                          design = ~ Group)
  } else {
    coldata = data.frame(Group = factor(groupVector, levels = c(group1, group2)),
                         batch = batch)
    colnames(coldata) = c("Group", "batch")

    dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts,
                                          colData = coldata,
                                          design = ~ batch + Group)
  }

  if(parallel_TF) {
    suppressPackageStartupMessages(library(BiocParallel, quietly = T))
    BiocParallel::register(MulticoreParam(ncores))
    dds = DESeq2::DESeq(dds, parallel = parallel_TF)
    res = DESeq2::results(dds)
  } else {
    dds = DESeq2::DESeq(dds)
    res = DESeq2::results(dds)
  }

  res = data.frame(res@listData, row.names = rownames(res))
  res = res[!is.na(res[,'padj']),]; res = res[res[,'padj'] < 0.05,]
  res[,'log2FoldChange'] = res[,'log2FoldChange'] * (-1); res = res[order(res[,'log2FoldChange'], decreasing = T), c(2,6)]
  colnames(res) = c('log2FC','padj')
  #myShowDEGBoxplot(DEGlist = res, DEGlistFC = res$log2FC, ReadCount = edata, GroupVector = groupVector, TopTailNum = 10, outIndex = '~/')

  edata.group1 = edata[rownames(res),groupVector == group1] #dim(edata.group1)
  ratioG1 = apply(edata.group1 > 0, 1, sum) / ncol(edata.group1)
  edata.group2 = edata[rownames(res),groupVector == group2] #dim(edata.group2)
  ratioG2 = apply(edata.group2 > 0, 1, sum) / ncol(edata.group2)

  res$PTG1 = ratioG1
  res$PTG2 = ratioG2

  return(res)
}


DEG_mast = function(edata, groupVector = NULL, batch = NULL, group1, group2 = NULL, expt = NULL, parallel_TF = FALSE, ncores = 1) {
  suppressPackageStartupMessages(library(MAST, quietly = T))

  if(is.null(expt)) {
    edata = edata[apply(edata > 0, 1, sum) > 2,]
  } else {
    edata = edata[apply(edata > 0, 1, sum) > (expt * ncol(edata)),]
  }

  fdata = data.frame(rownames(edata))
  rownames(fdata) = rownames(edata)
  M = data.frame(colnames(edata))

  sca = FromMatrix(log2(as.matrix(edata)+1), M, fdata)

  cdr2 = colSums(assay(sca)>0)
  colData(sca)$cngeneson = scale(cdr2)

  cond = factor(groupVector, levels = c(group1, group2))
  colData(sca)$condition = cond


  if(parallel_TF) {
    suppressPackageStartupMessages(library(BiocParallel, quietly = T))
    BiocParallel::register(MulticoreParam(ncores))
  }

  if(is.null(batch)) {
    zlmCond = zlm(~condition + cngeneson, sca, parallel = parallel_TF)
  } else {
    colData(sca)$batch = batch
    zlmCond = zlm(~condition + batch + cngeneson, sca, parallel = parallel_TF)
  }

  summaryCond = summary(zlmCond, doLRT=paste0('condition', group2))
  summaryDt = summaryCond$datatable

  fcHurdle = merge(summaryDt[contrast==paste0('condition', group2) & component=='H',.(primerid, `Pr(>Chisq)`)],
                   summaryDt[contrast==paste0('condition', group2) & component=='logFC',
                             .(primerid, coef, ci.hi, ci.lo)], by='primerid')

  fcHurdle[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  fcHurdle = fcHurdle[fcHurdle$FDR < 0.05,]

  res = data.frame(log2FC = fcHurdle$coef,
                   padj = fcHurdle$FDR,
                   row.names = fcHurdle$primerid)
  res$log2FC = res$log2FC * (-1); res = res[order(res$log2FC, decreasing = T),]

  edata.group1 = edata[rownames(res),groupVector == group1] #dim(edata.group1)
  ratioG1 = apply(edata.group1 > 0, 1, sum) / ncol(edata.group1)
  edata.group2 = edata[rownames(res),groupVector == group2] #dim(edata.group2)
  ratioG2 = apply(edata.group2 > 0, 1, sum) / ncol(edata.group2)

  res$PTG1 = ratioG1
  res$PTG2 = ratioG2
  #myShowDEGBoxplot(DEGlist = res, DEGlistFC = res$log2FC, ReadCount = edata, GroupVector = groupVector, TopTailNum = 40, outIndex = '~/')

  return(res)
}


DEG_monocle = function(edata, groupVector = NULL, batch = NULL, group1, group2 = NULL, expt = NULL, parallel_TF = FALSE, ncores = 1) {
  suppressPackageStartupMessages(library(monocle3, quietly = T))
  suppressPackageStartupMessages(library(dplyr, quietly = T))

  if(is.null(expt)) {
    edata = edata[apply(edata > 0, 1, sum) > 2,]
  } else {
    edata = edata[apply(edata > 0, 1, sum) > (expt * ncol(edata)),]
  }

  if(is.null(batch)) {
    cell_metadata= data.frame(group = factor(groupVector, levels = c(group1, group2)), row.names = colnames(edata))
    gene_annotation = data.frame(gene_short_name = rownames(edata), row.names = rownames(edata))

    cds = new_cell_data_set(edata,
                            cell_metadata = cell_metadata,
                            gene_metadata = gene_annotation)

    gene_fits = fit_models(cds, model_formula_str = "~group")
    fit_coefs = coefficient_table(gene_fits)
    emb_time_terms = fit_coefs %>% filter(term == paste0('group', group2))
    res = emb_time_terms %>% filter (q_value < 0.05) %>% select(gene_short_name, term, q_value, estimate)
    res = as.data.frame(res)
    row.names(res) = res[,1]
    res = data.frame(log2FC = res$estimate * (-1),
                     padj = res$q_value,
                     row.names = res$gene_short_name)
    res = res[order(res$log2FC, decreasing = T),]
    #myShowDEGBoxplot(DEGlist = res, DEGlistFC = res$log2FC, ReadCount = edata, GroupVector = groupVector, TopTailNum = 20, outIndex = '~/monocle_')
    edata.group1 = edata[rownames(res),groupVector == group1] #dim(edata.group1)
    ratioG1 = apply(edata.group1 > 0, 1, sum) / ncol(edata.group1)
    edata.group2 = edata[rownames(res),groupVector == group2] #dim(edata.group2)
    ratioG2 = apply(edata.group2 > 0, 1, sum) / ncol(edata.group2)

    res$PTG1 = ratioG1
    res$PTG2 = ratioG2
    #return(res)

  } else {
    cell_metadata= data.frame(group = factor(groupVector, levels = c(group1, group2)),
                              batch = batch,
                              row.names = colnames(edata))
    gene_annotation = data.frame(gene_short_name = rownames(edata), row.names = rownames(edata))

    cds = new_cell_data_set(edata,
                            cell_metadata = cell_metadata,
                            gene_metadata = gene_annotation)

    time_batch_models = fit_models(cds, model_formula_str = "~group + batch")#, expression_family="negbinomial")
    time_batch_models_coefs = coefficient_table(time_batch_models)
    #time_batch_models_coefs %>% filter(term != "(Intercept)") %>% select(gene_short_name, term, q_value, estimate)
    #evaluate_fits(time_batch_models)
    time_models = fit_models(cds, model_formula_str = "~group")#, expression_family="negbinomial")

    res = compare_models(time_batch_models, time_models) %>% select(gene_short_name, q_value)
    #return(res)
  }

  return(res)
}


DEG_desingle = function(edata, groupVector = NULL, group1, group2 = NULL, expt = NULL, parallel_TF = FALSE, ncores = 1) {
  suppressPackageStartupMessages(library(DEsingle, quietly = T))

  if(is.null(expt)) {
    edata = edata[apply(edata > 0, 1, sum) > 2,]
  } else {
    edata = edata[apply(edata > 0, 1, sum) > (expt * ncol(edata)),]
  }

  groupVector = factor(groupVector, levels = c(group1, group2))

  if(parallel_TF) {
    library(BiocParallel)

    # Set the parameters and register the back-end to be used
    param <- BiocParallel::MulticoreParam(workers = ncores, progressbar = TRUE)
    register(param)

    results <- DEsingle::DEsingle(counts = edata, group = groupVector, parallel = parallel_TF, BPPARAM = param)

  } else {
    results <- DEsingle::DEsingle(counts = edata, group = groupVector)

  }

  res = data.frame(log2FC = log2(results$foldChange),
                   padj = results$pvalue.adj.FDR,
                   row.names = rownames(results))
  res = res[res$padj < 0.05,]
  res = res[order(res$log2FC, decreasing = T),]

  #myShowDEGBoxplot(DEGlist = res, DEGlistFC = res$log2FC, ReadCount = edata, GroupVector = groupVector, TopTailNum = 20, outIndex = '~/monocle_')
  edata.group1 = edata[rownames(res), groupVector == group1] #dim(edata.group1)
  ratioG1 = apply(edata.group1 > 0, 1, sum) / ncol(edata.group1)
  edata.group2 = edata[rownames(res), groupVector == group2] #dim(edata.group2)
  ratioG2 = apply(edata.group2 > 0, 1, sum) / ncol(edata.group2)

  res$PTG1 = ratioG1
  res$PTG2 = ratioG2
  return(res)
}


#================#
# Cell development
#================#
Plot_Development_old = function(CS.data, method = 'CT', runWith = 'PCA', batch = F, colorBy = 'default', plot = TRUE, rev = F) {

  dimSe = c(1,2)
  CS.data = switch(EXPR       = method,
                   'PC'       = CellDev_Princurve(CS.data = CS.data, runWith = runWith),#, dimSe = dimSe),
                   'CT'       = CellDev_CytoTrace(CS.data = CS.data, used = 'counts', batch = batch),
                   'ICA'      = ,
                   'tSNE'     = ,
                   'Umap'     = )

  if (colorBy == 'default') {
    colorBy = colnames(colData(CS.data))[1]
  }

  cell.order = CS.data$cellOd
  if(rev) {
    cell.order = max(cell.order) - cell.order
  }

  df.plot = data.frame(reducedDim(CS.data, runWith)[,dimSe],
                       colorBy = as.factor(colData(CS.data)[,colorBy]),
                       color = as.character(colData(CS.data)[,colorBy]),
                       cellOd = cell.order,

                       stringsAsFactors = F)

  df.plot.od = df.plot[order(df.plot$cellOd, decreasing = F),]
  df.plot.od.tmp = data.frame()

  #### method 1
  #for(i in c(1:nrow(df.plot.od)-1)) {
  #  #i = 1
  #  df.plot.od.tmp[i,1] = 20*(df.plot.od[nrow(df.plot.od), 1] - df.plot.od[i, 1])/(nrow(df.plot.od) - 1) + df.plot.od[i, 1]
  #  df.plot.od.tmp[i,2] = 20*(df.plot.od[nrow(df.plot.od), 2] - df.plot.od[i, 2])/(nrow(df.plot.od) - 1) + df.plot.od[i, 2]
    #df.plot.od.tmp[i,1] = (df.plot.od[i+1, 1] - df.plot.od[i, 1])/(nrow(df.plot.od) - 1) + df.plot.od[i, 1]
    #df.plot.od.tmp[i,2] = (df.plot.od[i+1, 2] - df.plot.od[i, 2])/(nrow(df.plot.od) - 1) + df.plot.od[i, 2]
  #}
  #df.plot.od.tmp[nrow(df.plot.od),] = df.plot.od[nrow(df.plot.od),c(1,2)]
  #### method 1

  #### method 2
  x.ave = .myAverageWindow(VectorX = df.plot.od[,1], windowNum = round(nrow(df.plot.od)/5))
  y.ave = .myAverageWindow(VectorX = df.plot.od[,2], windowNum = round(nrow(df.plot.od)/5))
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
  #g.arrow


  library(colorspace, quietly = T)
  library(RColorBrewer, quietly = T)

  for(i in 1:nrow(df.plot)) {
    df.plot$color[i] = qualitative_hcl(length(unique(df.plot$colorBy)))[as.vector(as.integer(df.plot$colorBy))[i]]
  }

  library(ggplot2, quietly = T)
  #g = ggplot() +
  #  geom_point(data = df.plot, mapping = aes(x = df.plot[,1], y = df.plot[,2], color = colorBy)) +
  #  geom_point(data = df.fit.plot, mapping = aes(x = df.fit.plot[,1], y = df.fit.plot[,2]), shape = 20, color = 'black') +
  #  #ggtitle() +
  #  theme_bw() +
  #  theme(panel.grid.major = element_blank(),
  #        panel.grid.minor = element_blank(),
  #        axis.line = element_line(colour = "black")) +
  #  theme(legend.position="bottom",
  #        legend.title = element_blank()) +
  #  xlab(colnames(df.plot)[1]) +
  #  ylab(colnames(df.plot)[2])
  #g

  #df.plot.od = df.plot[order(),]
  #delta_long = 0.01
  #delta_lat = 0.01
  #g = ggplot() +
  #  geom_point(data = df.plot, mapping = aes(x = df.plot[,1], y = df.plot[,2], color = colorBy)) +
  #  geom_segment(data = df.plot.arrow, aes(#x = df.plot.arrow[,1],
  #                                         xend = df.plot.arrow[,2] + delta_long,
  #                                         #y = df.plot.arrow[,3],
  #                                         yend = df.plot.arrow[,4] + delta_lat),
  #               arrow = arrow(length = unit(0.1,"cm")))
  #g

  g.density.x = ggplot(data = df.plot, mapping = aes(x = cellOd, y = ..density.., color = colorBy)) +
    geom_density() +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black")) +
    theme(#axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()) +
    theme(legend.position = "none") +
    xlab('Cell order') +
    ylab('Density') #+
  #scale_y_reverse()
  #g.density.x


  g.legend = ggplot(data = df.plot, mapping = aes(x = 0, y = 0, color = colorBy)) +
    geom_point() +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black")) +
    theme(legend.position="bottom",
          legend.title = element_blank())
  #g.legend

  figures_without_legend = cowplot::plot_grid(
    g.arrow + theme(legend.position = "none"),
    g.density.x,

    align = "hv",
    nrow = 2,
    #rel_widths = c(0.5, 1),
    rel_heights = c(1, 0.5)
  )

  legends = cowplot::get_legend(g.legend)

  g.all = cowplot::plot_grid(figures_without_legend,
                             legends,
                             rel_heights = c(10,2), nrow = 2)
  print(g.all)

  return(CS.data)
}

Plot_Development = function(CS.data, runWith = 'Umap', dimSe = c(1,2), colorBy = 'default', rev = F, devWindowNum = 5) {

  if (colorBy == 'default') {
    colorBy = colnames(colData(CS.data))[1]
  }

  cell.order = CS.data$cellOd
  if(rev) {
    cell.order = max(cell.order) - cell.order
  }

  df.plot = data.frame(reducedDim(CS.data, runWith)[,dimSe],
                       colorBy = as.factor(colData(CS.data)[,colorBy]),
                       color = as.character(colData(CS.data)[,colorBy]),
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

Plot_Development_Princurve = function(CS.data, runWith = 'Umap', dimSe = c(1,2), colorBy = 'default', rev = F) {
  if (colorBy == 'default') {
    colorBy = colnames(colData(CS.data))[1]
  }

  cell.order = CS.data$cellOd
  if(rev) {
    cell.order = max(cell.order) - cell.order
  }

  suppressPackageStartupMessages(library(ggplot2, quietly = T))
  df.plot = data.frame(reducedDim(CS.data, runWith)[,dimSe],
                       fit_x = reducedDim(CS.data, 'PrinFit')[,1],
                       fit_y = reducedDim(CS.data, 'PrinFit')[,2],
                       colorBy = as.factor(colData(CS.data)[,colorBy]),
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

### Find pseudo trajectory ----

CellDevelopment <- function(CS.data,  method='Princurve_fit', runWith = 'PCA', dimSe = c(1,2),  used = 'counts', batch = F , clusterSe= 'cluster', seedBy){

  CS.data <- switch(EXPR = method,
                    "Princurve_fit" = CellDev_Princurve_fit(CS.data, runWith = runWith, dimSe = dimSe),
                    "CytoTrace" = CellDev_CytoTrace(CS.data, used = used, batch = batch),
                    "Slingshot" = CellDev_Slingshot(CS.data, clusterSe= clusterSe, runWith = runWith),
                    "Monocle" = CellDev_Monocle(CS.data, seedBy))
  return(CS.data)
}



CellDev_Princurve_fit = function(CS.data, runWith = 'PCA', dimSe = c(1,2)) {

  fit1 = princurve::principal_curve(SingleCellExperiment::reducedDim(CS.data, runWith)[,dimSe], plot = F, smoother = "smooth_spline")
  cell.order = fit1$lambda
  SummarizedExperiment::colData(CS.data)$cellOd = cell.order
  SingleCellExperiment::reducedDim(CS.data, 'PrinFit') = fit1$s
  return(CS.data)
}

CellDev_CytoTrace = function(CS.data, used = 'counts', batch = F) {
  edata = SummarizedExperiment::assay(CS.data, used)
  #suppressPackageStartupMessages(library(CytoTRACE, quietly = T))
  if(!batch) {
    results2 = CytoTRACE::CytoTRACE(mat = log2(edata + 1))
    SummarizedExperiment::colData(CS.data)$cellOd = results2$CytoTRACE
  } else {
    batch = SingleCellExperiment::colData(CS.data)[,'batch']
    results2 = CytoTRACE::CytoTRACE(mat = log2(edata + 1), batch = batch)
    SummarizedExperiment::colData(CS.data)$cellOd = results2$CytoTRACE
  }
  return(CS.data)
}

CellDev_Slingshot = function(CS.data, clusterSe= 'label', runWith = 'PCA') {
  #clusterSe = 'ClusterGraph'
  sim = slingshot::slingshot(CS.data, clusterLabels = clusterSe, reducedDim = runWith)
  SummarizedExperiment::colData(CS.data)$cellOd = sim@colData$slingPseudotime_1
  #plot(SingleCellExperiment::reducedDim(CS.data, runWith))
  #lines(slingshot::SlingshotDataSet(sim), lwd=2, col='black')
  return(CS.data)

}

CellDev_Monocle = function(CS.data, seedBy) {
  #seedBy = Cell Type or cluster:(Single or multiple)
  #E.g. cluster: 1:5 or 6
  cds = monocle3::new_cell_data_set(expression_data = SummarizedExperiment::assay(CS.data, 'counts'),
                                    cell_metadata = SingleCellExperiment::colData(CS.data),
                                    gene_metadata = data.frame(gene_short_name = rownames(SummarizedExperiment::assay(CS.data, 'counts')),
                                                               row.names = rownames(SummarizedExperiment::assay(CS.data, 'counts'))))
  col.name = colnames(SingleCellExperiment::colData(CS.data))[apply(SingleCellExperiment::colData(CS.data), 2, function(col) any(col == seedBy))]
  cds = monocle3::preprocess_cds(cds, num_dim = 50)
  if('batch' %in% colnames(SingleCellExperiment::colData(CS.data))) { cds = monocle3::align_cds(cds, alignment_group = "batch") }
  cds = monocle3::reduce_dimension(cds)
  cds = monocle3::cluster_cells(cds)
  cds = monocle3::learn_graph(cds)
  cell_ids = which(SingleCellExperiment::colData(cds)[, col.name] == seedBy)
  closest_vertex = cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex = as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes = igraph::V(monocle3::principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  cds = monocle3::order_cells(cds, root_pr_nodes=root_pr_nodes)
  SummarizedExperiment::colData(CS.data)$cellOd = monocle3::pseudotime(cds)[colnames(SummarizedExperiment::assay(CS.data, 'counts'))]
  return(CS.data)
}

#================#
# Plot markers
#================#
mks = c('T', 'KRT15')
Plot_MarkersScatter = function(CS.data, mks, used = 'logcounts', runWith = 'tSNE', dim.se = c(1,2), shapeBy = NULL, splitBy = NULL) {
  matrix.run = switch(EXPR       = runWith,
                      'PCA'      = reducedDim(x = CS.data, type = 'PCA')[,dim.se],
                      'ICA'      = reducedDim(x = CS.data, type = 'ICA')[,dim.se],
                      'tSNE'     = reducedDim(x = CS.data, type = 'tSNE'),
                      'Umap'     = reducedDim(x = CS.data, type = 'Umap'))

  edata = assay(CS.data, used)

  batch = CS.data@colData$batch
  pheno = colData(CS.data)

  library(ggplot2)
  library(RColorBrewer)
  library(scales)
  library(ggpubr)

  edata.mks = edata[mks,]


  if(!is.null(shapeBy)) {
    pheno = pheno[,shapeBy]
    edata.mks = cbind(data.frame(shape = pheno),
                      matrix.run,
                      as.data.frame(t(edata.mks)))

    df = tidyr::gather(data = edata.mks, key = 'gene', value = 'exp', -colnames(edata.mks)[1:3])

  } else {
    edata.mks = cbind(matrix.run,
                      as.data.frame(t(edata.mks)))

    df = tidyr::gather(data = edata.mks, key = 'gene', value = 'exp', -colnames(edata.mks)[1:2])
  }

  #df.label.x = summarySE(data = df, measurevar = colnames(df)[1], groupvars = 'shape')
  #df.label.y = summarySE(data = df, measurevar = colnames(df)[1], groupvars = 'shape')
  #df.label = data.frame(shape = df.label.x$shape,
  #                      x = df.label.x$UMAP_1,
  #                      y = df.label.y$UMAP_2)

  g = ggplot(data = df, mapping = aes(x = df[,1], y = df[,2])) +
    geom_point(mapping = aes (colour = exp), size = 1, data = df) +
    scale_colour_gradient2(low = "blue", mid = "white" ,high = "red",
                           midpoint = ((max(df$exp)+min(df$exp))/2),
                           space = "Lab", na.value = "midnightblue", guide = "colourbar",
                           limits=c(min(df$exp), max(df$exp)),
                           name = used) +
    #scale_shape_manual(values = 1:(nlevels((shape)))) +
    xlab(paste0(colnames(df)[1])) +
    ylab(paste0(colnames(df)[2])) +
    ggtitle(mks) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(legend.position = 'bottom')
  #g = g + guides(color = guide_legend(title="VST exp"), shape = guide_legend(title = 'Cluster'))
  #g = g + annotate(geom = 'text', x = df.label$x, y = df.label$y, #font = 2,
  #                 size = 5,
  #                 label = df.label$shape,
  #                 colour = 'cyan')

  return(g)
  #ggsave(filename = outName, plot = g, height = height, width = width, dpi = 600, useDingbats=FALSE)
}


#================#
# Other functions
#================#
.mySelectTopTailDEGs = function(degTable, FCVector, TopTailNum) {
  dif.gene = degTable
  difNum = TopTailNum

  dif.gene.pos = dif.gene[FCVector > 0,]
  dif.gene.neg = dif.gene[FCVector < 0,]

  if ((nrow(dif.gene.pos) < TopTailNum) & (nrow(dif.gene.pos) > 0)) {
    TopTailNum.pos = nrow(dif.gene.pos)
  } else if (nrow(dif.gene.pos) == 0) {
    TopTailNum.pos = 0
  } else {
    TopTailNum.pos = TopTailNum
  }

  if ((nrow(dif.gene.neg) < TopTailNum) & (nrow(dif.gene.neg) > 0)) {
    TopTailNum.neg = nrow(dif.gene.neg)
  } else if (nrow(dif.gene.neg) == 0) {
    TopTailNum.neg = 0
  } else {
    TopTailNum.neg = TopTailNum
  }

  dif.gene = dif.gene[order(FCVector, decreasing = T),]
  if (TopTailNum.pos == 0) {
    dif.gene.list = dif.gene[c(c((nrow(dif.gene) - TopTailNum.neg+1): nrow(dif.gene))),]
  } else if (TopTailNum.neg == 0) {
    dif.gene.list = dif.gene[c(c(1:TopTailNum.pos)),]
  } else {
    dif.gene.list = dif.gene[c(c(1:TopTailNum.pos), c((nrow(dif.gene) - TopTailNum.neg+1): nrow(dif.gene))),]
  }

  return(dif.gene.list)
}

.myAverageWindow = function(VectorX, windowNum) {
  VectorX_new = NULL
  for (i in (1:length(VectorX))){
    #i = 1
    tmp = VectorX[i:(i+windowNum-1)]
    tmp = tmp[!is.na(tmp)]
    tmp = tmp[!is.nan(tmp)]
    tmp = tmp[!is.infinite(tmp)]
    VectorX_new[i] = mean(tmp)
  }
  return(VectorX_new)
}

.mySelect_PCNum = function(CS.data, used = 'VGcounts', method = 'elbow') {
  if (used == 'VGcounts') {
    matrix.PJ = assay(altExp(CS.data), used)
    matrix.PJ = log2(matrix.PJ + 1)
  }

  if (used == 'PCA') {
    matrix.PJ = t(reducedDim(CS.data, used))
  }

  if (used == 'BEPCA') {
    matrix.PJ = t(reducedDim(CS.data, used))
  }


  suppressPackageStartupMessages(library(PCAtools, quietly = T))
  p = PCAtools::pca(matrix.PJ)
  #PCAtools::screeplot(p, axisLabSize = 18, titleLabSize = 22)
  #PCAtools::biplot(p)
  #PCAtools::plotloadings(p, labSize = 3)
  elbow = PCAtools::findElbowPoint(p$variance)

  horn = PCAtools::parallelPCA(matrix.PJ)

  ymax = max(elbow, horn$n)

  if (method == 'elbow') {
    pc.num = elbow
  }

  if (method == 'horn') {
    pc.num = elbow
  }

  return(pc.num)
}

.myTestMarker = function(CS.data, markers, used = 'counts', runWith = 'PCA') {

  if(used == 'counts') {
    edata = assay(CS.data, used)
    edata = log2(edata + 1)
  } else {
    edata = assay(CS.data, used)
  }

  emb = reducedDim(CS.data, runWith)

  matrix.exp = edata[markers,,drop = F]

  aql = reshape2::melt(t(matrix.exp))
  aql = aql[,-1]
  colnames(aql) = c('Gene', 'value')

  df_emb = emb
  for(i in c(1:(length(markers)-1))) {
    df_emb = rbind(df_emb, emb)
  }

  df_plot = cbind(df_emb,
                  aql)

  g = ggplot2::ggplot(mapping = aes(x = df_plot[,1], y = df_plot[,2]), data = aql) +
    ggplot2::geom_point(mapping = aes (colour = value), size = 1, data = df_plot) +
    ggplot2::scale_colour_gradient2(low = "blue", mid = "white" ,high = "red",
                                    midpoint = ((max(df_plot$value)+min(df_plot$value))/2),
                                    space = "Lab", na.value = "midnightblue", guide = "colourbar",
                                    limits=c(min(df_plot$value), max(df_plot$value))) +
    ggplot2::xlab(paste0(colnames(df_plot)[1])) +
    ggplot2::ylab(paste0(colnames(df_plot)[2])) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.line = element_line(colour = "black")) +
    ggplot2::facet_wrap(~Gene) +
    ggplot2::theme(legend.title = element_blank())
  return(g)
}

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
.summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}


## Norms the data within specified groups in a data frame; it normalizes each
## subject (identified by idvar) so that they have the same mean, within each group
## specified by betweenvars.
##   data: a data frame.
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   na.rm: a boolean that indicates whether to ignore NA's
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
    library(plyr)

    # Measure var on left, idvar + between vars on right of formula.
    data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
     .fun = function(xx, col, na.rm) {
        c(subjMean = mean(xx[,col], na.rm=na.rm))
      },
      measurevar,
      na.rm
    )

    # Put the subject means with original data
    data <- merge(data, data.subjMean)

    # Get the normalized data in a new column
    measureNormedVar <- paste(measurevar, "_norm", sep="")
    data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
                               mean(data[,measurevar], na.rm=na.rm)

    # Remove this subject mean column
    data$subjMean <- NULL

    return(data)
}

## Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, un-normed mean, normed mean (with same between-group mean),
##   standard deviation, standard error of the mean, and confidence interval.
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {

  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
    FUN=is.factor, FUN.VALUE=logical(1))

  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }

  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)

  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL

  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)

  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")

  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)

  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                           FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )

  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor

  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}
