
################################################################################
###################### METHODS FOR DEG ANALYSIS ################################
################################################################################

#' A function to identify differential expressed genes one group vs rest of the groups
#'
#' @param sce a SingleCellExperiment object
#' @param used run with counts
#' @param cmp group 1 for comparison. it may be cluster/disease condition
#' @param useBatch true/false
#' @param method  method to identify DEG
#'
#' @return  a list of tables for DEG
#'
#' DEG_Analyzer_All(sce, cmp = 1)
#'

DEG_Analyzer_All = function(sce, used = 'counts', cmp, useBatch = FALSE,
                            method = 'limma', groupVector = NULL, expt = 0.2, parallel_TF = FALSE, ncores = 1) {
  edata = assay(sce, used)
  #batch = sce@colData$batch
  pheno = colData(sce)

  if(useBatch) {
    batch = sce@colData$batch
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
#' A function to identify differential expressed genes one group vs rest of the groups
#'
#' @param sce a SingleCellExperiment object
#' @param used run with counts
#' @param cmp group 1 for comparison. it may be cluster/disease condition
#' @param useBatch true/false
#' @param method  method to identify DEG
#' @param cmp column name of the comparision vector
#' @return  a list of tables for DEG
#'
#' DEG_Analyzer_All(sce, cmp = 'cluster')
#'

DEG_Analyzer_c2c = function(sce, used = 'counts', cmp, useBatch = FALSE, method = 'limma', groupVector = NULL, expt = 0.2, parallel_TF = FALSE, ncores = 1) {
  #library(SingleCellExperiment)
  edata = SummarizedExperiment::assay(sce, used)
  #batch = SingleCellExperiment::colData(sce)[,'batch']
  pheno = SummarizedExperiment::colData(sce)

  if(useBatch) {
    batch = SingleCellExperiment::colData(sce)[,'batch']
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

#' A function to identify differential expressed genes one group vs another group
#'
#' @param sce a SingleCellExperiment object
#' @param used run with counts
#' @param group1 1st group  for comparison. it may be cluster/disease condition
#' @param group2  2nd group  for comparison. it may be cluster/disease condition
#' @param useBatch use batch T/F
#' @param method  method to identify DEG
#' @param cmp column name of the comparison vector
#' @return  a a DEG table
#'
#' DEG_Analyzer(sce, cmp = 'cluster', group1 =1, group2=2 , method = 'limma' )
#'

DEG_Analyzer = function(sce, used = 'counts', cmp, useBatch = FALSE, method = 'limma',
                        groupVector = NULL, group1, group2 = NULL, expt = 0.2, parallel_TF = FALSE, ncores = 1) {
  edata = assay(sce, used)
  #batch = sce@colData$batch
  pheno = colData(sce)

  if(useBatch) {
    batch = sce@colData$batch
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
        batch = sce@colData$batch[!is.na(groupVector)]
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
        batch = sce@colData$batch[!is.na(groupVector)]
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


#DEG_limma = function(sce, used = 'counts', cmp, groupVector = NULL, group1, group2 = NULL, expt = NULL) {
# A function to identify DEG using limma package
DEG_limma = function(edata, groupVector = NULL, group1, group2 = NULL, expt = NULL) {
  #edata = assay(sce, used)
  ##batch = sce@colData$batch
  #pheno = colData(sce)

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

# A function to identify DEG using edger package
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

# A function to identify DEG using deseq2 package
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

# A function to identify DEG using mast package
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

# A function to identify DEG using monocle package
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

# A function to identify DEG using desingle package
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




#' A function to make a DEG table
#'
#' @param DEG a Differentilly expressed gene Dataframe
#' @param fc Fold change cutoff
#' @param reduction PCA
#' @inputs pval p value cutoff
#' @return a dataframe
#' @examples
#'
#' DEG_table(DEG, fc=2, pval= 0.005)
#'


#DEG Table
DEG_table<- function(DEG, fc=2, pval= 0.005){
  library(ggpubr)
  library(ggrepel)
  d=DEG
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

################################################################################
########################## VISUALIZATION #######################################
################################################################################

#' A function to plot heatmap of differential expressed genes
#'
#' @param sce a SingleCellExperiment object
#' @param used run with counts
#' @param group1 group 1 for comparison. it may be cluster/disease condition
#' @param group2 group 2 for comparison. it may be cluster/disease condition
#'
#' @return  a heatmap
#'
#' Plot_HeatmapDEG(sce, group1 =1, group2 = 12)
#'
Plot_HeatmapDEG = function(sce, used = 'counts', by = NULL, group1 = NULL, group2 = NULL,
                           degTable, degTop = 10, expScale = F, expLog = T, outIndex = getwd(), savePlot = T) {
  #by: will take the first one as major to present the data
  edata = assay(sce, used)
  batch = sce@colData$batch
  pheno = colData(sce)

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


#' A function to make volcano plot from a DEG table
#'
#' @param deg_table a Differentilly expressed gene table (Dataframe)
#' @param fc Fold change cutoff
#' @param reduction PCA
#' @param pval p value cutoff
#' @return a a ggplot
#' @examples
#'
#' volcano_plot(deg_table, fc=2, pval= 0.005)
#'
volcano_plot<- function(deg_table, fc=2, pval= 0.005){
  library(ggpubr)
  library(ggrepel)
  d=deg_table
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
    ylab(paste0("-log10 ",names(d)[2]))+
    xlab("Log2FoldChange")+
    theme_pubr(base_size = 12,border=TRUE)+
    geom_hline(yintercept=-log10(pval), linetype="dashed",color = "black", linewidth=0.5)+
    geom_vline(xintercept=c(-fc,fc), linetype="dashed",color = "black", linewidth=0.5)+
    geom_label_repel(data=top,fontface="bold",color="purple",box.padding=unit(1, "lines"),
                     point.padding=unit(0.5, "lines"),segment.colour = "purple",segment.size = 0.5,segment.alpha = 0.5,max.overlaps = Inf)+
    annotate(geom = 'text', label = paste0('UP_Number: ', up_num), x = Inf, y = Inf, hjust = 1.1, vjust = 1.5)+
    annotate(geom = 'text', label = paste0('DOWN_Number: ', down_num), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5)+
    labs(color = "class")
  return(p)
}

#================#
# Plot markers
#================#

#' A function to Plot marker genes of interest
#'
#' @param sce a SingleCellExperiment object
#' @param method plot on expression value or on reduction
#' @param colour_by column name in colData
#' @param features gene symbols one or multiple
#' @return  plots
#' @examples
#'
#' plot_filterd_data(sce, colour_by="All")
#'
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

mks = c('T', 'KRT15')
Plot_MarkersScatter = function(sce, mks, used = 'logcounts', runWith = 'tSNE', dim.se = c(1,2), shapeBy = NULL, splitBy = NULL) {
  matrix.run = switch(EXPR       = runWith,
                      'PCA'      = reducedDim(x = sce, type = 'PCA')[,dim.se],
                      'ICA'      = reducedDim(x = sce, type = 'ICA')[,dim.se],
                      'tSNE'     = reducedDim(x = sce, type = 'tSNE'),
                      'Umap'     = reducedDim(x = sce, type = 'Umap'))

  edata = assay(sce, used)

  batch = sce@colData$batch
  pheno = colData(sce)

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

