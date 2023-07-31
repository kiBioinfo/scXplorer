
################################################################################
##################### METHODS FOR CELL TYPE ANNOTATION #########################
################################################################################


#' A function to automatically annotate Cell Type
#'
#' @param sce a SingleCellExperiment object
#' @param method SingleR/sctype
#' @param dataset for singleR method: celldex::BlueprintEncodeData()
#' @param TissueType Tissue Type
#' @param used count matrix raw or normalized
#' @return a cell type annotated SingleCellExperiment object
#' TissueType information is only required for sctype method
#' @examples
#'
#' Cell_Type_Annotation(data=sce,method='SingleR',dataset="BlueprintEncodeData", TissueType, used='counts')
#'

Cell_Type_Annotation<- function(data=sce,method='SingleR',
                                dataset="BlueprintEncodeData", TissueType,  clusters, used='counts')
{
  clusters=data@colData$cluster
  # sce= SingleCellExperiment object
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

################################################################################
########################## VISUALIZATION #######################################
################################################################################
