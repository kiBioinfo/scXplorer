# Read in file and perform validation.
library(xfun)
load_sce_obj <- function(path){
  errors <- c()


  # try to read in file


      if(tolower(tools::file_ext(path)) == "rds") {
        obj <- readRDS(path)
      }
      else if(tolower(tools::file_ext(path)) == "txt" || tolower(xfun::file_ext(path)) == "gz" || tolower(tools::file_ext(path)) == "csv"){
        obj <- vroom::vroom(path)
        colnames(obj)[1]<-"GeneID"
        if( any(str_detect(toupper(obj$GeneID), "ENSG")) | any(str_detect(toupper(obj$GeneID), "ENSMU"))){
          obj$GeneID <- stringr::str_remove(obj$GeneID, ".*_")
        }

        if(any(str_detect(obj$GeneID, "Mar-"))){
          obj=obj[!str_detect(obj$GeneID, "Mar-"),]
        }

        obj<-obj %>% remove_rownames() %>% column_to_rownames(var="GeneID")
      }


      else{
        obj<- c(errors, "Invalid  file.")
      }

  # Validate obj is a SingleCellExperiment object
  if (inherits(obj, "character")){
    errors <- c(errors)
    return(errors)
  }
  else if(typeof(obj)=="S4")
  {
    if(!inherits(obj, "SingleCellExperiment"))
    {
      obj<-as.SingleCellExperiment(obj)

    }
    return(obj)
  }
  else{
    obj = SingleCellExperiment(assay = list(counts = obj))
    if( any(str_detect(rownames(obj), "ENSG")) | any(str_detect(rownames(obj), "ENSMU"))){
      obj <- convert_to_Symbol(obj)
    }
    else{
      obj<- obj
    }
    return(obj)
  }

}

convert_to_Symbol <- function(obj){
  library(biomaRt)
  if( any(str_detect(rownames(obj), "ENSG")) & any(str_detect(rownames(obj), "ENSMU")))
  {
    # symb_mm <- mapIds(org.Mm.eg.db, keys=rownames(obj), keytype="ENSEMBL", column="SYMBOL")
    # symb_mm <- as.data.frame(symb_mm)
    # symb_mm$GeneID <- rownames(symb_mm)
    # colnames(symb_mm)[1]<- 'SYMBOL'
    # symb_mm <-symb_mm %>%  dplyr::filter(grepl('ENSMU', GeneID, ignore.case =T))
    #
    # symb_hs <- mapIds(org.Hs.eg.db, keys=rownames(obj), keytype="ENSEMBL", column="SYMBOL")
    # symb_hs <- as.data.frame(symb_hs)
    # symb_hs$GeneID <- rownames(symb_hs)
    # colnames(symb_hs)[1]<- 'SYMBOL'
    # symb_hs <-symb_hs %>%  dplyr::filter(grepl('ENSG', GeneID, ignore.case =T))
    # symb <- rbind(symb_hs, symb_mm)
    ensemble_h <- useMart("ensembl", dataset = c( 'hsapiens_gene_ensembl'))
    gene_info_h <- getBM(attributes=c( "ensembl_gene_id", 'hgnc_symbol'), filters="ensembl_gene_id", values= rownames(obj), mart=ensemble_h)
    ensemble_m <- useMart("ensembl", dataset = c( 'mmusculus_gene_ensembl'))
    gene_info_m <- getBM(attributes=c( "ensembl_gene_id", 'external_gene_name'), filters="ensembl_gene_id", values= rownames(obj), mart=ensemble_m)
    colnames(gene_info_m)[2] <- 'hgnc_symbol'
    names<- as.data.frame(rownames(obj))
    colnames(names)[1] <- 'ensembl_gene_id'
    symbol=rbind(gene_info_m, gene_info_h)
    dt<- merge(names, symbol, all= TRUE)
  }
  else if(any(str_detect(rownames(obj), "ENSG"))==TRUE & any(str_detect(rownames(obj), "ENSMU")) == FALSE) {
    ensemble_h <- useMart("ensembl", dataset = c( 'hsapiens_gene_ensembl'))
    gene_info_h <- getBM(attributes=c( "ensembl_gene_id", 'hgnc_symbol'), filters="ensembl_gene_id", values= rownames(obj), mart=ensemble_h)
    names<- as.data.frame(rownames(obj))
    colnames(names)[1] <- 'ensembl_gene_id'
    dt<- merge(names, gene_info_h, all= TRUE)
  }
  else if(any(str_detect(rownames(obj), "ENSG"))==FALSE & any(str_detect(rownames(obj), "ENSMU")) == TRUE) {
    ensemble_m <- useMart("ensembl", dataset = c( 'mmusculus_gene_ensembl'))
    gene_info_m <- getBM(attributes=c( "ensembl_gene_id", 'external_gene_name'), filters="ensembl_gene_id", values= rownames(obj), mart=ensemble_m)
    colnames(gene_info_m)[2] <- 'hgnc_symbol'
    names<- as.data.frame(rownames(obj))
    colnames(names)[1] <- 'ensembl_gene_id'
    dt<- merge(names, gene_info_m, all= TRUE)

  }
  rowData(obj)$ENSEMBL <- rownames(obj)
  rowData(obj)$SYMBOL <- dt$hgnc_symbol
  new.names <- rowData(obj)$SYMBOL
  missing.name <- is.na(new.names)
  new.names[missing.name] <- rowData(obj)$ENSEMBL[missing.name]
  dup.name <- new.names %in% new.names[duplicated(new.names)]
  new.names[dup.name] <- paste0(new.names, "_", rowData(obj)$ENSEMBL)[dup.name]
  rownames(obj) <- new.names

  return(obj)
}



metadata_modal<- function(current_data){
  if(is.null(current_data)){
    d=read.table("www/annotation_lab.txt",
                 header = TRUE)


   # d[is.na(d)] <- ""
  }

  else{
    d<-current_data
    anno<-c()
    if(class(d)=="SingleCellExperiment") anno$Cell<-rownames(d@colData)
    if(class(d)=="data.frame") anno$Cell<-rownames(d)

    d<-as.data.frame(anno)
    if(ncol(d)==1){
      d$Patients=NA
      d$Type=NA
      d$Batch=NA
    }
  }
  d<-as.data.frame(d)
  return(d)
}

add_QC_matrix<- function(sce, metadata=NULL){
  if(is.null(metadata)){
    count<- assay(sce, "counts")
    sce<- SingleCellExperiment(assay = list(counts = count))
  }
  else{
    metadata<- read.table(metadata, header = TRUE, sep = "\t")
    colnames(metadata)[1]<-"Cell"
    if(ncol(metadata)==1){
      metadata$Patients=NA
      metadata$Type=NA
      metadata$Batch=NA
    }
    #metadata<-metadata %>% remove_rownames() %>% column_to_rownames(var="Cell")
    count<- assay(sce, "counts")
    count<-count[,colnames(count) %in% metadata$Cell]
    metadata<- metadata[metadata$Cell %in% colnames(count),]
    sce<- SingleCellExperiment(assay = list(counts = count),
                               colData = data.frame(metadata))
    if(any(colnames(sce@colData)=='Batch')){
      colnames(colData(sce))[colnames(colData(sce))=='Batch'] = 'batch'
    }



     }



  #calculate %MT genes
  mito_genes <- rownames(sce)[grep("(?i)^MT-", rownames(sce))]
  # Ribosomal genes
  ribo_genes <- rownames(sce)[grep("(?i)^RP[SL]", rownames(sce))]

  # Hemoglobin genes - includes all genes starting with HB except HBP.
  hb_genes <- rownames(sce)[grep("(?i)^HB[^(P)]", rownames(sce))]

  #add QC
  sce <- scuttle::addPerCellQC(sce, flatten = T, subsets = list(mt = mito_genes, hb = hb_genes, ribo = ribo_genes))

  #rename colnames
  if( any(colnames(colData(sce))=="Cell")){
    sce@colData <- subset(sce@colData , select= -c(Cell))
  }

  colTable=colData(sce) %>% as.data.frame()  %>%
  dplyr::select(-c("total",  "subsets_mt_sum","subsets_mt_detected","subsets_hb_sum","subsets_hb_detected", "subsets_ribo_sum","subsets_ribo_detected" ))  %>%
  rename(nCount_RNA=sum, nFeature_RNA=detected, Mito_gene_percent=subsets_mt_percent, Hemoglobin_gene_percent=subsets_hb_percent, Ribosomal_gene_percent= subsets_ribo_percent)

  sce<- SingleCellExperiment(assay = list(counts = assay(sce,'counts')),
                             colData = colTable)
  sce@colData$Hemoglobin_gene_percent <- replace(sce@colData$Hemoglobin_gene_percent, is.na(sce@colData$Hemoglobin_gene_percent), 0)
  sce@colData$Mito_gene_percent <- replace(sce@colData$Mito_gene_percent, is.na(sce@colData$Mito_gene_percent), 0)
  sce@colData$Ribosomal_gene_percent <- replace(sce@colData$Ribosomal_gene_percent, is.na(sce@colData$Ribosomal_gene_percent), 0)
return(sce)

}

QC_Stats <- function(sce){

  min_UMI <- min(sce@colData$nCount_RNA)
  max_UMI <- max(sce@colData$nCount_RNA)
  mean_UMI <- mean(sce@colData$nCount_RNA)


  mean_gene <- mean(sce@colData$nFeature_RNA)
  min_gene <- min(sce@colData$nFeature_RNA)
  max_gene <- max(sce@colData$nFeature_RNA)

  min_Hemoglobin <- round(min(sce@colData$Hemoglobin_gene_percent), 2)
  max_Hemoglobin <- round(max(sce@colData$Hemoglobin_gene_percent), 2)
  mean_Hemoglobin <- round(mean(sce@colData$Hemoglobin_gene_percent), 2)

  min_mito <- round(min(sce@colData$Mito_gene_percent), 2)
  max_mito <- round(max(sce@colData$Mito_gene_percent), 2)
  mean_mito <- round(mean(sce@colData$Mito_gene_percent), 2)

  min_Ribosomal <- round(min(sce@colData$Ribosomal_gene_percent), 2)
  max_Ribosomal <- round(max(sce@colData$Ribosomal_gene_percent), 2)
  mean_Ribosomal <- round(mean(sce@colData$Ribosomal_gene_percent), 2)


  nGene_summary <- data.frame(c(min_gene, mean_gene, max_gene),
                              c(min_UMI, mean_UMI, max_UMI),
                              c(min_mito,mean_mito,max_mito),
                              c(min_Ribosomal,mean_Ribosomal, max_Ribosomal),
                              c(min_Hemoglobin, mean_Hemoglobin, max_Hemoglobin),
                              row.names = c("Min", "mean", "Max"))
  nGene_summary <- t(nGene_summary)
  row.names(nGene_summary) <- c("nGene", "nUMI","%_Mito","%_Ribosomal","%_Hemoglogin")
  return(nGene_summary)
}

#Filter raw data
filter_raw_data <- function(
    sce, #SingleCellExperiment object
    min_count_gene=3, #min times a gene is expressed
    max_count_gene=10000000, # max times a gene is expressed
    min_genes_per_cell=NULL, # min expressed genes per cell
    max_genes_per_cell=NULL, # max expressed genes per cell
    min_count_cell=NULL,  # min counts per cell
    max_count_cell=NULL,  # max counts per cell
    mt_percnt=NULL, #Mitochondrial gene percentage
    ribsoml_percnt=NULL #Ribosomal gene percentage
    #hemglbn_percnt=NULL #Hemoglobin gene percentage

)

{
  if(any(is.na(sce@colData$Patients)))
  {
    sce@colData <- subset(sce@colData, select = -c(Patients,Type,Batch, batch))
  }
  selected_f <- rownames(sce)[Matrix::rowSums(counts(sce)) >= min_count_gene & Matrix::rowSums(counts(sce)) < max_count_gene ]
  sce<-sce[selected_f]
  selected_ribo <- sce$Ribosomal_gene_percent >= ribsoml_percnt
  sce <- sce[,  selected_ribo]
  metadata<-sce %>%
    colData() %>% as.data.frame() %>% filter(nFeature_RNA>= min_genes_per_cell & nFeature_RNA< max_genes_per_cell & nCount_RNA>= min_count_cell
                                             & nCount_RNA< max_count_cell & Ribosomal_gene_percent>= ribsoml_percnt & Mito_gene_percent <= mt_percnt )
  count<- assay(sce, "counts")
  count<-count[,colnames(count) %in% rownames(metadata)]
  sce<- SingleCellExperiment(assay = list(counts = count),
                             colData = data.frame(metadata))
  sce@colData$All<-"sample"
  return(sce)

}
























