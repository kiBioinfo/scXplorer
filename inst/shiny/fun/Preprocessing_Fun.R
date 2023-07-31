
################################################################################
################ METHODS FOR DATA INPUTS AND FILTERING##########################
################################################################################

#' A function to read count matrix and convert to a SingleCellExperiment object
#' and will add official gene symbol if not present
#' Only support Human and mouse data
#'
#' @param path an absolute path to the file
#' @inputs txt/csv || txt.gz/csv.gz || .rds file
#' @return a SingleCellExperiment object with gene symbol
#' @examples
#'
#' load_sce_obj("~/Desktop/files/counts.csv")
#'
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

#' A function to convert Ensemble gene symbol to HUGO gene symbol
#'
#' Only support Human and mouse data
#'
#' @param obj an S4 object or a dataframe containing row names Ensemble gene symbol
#'
#' @return annotated file
#' @examples
#'
#' convert_to_Symbol(obj = obj)
#'
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

#' A function to create an annotation dataframe
#'
#'
#' @param current_data a dataframe containing annotation about the experiment
#' if current_data==NULL then the function will use annotation_lab.txt file from www folder
#'
#' @return formatted dataframe
#' @examples
#'
#' metadata_modal(current_data = dataframe) || metadata_modal(current_data==NULL)
#'

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
    #if no information is available then fill the dataframe with NA
    if(ncol(d)==1){
      d$Patients=NA
      d$Type=NA
      d$Batch=NA
    }
  }
  d<-as.data.frame(d)
  return(d)
}


#' A function to add QC matrix to the SingleCellExperiment object
#'
#'
#' @param sce & metadata  a SingleCellExperiment object and the path to the metadata.txt
#' file containing annotation about the experiment
#' if metadat file not available pass only the SingleCellExperiment object
#'
#'
#' @return SingleCellExperiment object with qc matrix information
#' @examples
#'
#' add_QC_matrix(sce, metadata= '~/Desktop/files/metadata.txt') || add_QC_matrix(sce, metadata=NULL)
#'

add_QC_matrix<- function(sce, metadata=NULL){
  if(is.null(metadata)){
    count<- assay(sce, "counts")
    sce<- SingleCellExperiment(assay = list(counts = count))
  }
  else{
    #Read the metadata file which is in .txt format
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
  #Replace NA with 0
  sce@colData$Hemoglobin_gene_percent <- replace(sce@colData$Hemoglobin_gene_percent, is.na(sce@colData$Hemoglobin_gene_percent), 0)
  sce@colData$Mito_gene_percent <- replace(sce@colData$Mito_gene_percent, is.na(sce@colData$Mito_gene_percent), 0)
  sce@colData$Ribosomal_gene_percent <- replace(sce@colData$Ribosomal_gene_percent, is.na(sce@colData$Ribosomal_gene_percent), 0)
  return(sce)

}


#' A function to create a table of QC matrix from a SingleCellExperiment object
#'
#'
#' @param sce  a SingleCellExperiment object
#' object containing  QC matrix information
#'
#'
#' @return a  table with qc matrix information
#' @examples
#'
#' QC_Stats(sce)
#'
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
#' A function to Filter raw data from a SingleCellExperiment object
#'
#'
#' @param sce  a SingleCellExperiment object as input and other parameters to filter raw data
#'
#'
#' @return a QC filtered SingleCellExperiment object
#' @examples
#'
#' filter_raw_data(sce,min_count_gene=3, max_count_gene=10000,min_genes_per_cell=100,
#max_genes_per_cell=10000,min_count_cell=3,max_count_cell=100000,mt_percnt=0.5,ribsoml_percnt=0.5 #Ribosomal gene percentage
#hemglbn_percnt=NULL #Hemoglobin gene percentage
#)
#'
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

#QC filter based on number of genes, count, mitochondrial gene%
#' A function to filter raw data
#'
#' @param sce a SingleCellExperiment object
#' @param totalCount minimum counts
#' @param totalGene minimum genes
#' @param MT mitochondrial %
#' @return  a SingleCellExperiment object with filtered data
#' @examples
#'
#' Filter_Matrix(sce)
#'
Filter_Matrix = function(sce = NULL, matrix = NULL, totalCount = 200, totalGene = 200, geneCapt = 2, MT = 0.1, cutOff = 1) {
  if(!is.null(sce)) {
    print("Filtering by S4!")
    edata = assay(sce, 'counts')
    batch = sce@colData$batch
    pheno = colData(sce)

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

    sce = SingleCellExperiment(assay = list(counts = edata, logcounts = log2(edata + 1)),
                               colData = pheno,
                               metadata = list(study = Study.Name))

    return(sce)

  } else {

    matrix = as.matrix(matrix)
    totalCount.num = apply(matrix, 2, sum)
    totalGene.num  = apply(matrix >= cutOff, 2, sum)
    geneCapt.num   = apply(matrix >= 1, 1, sum)

    matrix = matrix[geneCapt.num >= geneCapt, ((totalCount.num >= totalCount) & (totalGene.num >= totalGene))]

    return(matrix)
  }
}


################################################################################
########################## VISUALIZATION #######################################
################################################################################
#' A function Plot raw data
#'
#' @param sce a SingleCellExperiment object
#' @param colour_by column name in colData
#' @return a plots
#' @examples
#'
#' Plot_raw_data(sce, colour_by="All")
#'
#Plot raw data
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

#' A function Plot filtered data
#'
#' @param sce a SingleCellExperiment object
#' @param colour_by column name in colData
#' @return  plots
#' @examples
#'
#' plot_filterd_data(sce, colour_by="All")
#'
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




#plot total counts in all the samples
#' A function to Plot total counts in all the samples
#'
#' @param sce a SingleCellExperiment object
#' @param used count matrix
#' @param by column name in colData
#'
#' @return  plots
#' @examples
#'
#' Plot_TotalCount(sce, used = 'counts',by = colnames(pheno)[1] )
#'
Plot_TotalCount = function(sce = NULL, used = 'counts', by = NULL, outIndex = getwd(), width = 4, height = 4.2, savePlot = T) {

  edata = assay(sce, used)
  batch = sce@colData$batch
  pheno = colData(sce)

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

#plot total Genes in all the samples

#' A function to Plot total Genes in all the samples
#'
#' @param sce a SingleCellExperiment object
#' @param used count matrix
#' @param by column name in colData
#'
#' @return  plots
#' @examples
#'
#' Plot_TotalGene(sce, used = 'counts',by = colnames(pheno)[1] )
#'
Plot_TotalGene = function(sce = NULL, used = 'counts', by = NULL, cutOff = 1, outIndex = getwd(), width = 4, height = 4.2, savePlot = T) {

  edata = assay(sce, used)
  batch = sce@colData$batch
  pheno = colData(sce)

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

#plot mean counts in all the samples
#' A function to Plot mean counts in all the samples
#'
#' @param sce a SingleCellExperiment object
#' @param used count matrix
#' @param by column name in colData
#'
#' @return  plots
#' @examples
#'
#' Plot_MeanCount(sce, used = 'counts',by = colnames(pheno)[1] )
#'
Plot_MeanCount = function(sce = NULL, used = 'counts', by = NULL, outIndex = getwd(), width = 4, height = 4.2, savePlot = T) {

  edata = assay(sce, used)
  batch = sce@colData$batch
  pheno = colData(sce)

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




















