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

.mySelect_PCNum = function(sce, used = 'VGcounts', method = 'elbow') {
  if (used == 'VGcounts') {
    matrix.PJ = assay(altExp(sce), used)
    matrix.PJ = log2(matrix.PJ + 1)
  }

  if (used == 'PCA') {
    matrix.PJ = t(reducedDim(sce, used))
  }

  if (used == 'BEPCA') {
    matrix.PJ = t(reducedDim(sce, used))
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

.myTestMarker = function(sce, markers, used = 'counts', runWith = 'PCA') {

  if(used == 'counts') {
    edata = assay(sce, used)
    edata = log2(edata + 1)
  } else {
    edata = assay(sce, used)
  }

  emb = reducedDim(sce, runWith)

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
