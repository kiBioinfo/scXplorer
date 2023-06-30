library(htmltools)
withBusyIndicatorUI <- function(icon_name, button) {
  id <- button[['attribs']][['id']]
  n_s <- strsplit(id, split = "-")[[1]][1]
  div(
    `data-for-btn` = id,
    button,
    span(
      class = "btn-loading-container",
      shinyjs::hidden(
        img(id = paste0(n_s, "-load_", icon_name), src = "ajax-loader-bar.gif", class = "btn-loading-indicator"),
        img(id = paste0(n_s, "-check_", icon_name), src = "checkmark.jpg", class = "btn-loading-indicator")
      )
    )
  )
}
withBusyIndicator <- function(button) {
  id <- button[['attribs']][['id']]
  tagList(
    button,
    span(
      class = "btn-loading-container",
      `data-for-btn` = id,
      hidden(
        img(src = "ajax-loader-bar.gif", class = "btn-loading-indicator"),
        icon("check", class = "btn-done-indicator")
      )
    )
  )
}

convertGeneIDs <- function(inSCE, inSymbol, outSymbol, database="org.Hs.eg.db"){
  if (!(database %in% c("org.Ag.eg.db", "org.At.tair.db", "org.Bt.eg.db",
                        "org.Ce.eg.db", "org.Cf.eg.db", "org.Dm.eg.db", "org.Dr.eg.db",
                        "org.EcK12.eg.db", "org.EcSakai.eg.db", "org.Gg.eg.db", "org.Hs.eg.db",
                        "org.Mm.eg.db", "org.Mmu.eg.db", "org.Pf.plasmo.db", "org.Pt.eg.db",
                        "org.Rn.eg.db", "org.Sc.sgd.db", "org.Ss.eg.db", "org.Xl.eg.db"))){
    stop("The database you want to use, ", database, ", is not supported")
  }
  if (!(database %in% as.character(grep("^org\\.",
                                        utils::installed.packages()[, "Package"],
                                        value = TRUE)))){
    stop("The database you want to use, ", database, ", is not installed.")
  }
  if (!(database %in% (.packages()))){
    stop("You need to load the database to use it: library(", database, ")")
  }
  if (inSymbol == outSymbol){
    message("No conversion necessary.")
    return(inSCE)
  }
  indb <- get(paste(database))
  oldids <- rownames(inSCE)
  if (any(duplicated(oldids)) | any(is.na(oldids))){
    stop("problem with input IDs, duplicates or NAs found")
  }
  res <- AnnotationDbi::select(indb, keys = oldids, columns = c(outSymbol),
                               keytype = inSymbol)
  if (any(is.na(res[, outSymbol]))){
    message(sum(is.na(res[, outSymbol])), " ", outSymbol,
            " are NA after conversion. Removing")
    res <- res[!is.na(res[, outSymbol]), ]
    if (any(is.na(res[, outSymbol]))){
      stop("problem.")
    }
  }
  if (any(duplicated(res[, outSymbol]))){
    message(sum(duplicated(res[, outSymbol])), " ", outSymbol,
            " are duplicated after conversion. Removing additional copies")
    res <- res[!duplicated(res[, outSymbol]), ]
    if (any(duplicated(res[, outSymbol]))){
      stop("problem.")
    }
  }
  if (any(duplicated(res[, inSymbol]))){
    message(sum(duplicated(res[, inSymbol])), " ", inSymbol,
            " are duplicated after conversion. Removing additional copies")
    res <- res[!duplicated(res[, inSymbol]), ]
    if (any(duplicated(res[, inSymbol]))){
      stop("problem.")
    }
  }
  message(length(oldids), " ", inSymbol, " originally, ", nrow(res), " ",
          outSymbol, "s")
  newsce <- inSCE[res[, inSymbol], ]
  rownames(newsce) <- res[, outSymbol]
  return(newsce)
} 


my_icon = function (name, class = NULL, lib = "font-awesome") {
  
  prefixes <- list(`font-awesome` = "fa", glyphicon = "glyphicon")
  prefix <- prefixes[[lib]]
  if (is.null(prefix)) {
    stop("Unknown font library '", lib, "' specified. Must be one of ", 
         paste0("\"", names(prefixes), "\"", collapse = ", "))
  }
  iconClass <- ""
  if (!is.null(name)) {
    prefix_class <- prefix
    #if (prefix_class == "fa" && name %in% font_awesome_brands) {
    #  prefix_class <- "fab"
    #}
    iconClass <- paste0(prefix_class, " ", prefix, "-", name)
  }
  if (!is.null(class)) 
    iconClass <- paste(iconClass, class)
  iconTag <- tags$i(class = iconClass)
  if (lib == "font-awesome") {
    htmlDependencies(iconTag) <- htmlDependency("font-awesome", 
                                                "6.2.1", "./www/fontawesome/V6.2.1", 
                                                stylesheet = c("css/all.css"))
  }
  htmltools::browsable(iconTag)
}

helpPopup <- function(content, title = NULL) {
  a(href = "#",
    class = "popover-link",
    `data-toggle` = "popover",
    `data-title` = title,
    `data-content` = content,
    `data-html` = "true",
    `data-trigger` = "hover",
    icon("question-circle")
  )
}

valfn <- function(x) if(is.null(x) | is.na(x) | x=="") return(warning("Input data is incorrect."))



useAutoColor <- function(input, output, session = shiny::getDefaultReactiveDomain()) {
  input <- get("input", envir = parent.frame())
  # Now we need to set up a initial theme so that 
  # session$setCurrentTheme does not complain about 
  # changing the bootstrap version
  theme <- bslib::bs_theme(version = 4)
  shiny::shinyOptions(bootstrapTheme = theme)
  # input$dark_mode is created on the client
  shiny::observeEvent(input$dark_mode, {
    # We actually don't do anything fancy. We send the same theme as
    # the initial one
    session$setCurrentTheme(theme)
  })
}

#select columns to to use as colour in plot
names=c("nCount_RNA","nFeature_RNA", "Mito_gene_percent", "Hemoglobin_gene_percent", "Ribosomal_gene_percent" , "sizeFactor")