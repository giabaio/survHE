write.surv <- function(object,file,sheet=NULL,what="surv") {
  # Writes the survival summary to an excel file (helpful to then call the values in the Markov model)
  # object = a summary.flexsurvreg object containing the survival curves (with times, estimates and interval limits)
  # file = a string with the full path to the file name to be saved
  # sheet = a string with the name of the sheet to be created
  # what = the object to be exported. Possible values are:
  #        'surv' = a matrix with nsim rows and ntimes columns with the survival curve (one such matix for each configuration of the covariates)
  #        'sim' = a matrix with nsim rows and simulations for the survival parameters (scale, shape, rate, etc)
  
  # If xlsx is not installed, then request installation
  if(!isTRUE(requireNamespace("xlsx",quietly=TRUE))) {
    stop("You need to install the R package 'xlsx'. Please run in your R terminal:\n install.packages('xlsx')")
  }
  # But if it is installed, check if it's loaded and if not make its Namespace available so that all its functions are available
  if (isTRUE(requireNamespace("xlsx",quietly=TRUE))) {
    if (!is.element("xlsx", (.packages()))) {
      attachNamespace("xlsx")
    }
    # Extracts the relevant component of the make.surv output
    if(what=="surv") {
      export <- object$mat
      export.lab <- paste0(object$nsim," simulation(s) for the survival curve:")
      nobjs <- length(export)
      profile.lab <- unlist(lapply(1:nobjs,function(i) paste0("[[",i,"]] = ",rownames(object$des.mat)[i],"\n")))
      dims <- dim(export[[1]])
      # Finds the total number of rows necessary to write the simulations to the output file
      tot.rows <- dims[1]*nobjs + nobjs
      cn <- as.character(object$S[[1]][[1]][,1])
      for (i in 1:nobjs) {
        colnames(export[[i]]) <- paste0("t_",cn)
      }
    }
    
    # # Gets the extension of the file --- then decides whether to do write.csv or xlsx (do we want more formats??)
    # exts <- tools::file_ext(file)
    # if(exts=="csv") {
    #   out=export[[1]]
    #   if (nobjs>1) {
    #     for (i in 2:nobjs) {
    #       out <- rbind(out,rep(NA,ncol(out)),export[[i]])
    #     }
    #   }
    #   write.csv(out,file=file)
    # }
    
    if(is.null(sheet)) {sheet <- "Sheet 1"}
    
    # If it already exists, we need to append the data to a different sheet
    if (file.exists(file)) {
      wb <- xlsx::loadWorkbook(file)
      # If worksheet already exists needs to replace it & overwrite it
      if (sheet %in% names(xlsx::getSheets(wb))) {xlsx::removeSheet(wb,sheetName=sheet)}
      sheet <- xlsx::createSheet(wb,sheet)
      sr <- seq(from=1,by=(dims[1]+2),to=tot.rows)
      ex <- lapply(1:nobjs,function(i) xlsx::addDataFrame(export[[i]],sheet=sheet,startRow=sr[i],startColumn=1,row.names=T,col.names=T))
      # sc <- seq(from=1,by=5,length.out=nobjs)
      # for (i in 1:nobjs) {
      #   xlsx::addDataFrame(export[[i]],sheet=sheet,startRow=1,startColumn=sc[i],row.names=F)
      # }
      xlsx::saveWorkbook(wb,file)
    }
    
    # But if file does not exist, then create it
    if (!file.exists(file)) {
      exts <- tools::file_ext(file)
      ## Should put some restriction as to what file extensions we want?
      wb <- xlsx::createWorkbook(type=exts)
      sheet <- xlsx::createSheet(wb,sheet)
      sr <- seq(from=1,by=(dims[1]+2),to=tot.rows)
      ex <- lapply(1:nobjs,function(i) xlsx::addDataFrame(export[[i]],sheet=sheet,startRow=sr[i],startColumn=1,row.names=T,col.names=T))
      # sc <- seq(from=1,by=5,length.out=nobjs)
      # for (i in 1:nobjs) {
      #   xlsx::addDataFrame(object[[i]],sheet=sheet,startRow=1,startColumn=sc[i],row.names=F)
      # }
      xlsx::saveWorkbook(wb,file)
    }
    
    msg <- paste0("Written to file: ",file)
    cat(export.lab,"\n",profile.lab,"\n",msg,"\n")
  }
}
