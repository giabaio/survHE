#' write.surv
#' 
#' Writes the survival summary to an excel file (helpful to then call the
#' values in the Markov model)
#' 
#' Something will go here
#' 
#' @param object a summary.flexsurvreg object containing the survival curves
#' (with times, estimates and interval limits)
#' @param file a string with the full path to the file name to be saved
#' @param sheet a string with the name of the sheet to be created
#' @param what a string to describe what to be exported. Can either be 
#' 'surv' (default), which outputs the simulation(s) for the survival curves
#' or 'sim', which outputs the simulation(s) for the underlying model 
#' parameters. If there are several 'profiles', they get written in 
#' separate spreadsheets and a clear indication is given as the name of the
#' spreadsheet
#' @return A spreadsheet file with the simulation(s) of the relevant quantity
#' @author Gianluca Baio
#' @seealso \code{make.surv}
#' @template refs
#' @keywords Excel PSA
#' @examples
#' \dontrun{
#' # Loads an example dataset from 'flexsurv'
#' data(bc)
#' 
#' # Fits the same model using the 3 inference methods
#' mle = fit.models(formula=Surv(recyrs,censrec)~group,data=bc,
#'     distr="exp",method="mle")
#' p.mle = make.surv(mle)
#' write.surv(p.mle,file="test.xlsx")
#' }
#' @export write.surv
write.surv <- function(object,file,sheet=NULL,what="surv") {
  # Writes the survival summary to an excel file (helpful to then call the values in the Markov model)
  # object = a summary.flexsurvreg object containing the survival curves (with times, estimates and interval limits) - that's actually
  # the result of the PSA from 'make.surv'
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
      export.lab <- paste0(" ",object$nsim," simulation(s) for the survival curve")
    } else {
      export=object$sim
      export.lab <- paste0(" ",object$nsim," simulation(s) for the model parameters")
    }
    nobjs <- length(export)
    profile.lab=lapply(1:nrow(object$des.mat),function(x){
      object$des.mat %>% as_tibble() %>% select(!contains("(Intercept)")) %>% slice(x) %>% 
        round(digits=2) %>% mutate(strata=paste0(names(.),"=",.,collapse=","))
    }) %>% bind_rows(.) %>% pull(strata)
    dims <- dim(export[[1]])
    # Finds the total number of rows necessary to write the simulations to the output file
    tot.rows <- dims[1]*nobjs + nobjs

    if(is.null(sheet)) {
      sheet = profile.lab
    }
    
    # If it already exists, we need to append the data to a different sheet
    if (file.exists(file)) {
      wb <- xlsx::loadWorkbook(file)
      # If worksheet already exists needs to replace it & overwrite it
      for (i in 1:length(sheet)){
        if(sheet[i] %in% names(xlsx::getSheets(wb))) {xlsx::removeSheet(wb,sheetName=sheet[i])}
      }
      sheet <- lapply(sheet,function(i) xlsx::createSheet(wb,i))
      ex <- lapply(1:nobjs,function(i) xlsx::addDataFrame(export[[i]],sheet=sheet[[i]],startRow=1,startColumn=1,row.names=T,col.names=T))
      xlsx::saveWorkbook(wb,file)
    }
    
    # But if file does not exist, then create it
    if (!file.exists(file)) {
      exts <- tools::file_ext(file)
      ## Should put some restriction as to what file extensions we want?
      wb <- xlsx::createWorkbook(type=exts)
      sheet <- lapply(sheet,function(i) xlsx::createSheet(wb,i))
      ex <- lapply(1:nobjs,function(i) xlsx::addDataFrame(export[[i]],sheet=sheet[[i]],startRow=1,startColumn=1,row.names=T,col.names=T))
      xlsx::saveWorkbook(wb,file)
    }
    
    msg <- paste0("written to file: '",file,"'")
    cat(export.lab,msg,"\n","Profile(s):",paste("\n  ",profile.lab),"\n")
  }
}
