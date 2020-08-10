#' Helper function to run the survival models using MLE and flexsurv
#' 
#' @param x a (vector of) string(s) containing the name(s) of the model(s)
#' to be fitted
#' #' @param exArgs a list of extra arguments passed from the main 'fit.models' 
#' function
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso fit.models
#' @references Baio (2020). survHE
#' @keywords Parametric survival models Maximum likelihood estimation
runMLE <- function(x,exArgs) {
  ##### PROBABLY CAN REMOVE THIS ######
  # Checks that 'flexsurv' is loaded up. NB: ***Probably*** not needed, as 'flexsurv' is a primary dependency???
  if(!isTRUE(requireNamespace("flexsurv",quietly=TRUE))) {
    stop("You need to install the R package 'flexsurv'. Please run in your R terminal:\n install.packages('flexsurv')")
  }
  # Loads in the available models in each method
  availables <- load_availables()
  
  tic <- proc.time()
  # If user selects RPS model, then could also provide some optional arguments - uses flexsurv defaults
  if(x=="rps") {
    if(exists("bhazard",where=exArgs)) {bhazard <- exArgs$bhazard} else {bhazard <-NULL}
    if(exists("weights",where=exArgs)) {weights <- exArgs$weights} else {weights <- NULL}
    if(exists("subset",where=exArgs)) {subset <- exArgs$subset} else {subset <- NULL}
    if(exists("knots",where=exArgs)) {knots <- exArgs$knots} else {knots <- NULL}
    if(exists("k",where=exArgs)) {k <- exArgs$k} else {k <- 0}
    if(exists("bknots",where=exArgs)) {bknots <- exArgs$bknots} else {bknots <- NULL}
    if(exists("scale",where=exArgs)) {scale <- exArgs$scale} else {scale <- "hazard"}
    if(exists("timescale",where=exArgs)) {timescale <- exArgs$scale} else {timescale <- "log"}
    model <- flexsurv::flexsurvspline(formula=formula,data=data,k=k,knots=knots,bknots=bknots,scale=scale,
                                      timescale=timescale)
  } else {
  # If it's one of the other available models under MLE, then simply runs flexsurv::flexsurvreg
    model <- flexsurv::flexsurvreg(formula=formula,data=data,dist=x)
  }
  toc <- proc.time()-tic

  # Finally returns the output
  list(
    model=model,
    aic=model$AIC,
    bic=-2*model$loglik+model$npars*log(model$N),
    dic=NULL,
    time2run=toc[3]
  )
}
