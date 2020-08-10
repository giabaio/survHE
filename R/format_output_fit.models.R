#' Helper function to format the output of the modelling (produced either
#' by running 'runMLE', or 'runINLA', 'runHMC'), in a way that is consistent
#' with the architecture of 'survHE'
#' 
#' @param output The output of one of the helper functions used to run the
#' models.
#' @return \item{res}{A 'survHE' object containing all the relevant output
#' conveniently formatted}.
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso fit.models
#' @references Baio (2020). survHE
#' @keywords Parametric survival models Bayesian inference via Hamiltonian
#' Monte Carlo Bayesian inference via Integrated Nested Laplace Approximation
format_output_fit.models <- function(output) {
  # Model output
  models <- lapply(output, function(x) x$model)
  # Model fitting statistics
  model.fitting <- list(
    aic=unlist(lapply(output,function(x) x$aic)),
    bic=unlist(lapply(output,function(x) x$bic)),
    dic=unlist(lapply(output,function(x) x$dic))
  )
  # Miscellanea
  misc <- list(
    time2run= unlist(lapply(output, function(x) x$time2run)),
    formula=formula,
    km=km,
    data=data
  )
  
  # HMC-specific extra output
  if(method=="hmc"){
    # Completes the 'misc' and 'model.fitting' lists with additional output
    misc$data.stan <- lapply(output,function(x) x$data.stan)
    model.fitting$dic2 <- unlist(lapply(output,function(x) x$dic2))
    # Also, if 'save.stan' is set to TRUE, then saves the Stan model file(s) & data
    if (save.stan == TRUE) {
      write_model <- lapply(output,function(i) {
        model_code <- attr(model[[i]]$out@stanmodel,"model_code")
        con <- paste0(distr[i],".stan")
        writeLines(model_code, con = con)
        txt <- paste0("Model code saved to the file: ",con,"\n")
        cat(txt)
      })
    }
  }
  
  # Names the elements of the list
  names(models) <- labs

  # Formats all output in a list
  res <- list(models=models,model.fitting=model.fitting,method=method,misc=misc)
  # And sets its class attribute to "survHE"
  class(res) <- "survHE"
  return(res)
}