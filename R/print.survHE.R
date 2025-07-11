#' Print a summary of the survival model(s) fitted by \code{fit.models}
#' 
#' Prints the summary table for the model(s) fitted, with the estimate of the
#' parameters
#' 
#' 
#' @param x the \code{survHE} object (the output of the call to
#' \code{fit.models})
#' @param mod is the index of the model. Default value is 1, but the user can
#' choose which model fit to visualise, if the call to fit.models has a vector
#' argument for distr (so many models are fitted & stored in the same object)
#' @param \dots additional options, including: \code{digits} = number of
#' significant digits to be shown in the summary table (default = 6);
#' \code{original} = a flag to say whether the *original* table
#' from either \code{flexsurv} or \code{INLA} or \code{rstan} should be printed;
#' \code{print_priors} = a flag to sy whether the Stan table should also 
#' include a summary of the prior used in the HMC model
#' 
#' @author Gianluca Baio
#' @template refs
#' @keywords Parametric survival models
#' @examples
#' \dontrun{ 
#' data(bc)
#' 
#' mle = fit.models(formula=Surv(recyrs,censrec)~group,data=bc,
#'     distr="exp",method="mle")
#' print(mle)
#' }
#' 
#' @export print.survHE
print.survHE <- function(x,mod=1,...) {
  # Creates a print method for the objects in the class survHE
  # x is the survHE object (the output of the call to fit.models)
  # mod is the index of the model. Default value is 1, but the user can choose which model fit to visualise, 
  #     if the call to fit.models has a vector argument for distr (so many models are fitted & stored in the same object)
  # ... optional arguments
  # digits = number of *significant* digits to be shown in the summary table (default = 6)
  # original = a flag to say whether the *original* table from either INLA or MCMC should be printed
  # print_priors = a flag to say whether the Stan table should also print a summary of the prior distributions
  
  exArgs <- list(...)

  # Loads available models
  availables <- load_availables()
  
  # Can select the number of digits to be printed in the output table
  if(!exists("digits",where=exArgs)){digits=6} else {digits=exArgs$digits}
  if(!exists("original",where=exArgs)){original=FALSE} else {original=exArgs$original}
  # Aliases for 'original'
  if(exists("orig",exArgs)){original=exArgs$orig}
  
  # Now computes the stats, using different helpers depending on the underlying method
  # Can ask for the original output from either 'flexsurv', 'inla' or 'rstan'
  if(original==TRUE) {
    do.call(
      paste0("original_table_",x$method),
      args=list(x,mod,digits,...)
    )
  } # If not, go with the default formatting using the standardised 'survHE' output
  else {
    # First make the results table using the helper functions 
    res=do.call(
      paste0("get_stats_",x$method),
      args=list(x,mod)
    )
    # Now formats the table
    format_table(x,mod,res,digits)
  }
}

