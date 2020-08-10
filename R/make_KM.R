#' Helper function to make the Kaplan-Meier analysis of the underlying data
#' for a given formula and dataset
#' 
#' @param formula a formula specifying the model to be used, in the form
#' \code{Surv(time,event)~treatment[+covariates]} in flexsurv terms, or
#' \code{inla.surv(time,event)~treatment[+covariates]} in INLA terms.
#' @param method A string specifying the inferential method (\code{'mle'},
#' \code{'inla'} or \code{'hmc'}). If \code{method} is set to \code{'hmc'},
#' then \code{survHE} will write suitable model code in the Stan language
#' (according to the specified distribution), prepare data and initial values
#' and then run the model.
#' @param data A data frame containing the data to be used for the analysis.
#' This must contain data for the 'event' variable. In case there is no
#' censoring, then \code{event} is a column of 1s.
#' @return \item{ObjSurvfit}{A 'rms::npsurv' estimate of the KM curves}.
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso fit.models
#' @references Baio (2020). survHE
#' @keywords Kaplan-Meier estimate 
make_KM <- function(formula,data) {
  km.formula <- as.formula(gsub("inla.surv","Surv",deparse(formula)))
  # Computes the Kaplan Meier curve using the package "rms"
  ObjSurvfit <- rms::npsurv(      # Uses the function "npsurv" from the package "rms"
    formula = km.formula,         # to fit the model specified in the "formula" object
    data = data                   # to the dataset named "data"
  )
  return(ObjSurvfit)
}



# # Need to create a formula for the KM that is in the right format (flexsurv-like), 
# # so checks the current format
# if(is.na(pmatch("Surv",attributes(terms(formula))$variables[[2]][1]))) {
#   # If TRUE, then it's in inla.surv() terms. If FALSE, then formula is in Surv() terms
#   # If it's in inla.surv() needs to create a Surv() formula for the KM + check that the method isn't mle
#   km.formula <- as.formula(gsub("inla.surv","Surv",deparse(formula)))
#   # If the method was originally INLA but has been modified or if the method is mle but the formula is 
#   # in the wrong format, then need to change the formula to be consistent with the Surv() notation
#   if (method=="mle") {
#     formula <- km.formula
#   }
# } else {
#   km.formula <- formula
#   # If the method is inla but the formula is in Surv() terms, then change the formula
#   if (method=="inla") {
#     formula <- as.formula(gsub("Surv","INLA::inla.surv",deparse(formula)))
#   }
# }