#' Plot survival curves for the models fitted using \code{fit.models}
#' 
#' Plots the results of model fit.
#' 
#' 
#' @param ...  Must include at least one result object saved as the call to the
#' \code{fit.models} function.  Other possibilities are additional (mainly
#' graphical) options. These are: \code{xlab} = a string with the label for the
#' x-axis (default = "time") \code{ylab} = a string with the label for the
#' y-axis (default = "Survival") \code{lab.trt} = a (vector of) string(s)
#' indicating the labels associated with the strata defining the different
#' survival curves to plot. Default to the value used by the Kaplan Meier
#' estimate given in \code{fit.models} \code{cex.trt} = factor by which the
#' size of the font used to write the strata is resized (default = 0.8)
#' \code{n.risk} = logical. If TRUE (defaults) writes the number at risk at
#' different time points (as determined by the Kaplan Meier estimate)
#' \code{newdata} = a list (of lists) providing the values for the relevant
#' covariates If NULL, then will use the mean values for the covariates if at
#' least one is a continuous variable, or the combination of the categorical
#' covariates. \code{xlim} = a vector determining the limits for the x-axis
#' \code{colors} = a vector of characters defining the colours in which to plot
#' the different survival curves \code{labs} = a vector of characters defining
#' the names of the models fitted \code{add.km} = TRUE (whether to also add the
#' Kaplan Meier estimates of the data) \code{legend} = TRUE (whether to also 
#' add the legend to the graph)
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso Something will go here
#' @references Something will go here
#' @keywords Parametric survival models
#' @examples
#' 
#' data(bc)
#' 
#' mle = fit.models(formula=Surv(recyrs,censrec)~group,data=bc,
#'     distr="exp",method="mle")
#' plot(mle)
#' 
#' @export plot.survHE
plot.survHE <- function(...) {
  exArgs=list(...)
  if(exists("graph",where=exArgs)){graph=exArgs$graph} else {graph="base"}
  if(graph=="base") {
    do.call(survHE:::plot_base_survHE,exArgs)
  }
}
