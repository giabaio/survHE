#' Plot survival curves for the models fitted using \code{fit.models}
#' 
#' Plots the results of model fit.
#' 
#' @param ... Must include at least one result object saved as 
#' the call to the \code{fit.models} function. Nay include other 
#' optional parameters. These include whether the KM curve should be 
#' added \code{add.km} and whether the user specifies a profile of covariates 
#' (in the list \code{newdata}). Other possibilities are additional 
#' (mainly graphical) options. These are: \code{xlab} = a string with the label 
#' for the x-axis (default = "time") \code{ylab} = a string with the label for the
#' y-axis (default = "Survival") \code{lab.profile} = a (vector of) string(s)
#' indicating the labels associated with the strata defining the different
#' survival curves to plot. Default to the value used by the Kaplan Meier
#' estimate given in \code{fit.models}. \code{newdata} = a list (of lists) 
#' providing the values for the relevant covariates If NULL, then will use 
#' the mean values for the covariates if at least one is a continuous variable, 
#' or the combination of the categorical covariates. \code{xlim} = a vector 
#' determining the limits for the x-axis \code{colors} = a vector of characters 
#' defining the colours in which to plot the different survival curves 
#' \code{lab.profile} = a vector of characters defining the names of the models fitted 
#' \code{add.km} = TRUE (whether to also add the Kaplan Meier estimates of the data) 
#' \code{annotate} = FALSE (whether to also add text to highlight the observed vs
#' extrapolated data)
#' \code{legend.position} = a vector of proportions to place the legend. Default
#' to 'c(.75,.9)', which means 75% across the x-axis and 90% across the y-axis
#' \code{legend.title} = suitable instructions to format the title of the legend;
#' defaults to 'element_text(size=15,face="bold")' but there may be other 
#' arguments that can be added (using 'ggplot' facilities)
#' \code{legend.text} = suitable instructions to format the text of the legend;
#' defaults to 'element_text(colour="black", size=14, face="plain")' but there 
#' may be other arguments that can be added (using 'ggplot' facilities)
#' @author Gianluca Baio
#' @seealso \code{fit.models}, \code{write.surv}
#' @template refs
#' @keywords Parametric survival models
#' @examples
#' \dontrun{
#' data(bc)
#' 
#' mle = fit.models(formula=Surv(recyrs,censrec)~group,data=bc,
#'     distr="exp",method="mle")
#' inla = fit.models(formula=Surv(recyrs,censrec)~group,data=bc,
#'     distr="exp",method="inla")
#' plot(MLE=mle,INLA=inla)
#' }
#' 
#' @export plot.survHE
plot.survHE <- function(...) {
  
  # Collects all the extra arguments
  exArgs=list(...)
  
  # Finds out whether there are objects with no name (if so, they will be 'survHE' objects!)
  # If there are any, then needs to rename them to make the rest of the function work
  if(length(names(exArgs))==0) {
    # This is the case where the only argument(s) is/are unnamed 'survHE' object(s)
    names(exArgs)=paste0("Object",1:length(exArgs))
  }
  if(length(which(names(exArgs)==""))>0){
    names(exArgs)[which(names(exArgs)=="")] = paste0("Object",1:length(which(names(exArgs)=="")))
  }
  
  # The default is to go with the 'ggplot' version of the graph. 
  if (exists("graph",exArgs)) {graph=exArgs$graph} else {graph="ggplot"}

  # If so, then call the function 'plot_ggplot_survHE
  if(graph=="ggplot") {
    return(plot_ggplot_survHE(exArgs))
  }

  # If the user selects 'base' (only for back-compatibility), then runs the old code
  ### NB: Do I want this? (probably not...)
  if(graph=="base") {
    do.call(plot_base_survHE,exArgs)
  }
}
