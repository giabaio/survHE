#' Prints a summary table for the distribution the mean survival time for a
#' given model and data
#' 
#' Calculates the mean survival time as the area under the survival curve
#' 
#' A list comprising of the following elements
#' 
#' @param object a \code{survHE} object (resulting from the call to
#' \code{fit.models}
#' @param mod the model to be analysed (default = 1)
#' @param t the vector of times to be used in the computation. Default = NULL,
#' which means the observed times will be used. NB: the vector of times should
#' be: i) long enough so that S(t) goes to 0; and ii) dense enough so that the
#' approximation to the AUC is sufficiently precise
#' @param nsim the number of simulations from the survival curve distributions
#' to be used (to compute interval estimates)
#' @param \dots Additional options
#' @return \item{mean.surv}{ A matrix with the simulated values for the mean
#' survival times } \item{tab}{ A summary table }
#' @author Gianluca Baio
#' @seealso \code{fit.models}, \code{make.surv}
#' @template refs
#' @keywords Parametric survival models Mean survival time
#' @examples
#' 
#' data(bc)
#' 
#' mle = fit.models(formula=Surv(recyrs,censrec)~group,data=bc,
#'     distr="exp",method="mle")
#' summary(mle,nsim=100)
#' 
#' @export summary.survHE
summary.survHE <- function(object,mod=1,t=NULL,nsim=1000,...) {
  # Computes the estimated mean survival as the area under the survival curve
  # This is obtained using the trapezoidal method by taking the average of the "left" and "right" y-values.
  # object: is the output from a fit.models call
  # mod: the model to be analysed (default = 1)
  # t: the vector of times to be used in the computation. Default = NULL, which means the observed times will be used.
  #     NB: the vector of times should be: i) long enough so that S(t) goes to 0; and ii) dense enough so that
  #         the approximation to the AUC is sufficiently precise
  # nsim: number of simulations from the survival curve distributions to be used (to compute interval estimates)
  # stats: a logical value. If TRUE, also shows a table 
  # ...: optional arguments
  # newdata = a list (of lists), specifiying the values of the covariates at which the computation is performed. For example
  #           'list(list(arm=0),list(arm=1))' will create two survival curves, one obtained by setting the covariate 'arm'
  #           to the value 0 and the other by setting it to the value 1. In line with 'flexsurv' notation, the user needs
  #           to either specify the value for *all* the covariates or for none (in which case, 'newdata=NULL', which is the
  #           default). If some value is specified and at least one of the covariates is continuous, then a single survival
  #           curve will be computed in correspondence of the average values of all the covariates (including the factors, 
  #           which in this case are expanded into indicators). The order of the variables in the list *must* be the same
  #           as in the formula used for the model
  # labs: a vector of strings giving the names of the "profile" of covariates for which the mean survival times are computed
  #
  # NB: NEED TO FIX THIS FOR THE POLY-WEIBULL
  #
  # Defines the utility function to compute the stats table
  make.stats <- function(x, dim = 2) {
    bugs.stats <- function(x) {
      c(mean(x), sd(x), quantile(x, 0.025), median(x), quantile(x, 0.975))
    }
    if (is.null(dim(x)) == TRUE) {
      tab <- bugs.stats(x)
      names(tab) <- c("mean", "sd", "2.5%", "median", "97.5%")
    }
    if (is.null(dim(x)) == FALSE) {
      tab <- t(apply(x, dim, function(x) bugs.stats(x)))
      colnames(tab) <- c("mean", "sd", "2.5%", "median", "97.5%")
    }
    return(tab)
  }
  
  exArgs <- list(...)
  if (!exists("newdata",where=exArgs)) {newdata <- NULL} else {newdata <- exArgs$newdata}
  if (!exists("labs",where=exArgs)) {labs <- NULL} else {labs <- exArgs$labs}
  if(is.null(t)) {
    if(object$misc$model_name[mod]=="pow") {
      t <- sort(unique(object$misc$km[[mod]]$time))
    } else {
      t <- sort(unique(object$misc$km$time))
    }
  }
  
  psa <- make.surv(object,mod=mod,t=t,nsim=nsim,newdata=newdata)
  rlabs <- rownames(psa$des.mat)
  if (!is.null(rlabs)) {
    rlabs <- gsub("^1,","",rlabs)
  } else {
    rlabs <- rep("",length(psa$sim))
  }
  if(!is.null(labs) & length(labs)==length(rlabs)) {rlabs <- labs}
  
  mean.surv=matrix(unlist(
    lapply(psa$mat,function(i) {
      lapply(1:psa$nsim,function(j) {
        xvar=i$t
        yvar=i[,(j+1)]
        sum(diff(xvar) * (head(yvar,-1)+tail(yvar,-1)), na.rm=T)/2
      })
    })
  ),nrow=psa$nsim,byrow=F)

  if (ncol(mean.surv)==length(names(object$misc$km$strata))) {
    colnames(mean.surv) <- names(object$misc$km$strata)
  }
  
  tab <- NULL
  if(psa$nsim>1) {
    tab <- make.stats(mean.surv)
    rownames(tab) <- rlabs
    if(!is.null(names(object$misc$km$strata))) {
      if (ncol(mean.surv)==length(names(object$misc$km$strata))) {
        rownames(tab) <- names(object$misc$km$strata)
      } else {
        rownames(tab) <- rlabs
      }
    } else {
      rownames(tab) <- rlabs
    }
    cat("\nEstimated average survival time distribution* \n")
    print(tab)
    cat(paste0("\n*Computed over the range: [",paste(format(range(t),digits=4,nsmall=3),collapse="-"),"] using ",psa$nsim," simulations.\nNB: Check that the survival curves tend to 0 over this range!\n"))
  }
  return(invisible(list(mean.surv=mean.surv,tab=tab)))
}
