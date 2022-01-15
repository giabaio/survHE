## SET OF UTILITY FUNCTIONS TO INCLUDE SURVIVAL ANALYSIS RESULTS INTO A HEALTH ECONOMIC MODEL
## Gianluca Baio + Will Browne + Peter Konings (10 Jan 2017)
#' Fit parametric survival analysis for health economic evaluations
#' 
#' Runs the survival analysis with several useful options, using either MLE
#' (via flexsurv) or a Bayesian approach (via R-INLA or rstan)
#' 
#' On object in the class \code{survHE} containing the following elements
#' 
#' @param formula a formula specifying the model to be used, in the form
#' \code{Surv(time,event)~treatment[+covariates]} for flexsurv, or
#' \code{inla.surv(time,event)~treatment[+covariates]} for INLA
#' @param data A data frame containing the data to be used for the analysis.
#' This must contain data for the 'event' variable. In case there is no
#' censoring, then \code{event} is a column of 1s.
#' @param distr a (vector of) string(s) containing the name(s) of the model(s)
#' to be fitted.  Available options are:
#' 
#' \code{flexsurv}:
#' "exponential","gamma","genf","gengamma","gompertz","weibull",
#' "weibullPH","loglogistic","lognormal" \code{INLA}:
#' "exponential","weibull","lognormal","loglogistic" \code{hmc}:
#' "exponential","gamma","genf","gengamma","gompertz","weibull","weibullPH",
#' "loglogistic","lognormal"
#' @param method A string specifying the inferential method (\code{'mle'},
#' \code{'inla'} or \code{'hmc'}). If \code{method} is set to \code{'hmc'},
#' then \code{survHE} will write suitable model code in the Stan language
#' (according to the specified distribution), prepare data and initial values
#' and then run the model.
#' @param \dots Additional options (for INLA or HMC).
#' 
#' **INLA** specific options \code{dz} = defines the step length for the grid
#' search over the hyperparameters space (default = 0.1) \code{diff.logdens} =
#' defines the difference in the log-density for the hyperparameters to stop
#' integration (default = 5) \code{control.fixed} = defines the default for the
#' priors, unless specified by the user.  Default values are prior mean = 0 for
#' *all* fixed effects prior var = 1000 for *all* fixed effects prior mean = 0
#' for the intercept prior prec -> 0 for the intercept \code{control.family} =
#' a list of options. If distr is a vector, then can be provided as a named
#' list of options, for example something like this:
#' \code{control.family=list(weibull=list(param=c(.1,.1)),lognormal=list(initial=2))}
#' the names of the elements of the list need to be the same as those given in
#' the vector \code{distr}
#' 
#' **HMC** specific options \code{chains} = number of chains to run in the HMC
#' (default = 2) \code{iter} = total number of iterations (default = 2000)
#' \code{warmup} = number of warmup iterations (default = iter/2) \code{thin} =
#' number of thinning (default = 1) \code{control} = a list specifying
#' Stan-related options, eg \code{control=list(adapt_delta=0.85)} (default =
#' NULL) \code{seed} = the random seed (to make things replicable) \code{pars}
#' = a vector of parameters (string, default = NA) \code{include} = a logical
#' indicator (if FALSE, then the pars are not saved; default = TRUE)
#' \code{priors} = a list (of lists) specifying the values for the parameters
#' of the prior distributions in the models \code{save.stan} = a logical
#' indicator (default = FALSE). If TRUE, then saves the data list for Stan and
#' the model file(s)
#' @return \item{models}{ A list containing the fitted models. These contain
#' the output from the original inference engine (\code{flexsurv}, \code{INLA}
#' or \code{rstan}). Can be processed using the methods specific to the
#' original packages, or via \code{survHE}-specific methods (such as
#' \code{plot}, \code{print}) or other specialised functions (eg to extrapolate
#' the survival curves, etc). } \item{model.fitting}{ A list containing the
#' output of the model-fit statistics (AIC, BIC, DIC). The AIC and BIC are
#' estimated for all methods, while the DIC is only estimated when using
#' Bayesian inference. } \item{method}{ A string indicating the method used to
#' fit the model, ie \code{'mle'}, \code{'inla'} or \code{'hmc'}.  }
#' \item{misc}{ A list containing the time needed to run the model(s) (in
#' seconds), the formula used, the results of the Kaplan-Meier analysis (which
#' is automatically performed using \code{npsurv}) and the original data frame.
#' }
#' @author Gianluca Baio
#' @seealso \code{make.surv}
#' @template refs
#' @keywords Parametric survival models Bayesian inference via Hamiltonian
#' Monte Carlo Bayesian inference via Integrated Nested Laplace Approximation
#' @examples
#' \dontrun{
#' # Loads an example dataset from 'flexsurv'
#' data(bc)
#' 
#' # Fits the same model using the 3 inference methods
#' mle = fit.models(formula=Surv(recyrs,censrec)~group,data=bc,
#'     distr="exp",method="mle")
#' inla = fit.models(formula=Surv(recyrs,censrec)~group,data=bc,
#'     distr="exp",method="inla")
#' hmc = fit.models(formula=Surv(recyrs,censrec)~group,data=bc,
#'     distr="exp",method="hmc")
#'     
#' # Prints the results in comparable fashion using the survHE method
#' print(mle)
#' print(inla)
#' print(hmc)
#' 
#' # Or visualises the results using the original packages methods
#' print(mle,original=TRUE)
#' print(inla,original=TRUE)
#' print(hmc,original=TRUE)
#' 
#' # Plots the survival curves and estimates
#' plot(mle)
#' plot(mle,inla,hmc,labs=c("MLE","INLA","HMC"),colors=c("black","red","blue"))
#' }
#' 
#' @export fit.models
fit.models <- function(formula = NULL, data , distr = NULL, method = "mle", ...) {
  # Captures the call
  call=match.call()
  
  # Lists all the additional inputs
  exArgs <- list(...)
  # Adds the 'formula' to exArgs, so it can be used by 'runHMC' and 'runINLA'
  exArgs$formula <- formula
  # Adds the 'data' to exArgs so it can be used by 'runHMC', 'runMLE' and 'runINLA'
  exArgs$data=data
  # Adds the 'call' to exArgs
  exArgs$call=call
  
  # Avoids the 'no visible binding for global variable' error, when compiling
  #model <- NULL
  
  # Needs to specify either the formula or the list of variables!
  if(is.null(formula)) {
    stop("You need to specify a model 'formula', e.g. 'formula=Surv(time,event)~treat'")
  }
  # ensures method is lower case
  method <- tolower(method)
  # ensures method is one of "mle","inla", "mcmc"
  if(!method %in% c("hmc","inla","mle")) {
    stop("Methods available for use are 'mle', 'hmc' or 'inla'")
  }
  
  # Check whether the selected distribution(s) can be implemented with the selected method
  # (and if not, falls back to 'mle')
  method=check_distributions(method,distr)

  # MLE -----
  # If method = MLE, then fits the model(s) using flexsurvreg
  if (method=="mle") {
    # Runs the models using the helper 'runMLE' and use the helper 'format_output_fit.models 
    res <- format_output_fit.models(lapply(distr,function(x) runMLE(x,exArgs)),method,distr,formula,data)
  }
  
  # INLA -----
  # If method = INLA, then fits model(s) using inla
  if (method=="inla") {
    if (!isTRUE(requireNamespace("survHEinla", quietly = TRUE))) {
      stop("You need to install the packages 'survHEinla'. Please run in your R terminal:\n remotes::install.github('giabaio/survHE', ref='inla')")
    }
    # If survHEinla is installed but not loaded then attach the Namespace (so that all the relevant functions are available)
    if (isTRUE(requireNamespace("survHEinla", quietly = TRUE))) {
      if (!is.element("survHEinla", (.packages()))) {
        attachNamespace("survHEinla")
      }
      res <- format_output_fit.models(lapply(distr,function(x) runINLA(x,exArgs)),method,distr,formula,data)
    }
  }
  
  # HMC -----
  if (method == "hmc") {
    if (!isTRUE(requireNamespace("survHEhmc", quietly = TRUE))) {
      stop("You need to install the packages 'survHEhmc'. Please run in your R terminal:\n remotes::install.github('giabaio/survHE', ref='hmc')")
    }
    # If survHEhmc is installed but not loaded then attach the Namespace (so that all the relevant functions are available)
    if (isTRUE(requireNamespace("survHEhmc", quietly = TRUE))) {
      if (!is.element("survHEhmc", (.packages()))) {
        attachNamespace("survHEhmc")
      }
      res <- format_output_fit.models(lapply(distr,function(x) runHMC(x,exArgs)),method,distr,formula,data)
    }
  }

  return(res)
}
