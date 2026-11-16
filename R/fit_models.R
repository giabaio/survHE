## SET OF UTILITY FUNCTIONS TO INCLUDE SURVIVAL ANALYSIS RESULTS INTO A HEALTH ECONOMIC MODEL
## Gianluca Baio + Will Browne + Peter Konings (10 Jan 2017)
#' Fit parametric survival analysis for health economic evaluations
#' 
#' Runs the survival analysis with several useful options, using either MLE
#' (via flexsurv) or a Bayesian approach (via R-INLA or rstan)
#' 
#' Returns an object in the class \code{survHE}, containing the elements
#' described in the Value section below.
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
#' \code{mle} (via \code{flexsurv}):
#' "exponential","gamma","genf","gengamma","gompertz","weibull",
#' "weibullPH","loglogistic","lognormal","rps" (Royston-Parmar splines).
#' \code{inla}: "exponential","gamma","gompertz","weibull","weibullPH",
#' "loglogistic","lognormal".
#' \code{hmc}: "exponential","gamma","genf","gengamma","gompertz","weibull",
#' "weibullPH","loglogistic","lognormal","rps","polyweibull".
#' If a distribution is not available for the selected \code{method}, then
#' \code{survHE} falls back to \code{method='mle'}.
#' @param method A string specifying the inferential method (\code{'mle'},
#' \code{'inla'} or \code{'hmc'}). If \code{method} is set to \code{'hmc'},
#' then \code{survHE} will write suitable model code in the Stan language
#' (according to the specified distribution), prepare data and initial values
#' and then run the model.
#' @param \dots Additional options passed to the fitting method.
#'
#' **MLE** specific options: \code{inits} = user-defined initial values for
#' the optimisation in \code{flexsurv::flexsurvreg} (or
#' \code{flexsurv::flexsurvspline} for \code{distr="rps"}). When fitting a
#' single distribution, \code{inits} can be a numeric vector or a list
#' containing one numeric vector. When fitting several distributions,
#' \code{inits} must be a list with one element per distribution, in the same
#' order as \code{distr}; each element is either a numeric vector or
#' \code{NULL} (to use the \code{flexsurv} default for that model). Initial
#' values are given on the natural scale of the parameters, in the order used
#' by \code{flexsurv}: the baseline parameters (e.g. \code{shape} before
#' \code{scale} for the Weibull; \code{gamma0},...,\code{gamma(k+1)} for the
#' Royston-Parmar model) followed by the covariate effects. If only the baseline
#' parameters are given, the covariate effects are initialised at 0. If
#' \code{inits} is not specified, the \code{flexsurv} defaults are used for
#' all distributions. \code{weights} and \code{subset} are passed to
#' \code{flexsurv} for all distributions. For \code{distr="rps"},
#' \code{k} = number of internal knots (default = 0), \code{knots} = locations
#' of the internal knots, \code{bknots} = locations of the boundary knots,
#' \code{scale} = \code{"hazard"} (default), \code{"odds"} or \code{"normal"}
#' and \code{timescale} = \code{"log"} (default) or \code{"identity"}; see
#' \code{flexsurv::flexsurvspline} for details.
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
#' @keywords Parametric survival models Bayesian inference via Hamiltonian Monte Carlo Bayesian inference via Integrated Nested Laplace Approximation
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
#' # Uses user-defined initial values for the MLE (one set per distribution)
#' mle2 = fit.models(formula=Surv(recyrs,censrec)~group,data=bc,
#'     distr=c("exp","weibull"),method="mle",inits=list(0.05,c(1.5,10)))
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
  # If method = MLE, then fits the model(s) using flexsurvreg (or flexsurvspline for the RPS)
  if(method=="mle") {
    # Retrieves (and checks) the user-specified initial values, if any
    inits <- exArgs$inits
    if(!is.null(inits)) {
      # Allows a numeric vector when fitting a single distribution
      if(is.numeric(inits) && length(distr)==1) {inits <- list(inits)}
      # Requires one list element per distribution
      if(!is.list(inits) || length(inits)!=length(distr)) {
        stop("'inits' must be a list with one element per distribution, or a numeric vector for a single distribution")
      }
      # Each element must be a numeric vector or NULL
      if(!all(vapply(inits,function(x) is.null(x) || is.numeric(x),logical(1)))) {
        stop("Each element of 'inits' must be numeric or NULL")
      }
    }
    # Runs the models using the helper 'runMLE' (each with its own initial values) and 
    # then uses the helper 'format_output_fit.models'
    output <- lapply(seq_along(distr),function(i) {
      modelArgs <- exArgs
      if(!is.null(inits)) {modelArgs$inits <- inits[[i]]}
      runMLE(distr[[i]],modelArgs)
    })
    res <- format_output_fit.models(output,method,distr,formula,data)
  }
  
  # INLA -----
  # If method = INLA, then fits model(s) using inla
  if (method=="inla") {
    if (!isTRUE(requireNamespace("survHEinla", quietly = TRUE))) {
      stop("You need to install the packages 'survHEinla'. Please run in your R terminal:\n remotes::install_github('giabaio/survHEinla')")
    }
    # If survHEinla is installed but not loaded then attach the Namespace (so that all the relevant functions are available)
    if (isTRUE(requireNamespace("survHEinla", quietly = TRUE))) {
      if (!is.element("survHEinla", (.packages()))) {
        attachNamespace("survHEinla")
      }
      res <- format_output_fit.models(lapply(distr,function(x) survHEinla::runINLA(x,exArgs)),method,distr,formula,data)
    }
  }
  
  # HMC -----
  if (method == "hmc") {
    if (!isTRUE(requireNamespace("survHEhmc", quietly = TRUE))) {
      stop("You need to install the packages 'survHEhmc'. Please run in your R terminal:\n remotes::install_github('giabaio/survHEhmc')")
    }
    # If survHEhmc is installed but not loaded then attach the Namespace (so that all the relevant functions are available)
    if (isTRUE(requireNamespace("survHEhmc", quietly = TRUE))) {
      if (!is.element("survHEhmc", (.packages()))) {
        attachNamespace("survHEhmc")
      }
      res <- format_output_fit.models(lapply(distr,function(x) survHEhmc::runHMC(x,exArgs)),method,distr,formula,data)
    }
  }

  return(res)
}
