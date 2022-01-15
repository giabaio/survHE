#' Helper function to run the survival models using MLE and flexsurv
#' 
#' @param x a string containing the name of the model
#' to be fitted
#' #' @param exArgs a list of extra arguments passed from the main 'fit.models' 
#' function
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso fit.models
#' @references Baio (2020). survHE
#' @keywords Parametric survival models Maximum likelihood estimation
#' @noRd 
runMLE <- function(x,exArgs) {
  ##### PROBABLY CAN REMOVE THIS ######
  # Checks that 'flexsurv' is loaded up. NB: ***Probably*** not needed, as 'flexsurv' is a primary dependency???
  #if(!isTRUE(requireNamespace("flexsurv",quietly=TRUE))) {
  #  stop("You need to install the R package 'flexsurv'. Please run in your R terminal:\n install.packages('flexsurv')")
  #}
  # Loads the model formula & data
  formula <- exArgs$formula
  data=exArgs$data
  
  # Loads in the available models in each method
  availables <- load_availables()
  # Uses the helper 'manipulated_distributions' to create the vectors distr, distr3 and labs
  d3 <- manipulate_distributions(x)$distr3
  x <- manipulate_distributions(x)$distr
  
  tic <- proc.time()
  # If user selects RPS model, then could also provide some optional arguments - uses flexsurv defaults
  if(x=="survspline") {
    if(exists("bhazard",where=exArgs)) {bhazard <- exArgs$bhazard} else {bhazard <-NULL}
    if(exists("weights",where=exArgs)) {weights <- exArgs$weights} else {weights <- NULL}
    if(exists("subset",where=exArgs)) {subset <- exArgs$subset} else {subset <- NULL}
    if(exists("knots",where=exArgs)) {knots <- exArgs$knots} else {knots <- NULL}
    if(exists("k",where=exArgs)) {k <- exArgs$k} else {k <- 0}
    if(exists("bknots",where=exArgs)) {bknots <- exArgs$bknots} else {bknots <- NULL}
    if(exists("scale",where=exArgs)) {scale <- exArgs$scale} else {scale <- "hazard"}
    if(exists("timescale",where=exArgs)) {timescale <- exArgs$scale} else {timescale <- "log"}
    model <- do.call(
      flexsurv::flexsurvspline,
      args=list(
        formula=formula,
        data=quote(data),
        k=k,
        knots=knots,
        bknots=bknots,
        scale=scale,
        timescale=timescale,
        weights=weights,
        subset=subset
      )
    )
    ###model <- flexsurv::flexsurvspline(formula=formula,data=data,k=k,knots=knots,bknots=bknots,scale=scale,
    ###                                  timescale=timescale)
    if(!is.null(weights)) {model$call$weights=exArgs$call$weights} else {model$call$weights=NULL}
    if(!is.null(subset)) {model$call$subset=exArgs$call$subset} else {model$call$subset=NULL}
  } else {
    # If it's one of the other available models under MLE, then simply runs flexsurv::flexsurvreg
    # But allows to use weight and subset
    if(exists("weights",where=exArgs)) {weights <- exArgs$weights} else {weights=NULL}
    if(exists("subset",where=exArgs)) {subset <- exArgs$subset} else {subset <- NULL}

    ###model <- flexsurv::flexsurvreg(formula=formula,data=data,dist=x,weights=weights)
    model <- do.call("flexsurvreg",args=list(
      formula=formula,
      data=quote(data),
      dist=x,
      weights=weights,
      subset=subset
    ))
    # Fix the 'call' by using parts of the original call to fit.models
    if(!is.null(weights)) {model$call$weights=exArgs$call$weights} else {model$call$weights=NULL}
    if(!is.null(subset)) {model$call$subset=exArgs$call$subset} else {model$call$subset=NULL}
  }
  toc <- proc.time()-tic
  
  # Replaces a field used in 'make.surv' to indicate the model used (standardised across models)
  model_name <- d3
  
  # Finally returns the output
  list(
    model=model,
    aic=model$AIC,
    bic=-2*model$loglik+model$npars*log(model$N),
    dic=NULL,
    time2run=toc[3],
    model_name=model_name
  )
}



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
#' @noRd 
make_KM <- function(formula,data) {
  km.formula <- as.formula(gsub("inla.surv","Surv",deparse(formula)))
  # Computes the Kaplan Meier curve using the package "rms"
  ObjSurvfit <- rms::npsurv(      # Uses the function "npsurv" from the package "rms"
    formula = km.formula,         # to fit the model specified in the "formula" object
    data = data                   # to the dataset named "data"
  )
  return(ObjSurvfit)
}




#' Helper function to format the output of the modelling (produced either
#' by running 'runMLE', or 'runINLA', 'runHMC'), in a way that is consistent
#' with the architecture of 'survHE'
#' 
#' @param output The output of one of the helper functions used to run the
#' models.
#' @param method The method used to do the estimation
#' @param distr The abbreviated name for the distribution to be used
#' @param formula The model formula
#' @param data The dataset used
#' @return \item{res}{A 'survHE' object containing all the relevant output
#' conveniently formatted}.
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso fit.models
#' @references Baio (2020). survHE
#' @keywords Parametric survival models Bayesian inference via Hamiltonian
#' Monte Carlo Bayesian inference via Integrated Nested Laplace Approximation
#' @noRd 
format_output_fit.models <- function(output,method,distr,formula,data) {
  
  # Uses the helper 'manipulated_distributions' to create the vector labs
  labs <- manipulate_distributions(distr)$labs
  
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
    data=data,
    model_name=unlist(lapply(output,function(x) x$model_name))
  )
  if(any(distr=="polyweibull")) {
    misc$km=lapply(formula,function(f) make_KM(f,data))
  } else {
    misc$km=make_KM(formula,data)
  }
  
  # HMC-specific extra output
  if(method=="hmc"){
    # Completes the 'misc' and 'model.fitting' lists with additional output
    misc$data.stan <- lapply(output,function(x) x$data.stan)
    model.fitting$dic2 <- unlist(lapply(output,function(x) x$dic2))
  }
  if(method=="inla") {
    model.fitting$dic2 = unlist(lapply(output,function(x) x$dic2))
  }
  
  # Names the elements of the list
  names(models) <- labs
  
  # Formats all output in a list
  res <- list(models=models,model.fitting=model.fitting,method=method,misc=misc)
  # And sets its class attribute to "survHE"
  class(res) <- "survHE"
  return(res)
}



#' Helper function to provide a list of models available in each method
#' 
#' @return \item{availables}{A list of models available in each method}.
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso fit.models
#' @references Baio (2020). survHE
#' @keywords Parametric survival models Bayesian inference via Hamiltonian
#' Monte Carlo Bayesian inference via Integrated Nested Laplace Approximation
#' @noRd 
load_availables <- function() {
  # INLA can only do a limited set of models (for now) so if user has selected
  # one that is not available, then falls back on MLE analysis
  availables=list(
    mle=c("genf" = "gef",
          "genf.orig" = "gof",
          "gengamma" = "gga",
          "gengamma.orig" = "ggo",
          "exp" = "exp",
          "weibull" = "wei",
          "weibullPH" = "wph",
          "lnorm" = "lno",
          "gamma" = "gam",
          "gompertz" = "gom",
          "llogis" = "llo",
          "lognormal" = "lno",
          "rps" = "rps"
    ),
    inla=c("exponential" = "exp",
           "weibull" = "wei",
           "weibullPH" = "wph",
           "lognormal" = "lno",
           "loglogistic" = "llo",
           "rps" = "rps",
           "gompertz" = "gom"      # added Mar 19, 2021
    ),
    hmc=c("Exponential" = "exp",
          "Gamma" = "gam",
          "GenF" = "gef",
          "GenGamma" = "gga",
          "Gompertz" = "gom",
          "PolyWeibull" = "pow",
          "RP" = "rps",
          "WeibullAF" = "wei",
          "WeibullPH" = "wph",
          "logLogistic" = "llo",
          "logNormal" = "lno"
    )
  )
  return(availables)
}



#' Helper function to manipulate the strings of text defining the 
#' distributions selected by the user so they are consistent with the
#' various methods
#' 
#' @param x A string with the distribution name selected by the user.
#' @return \item{list}{A list containing the modified name of the 
#' distribution, the acronym (3-letters abbreviation), or the
#' labels (humane-readable name)}.
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso fit.models
#' @references Baio (2020). survHE
#' @keywords Parametric survival models Bayesian inference via Hamiltonian
#' Monte Carlo Bayesian inference via Integrated Nested Laplace Approximation
#' @noRd 
manipulate_distributions <- function(x){
  # selected model checks -----
  matchTable = list(
    "exp" = c("exponential", "exp"),
    "wei" = c("weibull", "weibullaft", "weiaft", "waft", "weibullaf", "weiaf", "waf", "wei"),
    "wph" = c("weibullph", "weiph", "wph"),
    "gam" = c("gamma", "gam", "gma"),
    "lno" = c("lognormal", "lnormal", "lnorm", "lognorm", "lno"),
    "llo" = c("loglogistic", "loglog", "llogistic", "llogis", "llo", "llogist"),
    "gga" = c("generalisedgamma", "generalizedgamma", "ggamma", "gengamma", "gga", "ggam"),
    "ggo" = c("gengamma.orig", "ggo"),
    "gef" = c("generalisedf", "generalizedf", "genf", "gef"),
    "gof" = c("genf.orig", "gof"),
    "gom" = c("gompertz", "gpz", "gomp", "gompz", "gom"),
    "rps" = c("roystonparmar", "roystonparmarsplines", "roystonparmarspline", "spline", "splines", "rps"),
    "pow" = c("polyweibull","pow","PolyWeibull")
  )
  # Human readable label
  labelTable = c(
    "exp" = "Exponential",
    "wei" = "Weibull (AFT)",
    "wph" = "Weibull (PH)",
    "gam" = "Gamma",
    "lno" = "log-Normal", 
    "llo" = "log-Logistic",
    "gga" = "Gen. Gamma", "ggo" = "Gen. Gamma (orig parametrisation)",
    "gef" = "Gen. F", "gof" = "Gen. F (orig parametrisation)",
    "gom" = "Gompertz",
    "rps" = "Royston-Parmar",
    "pow" = "Poly-Weibull")
  # Labels used by R to define p..., r... and d... commands
  labelR = c(
    "exp" = "exp",
    "wei" = "weibull", 
    "wph" = "weibullPH", 
    "gam" = "gamma",
    "lno" = "lnorm", 
    "llo" = "llogis",
    "gga" = "gengamma", 
    "ggo" = "gengamma.orig",
    "gef" = "genf",
    "gof" = "genf.orig",
    "gom" = "gompertz",
    "rps" = "survspline",
    "pow" = "polyweibull"
  )
  
  distr = gsub("[ ]*[-]*", "", tolower(x))
  isDistrUnmatched = which(!sapply(
    1:length(distr),
    '%in%',
    unname(unlist(sapply(matchTable, match, distr)))))
  if (length(isDistrUnmatched) > 0) {
    stop(paste0("Distribution ", paste(distr[isDistrUnmatched], collapse = ", "), " could not be matched."))
  }
  
  distr3 <- numeric()
  for (i in 1:length(distr)) {
    distr3[i] <- names(which(unlist(lapply(matchTable,function(x) distr[i]%in%x))))  
  }
  labs <- unname(labelTable[distr3])
  distr <- unname(labelR[distr3])
  
  list(distr=distr,distr3=distr3,labs=labs)
}


#' Helper function to check that the distribution(s) provided by the user are
#' consistent with the method chosen for inference.
#' 
#' \code{'inla'} or \code{'hmc'}). If \code{method} is set to \code{'hmc'},
#' then \code{survHE} will write suitable model code in the Stan language
#' (according to the specified distribution), prepare data and initial values
#' and then run the model.
#' @param distr3 A vector of distribution labels (as created by 'fit.models'). 
#' It's a 3-letters label to identify the distributions
#' @param availables A list with the distributions available for each method.
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso fit.models
#' @references Baio (2020). survHE
#' @keywords Parametric survival models 
#' @noRd 
check_distributions <- function(method,distr) {
  # Loads in the available models in each method
  availables <- load_availables()
  # Uses the helper 'manipulated_distributions' to create the vectors distr, distr3 and labs
  distr3 <- manipulate_distributions(distr)$distr3

  # If 'method' is either 'inla' or 'hmc but we're trying to run a model that is not available, then
  # falls back to 'mle'
  if(method %in% c("inla","hmc")) {
    if(!all(distr3 %in% availables[[method]])) {
      ####modelsString <- unname(labelTable[availables[[method]]])
      modelsString <- unname(manipulate_distributions(availables[[method]])$labs)
      modelsString[length(modelsString)] = paste0("or ", modelsString[length(modelsString)])
      message(paste0(
        "NB: ",toupper(method)," can only fit ",
        paste(modelsString, collapse = ", "),
        " parametric survival models. Falling back on MLE analysis")
      )
      method <- "mle"
    }
  }
  
  # 'mle' can implement all the possible models, except the PolyWeibull
  # In this case, I choose to *stop* execution, rather than falling back to 'hmc'!
  if (method == "mle") {
    if(!all(distr3 %in% availables[[method]])) {
      stop(paste0("The Poly-Weibull model is only implemented under method='hmc'.
       Please set this option in your call to 'fit.models'"))
    }
  }
  return(method)
}
