#' Helper function to run the survival models using INLA
#' for a given formula and dataset
#' 
#' @param x a (vector of) string(s) containing the name(s) of the model(s)
#' to be fitted
#' @param exArgs a list of extra arguments passed from the main 'fit.models' 
#' function
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso fit.models
#' @references Baio (2020). survHE
#' @keywords Parametric survival models Integrated Nested Laplace Approximation
#' @noRd 
runINLA <- function(x,exArgs) {
  # First checks whether INLA is installed (it's only a suggestion, not a full dependency)
  if (!isTRUE(requireNamespace("INLA", quietly = TRUE))) {
    stop("You need to install the packages 'INLA'. Please run in your R terminal:\n install.packages('INLA', repos='https://inla.r-inla-download.org/R/stable')")
  }
  # If INLA is installed but not loaded then attach the Namespace (so that all the relevant functions are available)
  if (isTRUE(requireNamespace("INLA", quietly = TRUE))) {
    if (!is.element("INLA", (.packages()))) {
      attachNamespace("INLA")
    }
  }
  
  # Loads the model formula & data
  formula <- exArgs$formula
  data=exArgs$data
  
  # Loads in the available models in each method
  availables <- load_availables()
  # Uses the helper 'manipulated_distributions' to create the vectors distr, distr3 and labs
  d3 <- manipulate_distributions(x)$distr3
  method <- "inla"
  
  # Set up optional parameters to default values if the user hasn't done it themselves
  # 1. defines the step length for the grid search over the hyperparameters space
  if(exists("dz",where=exArgs)) {dz <- exArgs$dz} else {dz <- 0.1}
  # 2. defines the difference in the log-density for the hyperparameters to stop integration
  if(exists("diff.logdens",where=exArgs)) {diff.logdens <- exArgs$diff.logdens} else {diff.logdens <- 5}
  # 3. defines the default for the priors, unless specified by the user
  if(exists("control.fixed",where=exArgs)) {
    control.fixed <- exArgs$control.fixed
  } else {
    control.fixed <- INLA::inla.set.control.fixed.default()
    # prior mean = 0 for *all* fixed effects
    # prior var = 1000 for *all* fixed effects
    # prior mean = 0 for the intercept
    # prior prec -> 0 for the intercept 
    ## This makes the priors consistent with the defaults in HMC
    ## The available models all have sd=5 in HMC, which translates to a precision of 1/25!
    control.fixed$prec <- control.fixed$prec.intercept <- 1/(5^2)
  }
  # Recomputes the three-letters code for the distributions and the INLA-specific name
  d <- names(availables[[method]][match(d3, availables[[method]])])
  # Needs to create a copy of the model name because if it's 'WeibullPH' INLA won't find it...
  if(d=="weibullPH") {d="weibull"} 
  # If 'control.family' is specified for the distribution currently used, then use the values in
  # 'exArgs'. If not, or if specified but for another distribution, use INLA defaults
  cf <- INLA::inla.set.control.family.default()
  if(exists("control.family",where=exArgs)) {
    abbrs=manipulate_distributions(names(exArgs$control.family))$distr3
    pos=grep(d3,abbrs)
    if(length(pos)>0) {
      cf <- exArgs$control.family[[pos]]
    }
  }
  if(exists("verbose",where=exArgs)) {verbose <- exArgs$verbose} else {verbose <- FALSE}
  
  # 4. Finally runs INLA
  # Ensures that the formula is in INLA terms. If not, make it 
  if(!grepl("inla.surv",deparse(formula))) {formula <- as.formula(gsub("Surv","INLA::inla.surv",deparse(formula)))}
  # As of 9 Jan 2017, INLA is creating new distribution names for survival models, so needs to update the name.
  # Also, there are two variants of the Weibull model (PH vs AFT) so need to identify that too
  if(d3=="wph") {
    cf$variant <- 0
  } else if (d3=="wei") {
    cf$variant <- 1
  }
  if(d3=="llo") {
    cf$variant=1
  }
  model <- INLA::inla(formula,family=paste0(d,"surv"),data=data,control.compute=list(config=TRUE,dic=TRUE),
                      control.inla=list(int.strategy="grid",dz=dz,diff.logdens=diff.logdens),
                      control.fixed=control.fixed,control.family=cf,verbose=verbose
  )
  
  # Now re-writes the formula in general terms (without linking to INLA::inla.surv)
  formula <- as.formula(gsub("INLA::inla.surv","Surv",deparse(formula)))
  # Adds a field used in 'make.surv' to indicate the model used
  model_name <- d3
  # And reset if original name *was* 'WeibullPH'
  if(d3=="wph") {d="weibullPH"} 
  
  # Finally returns the output
  list(
    model=model,
    aic=2*model$dic$p.eff+model$dic$deviance.mean,
    bic=model$dic$deviance.mean+model$dic$p.eff*log(model$size.linear.predictor$n),
    dic=model$dic$dic,
    time2run=model$cpu.used["Total"],
    model_name=model_name
  )
}



#' Helper function to run the survival models using HMC (rstan)
#' for a given formula and dataset
#' 
#' @param x a (vector of) string(s) containing the name(s) of the model(s)
#' to be fitted
#' @param exArgs a list of extra arguments passed from the main 'fit.models' 
#' function
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso fit.models
#' @references Baio (2020). survHE
#' @keywords Parametric survival models Hamiltonian Monte Carlo
#' @noRd 
runHMC <- function(x,exArgs) {
  # First checks whether INLA is installed (it's only a suggestion, not a full dependency)
  if (!isTRUE(requireNamespace("rstan", quietly = TRUE))) {
    stop("You need to install the R package 'rstan'. Please run in your R terminal:\n install.packages('rstan')")
  }
  
  # Loads the model formula & data
  formula <- exArgs$formula
  data=exArgs$data
  
  # Loads in the available models in each method
  availables <- load_availables()
  # Uses the helper 'manipulated_distributions' to create the vectors distr, distr3 and labs
  d3 <- manipulate_distributions(x)$distr3
  method <- "hmc"
  
  # Now runs the model
  # Set up optional parameters to default values if the user hasn't done it themselves
  if (exists("chains", where = exArgs)) {chains <- exArgs$chains} else {chains <- 2} # DO WE WANT 4???
  if (exists("iter", where = exArgs)) {iter <- exArgs$iter} else {iter <- 2000}
  if (exists("warmup", where = exArgs)) {warmup <- exArgs$warmup} else {warmup <- floor(iter/2)}
  if (exists("thin", where = exArgs)) {thin <- exArgs$thin} else {thin <- 1}
  if (exists("control", where = exArgs)) {
    control <- exArgs$control 
  } else {
    control <- list(NULL)
  }
  if (exists("seed",where = exArgs)) {seed <- exArgs$seed} else {seed <- sample.int(.Machine$integer.max, 1)}
  if (exists("pars",where = exArgs)) {pars <- exArgs$pars} else {
    pars <- c("lambda_cens","lambda_obs","cens","d","lp__","loglambda_cens","loglambda_obs","mu","logP","linpred")
  }
  if (exists("include",where = exArgs)) {include <- exArgs$include} else {include <- FALSE}
  if (exists("k",where = exArgs)) {k <- exArgs$k} else {k <- 0}
  if (exists("cores",where = exArgs)) {cores <- exArgs$cores} else {cores <- 1}
  if (exists("init",where = exArgs)) {init <- exArgs$init} else {init="random"}
  if (exists("save.stan",where=exArgs)) {save.stan <- exArgs$save.stan} else {save.stan=FALSE}
  if (exists("refresh",where = exArgs)) {refresh=exArgs$refresh} else {refresh=max(iter/10, 1)}
  
  # Recomputes the three-letters code for the distributions and the HMC-specific name
  d <- names(availables[[method]][match(d3, availables[[method]])])
  
  # Loads the pre-compiled models
  dso <- stanmodels[[d]]
  
  # Create the data list
  data.stan <- make_data_stan(formula,data,d3,exArgs)
  
  # Now runs Stan to sample from the posterior distributions
  tic <- proc.time()
  model <- rstan::sampling(dso,data.stan,chains=chains,iter=iter,warmup=warmup,thin=thin,seed=seed,
                           control=control,pars=pars,include=include,cores=cores,init=init,refresh=refresh)
  toc <- proc.time()-tic
  time_survHE <- toc[3]
  # rstan does have its way of computing the running time, but it may not be the actual one when running multiple
  # chains. 
  time_stan <- sum(rstan::get_elapsed_time(model))
  
  # Uses the helper function to compute the *IC
  ics <- compute_ICs_stan(model,d3,data.stan)
  
  # If 'save.stan' is set to TRUE, then saves the Stan model file(s) & data
  if(save.stan) {
    model_code <- attr(model@stanmodel,"model_code")
    con <- paste0(d,".stan")
    writeLines(model_code,con=con)
    cat(paste0("Model code saved to the file: ",con,"\n"))
  }
  
  # Adds a field used in 'make.surv' to indicate the model used
  model_name <- d3
  
  # Finally returns the output
  list(
    model=model,
    aic=ics$aic,
    bic=ics$bic,
    dic=ics$dic,
    dic=ics$dic2,
    time2run=pmin(time_survHE,time_stan),
    data.stan=data.stan,
    save.stan=save.stan,
    model_name=model_name
  )
}



#' Helper function to create data in the correct format for rstan
#' 
#' @param formula a formula specifying the model to be used, in the form
#' \code{Surv(time,event)~treatment[+covariates]} in flexsurv terms, or
#' \code{inla.surv(time,event)~treatment[+covariates]} in INLA terms.
#' @param data A data frame containing the data to be used for the analysis.
#' This must contain data for the 'event' variable. In case there is no
#' censoring, then \code{event} is a column of 1s.
#' @return \item{data.stan}{A list containing the variables needed to pass
#' to 'stan' when calling \code{fit.models} with \code{method="hmc"}}.
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso fit.models
#' @references Baio (2020). survHE
#' @keywords Parametric survival models Bayesian inference via Hamiltonian
#' Monte Carlo 
#' @noRd 
make_data_stan=function(formula,data,distr3,exArgs) {
  # Loads in the available models in each method
  availables <- load_availables()
  method <- "hmc"
  
  # Modifies the original formula to separate 'time' and 'event'
  formula_temp <- update(formula,paste(all.vars(formula,data)[1],"~",all.vars(formula,data)[2],"+."))
  # Creates a model.frame + renames the variables to conform with stan's expectations later
  mf <- as_tibble(model.frame(formula_temp,data)) %>% 
    rename(time=1, event=2) %>% rename_if(is.factor,.funs=~gsub("as.factor[( )]","",.x)) %>% 
    rename_if(is.factor,.funs=~gsub("[( )]","",.x)) %>% 
    bind_cols(as_tibble(model.matrix(formula_temp,data)) %>% select(contains("Intercept"))) %>% 
    select(time,event,contains("Intercept"),everything()) %>% tibble::rownames_to_column("ID")
  
  # Now arrange data in list to pass to 'stan'
  # NB: Need different formatting depending on the underlying sampling distribution
  if (distr3 %in% c("gam","gga","gef")) {
    # If model is Gamma, GenGamma or GenF, then use the "obs vs" censored format
    data.stan <- list(t=(mf %>% filter(event==1))$time,
                      d=(mf %>% filter(event==0))$time,
                      n_obs=mf %>% filter(event==1) %>% with(nrow(.)),
                      n_cens=mf %>% filter(event==0) %>% with(nrow(.))
    )
    data.stan$X_obs <- matrix(model.matrix(formula,data)[(mf %>% filter(event==1))$ID,],nrow=data.stan$n_obs)
    data.stan$X_cens <- matrix(model.matrix(formula,data)[(mf %>% filter(event==0))$ID,],nrow=data.stan$n_cens)
    data.stan$H=ncol(data.stan$X_obs)
    
    # NB: Stan doesn't allow vectors of size 1, so if there's only one covariate (eg intercept only), needs a little trick
    if (data.stan$H==1) {
      data.stan$X_obs <- cbind(data.stan$X_obs,rep(0,data.stan$n_obs))
      data.stan$X_cens <- cbind(data.stan$X_cens,rep(0,data.stan$n_cens))
      data.stan$H <- ncol(data.stan$X_obs)
    }
  }
  
  if (distr3 %in% c("exp", "gom", "wei", "wph", "llo", "lno","pow")) {
    # If it's one of the others (except polyweibull), use the "h,S" format
    data.stan <- list(t=(mf$time),
                      d=mf$event,
                      n=nrow(mf),
                      X=matrix(model.matrix(formula,data),nrow=nrow(mf)),
                      H=ncol(model.matrix(formula,data))
    )
    # NB: Stan doesn't allow vectors of size 1, so if there's only one covariate (eg intercept only), needs a little trick
    if (data.stan$H==1) {
      data.stan$X <- cbind(data.stan$X,rep(0,data.stan$n))
      data.stan$H <- ncol(data.stan$X)
    }
  }

  # RPS needs specific data to be created and stored in the 'data.stan' list
  if (distr3=="rps"){
    # If it's Royston-Parmar splines, then gets the correct data 
    if (exists("k",where = exArgs)) {k <- exArgs$k} else {k <- 0}
    knots <- quantile(log((mf %>% filter(event==1))$time), seq(0,1,length=k+2))
    # Uses flexsurv to compute the basis and derivatives of the basis
    B <- flexsurv::basis(knots,log(mf$time))
    DB <- flexsurv::dbasis(knots,log(mf$time))
    # Now checks to see whether the user wants to specify covariates and removes the intercept from the formula (for identifiability)
    mm <- model.matrix(formula,data)[,-1]
    # a. if the formula is ~ 1, then adds two fictional covariates of all 0s
    if (length(mm)<1) {
      mm <- matrix(rep(0,nrow(mf)),nrow=nrow(mf),ncol=2)
    }
    # b. in case there's only one covariate, then adds another fake covariate of all 0s
    if (is.null(dim(mm))) {
      mm <- cbind(mm,rep(0,length(mm)))
    }
    data.stan <- list(t=mf$time,
                      d=mf$event,
                      n=nrow(mf),
                      M=k,
                      X=mm,
                      H=ncol(mm),
                      B=B,
                      DB=DB,
                      mu_gamma=rep(0,k+2),
                      sigma_gamma=rep(5,k+2),
                      knots=knots
    )
  }
  
  ## SETS THE PRIORS
  # Adds the values for the parameters of the prior distributions
  data.stan$mu_beta = rep(0, data.stan$H)
  # For the models *not* on the log scale, needs higher values for the sd of the coefficients
  if (distr3 %in% c("gef","gga","lno")) {
    data.stan$sigma_beta <- rep(100, data.stan$H)
  } else {
    data.stan$sigma_beta <- rep(5, data.stan$H)
  }
  
  # Ancillary parameters
  if (distr3 == "gef") {
    data.stan$a_sigma = data.stan$b_sigma = 0.1
    data.stan$mu_P = 0
    data.stan$sigma_P = 0.5
    data.stan$mu_Q = 0
    data.stan$sigma_Q = 2.5
  } else if (distr3 == "gga") {
    data.stan$a_sigma = data.stan$b_sigma = 0.1
    data.stan$mu_Q = 0
    data.stan$sigma_Q = 100
  } else if (distr3 %in% c("gam", "llo", "wei", "wph")) {
    data.stan$a_alpha = data.stan$b_alpha = 0.1
  } else if (distr3 %in% c("lno", "gom")) {
    data.stan$a_alpha = 0
    data.stan$b_alpha = 5
  }
  
  # Sets the default priors to an empty list, which gets modified if the user has specified some values
  # NB: the user needs to specify a **named** list, where the names have to match the distributions in the
  # correct order!
  # recodes 
  d <- names(availables[[method]][match(distr3, availables[[method]])])
  priors <- list()
  if(exists("priors",where=exArgs)) {
    # Checks whether the user has specified the current model in the named list 'priors'
    abbrs=manipulate_distributions(names(exArgs$priors))$distr3
    pos=grep(distr3,abbrs)
    if(length(pos)>0) {
      priors = exArgs$priors[[pos]]
    }
  }
  # Linear predictor coefficients
  if(!is.null(priors$mu_beta)) {
    data.stan$mu_beta=priors$mu_beta
  }
  if(!is.null(priors$sigma_beta)) {
    data.stan$sigma_beta <- priors$sigma_beta
  }
  if(!is.null(priors$mu_gamma) & distr3=="rps") {
    data.stan$mu_gamma <- priors$mu_gamma
  }
  if(!is.null(priors$sigma_gamma) & distr3=="rps") {
    data.stan$sigma_gamma <- priors$sigma_gamma
  }
  
  # Ancillary parameters
  if(!is.null(priors$a_sigma)) {data.stan$a_sigma=priors$a_sigma}
  if(!is.null(priors$b_sigma)) {data.stan$b_sigma=priors$b_sigma}
  if(!is.null(priors$mu_P)) {data.stan$mu_P=priors$mu_P}
  if(!is.null(priors$sigma_P)) {data.stan$sigma_P=priors$sigma_P}
  if(!is.null(priors$mu_Q)) {data.stan$mu_Q=priors$mu_Q}
  if(!is.null(priors$sigma_Q)) {data.stan$sigma_Q=priors$sigma_Q}
  if(!is.null(priors$a_alpha)) {data.stan$a_alpha=priors$a_alpha}
  if(!is.null(priors$b_alpha)) {data.stan$b_alpha=priors$b_alpha}
  
  # Returns the list of data
  data.stan
}



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
    model <- flexsurv::flexsurvspline(formula=formula,data=data,k=k,knots=knots,bknots=bknots,scale=scale,
                                      timescale=timescale)
  } else {
    # If it's one of the other available models under MLE, then simply runs flexsurv::flexsurvreg
    model <- flexsurv::flexsurvreg(formula=formula,data=data,dist=x)
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
           "rps" = "rps"
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


#' Helper function to compute the information criteria statistics
#' when using hmc as the inferential engine. 'rstan' does not do
#' DIC automatically and AIC/BIC are also not standard for Bayesian
#' models, so can compute them post-hoc by manipulating the 
#' likelihood functions.
#' 
#' @param model The 'rstan' object with the model fit
#' @param distr3 The 'rstan' object with the model fit
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
compute_ICs_stan <- function(model,distr3,data.stan) {
  # Computes the log-likelihood 
  beta <- rstan::extract(model)$beta
  # To safeguard against very asymmetric densities use the median (instead of the mean)
  beta.hat <- apply(beta,2,median)
  if (distr3 %in% c("exp", "wei", "wph", "gom", "lno", "llo","rps")) {
    linpred <- beta%*%t(data.stan$X)
    linpred.hat <- beta.hat%*%t(data.stan$X)
  } else if(distr3=="pow") {
    alpha <- rstan::extract(model)$alpha
    beta.hat <- apply(beta,c(2,3),median)
    alpha.hat <- apply(alpha,2,median)
    linpred <- lapply(1:data.stan$M,function(m) {
      beta[,m,]%*%t(data.stan$X[m,,])
    })
    linpred.hat <- lapply(1:data.stan$M,function(m) {
      beta.hat[m,]%*%t(data.stan$X[m,,])
    })
  } else {
    linpred <- list(
      lo=beta%*%t(data.stan$X_obs),
      lc=beta%*%t(data.stan$X_cens)
    )
    linpred.hat <- list(
      lo.bar=beta.hat%*%t(data.stan$X_obs),
      lc.bar=beta.hat%*%t(data.stan$X_cens)
    )
  }

  # Now computes the densities using the helper functions
  out=do.call(what=paste0("lik_",distr3),args=list(distr3,linpred,linpred.hat,model,data.stan))
  # Extracts relevant variables from the list (See if there's a better way to do it!)
  logf=out$logf
  logf.hat=out$logf.hat
  npars=out$npars
  f=out$f
  f.bar=out$f.bar
  s=out$s
  s.bar=out$s.bar
  
  # Now computes the log-likelihood and then deviance and DIC, AIC, BIC
  if (distr3 %in% c("gam","gga","gef")) {
    loglik <- compute.loglik(f,s)
    loglik.bar <- compute.loglik(f.bar,s.bar)
    data.stan$n <- data.stan$n_obs+data.stan$n_cens
  } else if(distr3=="pow") {
    # If the model is Poly-Weibull, then computes these quantities directly
    h <- log_s <- array(NA,c(nrow(alpha),data.stan$n,data.stan$M))
    h_bar <- log_s_bar <- matrix(NA,data.stan$n,data.stan$M)
    for (m in 1:data.stan$M) {
      h_bar[,m] <- alpha.hat[m]*exp(linpred.hat[[m]])*data.stan$t^(alpha.hat[m]-1)
      log_s_bar[,m] <- exp(linpred.hat[[m]])*data.stan$t^(alpha.hat[m])
      for (i in 1:nrow(linpred[[m]])) {
        h[i,,m] <- alpha[i,m]*exp(linpred[[m]][i,])*data.stan$t^(alpha[i,m]-1)
        log_s[i,,m] <- exp(linpred[[m]][i,])*data.stan$t^(alpha[i,m])
      }
    }
    d_log_sum_h <- matrix(NA,nrow(alpha),data.stan$n)
    for (i in 1:nrow(alpha)) {
      d_log_sum_h[i,] <- data.stan$d * log(rowSums(h[i,,]))
    }
    loglik.bar <- sum(data.stan$d*log(rowSums(h_bar))-rowSums(log_s_bar))
    loglik <- rowSums(d_log_sum_h) - rowSums(log_s,2)
  } else {
    loglik <- apply(logf,1,sum)
    loglik.bar <- apply(logf.hat,1,sum)
  }
  D.theta <- -2*loglik
  D.bar <- -2*loglik.bar
  pD <- mean(D.theta) - D.bar
  pV <- 0.5*var(D.theta)
  dic <- mean(D.theta)+pD
  dic2 <- mean(D.theta)+pV
  # Approximates AIC & BIC using the mean deviance and the number of nominal parameters
  aic <- D.bar+2*npars                   #mean(D.theta)+2*pD
  bic <- D.bar+npars*log(data.stan$n)    #mean(D.theta)+pD*log(data.stan$n)
    
  list(aic=aic,bic=bic,dic=dic,dic2=dic2)
}

### Little function to compute the log-likelihood (for the obs vs censored cases)
compute.loglik <- function(f,s) {
  loglik <- (apply(log(f),1,sum) + apply(log(s),1,sum))
  return(loglik)
}

### These are utility functions to compute the log-density for the models fitted by hmc
lik_pow <- function(x,linpred,linpred.hat,model,data.stan) {
  # Little trick for the Poly-Weibull. Because the distribution is more complex, 
  # simply computes the relevant variables in the main code and this returns mostly NULL
  npars <- data.stan$M + sum(unlist(lapply(1:data.stan$M,function(m) {sum(1-apply(data.stan$X[m,,],2,function(x) all(x==0)))})))
  list(logf=NULL,logf.hat=NULL,npars=npars,f=NULL,f.bar=NULL,s=NULL,s.bar=NULL)
}

lik_exp <- function(x,linpred,linpred.hat,model,data.stan) {
  logf <- matrix(
    unlist(lapply(1:nrow(linpred), function(i) {
      data.stan$d * log(hexp(data.stan$t, exp(linpred[i,]))) +
        log(1 - pexp(data.stan$t, exp(linpred[i,])))
    })),
    nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(
    data.stan$d * log(hexp(data.stan$t, exp(linpred.hat))) +
      log(1 - pexp(data.stan$t, exp(linpred.hat))),
    nrow = 1)
  # Number of parameters (for AIC): rate + covariates
  npars <- 1 + sum(1 - apply(data.stan$X, 2, function(x) all(x == 0)))
  list(logf=logf,logf.hat=logf.hat,npars=npars,f=NULL,f.bar=NULL,s=NULL,s.bar=NULL)
}

lik_wei <- function(x,linpred,linpred.hat,model,data.stan) {
  shape <- alpha <- as.numeric(rstan::extract(model)$alpha)
  shape.hat <- median(shape)
  logf <- matrix(
    unlist(lapply(1:nrow(linpred), function(i) {
      data.stan$d * log(hweibull(data.stan$t, shape[i], exp(linpred[i, ]))) + 
        log(1 - pweibull(data.stan$t, shape[i], exp(linpred[i, ])))
    })), 
    nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(
    data.stan$d * log(hweibull(data.stan$t, shape.hat, exp(linpred.hat))) + 
      log(1 - pweibull(data.stan$t, shape.hat, exp(linpred.hat))), 
    nrow = 1)
  # Number of parameters (for AIC): shape, scale + covariates
  npars <- 2 + sum(1 - apply(data.stan$X, 2, function(x) all(x == 0)))
  list(logf=logf,logf.hat=logf.hat,npars=npars,f=NULL,f.bar=NULL,s=NULL,s.bar=NULL)
}

lik_wph <- function(x,linpred,linpred.hat,model,data.stan) {
  shape <- alpha <- as.numeric(rstan::extract(model)$alpha)
  shape.hat = median(shape)
  logf <- matrix(
    unlist(lapply(1:nrow(linpred), function(i) {
      data.stan$d * log(hweibullPH(data.stan$t, shape[i], exp(linpred[i, ]))) +
        log(1 - pweibullPH(data.stan$t, shape[i], exp(linpred[i, ])))
    })), nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(
    data.stan$d * log(hweibullPH(data.stan$t, shape.hat, exp(linpred.hat))) +
      log(1 - pweibullPH(data.stan$t, shape.hat, exp(linpred.hat))), 
    nrow = 1)
  # Number of parameters (for AIC): shape, scale + covariates
  npars <- 2 + sum(1 - apply(data.stan$X, 2, function(x) all(x == 0)))
  list(logf=logf,logf.hat=logf.hat,npars=npars,f=NULL,f.bar=NULL,s=NULL,s.bar=NULL)
}

lik_gom <- function(x,linpred,linpred.hat,model,data.stan) {
  shape <- alpha <- as.numeric(rstan::extract(model)$alpha)
  shape.hat = median(shape)
  logf <- matrix(
    unlist(lapply(1:nrow(linpred), function(i) {
      data.stan$d * log(hgompertz(data.stan$t, shape = shape[i], rate = exp(linpred[i, ]))) +
        log(1 - pgompertz(shape = data.stan$t, shape[i], rate = exp(linpred[i, ])))
    })), 
    nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(
    data.stan$d * log(hgompertz(data.stan$t, shape.hat, exp(linpred.hat))) +
      log(1 - pgompertz(data.stan$t, shape.hat, exp(linpred.hat))),
    nrow = 1)
  # Number of parameters (for AIC): shape, rate + covariates
  npars <- 2 + sum(1 - apply(data.stan$X, 2, function(x) all(x == 0)))
  list(logf=logf,logf.hat=logf.hat,npars=npars,f=NULL,f.bar=NULL,s=NULL,s.bar=NULL)
}

lik_gam <- function(x,linpred,linpred.hat,model,data.stan) {
  shape <- alpha <- as.numeric(rstan::extract(model)$alpha)
  shape.bar <- median(shape)
  lo <- exp(linpred$lo)
  lc <- exp(linpred$lc)
  lo.bar <- exp(linpred.hat$lo.bar)
  lc.bar <- exp(linpred.hat$lc.bar)
  f = matrix(
    unlist(lapply(1:nrow(lo), function(i) 
      dgamma(data.stan$t, shape[i], lo[i, ]))),
    nrow = nrow(lo), byrow = T)
  f.bar = matrix(
    unlist(lapply(1:nrow(lo.bar), function(i)
      dgamma(data.stan$t, shape.bar, lo.bar[i, ]))), 
    nrow = 1, byrow = T)
  s = matrix(
    unlist(lapply(1:nrow(lc), function(i) 
      1 - pgamma(data.stan$d, shape[i], lc[i, ]))),
    nrow = nrow(lc), byrow = T)
  s.bar = matrix(
    unlist(lapply(1:nrow(lc.bar), function(i) 
      1 - pgamma(data.stan$d, shape.bar, lc.bar[i, ]))),
    nrow = 1, byrow = T)
  # Number of parameters (for AIC): shape, rate + covariates
  npars <- 2 + sum(1 - apply(data.stan$X_obs, 2, function(x) all(x == 0)))
  list(logf=NULL,logf.hat=NULL,npars=npars,f=f,f.bar=f.bar,s=s,s.bar=s.bar)
}

lik_gga <- function(x,linpred,linpred.hat,model,data.stan) {
  q = as.numeric(rstan::extract(model)$Q)
  q.bar = median(q)
  scale = as.numeric(rstan::extract(model)$sigma)
  scale.bar = median(scale)
  lo <- (linpred$lo)
  lc <- (linpred$lc)
  lo.bar <- (linpred.hat$lo.bar)
  lc.bar <- (linpred.hat$lc.bar)
  f = matrix(
    unlist(lapply(1:nrow(lo), function(i)
      dgengamma(data.stan$t, lo[i, ], scale[i], q[i]))), 
    nrow = nrow(lo), byrow = T)
  f.bar = matrix(
    unlist(lapply(1:nrow(lo.bar), function(i)
      dgengamma(data.stan$t, lo.bar[i, ], scale.bar, q.bar))), 
    nrow = 1, byrow = T)
  s = matrix(
    unlist(lapply(1:nrow(lc), function(i)
      1 - pgengamma(data.stan$d, lc[i, ], scale[i], q[i]))), 
    nrow = nrow(lc), byrow = T)
  s.bar = matrix(
    unlist(lapply(1:nrow(lc.bar), function(i) 
      1 - pgengamma(data.stan$d, lc.bar[i, ], scale.bar, q.bar))), 
    nrow = 1, byrow = T)
  # Number of parameters (for AIC): mu, sigma, Q + covariates
  npars <- 3 + sum(1 - apply(data.stan$X_obs, 2, function(x) all(x == 0)))
  list(logf=NULL,logf.hat=NULL,npars=npars,f=f,f.bar=f.bar,s=s,s.bar=s.bar)
}

lik_gef <- function(x,linpred,linpred.hat,model,data.stan) {
  Q = as.numeric(rstan::extract(model)$Q)
  Q.bar = median(Q)
  P = as.numeric(rstan::extract(model)$P)
  P.bar = median(P)
  sigma = as.numeric(rstan::extract(model)$sigma)
  sigma.bar = mean(sigma)
  lo <- (linpred$lo)
  lc <- (linpred$lc)
  lo.bar <- (linpred.hat$lo.bar)
  lc.bar <- (linpred.hat$lc.bar)
  f = matrix(
    unlist(lapply(1:nrow(lo), function(i) 
      dgenf(data.stan$t, lo[i, ], sigma[i], Q[i], P[i]))), 
    nrow = nrow(lo), byrow = T)
  f.bar = matrix(
    unlist(lapply(1:nrow(lo.bar), function(i) 
      dgenf(data.stan$t, lo.bar[i, ], sigma.bar, Q.bar, P.bar))), 
    nrow = 1, byrow = T)
  s = matrix(
    unlist(lapply(1:nrow(lc), function(i)
      1 - pgenf(data.stan$d, lc[i, ], sigma[i], Q[i], P[i]))),
    nrow = nrow(lc), byrow = T)
  s.bar = matrix(
    unlist(lapply(1:nrow(lc.bar), function(i) 
      1 - pgenf(data.stan$d, lc.bar[i, ], sigma.bar, Q.bar, P.bar))),
    nrow = 1, byrow = T)
  # Number of parameters (for AIC): mu, sigma, Q, P + covariates
  npars <- 4 + sum(
    1 - apply(data.stan$X_obs, 2, function(x) all(x == 0)))
  list(logf=NULL,logf.hat=NULL,npars=npars,f=f,f.bar=f.bar,s=s,s.bar=s.bar)
}

lik_lno <- function(x,linpred,linpred.hat,model,data.stan) {
  sigma = as.numeric(rstan::extract(model)$alpha)
  sigma.hat = median(sigma)
  logf <- matrix(
    unlist(lapply(1:nrow(linpred), function(i) {
      data.stan$d * log(hlnorm(data.stan$t, (linpred[i, ]), sigma[i])) +
        log(1 - plnorm(data.stan$t, (linpred[i, ]), sigma[i]))
    })), nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(
    data.stan$d * log(hlnorm(data.stan$t, (linpred.hat), sigma.hat)) +
      log(1 - plnorm(data.stan$t, (linpred.hat), sigma.hat)),
    nrow = 1)
  # Number of parameters (for AIC): meanlog, sdlog + covariates
  npars <- 2 + sum(
    1 - apply(data.stan$X, 2, function(x) all(x == 0)))
  list(logf=logf,logf.hat=logf.hat,npars=npars,f=NULL,f.bar=NULL,s=NULL,s.bar=NULL)
}

lik_llo <- function(x,linpred,linpred.hat,model,data.stan) {
  sigma = as.numeric(rstan::extract(model)$alpha)
  sigma.hat = median(sigma)
  logf <- matrix(
    unlist(lapply(1:nrow(linpred), function(i) {
      data.stan$d * log(hllogis(data.stan$t, sigma[i], exp(linpred[i, ]))) +
        log(1 - pllogis(data.stan$t, sigma[i], exp(linpred[i, ])))
    })), nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(
    data.stan$d * log(hllogis(data.stan$t, sigma.hat, exp(linpred.hat))) +
      log(1 - pllogis(data.stan$t, sigma.hat, exp(linpred.hat))),
    nrow = 1)
  # Number of parameters (for AIC): shape, scale + covariates
  npars <- 2 + sum(
    1 - apply(data.stan$X, 2, function(x) all(x == 0)))
  list(logf=logf,logf.hat=logf.hat,npars=npars,f=NULL,f.bar=NULL,s=NULL,s.bar=NULL)
}

lik_rps <- function(x,linpred,linpred.hat,model,data.stan) {
  gamma <- rstan::extract(model)$gamma
  gamma.hat <- apply(gamma, 2, median)
  # Needs to reformat linpred.hat to a vector for consistency
  linpred.hat <- as.numeric(linpred.hat)
  logf <- data.stan$d * (-log(data.stan$t) + log(gamma %*% t(data.stan$DB)) + gamma %*% t(data.stan$B) + linpred) -
    exp(gamma %*% t(data.stan$B) + linpred)
  logf.hat <- t(
    data.stan$d * (-log(data.stan$t) + log(data.stan$DB %*% gamma.hat) + data.stan$B %*% gamma.hat + linpred.hat) -
      exp(data.stan$B %*% gamma.hat + linpred.hat))
  # Number of parameters (for AIC): gamma + covariates
  npars <- length(gamma.hat) + sum(
    apply(data.stan$X, 2, function(x) 1 - all(x == 0)))
  list(logf=logf,logf.hat=logf.hat,npars=npars,f=NULL,f.bar=NULL,s=NULL,s.bar=NULL)
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
    }
    method <- "mle"
  }
  
  # 'mle' can implement all the possible models, excpet the PolyWeibull
  # In this case, I choose to *stop* execution, rather than falling back to 'hmc'!
  if (method == "mle") {
    if(!all(distr3 %in% availables[[method]])) {
      stop(paste0("The Poly-Weibull model is only implemented under method='hmc'.
       Please set this option in your call to 'fit.models'"))
    }
  }
}
