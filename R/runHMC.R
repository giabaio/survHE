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
runHMC <- function(x,exArgs) {
  # First checks whether INLA is installed (it's only a suggestion, not a full dependency)
  if (!isTRUE(requireNamespace("rstan", quietly = TRUE))) {
    stop("You need to install the R package 'rstan'. Please run in your R terminal:\n install.packages('rstan')")
  }
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
  
  # Recomputes the three-letters code for the distributions and the HMC-specific name
  d <- names(availables[[method]][match(d3, availables[[method]])])

  # Loads the pre-compiled models
  dso <- stanmodels[[d]]

  # Create the data list
  data.stan <- make_data_stan(formula,data,d3,exArgs)
  
  # Now runs Stan to sample from the posterior distributions
  tic <- proc.time()
  model <- rstan::sampling(dso,data.stan,chains=chains,iter=iter,warmup=warmup,thin=thin,seed=seed,
                           control=control,pars=pars,include=include,cores=cores,init=init)
  toc <- proc.time()-tic
  time_survHE <- toc[3]
  # rstan does have its way of computing the running time, but it may not be the actual one when running multiple
  # chains. 
  time_stan <- sum(rstan::get_elapsed_time(model))
  
  # Uses the helper function to compute the *IC
  ics <- compute_ICs_stan(model,d3,data.stan)
  
  # Finally returns the output
  list(
    model=model,
    aic=ics$aic,
    bic=ics$bic,
    dic=ics$dic,
    dic=ics$dic2,
    time2run=pmin(time_survHE,time_stan),
    data.stan=data.stan,
    save.stan=save.stan
  )
}