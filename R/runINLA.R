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
runINLA <- function(x,exArgs) {
  # First checks whether INLA is installed (it's only a suggestion, not a full dependency)
  if (!isTRUE(requireNamespace("INLA", quietly = TRUE))) {
    stop("You need to install the packages 'INLA'. Please run in your R terminal:\n install.packages('INLA', repos='https://www.math.ntnu.no/inla/R/stable')")
  }
  # If INLA is installed but not loaded then attach the Namespace (so that all the relevant functions are available)
  # if (isTRUE(requireNamespace("INLA", quietly = TRUE))) {
  #   if (!is.element("INLA", (.packages()))) {
  #     attachNamespace("INLA")
  #   }
  # }
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
  # If 'control.family' is specified for the distribution currently used, then use the values in
  # 'exArgs'. If not, or if specified but for another distribution, use INLA defaults
  cf <- INLA::inla.set.control.family.default()
  if(exists("control.family",where=exArgs)) {
    if(d==names(exArgs$control.family)) {
      cf <- exArgs$control.family[[d]]
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
  model <- INLA::inla(formula,family=paste0(d,"surv"),data=data,control.compute=list(config=TRUE,dic=TRUE),
                      control.inla=list(int.strategy="grid",dz=dz,diff.logdens=diff.logdens),
                      control.fixed=control.fixed,control.family=cf,verbose=verbose
  )
  
  # Now re-writes the formula in general terms (without linking to INLA::inla.surv)
  formula <- as.formula(gsub("INLA::inla.surv","Surv",deparse(formula)))
  # Adds a field used in 'make.surv' to indicate the model used
  model$dlist$name <- d

    # Finally returns the output
  list(
    model=model,
    aic=2*model$dic$p.eff+model$dic$deviance.mean,
    bic=model$dic$deviance.mean+model$dic$p.eff*log(model$size.linear.predictor$n),
    dic=model$dic$dic,
    time2run=model$cpu.used["Total"]
  )
}