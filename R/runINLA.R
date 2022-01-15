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
  availables <- survHE:::load_availables()
  # Uses the helper 'manipulated_distributions' to create the vectors distr, distr3 and labs
  d3 <- survHE:::manipulate_distributions(x)$distr3
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
    # Default priors
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
  
  # Gets the name of the time variable so that the times can be rescaled in [0-1]
  time_name=all.vars(formula)[1]
  d_name=all.vars(formula)[2]
  time_max=data %>% select(!!sym(time_name)) %>% with(max(.))
  
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
  model <- INLA::inla(formula,family=paste0(d,"surv"),#data=data,
                      # Rescales the times in [0-1] on the fly
                      data=data %>% mutate(!!sym(time_name) := !!sym(time_name) / max(!!sym(time_name))),
                      control.compute=list(config=TRUE,dic=TRUE),
                      control.inla=list(int.strategy="grid",dz=dz,diff.logdens=diff.logdens),
                      control.fixed=control.fixed,control.family=cf,verbose=verbose
  )
  
  # Now re-writes the formula in general terms (without linking to INLA::inla.surv)
  formula <- as.formula(gsub("INLA::inla.surv","Surv",deparse(formula)))
  # Adds a field used in 'make.surv' to indicate the model used
  model_name <- d3
  # And reset if original name *was* 'WeibullPH'
  if(d3=="wph") {d="weibullPH"} 
  
  # Now computes the densities using the helper functions
  t=data %>% select(!!sym(time_name)) %>% pull
  d=data %>% select(!!sym(d_name)) %>% pull
  out=do.call(what=paste0("lik_",d3,"_inla"),args=list(model,nsim=1000,time_max,t,d,formula,data))
  # Extracts relevant variables from the list (See if there's a better way to do it!)
  logf=out$logf
  logf.hat=out$logf.hat
  npars=out$npars
  # Now computes the log-likelihood
  loglik <- apply(logf,1,sum)
  loglik.bar <- apply(logf.hat,1,sum)
  
  D.theta <- -2*loglik
  D.bar <- -2*loglik.bar
  pD <- mean(D.theta) - D.bar
  pV <- 0.5*var(D.theta)
  dic <- mean(D.theta)+pD
  dic2 <- mean(D.theta)+pV
  # Approximates AIC & BIC using the mean deviance and the number of nominal parameters
  aic <- D.bar+2*npars
  bic <- D.bar+npars*log(nrow(data))
  
  # Finally returns the output
  list(
    model=model,
# These would be the values automatically computed by INLA - but they are not OK because the model is on the [0-1]
# rather than the original time scales!
#    aic=2*model$dic$p.eff+model$dic$deviance.mean,
#    bic=model$dic$deviance.mean+model$dic$p.eff*log(model$size.linear.predictor$n),
#    dic=model$dic$dic,
    aic=aic,
    bic=bic,
    dic=dic,
    dic2=dic2,
    time2run=model$cpu.used["Total"],
    model_name=model_name
  )
}

# Little utility function to compute the mode of a vector
Mode <- function(x) {
  if (is.numeric(x)) {
    x_table <- table(x)
    return(as.numeric(names(x_table)[which.max(x_table)]))
  }
}

### These are utility functions to compute the log-density for the models fitted by INLA
lik_wei_inla <- function(model,nsim,time_max,t,d,formula,data) {
  shape=INLA::inla.rmarginal(nsim,model$marginals.hyperpar[[1]])
  beta=lapply(1:nrow(model$summary.fixed),function(i) {
    INLA::inla.rmarginal(nsim,model$marginals.fixed[[i]])
  })
  names(beta)=colnames(model$model.matrix)
  beta=beta %>% bind_cols()
  # If there is an intercept, then rescale it using the suitable back-transformation to account for the fact that the
  # original INLA model is run on times in [0-1], instead of the original ones
  if(grep("(Intercept)",colnames(model$model.matrix))>0) {
    beta[,grep("(Intercept)",colnames(model$model.matrix))]=beta[,grep("(Intercept)",colnames(model$model.matrix))]-log(time_max)
  }
  shape.hat = mean(shape)
  beta.hat=beta %>% summarise_all(mean) %>% as.numeric()
  linpred=as.matrix(beta)%*%t(model.matrix(formula,data))
  linpred.hat=beta.hat%*%t(model.matrix(formula,data))
  
  logf <- matrix(
    unlist(lapply(1:nrow(linpred), function(i) {
      d * log(hweibull(t, shape = shape[i], scale = exp(-linpred[i, ]))) +
        log(1 - pweibull(t, shape[i], scale = exp(-linpred[i, ])))
    })), 
    nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(
    d * log(hweibull(t, shape.hat, exp(-linpred.hat))) +
      log(1 - pweibull(t, shape.hat, exp(-linpred.hat))),
    nrow = 1)
  # Number of parameters (for AIC): shape, rate + covariates
  npars <- 2 + (length(all.vars(formula))-1)
  list(logf=logf,logf.hat=logf.hat,npars=npars,f=NULL,f.bar=NULL)
}

lik_wph_inla <- function(model,nsim,time_max,t,d,formula,data) {
  shape=INLA::inla.rmarginal(nsim,model$marginals.hyperpar[[1]])
  beta=lapply(1:nrow(model$summary.fixed),function(i) {
    INLA::inla.rmarginal(nsim,model$marginals.fixed[[i]])
  })
  names(beta)=colnames(model$model.matrix)
  beta=beta %>% bind_cols()
  # If there is an intercept, then rescale it using the suitable back-transformation to account for the fact that the
  # original INLA model is run on times in [0-1], instead of the original ones
  if(grep("(Intercept)",colnames(model$model.matrix))>0) {
    beta[,grep("(Intercept)",colnames(model$model.matrix))]=(beta[,grep("(Intercept)",colnames(model$model.matrix))]+log(time_max))*(-shape)
  }
  shape.hat = mean(shape)
  beta.hat=beta %>% summarise_all(mean) %>% as.numeric()
  linpred=as.matrix(beta)%*%t(model.matrix(formula,data))
  linpred.hat=beta.hat%*%t(model.matrix(formula,data))
  
  logf <- matrix(
    unlist(lapply(1:nrow(linpred), function(i) {
      d * log(hweibullPH(t, shape = shape[i], scale = exp(linpred[i, ]))) +
        log(1 - pweibullPH(t, shape[i], scale = exp(linpred[i, ])))
    })), 
    nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(
    d * log(hweibullPH(t, shape.hat, exp(linpred.hat))) +
      log(1 - pweibullPH(t, shape.hat, exp(linpred.hat))),
    nrow = 1)
  # Number of parameters (for AIC): shape, rate + covariates
  npars <- 2 + (length(all.vars(formula))-1)
  list(logf=logf,logf.hat=logf.hat,npars=npars,f=NULL,f.bar=NULL)
}

lik_exp_inla <- function(model,nsim,time_max,t,d,formula,data) {
  beta=lapply(1:nrow(model$summary.fixed),function(i) {
    INLA::inla.rmarginal(nsim,model$marginals.fixed[[i]])
  })
  names(beta)=colnames(model$model.matrix)
  beta=beta %>% bind_cols()
  # If there is an intercept, then rescale it using the suitable back-transformation to account for the fact that the
  # original INLA model is run on times in [0-1], instead of the original ones
  if(grep("(Intercept)",colnames(model$model.matrix))>0) {
    beta[,grep("(Intercept)",colnames(model$model.matrix))]=(beta[,grep("(Intercept)",colnames(model$model.matrix))]-log(time_max))
  }
  beta.hat=beta %>% summarise_all(mean) %>% as.numeric()
  linpred=as.matrix(beta)%*%t(model.matrix(formula,data))
  linpred.hat=beta.hat%*%t(model.matrix(formula,data))
  
  logf <- matrix(
    unlist(lapply(1:nrow(linpred), function(i) {
      d * log(hexp(t, rate = exp(linpred[i, ]))) +
        log(1 - pexp(t, rate = exp(linpred[i, ])))
    })), 
    nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(
    d * log(hexp(t, exp(linpred.hat))) +
      log(1 - pexp(t, exp(linpred.hat))),
    nrow = 1)
  # Number of parameters (for AIC): shape, rate + covariates
  npars <- 1 + (length(all.vars(formula))-1)
  list(logf=logf,logf.hat=logf.hat,npars=npars,f=NULL,f.bar=NULL)
}

lik_lno_inla <- function(model,nsim,time_max,t,d,formula,data) {
  prec=INLA::inla.rmarginal(nsim,model$marginals.hyperpar[[1]])
  sdlog=sqrt(1/prec)
  beta=lapply(1:nrow(model$summary.fixed),function(i) {
    INLA::inla.rmarginal(nsim,model$marginals.fixed[[i]])
  })
  names(beta)=colnames(model$model.matrix)
  beta=beta %>% bind_cols()
  # If there is an intercept, then rescale it using the suitable back-transformation to account for the fact that the
  # original INLA model is run on times in [0-1], instead of the original ones
  if(grep("(Intercept)",colnames(model$model.matrix))>0) {
    beta[,grep("(Intercept)",colnames(model$model.matrix))]=(beta[,grep("(Intercept)",colnames(model$model.matrix))]+log(time_max))
  }
  sdlog.hat = mean(sdlog)
  beta.hat=beta %>% summarise_all(mean) %>% as.numeric()
  linpred=as.matrix(beta)%*%t(model.matrix(formula,data))
  linpred.hat=beta.hat%*%t(model.matrix(formula,data))
  
  logf <- matrix(
    unlist(lapply(1:nrow(linpred), function(i) {
      d * log(hlnorm(t, meanlog=(linpred[i, ]), sdlog=sdlog[i])) +
        log(1 - plnorm(t, meanlog=(linpred[i, ]), sdlog=sdlog[i]))
    })), 
    nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(
    d * log(hlnorm(t, (linpred.hat), sdlog.hat)) +
      log(1 - plnorm(t, (linpred.hat), sdlog.hat)),
    nrow = 1)
  # Number of parameters (for AIC): shape, rate + covariates
  npars <- 2 + (length(all.vars(formula))-1)
  list(logf=logf,logf.hat=logf.hat,npars=npars,f=NULL,f.bar=NULL)
}

lik_llo_inla <- function(model,nsim,time_max,t,d,formula,data) {
  shape=INLA::inla.rmarginal(nsim,model$marginals.hyperpar[[1]])
  beta=lapply(1:nrow(model$summary.fixed),function(i) {
    INLA::inla.rmarginal(nsim,model$marginals.fixed[[i]])
  })
  names(beta)=colnames(model$model.matrix)
  beta=beta %>% bind_cols()
  # If there is an intercept, then rescale it using the suitable back-transformation to account for the fact that the
  # original INLA model is run on times in [0-1], instead of the original ones
  if(grep("(Intercept)",colnames(model$model.matrix))>0) {
    beta[,grep("(Intercept)",colnames(model$model.matrix))]=beta[,grep("(Intercept)",colnames(model$model.matrix))]-log(time_max)
  }
  shape.hat = mean(shape)
  beta.hat=beta %>% summarise_all(mean) %>% as.numeric()
  linpred=as.matrix(beta)%*%t(model.matrix(formula,data))
  linpred.hat=beta.hat%*%t(model.matrix(formula,data))
  
  logf <- matrix(
    unlist(lapply(1:nrow(linpred), function(i) {
      d * log(hllogis(t, shape = shape[i], scale = exp(-linpred[i, ]))) +
        log(1 - pllogis(t, shape[i], scale = exp(-linpred[i, ])))
    })), 
    nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(
    d * log(hllogis(t, shape.hat, exp(-linpred.hat))) +
      log(1 - pllogis(t, shape.hat, exp(-linpred.hat))),
    nrow = 1)
  # Number of parameters (for AIC): shape, rate + covariates
  npars <- 2 + (length(all.vars(formula))-1)
  list(logf=logf,logf.hat=logf.hat,npars=npars,f=NULL,f.bar=NULL)
}

lik_gom_inla <- function(model,nsim,time_max,t,d,formula,data) {
  shape=INLA::inla.rmarginal(nsim,model$marginals.hyperpar[[1]])/time_max
  beta=lapply(1:nrow(model$summary.fixed),function(i) {
    INLA::inla.rmarginal(nsim,model$marginals.fixed[[i]])
  })
  names(beta)=colnames(model$model.matrix)
  beta=beta %>% bind_cols()
  if(grep("(Intercept)",colnames(model$model.matrix))>0) {
    beta[,grep("(Intercept)",colnames(model$model.matrix))]=beta[,grep("(Intercept)",colnames(model$model.matrix))]-log(time_max)
  }
  shape.hat = mean(shape)
  beta.hat=beta %>% summarise_all(mean) %>% as.numeric()
  linpred=as.matrix(beta)%*%t(model.matrix(formula,data))
  linpred.hat=beta.hat%*%t(model.matrix(formula,data))
  
  logf <- matrix(
    unlist(lapply(1:nrow(linpred), function(i) {
      d * log(hgompertz(t, shape = shape[i], rate = exp(linpred[i, ]))) +
        log(1 - pgompertz(t, shape[i], rate = exp(linpred[i, ])))
    })), 
    nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(
    d * log(hgompertz(t, shape.hat, exp(linpred.hat))) +
      log(1 - pgompertz(t, shape.hat, exp(linpred.hat))),
    nrow = 1)
  # Number of parameters (for AIC): shape, rate + covariates
  npars <- 2 + (length(all.vars(formula))-1)
  list(logf=logf,logf.hat=logf.hat,npars=npars,f=NULL,f.bar=NULL,s=NULL,s.bar=NULL)
}
