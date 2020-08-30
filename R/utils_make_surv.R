#' Helper function to manipulate the survival distributions and 
#' compute relevant quantities (*ICs, survival curves, etc)
#' 
#' @param x A string with the distribution abbreviation ('distr3')
#' @return \item{list}{A list containing the modified name of the 
#' distribution, the acronym (3-letters abbreviation), or the
#' labels (humane-readable name)}.
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso fit.models
#' @references Baio (2020). survHE
#' @keywords 
make_sim_mle <- function(m,t,X,nsim,newdata,dist,...) {
  # Simulates from the distribution of the model parameters - takes 100000 bootstrap samples
  nboot=100000
  B=ifelse(nsim<nboot,nboot,nsim)
  if(is.null(newdata)) {
    # NB: 'flexsurv' needs to exclude the intercept
    if(grep("Intercept",colnames(X))>0) {
      # If the intercept is part of the design matrix X then remove it
      X=matrix(X[,-grep("Intercept",colnames(X))],nrow=nrow(X))
    } 
    sim=flexsurv::normboot.flexsurvreg(m,B=B,X=as.matrix(X))
  } else {
    sim <- lapply(1:nrow(X),function(i) flexsurv::normboot.flexsurvreg(m,B=B,newdata=newdata[[i]]))
  }
  # Then if 'nsim'=1, then take the average over the bootstrap samples
  if(nsim==1) {
    if(nrow(X)==1) {
      sim=list(sim %>% as_tibble() %>% summarise_all(mean) %>% as.matrix(.,ncol=ncol(X)))
    } else {
      sim=lapply(sim,function(x) x %>% as_tibble() %>% summarise_all(mean) %>% as.matrix(.,ncol=ncol(X)))
    }
  } 
  # If nsim<=5000 (number of bootstrap samples), then samples only 'nsim' of them
  if(nsim>1 & nsim<nboot) {
    if (nrow(X)==1) {
      sim=list(sim %>% as_tibble %>% sample_n(ifelse(nsim<nboot,nsim,B),replace=FALSE) %>% 
        as.matrix(.,nrow=nsim,ncol=ncol(X)))
    } else {
      sim=lapply(sim,function(x) x %>% as_tibble %>% sample_n(ifelse(nsim<nboot,nsim,B),replace=FALSE) %>% 
                   as.matrix(.,nrow=nsim,ncol=ncol(X)))
    }
  }
  return(sim)
}

make_sim_inla <- function(m,t,X,nsim,newdata,dist,...) {
  # Simulates from the distribution of the model parameters
  
  exArgs_inla <- list(...)
  
  # If 'nsim'=1 then use the point estimates for the coefficients
  if(nsim==1) {
    beta=m$summary.fixed[,"mean"]
    linpred=beta%*% t(X)
    if(!is.null(m$summary.hyperpar)) {
      alpha=m$summary.hyperpar$mean
    } else {
      alpha=0
    }
    # Now uses the helper function to rescale the parameters and retrieve the correct simulations
    sim <- lapply(1:ncol(linpred),function(x) rescale.inla(linpred[,x],alpha,dist))
    
  } else {
    # Selects the variables to sample from the posterior distributions (excludes the linear predictor)
    selection <- eval(parse(text=paste0('list(Predictor=-c(1:',nrow(exArgs_inla$data),'),',
                                        paste0("`",colnames(model.matrix(exArgs_inla$formula,exArgs_inla$data)),"`",'=1',collapse=','),')')))
    # Runs INLA to get nsim samples from the posterior of the parameters 
    jpost <- INLA::inla.posterior.sample(nsim,m,selection=selection)
    # Then computes the linear predictor
    beta <- matrix(unlist(lapply(jpost,function(x){x$latent})),ncol=length(jpost[[1]]$latent),nrow=nsim,byrow=T)
    linpred <- beta %*% t(X)
    # And the ancillary parameters (which may or may not exist for a given model)
    if(!is.null(jpost[[1]]$hyperpar)) {
      alpha <- matrix(unlist(lapply(jpost,function(x){x$hyperpar})),ncol=length(jpost[[1]]$hyperpar),nrow=nsim,byrow=T)
    } else {
      alpha <- 0
    }
    
    # Now uses the helper function to rescale the parameters and retrieve the correct simulations
    sim <- lapply(1:ncol(linpred),function(x) rescale.inla(linpred[,x],alpha,dist))
  }
  return(sim)
}

make_sim_hmc <- function(m,t,X,nsim,newdata,dist,...) {

  # Extracts the model object from the survHE output
  iter_stan <- m@stan_args[[1]]$iter
  
  linpred <- rstan::extract(m)$beta %*% t(X)
  # Stores the values returned by rstan into the list 'sim'
  sim <- lapply(1:nrow(X),function(x) {
    do.call(paste0("rescale_hmc_",dist),
            args=list(m,X,linpred[,x]))
  })

  if(nsim>iter_stan) {
    stop("Please select a value for 'nsim' that is less than or equal to the value set in the call to 'fit.models'")
  }
  if(nsim==1) {
    # If the user requested only 1 simulation, then take the mean value
    sim <- lapply(sim,function(x) as.matrix(as_tibble(x) %>% summarise_all(mean),nrow=1,ncol=ncol(x)))
  }
  if(nsim>1 & nsim<iter_stan) {
    # If the user selected a number of simulation < the one from rstan, then select a random sample 
    sim <- lapply(sim,function(x) as.matrix(as_tibble(x) %>% sample_n(nsim,replace=FALSE),nrow=nsim,ncol=ncol(x)))
  }
  return(sim)
}

rescale_hmc_exp <- function(m,X,linpred){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # Exponential distribution
  sim <- as.matrix(exp(linpred)) # exp(rstan::extract(m)$beta %*% t(X)); 
  colnames(sim) <- "rate"
  return(sim)
}

rescale_hmc_wei <- function(m,X,linpred){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # Weibull distribution
  shape <- as.numeric(rstan::extract(m)$alpha)
  scale <- exp(linpred) #exp(rstan::extract(m)$beta %*% t(X))
  sim <- cbind(shape,scale) 
  colnames(sim) <- c("shape","scale")
  return(sim)
}

rescale_hmc_wph <- function(m,X,linpred){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # Weibull PH distribution
  shape <- as.numeric(rstan::extract(m)$alpha)
  scale <- exp(linpred) #exp(rstan::extract(m)$beta %*% t(X))
  sim <- cbind(shape,scale) 
  colnames(sim) <- c("shape","scale")
  return(sim)
}

rescale_hmc_gom <- function(m,X,linpred){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # Gompertz distribution
  shape <- as.numeric(rstan::extract(m)$alpha)
  rate <- exp(linpred) #exp(rstan::extract(m)$beta %*% t(X))
  sim <- cbind(shape,rate) 
  colnames(sim) <- c("shape","rate")
  return(sim)
}

rescale_hmc_gam <- function(m,X,linpred){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # Gamma distribution
  shape <- as.numeric(rstan::extract(m)$alpha)
  rate <- 1/exp(linpred) #exp(rstan::extract(m)$beta %*% t(X))
  sim <- cbind(shape,rate) 
  colnames(sim) <- c("shape","rate")
  return(sim)
}

rescale_hmc_gga <- function(m,X,linpred){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # Generalised Gamma distribution
  Q <- as.numeric(rstan::extract(m)$Q)
  sigma <- as.numeric(rstan(extract(m)$sigma))
  mu <- exp(linpred) #exp(rstan::extract(m)$beta %*% t(X))
  sim <- cbind(mu,sigma,Q) 
  colnames(sim) <- c("mu","sigma","Q")
  return(sim)
}

rescale_hmc_gef <- function(m,X,linpred){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # Generalised F distribution
  Q <- as.numeric(rstan::extract(m)$Q)
  P <- as.numeric(rstan::extract(m)$P)
  sigma <- as.numeric(rstan(extract(m)$sigma))
  mu <- exp(linpred) #exp(rstan::extract(m)$beta %*% t(X))
  sim <- cbind(mu,sigma,Q) 
  colnames(sim) <- c("mu","sigma","Q","P")
  return(sim)
}

rescale_hmc_lno <- function(m,X,linpred){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # logNormal distribution
  sdlog <- as.numeric(rstan::extract(m)$alpha)
  meanlog <- linpred #rstan::extract(m)$beta %*% t(X)
  sim <- cbind(meanlog,sigma) 
  colnames(sim) <- c("meanlog","sdlog")
  return(sim)
}

rescale_hmc_llo <- function(m,X){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # logNormal distribution
  shape <- as.numeric(rstan::extract(m)$alpha)
  scale <- exp(rstan::extract(m)$beta %*% t(X))
  sim <- cbind(shape,scale) 
  colnames(sim) <- c("shape","scale")
  return(sim)
}

rescale_hmc_rps <- function(m,X,linpred) {
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # RPS
  gamma <- rstan::extract(m)$gamma
  # If X has an intercept needs to remove it as RPS doesn't have one
  if(any(grepl("Intercept",colnames(X)))) {
    X <- as.matrix(as_tibble(X) %>% select(-`(Intercept)`))
  }
  offset <- linpred #rstan::extract(m)$beta %*% t(X) 
  sim <- cbind(offset,gamma)
  colnames(sim) <- c("offset",paste0("gamma",1:ncol(gamma)))
  return(sim)
}


rescale.inla <- function(linpred,alpha,distr) {
  if (distr=="wei") {
    shape <- alpha
    # NB: As of Jan 11 2017, there's a little mistake in INLA and so need to minus the linpred HERE
    scale <- exp(-linpred)
    sim <- cbind(shape,scale) 
    colnames(sim) <- c("shape","scale")
  }
  if (distr=="wph") {
    shape <- alpha
    scale <- exp(linpred)
    sim <- cbind(shape,scale) 
    colnames(sim) <- c("shape","scale")
  }
  if (distr=="exp") {
    sim <- as.matrix(exp(linpred))
    colnames(sim) <- "rate" 
  }
  if (distr=="llo") {
    shape <- alpha
    scale <- exp(linpred)
    sim <- cbind(shape,scale) 
    colnames(sim) <- c("shape","scale")
  }
  if (distr=="lno") {
    mulog <- linpred
    sdlog <- 1/sqrt(alpha)
    sim <- cbind(mulog,sdlog) 
    colnames(sim) <- c("meanlog","sdlog")
  }
  return(sim)
}

compute_surv_curve <- function(sim,exArgs,nsim,dist,t) {  
  # Computes the survival curves
  args <- args_surv()
  distr <- manipulate_distributions(dist)$distr 
  
  if(dist=="rps") {
    # RPS-related options
    if(exists("scale",where=exArgs)) {scale=exArgs$scale} else {scale="hazard"}
    if(exists("timescale",where=exArgs)) {timescale=exArgs$timescale} else {timescale="log"}
    #if(exists("offset",where=exArgs)) {offset=exArgs$offset} else {offset=0}
    if(exists("log",where=exArgs)) {log=exArgs$log} else {log=FALSE}
    knots <- exArgs$data.stan$knots
    mat <- lapply(sim,function(x) {
      gamma=as_tibble(x) %>% select(contains("gamma"))
      offset=as_tibble(x) %>% select(offset)
      matrix(
        unlist(
          lapply(1:nsim,function(i){
            1-do.call(psurvspline,args=list(
              q=t,
              gamma=as.numeric(gamma %>% slice(i)),
              beta=0,
              X=0,
              knots=knots,
              scale=scale,
              timescale=timescale,
              offset=as.numeric(offset%>% slice(i)),
              log=log
            ))
          })
        )
      )
    })
  } else {
    mat <- lapply(sim,function(x) {
      matrix(
        unlist(
          lapply(1:nsim,function(i) {
            1-do.call(paste0("p",distr),args=eval(parse(text=args[[distr]])))
          })
        ),nrow=length(t),ncol=nsim,byrow=FALSE
      )
    })
  }
  for (i in 1:length(mat)){colnames(mat[[i]])=paste0("S_",1:nsim)}
  mat <- mat %>% lapply(function(x) bind_cols(as_tibble(t),as_tibble(x)) %>% rename(t=value))

  return(mat)
}

# Utility function to define the arguments needed to compute the cumulative distribution,
## with which to derive the survival function
args_surv <- function() {
    list(
      exp='list(t,rate=x[,"rate"][i])',
      ##weibull='list(t,shape=(x %>% as_tibble %>% slice(i))$shape,scale=(x %>% as_tibble %>% slice(i))$scale)',
      weibull='list(t,shape=x[,"shape"][i],scale=x[,"scale"][i])',
      weibullPH='list(t,shape=x[,"shape"][i],scale=x[,"scale"][i])',
      gamma='list(t,shape=x[,"shape"][i],rate=x[,"rate"][i])',
      gengamma='list(t,mu=x[,"mu"][i],sigma=x[,"sigma"][i],Q=x[,"Q"][i])',
      gompertz='list(t,shape=x[,"shape"][i],rate=x[,"rate"][i])',
      lnorm='list(t,meanlog=x[,"meanlog"][i],sdlog=x[,"sdlog"][i])',
      llogis='list(t,shape=x[,"shape"][i],scale=x[,"scale"][i])',
      genf='list(t,mu=x[,"mu"][i],sigma=x[,"sigma"][i],Q=x[,"Q"][i],P=x[,"P][i])',
      rps='list(t,gamma=x[,"gamma"][i],knots=m$knots,scale=scale,timescale=timescale,offset=offset,log=log)',
      survspline='list(t,gamma=x[,"gamma"][i],knots=m$knots,scale=scale,timescale=timescale,offset=offset,log=log)'
  )
}



#' Helper function to create the covariates profile to use in the 
#' computation of the survival curves
#' 
#' @param formula a formula specifying the model to be used, in the form
#' \code{Surv(time,event)~treatment[+covariates]} in flexsurv terms, or
#' \code{inla.surv(time,event)~treatment[+covariates]} in INLA terms.
#' @param data A data frame containing the data to be used for the analysis.
#' This must contain data for the 'event' variable. In case there is no
#' censoring, then \code{event} is a column of 1s.
#' @param newdata a list (of lists), specifiying the values of the covariates
#' at which the computation is performed. For example
#' \code{list(list(arm=0),list(arm=1))} will create two survival curves, one
#' obtained by setting the covariate \code{arm} to the value 0 and the other by
#' setting it to the value 1. In line with \code{flexsurv} notation, the user
#' needs to either specify the value for *all* the covariates or for none (in
#' which case, \code{newdata=NULL}, which is the default). If some value is
#' specified and at least one of the covariates is continuous, then a single
#' survival curve will be computed in correspondence of the average values of
#' all the covariates (including the factors, which in this case are expanded
#' into indicators).
#' @return \item{X}{A matrix with the covariates profile selected to compute
#' the survival curves}.
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso make.surv
#' @references Baio (2020). survHE
#' @keywords Parametric survival models
make_profile_surv <- function(formula,data,newdata) {
  # Checks how many elements are given in 'newdata'
  n.elements <- ifelse(is.null(newdata),0,length(newdata))
  n.provided <- unlist(lapply(newdata,function(x) length(x)))
  
  # Temporarily re-writes the model formula to avoid issues with naming
  formula_temp <- update(formula,paste(all.vars(formula,data)[1],"~",all.vars(formula,data)[2],"+."))
  # Creates a tibble with all the covariates in their original format
  covs <- data %>% model.frame(formula_temp,.) %>% as_tibble(.) %>%  select(-c(1:2)) %>% 
    rename_if(is.factor,.funs=~gsub("as.factor[( )]","",.x)) %>% 
    rename_if(is.factor,.funs=~gsub("[( )]","",.x)) %>% 
    bind_cols(as_tibble(model.matrix(formula_temp,data)) %>% select(contains("Intercept"))) %>% 
    select(`(Intercept)`,everything())
  ncovs <- ncol(covs) - 1
  # Selects the subset of categorical covariates
  is.fac <- covs %>% select(where(is.factor))
  fac.levels=lapply(is.fac,levels)
  nfacts <- covs %>% select(where(is.factor)) %>% with(ncol(.))
  
  # If formula is in 'inla' terms now change it back to 'flexsurv' terms
  formula_temp <- as.formula(gsub("inla.surv","Surv",deparse(formula)))
  # Computes the "average" profile of the covariates
  X <- data %>% model.matrix(formula_temp,.) %>% as_tibble(.) %>% summarise_all(mean) 
  # If there's at least one factor with more than 2 levels, then do *not* rename the columns to the simpler version
  # which only has the name of the variable (rather than the combination, eg 'groupMedium' that R produces)
  if(all(unlist(lapply(fac.levels,length))<2)) {colnames(X)=colnames(covs)}
  
  # The way the object X *must* be formatted depends on which way it's been generated.
  # The point is that it *always* has to be a matrix for other functions to process
  if(n.elements==0){
    # If all the covariates are factors, then get survival curves for *all* the combinations
    if(nfacts==ncovs & nfacts>0) {
      X=unique(model.matrix(formula,data))
      # If there's at least one factor with more than 2 levels, then do *not* rename the columns to the simpler version
      # which only has the name of the variable (rather than the combination, eg 'groupMedium' that R produces)
      if(all(unlist(lapply(fac.levels,length))<2)) {colnames(X)=colnames(covs)}
      # X=apply(
      #   covs %>% unique %>%  mutate(across(.cols=everything(),.fns=factor)) %>% 
      #     mutate(across(.cols=everything(),.fns=as.numeric)),2,as.numeric
      # )
      #X <- apply(covs %>% unique(.),2,as.numeric)
    } else {
      X <- as.matrix(X,nrow=nrow(X),ncol=ncol(X))
    }
  }
  
  # If 'newdata' provides a specific (set of) profile(s), then use that
  if (n.elements>=1) {
    # This means that n.provided will also be > 0 but needs to check it's the right number and if not, stop
    if (!all(n.provided==ncovs)) {
      stop("You need to provide data for *all* the covariates specified in the model, in the list 'newdata'")
    } else {
      # Creates a design matrix containing the information provided in 'newdata'
      X <- bind_rows(lapply(newdata,function(x) as_tibble(x)))
      # But if there's an intercept in the model matrix, then add it
      if ("(Intercept)" %in% colnames(model.matrix(formula,data))) {
        X <- X %>% mutate(`(Intercept)`=1) %>% select(`(Intercept)`,everything())
      }
      X <- as.matrix(X)
    }
  }
  # Returns the output (the matrix with the covariates profile)
  return(X)
}
