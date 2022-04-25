#' Helper function to make the simulations for the survival curves using MLE
#' for a given formula and dataset
#' 
#' @param m The output of a 'survHE' object including the model
#' @param t A vector of times for which the survival curves are to be
#' computed
#' @param X The covariates profile for which to compute the survival curve
#' @param nsim The number of simulations to be computed
#' @param newdata A list with the "new data". This is a specific profile of
#' covariates in correspondence of which to compute the survival curves.
#' @param dist The string identifying the abbreviated name for the 
#' underlying model used 
#' @return \item{sim}{A list containing the generated simulations.}
#' @author Gianluca Baio
#' @seealso make.surv
#' @references Baio (2020). survHE
#' @keywords MLE
#' @noRd 
make_sim_mle <- function(m,t,X,nsim,newdata,dist,summary_stat,...) {
  # Simulates from the distribution of the model parameters - takes 100000 bootstrap samples
  nboot=100000
  B=ifelse(nsim<nboot,nboot,nsim)
  if(is.null(newdata)) {
    #####X=X %>% as_tibble()
    # NB: 'flexsurv' needs to exclude the intercept
    if(grep("Intercept",colnames(X))>0) {
      # If the intercept is part of the design matrix X then remove it (to make normboot work!)
      X=matrix(X[,-grep("Intercept",colnames(X))],nrow=nrow(X))
    } 
    # If X has only one row, needs to create a list, with length equal to the number of profiles (=nrow(X))
    if(nrow(X)==1) {
      sim=list(flexsurv::normboot.flexsurvreg(m,B=B,X=as.matrix(X)))
    } else {
      # Otherwise normboot will take care of it with the proper length for the automatically created list
      sim=flexsurv::normboot.flexsurvreg(m,B=B,X=as.matrix(X))
    }
  } else {
    # If there are newdata, then create the list of sims using it
    sim <- lapply(1:nrow(X),function(i) flexsurv::normboot.flexsurvreg(m,B=B,newdata=newdata[[i]]))
  }
  # Then if 'nsim'=1, then take the average over the bootstrap samples. 
  if(nsim==1) {
    sim=lapply(sim,function(x) x %>% as_tibble() %>% summarise_all(summary_stat) %>% as.matrix(.,ncol=ncol(X)))
  } 
  # If nsim<=5000 (number of bootstrap samples), then samples only 'nsim' of them
  if(nsim>1 & nsim<nboot) {
    sim=lapply(sim,function(x) x %>% as_tibble %>% sample_n(ifelse(nsim<nboot,nsim,B),replace=FALSE) %>% 
                 as.matrix(.,nrow=nsim,ncol=ncol(X)))
  }
  return(sim)
}


#' Helper function to make the simulations for the survival curves using INLA
#' for a given formula and dataset
#' 
#' @param m The output of a 'survHE' object including the model
#' @param t A vector of times for which the survival curves are to be
#' computed
#' @param X The covariates profile for which to compute the survival curve
#' @param nsim The number of simulations to be computed
#' @param newdata A list with the "new data". This is a specific profile of
#' covariates in correspondence of which to compute the survival curves.
#' @param dist The string identifying the abbreviated name for the 
#' underlying model used 
#' @return \item{sim}{A list containing the generated simulations.}
#' @author Gianluca Baio
#' @seealso make.surv
#' @keywords INLA
#' @noRd 
make_sim_inla <- function(m,t,X,nsim,newdata,dist,time_max,...) {
  # Simulates from the distribution of the model parameters
  
  exArgs_inla <- list(...)
  
  # If 'nsim'=1 then use the point estimates for the coefficients
  if(nsim==1) {
    beta=matrix(m$summary.fixed[,"mean"],nrow=1)
    if(!is.null(m$summary.hyperpar)) {
      alpha=m$summary.hyperpar$mean
    } else {
      alpha=0
    }
  } else {
    # Selects the variables to sample from the posterior distributions (excludes the linear predictor)
    selection <- eval(parse(text=paste0('list(Predictor=-c(1:',nrow(exArgs_inla$data),'),',
                                        paste0("`",colnames(model.matrix(exArgs_inla$formula,exArgs_inla$data)),"`",'=1',collapse=','),')')))
    # Runs INLA to get nsim samples from the posterior of the parameters 
    jpost <- INLA::inla.posterior.sample(nsim,m,selection=selection)
    # Then computes the linear predictor
    beta <- matrix(unlist(lapply(jpost,function(x){x$latent})),ncol=length(jpost[[1]]$latent),nrow=nsim,byrow=T)
    # And the ancillary parameters (which may or may not exist for a given model)
    if(!is.null(jpost[[1]]$hyperpar)) {
      alpha <- matrix(unlist(lapply(jpost,function(x){x$hyperpar})),ncol=length(jpost[[1]]$hyperpar),nrow=nsim,byrow=T)
    } else {
      alpha <- 0
    }
  }
  ###############################################################################################################
  ## Needs to rescale the parameters because INLA is run on a time range [0-1] for computational stability
  if(dist%in%c("wei","exp","llo")) {
    if(grep("Intercept",colnames(X))>0){
      beta[,grep("Intercept",colnames(X))]=beta[,grep("Intercept",colnames(X))]-log(time_max)
    }
  }
  if(dist=="wph") {
    if(grep("Intercept",colnames(X))>0){
      beta[,grep("Intercept",colnames(X))]=(beta[,grep("Intercept",colnames(X))]+log(time_max))*(-alpha)
    }
  }
  if(dist=="lno") {
    if(grep("Intercept",colnames(X))>0){
      beta[,grep("Intercept",colnames(X))]=(beta[,grep("Intercept",colnames(X))]+log(time_max))
    }
  }
  if(dist=="gom") {
    alpha=alpha/time_max
    if(grep("Intercept",colnames(X))>0){
      beta[,grep("Intercept",colnames(X))]=beta[,grep("Intercept",colnames(X))]-log(time_max)
    }
  }
  ###############################################################################################################
  linpred <- beta %*% t(X)
  
  # Now uses the helper function to rescale the parameters and retrieve the correct simulations
  sim <- lapply(1:ncol(linpred),function(x) rescale.inla(linpred[,x],alpha,dist))
  
  return(sim)
}


#' Helper function to make the simulations for the survival curves using HMC
#' for a given formula and dataset
#' 
#' @param m The output of a 'survHE' object including the model
#' @param t A vector of times for which the survival curves are to be
#' computed
#' @param X The covariates profile for which to compute the survival curve
#' @param nsim The number of simulations to be computed
#' @param newdata A list with the "new data". This is a specific profile of
#' covariates in correspondence of which to compute the survival curves.
#' @param dist The string identifying the abbreviated name for the 
#' underlying model used 
#' @return \item{sim}{A list containing the generated simulations.}
#' @author Gianluca Baio
#' @seealso make.surv
#' @references Baio (2020). survHE
#' @keywords HMC
#' @noRd 
make_sim_hmc <- function(m,t,X,nsim,newdata,dist,summary_stat,...) {

  # Extracts the model object from the survHE output
  iter_stan <- m@stan_args[[1]]$iter
  
  beta=rstan::extract(m)$beta
  # Takes care of a couple of exceptions...
  # 1. If the model is intercept only then needs to remove the extra column added as trick in 'make_data_stan'
  if (ncol(X)==1) {
    # Need to make this a matrix otherwise it breaks in the case of the rps when it tries to remove the 
    # remaining column!
    beta=beta[,1] %>% as.matrix()
  }
  # 2. RPS has a weird construction and needs to remove the intercept if it's present
  if(dist=="rps" & any(grepl("Intercept",colnames(X)))) {
    X <- as.matrix(as_tibble(X) %>% select(-`(Intercept)`))
    beta=beta[,-ncol(beta)]
  }
  
  linpred <-  beta %*% t(X)
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
    sim <- lapply(sim,function(x) as.matrix(as_tibble(x) %>% summarise_all(summary_stat),nrow=1,ncol=ncol(x)))
  }
  if(nsim>1 & nsim<iter_stan) {
    # If the user selected a number of simulation < the one from rstan, then select a random sample 
    sim <- lapply(sim,function(x) as.matrix(as_tibble(x) %>% sample_n(nsim,replace=FALSE),nrow=nsim,ncol=ncol(x)))
  }
  return(sim)
}

#' Helper function to make the rescale the original simulations to the 
#' list sim to be used by 'make.surv' Exponential distribution/HMC
#' 
#' @param m The output of a 'survHE' object including the model
#' @param X The covariates profile for which to compute the survival curve
#' @param linpred The linear predictor obtained by multiplying the 
#' simulated values for the model coefficients and the covariates profile
#' @return \item{sim}{A list containing the generated simulations.}
#' @author Gianluca Baio
#' @seealso make.surv
#' @references Baio (2020). survHE
#' @keywords HMC Exponential
#' @noRd 
rescale_hmc_exp <- function(m,X,linpred){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # Exponential distribution
  sim <- as.matrix(exp(linpred)) # exp(rstan::extract(m)$beta %*% t(X)); 
  colnames(sim) <- "rate"
  return(sim)
}

#' Helper function to make the rescale the original simulations to the 
#' list sim to be used by 'make.surv' WeibullAFT distribution/HMC
#' 
#' @param m The output of a 'survHE' object including the model
#' @param X The covariates profile for which to compute the survival curve
#' @param linpred The linear predictor obtained by multiplying the 
#' simulated values for the model coefficients and the covariates profile
#' @return \item{sim}{A list containing the generated simulations.}
#' @author Gianluca Baio
#' @seealso make.surv
#' @references Baio (2020). survHE
#' @keywords HMC Weibull AFT
#' @noRd 
rescale_hmc_wei <- function(m,X,linpred){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # Weibull distribution
  shape <- as.numeric(rstan::extract(m)$alpha)
  scale <- exp(linpred) #exp(rstan::extract(m)$beta %*% t(X))
  sim <- cbind(shape,scale) 
  colnames(sim) <- c("shape","scale")
  return(sim)
}

#' Helper function to make the rescale the original simulations to the 
#' list sim to be used by 'make.surv' WeibullPH distribution/HMC
#' 
#' @param m The output of a 'survHE' object including the model
#' @param X The covariates profile for which to compute the survival curve
#' @param linpred The linear predictor obtained by multiplying the 
#' simulated values for the model coefficients and the covariates profile
#' @return \item{sim}{A list containing the generated simulations.}
#' @author Gianluca Baio
#' @seealso make.surv
#' @references Baio (2020). survHE
#' @keywords HMC Weibull PH
#' @noRd 
rescale_hmc_wph <- function(m,X,linpred){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # Weibull PH distribution
  shape <- as.numeric(rstan::extract(m)$alpha)
  scale <- exp(linpred) #exp(rstan::extract(m)$beta %*% t(X))
  sim <- cbind(shape,scale) 
  colnames(sim) <- c("shape","scale")
  return(sim)
}

#' Helper function to make the rescale the original simulations to the 
#' list sim to be used by 'make.surv' Gompertz distribution/HMC
#' 
#' @param m The output of a 'survHE' object including the model
#' @param X The covariates profile for which to compute the survival curve
#' @param linpred The linear predictor obtained by multiplying the 
#' simulated values for the model coefficients and the covariates profile
#' @return \item{sim}{A list containing the generated simulations.}
#' @author Gianluca Baio
#' @seealso make.surv
#' @references Baio (2020). survHE
#' @keywords HMC Gompertz
#' @noRd 
rescale_hmc_gom <- function(m,X,linpred){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # Gompertz distribution
  shape <- as.numeric(rstan::extract(m)$alpha)
  rate <- exp(linpred) #exp(rstan::extract(m)$beta %*% t(X))
  sim <- cbind(shape,rate) 
  colnames(sim) <- c("shape","rate")
  return(sim)
}

#' Helper function to make the rescale the original simulations to the 
#' list sim to be used by 'make.surv' Gamma distribution/HMC
#' 
#' @param m The output of a 'survHE' object including the model
#' @param X The covariates profile for which to compute the survival curve
#' @param linpred The linear predictor obtained by multiplying the 
#' simulated values for the model coefficients and the covariates profile
#' @return \item{sim}{A list containing the generated simulations.}
#' @author Gianluca Baio
#' @seealso make.surv
#' @references Baio (2020). survHE
#' @keywords HMC Gamma
#' @noRd 
rescale_hmc_gam <- function(m,X,linpred){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # Gamma distribution
  shape <- as.numeric(rstan::extract(m)$alpha)
  # NB: for the Gamma model, the linpred is already the log-rate!
  rate <- exp(linpred)
  sim <- cbind(shape,rate) 
  colnames(sim) <- c("shape","rate")
  return(sim)
}

#' Helper function to make the rescale the original simulations to the 
#' list sim to be used by 'make.surv' GenGamma distribution/HMC
#' 
#' @param m The output of a 'survHE' object including the model
#' @param X The covariates profile for which to compute the survival curve
#' @param linpred The linear predictor obtained by multiplying the 
#' simulated values for the model coefficients and the covariates profile
#' @return \item{sim}{A list containing the generated simulations.}
#' @author Gianluca Baio
#' @seealso make.surv
#' @references Baio (2020). survHE
#' @keywords HMC Generalised Gamma
#' @noRd 
rescale_hmc_gga <- function(m,X,linpred){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # Generalised Gamma distribution
  Q <- as.numeric(rstan::extract(m)$Q)
  sigma <- as.numeric(rstan::extract(m)$sigma)
  # NB: The linpred is already the mean mu on the natural scale so no need to exponentiate!
  mu <- linpred
  sim <- cbind(mu,sigma,Q) 
  colnames(sim) <- c("mu","sigma","Q")
  return(sim)
}

#' Helper function to make the rescale the original simulations to the 
#' list sim to be used by 'make.surv' Gen F distribution/HMC
#' 
#' @param m The output of a 'survHE' object including the model
#' @param X The covariates profile for which to compute the survival curve
#' @param linpred The linear predictor obtained by multiplying the 
#' simulated values for the model coefficients and the covariates profile
#' @return \item{sim}{A list containing the generated simulations.}
#' @author Gianluca Baio
#' @seealso make.surv
#' @references Baio (2020). survHE
#' @keywords HMC Generalised F
#' @noRd 
rescale_hmc_gef <- function(m,X,linpred){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # Generalised F distribution
  Q <- as.numeric(rstan::extract(m)$Q)
  P <- as.numeric(rstan::extract(m)$P)
  sigma <- as.numeric(rstan::extract(m)$sigma)
  # NB: The linpred is already the mean mu on the natural scale so no need to exponentiate!
  mu <- (linpred) 
  sim <- cbind(mu,sigma,Q,P) 
  colnames(sim) <- c("mu","sigma","Q","P")
  return(sim)
}

#' Helper function to make the rescale the original simulations to the 
#' list sim to be used by 'make.surv' logNormal distribution/HMC
#' 
#' @param m The output of a 'survHE' object including the model
#' @param X The covariates profile for which to compute the survival curve
#' @param linpred The linear predictor obtained by multiplying the 
#' simulated values for the model coefficients and the covariates profile
#' @return \item{sim}{A list containing the generated simulations.}
#' @author Gianluca Baio
#' @seealso make.surv
#' @references Baio (2020). survHE
#' @keywords HMC logNormal
#' @noRd 
rescale_hmc_lno <- function(m,X,linpred){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # logNormal distribution
  sdlog <- as.numeric(rstan::extract(m)$alpha)
  meanlog <- linpred 
  sim <- cbind(meanlog,sdlog) 
  colnames(sim) <- c("meanlog","sdlog")
  return(sim)
}

#' Helper function to make the rescale the original simulations to the 
#' list sim to be used by 'make.surv' logLogistic distribution/HMC
#' 
#' @param m The output of a 'survHE' object including the model
#' @param X The covariates profile for which to compute the survival curve
#' @param linpred The linear predictor obtained by multiplying the 
#' simulated values for the model coefficients and the covariates profile
#' @return \item{sim}{A list containing the generated simulations.}
#' @author Gianluca Baio
#' @seealso make.surv
#' @references Baio (2020). survHE
#' @keywords HMC logLogistic
#' @noRd 
rescale_hmc_llo <- function(m,X,linpred){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # logNormal distribution
  shape <- as.numeric(rstan::extract(m)$alpha)
  scale <- exp(linpred)
  sim <- cbind(shape,scale) 
  colnames(sim) <- c("shape","scale")
  return(sim)
}

#' Helper function to make the rescale the original simulations to the 
#' list sim to be used by 'make.surv' RPS distribution/HMC
#' 
#' @param m The output of a 'survHE' object including the model
#' @param X The covariates profile for which to compute the survival curve
#' @param linpred The linear predictor obtained by multiplying the 
#' simulated values for the model coefficients and the covariates profile
#' @return \item{sim}{A list containing the generated simulations.}
#' @author Gianluca Baio
#' @seealso make.surv
#' @references Baio (2020). survHE
#' @keywords HMC Royston-Parmar splines
#' @noRd 
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


#' Helper function to make the rescale the original simulations to the 
#' list sim to be used by 'make.surv' INLA
#' 
#' @param linpred The linear predictor obtained by multiplying the 
#' simulated values for the model coefficients and the covariates profile
#' @param alpha The hyperparameter from the INLA model
#' @param distr The abbreviated name of the underlying distribution
#' @return \item{sim}{A list containing the generated simulations.}
#' @author Gianluca Baio
#' @seealso make.surv
#' @references Baio (2020). survHE
#' @keywords HMC Exponential
#' @noRd 
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
    scale <- exp(-linpred)
    sim <- cbind(shape,scale) 
    colnames(sim) <- c("shape","scale")
  }
  if (distr=="lno") {
    mulog <- linpred
    sdlog <- 1/sqrt(alpha)
    sim <- cbind(mulog,sdlog) 
    colnames(sim) <- c("meanlog","sdlog")
  }
  if (distr=="gom") {
    shape <- alpha
    rate <- exp(linpred)
    sim = cbind(shape,rate)
    colnames(sim) <- c("shape","rate")
  }
  return(sim)
}


#' Helper function to make the compute the survival curves
#' 
#' @param sim A list containing the simulations for the relevant parameters
#' @param exArgs A list of extra arguments, as specified in 'make.surv'
#' @param nsim The number of simulations included
#' @param dist The abbreviated name of the underlying distribution
#' @param t The vector of times to be used in the x-axis
#' @return \item{mat}{A matrix of simulated values for the survival curves}
#' @author Gianluca Baio
#' @seealso make.surv
#' @references Baio (2020). survHE
#' @keywords Survival curves
#' @noRd 
compute_surv_curve <- function(sim,exArgs,nsim,dist,t,method,X) {  
  # Computes the survival curves
  args <- args_surv()
  distr <- manipulate_distributions(dist)$distr 
  
  if(dist=="rps") {
    # RPS-related options
    if(exists("scale",where=exArgs)) {scale=exArgs$scale} else {scale="hazard"}
    if(exists("timescale",where=exArgs)) {timescale=exArgs$timescale} else {timescale="log"}
    if(exists("log",where=exArgs)) {log=exArgs$log} else {log=FALSE}
    
    if(method=="hmc") {
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
          ),nrow=length(t),ncol=nsim,byrow=FALSE
        )
      })
    } 
    if(method=="mle") {
      # First needs to fiddle with the matrix of covariates profile
      if("(Intercept)"%in%colnames(X)){
        X=X %>% as_tibble() %>% select(-"(Intercept)")
      } else {
        X=X %>% as_tibble()
      }
###      if(exists("offset",where=exArgs)) {offset=exArgs$offset} else {offset=0}
      mat=lapply(sim,function(x) {
        gamma=as_tibble(x) %>% select(contains("gamma"))
        matrix(unlist(
          lapply(1:nsim,function(i) {
            1-do.call(psurvspline,args=list(
              q=t,
              gamma=as.numeric(gamma %>% slice(i)),
              beta=0,
              X=0,
              knots=exArgs$knots,
              scale=scale,
              timescale=timescale,
              offset=0,
              log=log
            ))
          })
          ),nrow=length(t),ncol=nsim,byrow=FALSE
        )
      })
    }
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
  mat <- mat %>% lapply(function(x) bind_cols(tibble(t=t),as_tibble(x)))

  return(mat)
}


#' Utility function to define the arguments needed to compute the cumulative 
#' distribution, with which to derive the survival function
#' 
#' @author Gianluca Baio
#' @seealso make.surv
#' @references Baio (2020). survHE
#' @keywords Survival curves
#' @noRd 
args_surv <- function() {
    list(
      exp='list(t,rate=x[,"rate"][i])',
      weibull='list(t,shape=x[,"shape"][i],scale=x[,"scale"][i])',
      weibullPH='list(t,shape=x[,"shape"][i],scale=x[,"scale"][i])',
      gamma='list(t,shape=x[,"shape"][i],rate=x[,"rate"][i])',
      gengamma='list(t,mu=x[,"mu"][i],sigma=x[,"sigma"][i],Q=x[,"Q"][i])',
      gompertz='list(t,shape=x[,"shape"][i],rate=x[,"rate"][i])',
      lnorm='list(t,meanlog=x[,"meanlog"][i],sdlog=x[,"sdlog"][i])',
      llogis='list(t,shape=x[,"shape"][i],scale=x[,"scale"][i])',
      genf='list(t,mu=x[,"mu"][i],sigma=x[,"sigma"][i],Q=x[,"Q"][i],P=x[,"P"][i])',
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
#' @param newdata a list **of lists**, specifiying the values of the covariates
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
#' @noRd 
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
    bind_cols(as_tibble(model.matrix(formula_temp,data)) %>% select(contains("Intercept"))) 
  if("(Intercept)"%in% names(covs)) {
    covs=covs %>% select(`(Intercept)`,everything())
  }
  
  ncovs <- covs %>% select(-contains("Intercept")) %>% with(ncol(.))
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
  if(all(unlist(lapply(fac.levels,length))<=2)) {colnames(X)=colnames(covs)}
  
  # The way the object X *must* be formatted depends on which way it's been generated.
  # The point is that it *always* has to be a matrix for other functions to process
  if(n.elements==0){
    # If all the covariates are factors, then get survival curves for *all* the combinations
    if(nfacts==ncovs & nfacts>0) {
      X=unique(model.matrix(formula,data))
      # If there's at least one factor with more than 2 levels, then do *not* rename the columns to the simpler version
      # which only has the name of the variable (rather than the combination, eg 'groupMedium' that R produces)
      if(all(unlist(lapply(fac.levels,length))<2)) {colnames(X)=colnames(covs)}
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
      # Turns the list 'newdata' into a data.frame & ensures the factors are indeed factors
      nd=do.call(rbind.data.frame,newdata) %>% as_tibble() %>% 
        mutate(across(names(is.fac),factor))
      # Augments the original data with the values supplied in 'newdata'. To make things work, 
      # needs to change the type of variables that are 'factors' in the analysis as well as remove
      # the NAs and replace with 0s
      aug_data=data %>% mutate(across(names(is.fac),factor)) %>% add_row(nd) %>% replace(is.na(.),0)
      # Creates a model frame with IDs
      mf=model.frame(formula,aug_data) %>% as_tibble() %>% select(-1) %>% 
        mutate(id=row_number()) %>% rename_if(is.factor,.funs=~gsub("as.factor[( )]","",.x)) %>% 
        rename_if(is.factor,.funs=~gsub("[( )]","",.x))
      # Now creates the 'model matrix' with the combination of all the factors
      mm=model.matrix(formula,aug_data) %>% as_tibble() %>% mutate(id=row_number())
      mf=suppressMessages(mf %>% right_join(nd))
      # And selects only the rows that match with the profile selected in 'newdata'
      X=as.matrix(mm %>% filter(id %in% mf$id) %>% select(-id) %>% unique,drop=FALSE)
    }
  }
  
  # Returns the output (the matrix with the covariates profile)
  return(X)
}


# Specialised function to make the PSA values and survival curves for the Poly-Weibull/HMC model
make_surv_pw=function(fit,mod,t,newdata,nsim,exArgs) {
  
  # Extracts the model object and the data from the survHE output
  m <- fit$models[[mod]]
  data <- fit$misc$data
  # Create a vector of times, if the user hasn't provided one, based on the observed data
  if(is.null(t)) {
    t <- sort(unique(fit$misc$km[[mod]]$time))
  }
  
  # Makes sure the distribution name(s) vector is in a useable format
  dist <- fit$misc$model_name[mod]
  
  # Now creates the profile of covariates for which to compute the survival curves
  X <- lapply(1:length(fit$misc$formula),function(i) make_profile_surv(fit$misc$formula[[i]],data,newdata[[i]]))
  #X <- lapply(fit$misc$formula,function(f) make_profile_surv(f,data,newdata))
  
  # Extracts the model object from the survHE output
  iter_stan <- m@stan_args[[1]]$iter
  # Coefficients (NB array with size nsim, M=number of mixture components, H=maximum number of covariates)
  beta=rstan::extract(m)$beta
  alpha=rstan::extract(m)$alpha
  
  # Checks that the coefficients and covariate matrix are conformable
  # This is because the formula has different length for the components so needs to
  # create "fake" covariates for the components with fewer covariates...
  truecovs=unlist(lapply(X,function(i) ncol(i)))
  ncomponents=length(truecovs)
  linpred=lapply(1:ncomponents,function(i) beta[,i,1:truecovs[i]]%*%t(X[[i]]))
  
  # Creates a list of lists containing the simulations for the rescaled model parameters
  # 'sim' has M elements (as many as the components of the mixture). Each of these has
  # as many elements as there are in the profile matrix for that component
  sim=lapply(1:ncomponents,function(k) {
    lapply(1:nrow(X[[k]]),function(i) {
      cbind(
        shape=as.numeric(rstan::extract(m)$alpha[,k]),
        rate=exp(linpred[[k]][,i])
      ) %>% as_tibble()
    })
  })
  
  # Selects the relevant number of simulations depending on the value of 'nsim'
  if(nsim>iter_stan) {
    stop("Please select a value for 'nsim' that is less than or equal to the value set in the call to 'fit.models'")
  }
  if(nsim==1) {
    # If the user requested only 1 simulation, then take the mean value
    sim=lapply(sim,function(x){
      lapply(1:length(x),function(i) {
        as.matrix(as_tibble(x[[i]]) %>% summarise_all(mean),nrow=1,ncol=ncol(x[[i]]))
      })
    })
  }
  if(nsim>1 & nsim<iter_stan) {
    # If the user selected a number of simulation < the one from rstan, then select a random sample 
    sim=lapply(sim,function(x) {
      lapply(1:length(x),function(i) {
        as.matrix(as_tibble(x[[i]]) %>% sample_n(nsim,replace=FALSE),nrow=nsim,ncol=ncol(x[[i]]))
      })
    })
  }
  
  # Now makes some work to compute the survival curves. First rescale the rates by times and shape
  vals=lapply(sim,function(k) {
    lapply(1:length(k),function(i) {
      matrix(
        unlist(
          lapply(1:nsim,function(j) {
            as.numeric(k[[i]][j,"rate"])*t^as.numeric(k[[i]][j,"shape"])
          })
        ),nrow=length(t),ncol=nsim
      )
    })
  })
  dims=unlist(lapply(vals,length))
  # If 'vals' has 1 element only in each component, then obviously need to use these
  if(all(dims==1)) {
    combs=matrix(1,nrow=1,ncol=length(vals))
  } else {
    # Creates a matrix of indicators for the combination of the covariates
    combs=expand.grid(lapply(X,function(i) 1:ncol(i)))
  }
  
  
  # Creates the list of matrices with the simulations from the survival curves
  mat=lapply(1:nrow(combs),function(i) {
    eval(parse(text=paste0("exp(- (",paste0("vals[[",1:ncol(combs),"]][[",combs[i,],"]]",collapse="+"),") )")))
  })
  
  for (i in 1:length(mat)){colnames(mat[[i]])=paste0("S_",1:nsim)}
  mat <- mat %>% lapply(function(x) bind_cols(as_tibble(t),as_tibble(x)) %>% rename(t=value))
  
  # Updates the model formulae to a single one with *all* the covariates
  if(is.null(newdata)) {
    comb.formula=update(fit$misc$formula[[1]],paste("~.+",paste(lapply(fit$misc$formula,function(i) terms(i)[[3]]),collapse="+")))
    X=make_profile_surv(comb.formula,data,newdata)
  } else {
    X=matrix(0,nrow=length(mat),ncol=1);
    colnames(X)="Custom profile"
  }

  list(sim=sim,mat=mat,X=X,t=t)
}