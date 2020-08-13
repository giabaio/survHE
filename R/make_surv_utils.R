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
make_sim_mle <- function(m,t,X,nsim,newdata,...) {
  # Simulates from the distribution of the model parameters
  # If the user has not provided any new data, then *all* the possible covariates are averaged out and so there's
  # only one survival curve computed
  if(is.null(newdata)) {
    sim <- lapply(1:dim(X)[1],function(i) flexsurv::normboot.flexsurvreg(m,B=nsim,X=matrix(X[i,-1],nrow=1)))
  } 
  # If the user specifies a profile of covariates (in the list 'newdata') then use that. In this case
  if(!is.null(newdata)) {
    sim <- flexsurv::normboot.flexsurvreg(m,B=nsim,newdata=X)
  }
  return(sim)
}

make_sim_inla <- function(m,t,X,nsim,newdata,...) {
  # Simulates from the distribution of the model parameters
  # Selects the variables to sample from the posterior distributions (excludes the linear predictor)
  selection <- eval(parse(text=paste0('list(Predictor=-c(1:',nrow(data),'),',
                   paste0("`",colnames(model.matrix(formula,data)),"`",'=1',collapse=','),')')))
  # Runs INLA to get nsim samples from the posterior of the parameters 
  jpost <- INLA:::inla.posterior.sample(nsim,m,selection=selection)
  # Then computes the linear predictor
  linpred <- matrix(unlist(lapply(jpost,function(x){x$latent})),ncol=length(jpost[[1]]$latent),nrow=nsim,byrow=T)%*%t(X)
  # And the ancillary parameters (which may or may not exist for a given model)
  alpha <- matrix(unlist(lapply(jpost,function(x){x$hyperpar})),ncol=length(jpost[[1]]$hyperpar),nrow=nsim,byrow=T)
  
  # Now uses the helper function to rescale the parameters and retrieve the correct simulations
  sim <- lapply(1:ncol(linpred),function(x) rescale.inla(linpred[,x],alpha,m$dlist$name))
  return(sim)
}

make_sim_hmc <- function(m,t,X,nsim,newdata,...) {
  #### TO DO
  # need to filter the cases 
  # when nsim=1, use all simulations from stan, but then take mean value of parameters
  # when nsim=number of simulations from stan, use all
  # when nsim<number of simulations from stan, then use a sample of nsim from the original ones
  # when nsim>number of simulations from stan, return error
  #
  # Then need to figure out a clever way to derive the linpred and the other params
  
  # Extracts the model object from the survHE output
  exArgs <- list(...)
  dist <- exArgs$dist
  nsim <- exArgs$nsim
  
  # Stores the values returned by rstan into the list 'sim'
  sim <- lapply(1:nrow(X),function(x) do.call(paste0("rescale_hmc_",dist),
                                              args=list(m,X)))
  if(nsim>m@stan_args[[1]]$iter) {
    stop("Please select a value for 'nsim' that is less than or equal to the value set in the call to 'fit.models'")
  }
  if(nsim==1) {
    # If the user requested only 1 simulation, then take the mean value
    sim <- lapply(sim,function(x) apply(x,2,mean))
  }
  if(nsim<m@stan_args[[1]]$iter) {
    # If the user selected a number of simulation < the one from rstan, then select a random sample 
    sim <- lapply(sim,function(x) x[sample(1:nrow(sim[[1]]),nsim,replace=FALSE),])
  }
  return(sim)
}

rescale_hmc_exp <- function(m,X){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # Exponential distribution
  sim <- as.matrix(exp(rstan::extract(m)$beta %*% t(X))); 
  colnames(sim) <- "rate"
  return(sim)
}

rescale_hmc_wei <- function(m,X){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # Weibull distribution
  shape <- as.numeric(rstan::extract(m)$alpha)
  scale <- exp(rstan::extract(m)$beta %*% t(X))
  sim <- cbind(shape,scale) 
  colnames(sim) <- c("shape","scale")
  return(sim)
}

rescale_hmc_wph <- function(m,X){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # Weibull PH distribution
  shape <- as.numeric(rstan::extract(m)$alpha)
  scale <- exp(rstan::extract(m)$beta %*% t(X))
  sim <- cbind(shape,scale) 
  colnames(sim) <- c("shape","scale")
  return(sim)
}

rescale_hmc_gom <- function(m,X){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # Gompertz distribution
  shape <- as.numeric(rstan::extract(m)$alpha)
  rate <- exp(rstan::extract(m)$beta %*% t(X))
  sim <- cbind(shape,rate) 
  colnames(sim) <- c("shape","rate")
  return(sim)
}

rescale_hmc_gam <- function(m,X){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # Gamma distribution
  shape <- as.numeric(rstan::extract(m)$alpha)
  rate <- 1/exp(rstan::extract(m)$beta %*% t(X))
  sim <- cbind(shape,rate) 
  colnames(sim) <- c("shape","rate")
  return(sim)
}

rescale_hmc_gga <- function(m,X){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # Generalised Gamma distribution
  Q <- as.numeric(rstan::extract(m)$Q)
  sigma <- as.numeric(rstan(extract(m)$sigma))
  mu <- exp(rstan::extract(m)$beta %*% t(X))
  sim <- cbind(mu,sigma,Q) 
  colnames(sim) <- c("mu","sigma","Q")
  return(sim)
}

rescale_hmc_gef <- function(m,X){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # Generalised F distribution
  Q <- as.numeric(rstan::extract(m)$Q)
  P <- as.numeric(rstan::extract(m)$P)
  sigma <- as.numeric(rstan(extract(m)$sigma))
  mu <- exp(rstan::extract(m)$beta %*% t(X))
  sim <- cbind(mu,sigma,Q) 
  colnames(sim) <- c("mu","sigma","Q","P")
  return(sim)
}

rescale_hmc_lno <- function(m,X){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # logNormal distribution
  sdlog <- as.numeric(rstan::extract(m)$alpha)
  meanlog <- rstan::extract(m)$beta %*% t(X)
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
    sim <- exp(linpred)
    colnames(sim) <- "rate" 
  }
  if (distr=="llo") {
    shape <- alpha
    scale <- exp(linpred)
    sim <- cbind(shape,scale) 
    colnames(sim) <- c("shape","scale")
  }
  if (m$dlist$name=="lno") {
    mulog <- linpred
    sdlog <- 1/sqrt(alpha)
    sim <- cbind(mulog,sdlog) 
    colnames(sim) <- c("meanlog","sdlog")
  }
  return(sim)
}

compute_surv_curve <- function(sim,exArgs,nsim,dist,t) {  
  # Computes the survival curves
  args=args_surv()
  
  # RPS-related options
  if(exists("scale",where=exArgs)) {scale=exArgs$scale} else {scale="hazard"}
  if(exists("timescale",where=exArgs)) {timescale=exArgs$timescale} else {timescale="log"}
  if(exists("offset",where=exArgs)) {offset=exArgs$offset} else {offset=0}
  if(exists("log",where=exArgs)) {log=exArgs$log} else {log=FALSE}
  S <- lapply(sim,function(x) {
    cbind(t,matrix(
      unlist(
        lapply(1:nsim,function(i) {
          1-do.call(paste0("p",dist),args=eval(parse(text=args[[dist]])))
        })
      ),nrow=length(t),ncol=nsim,byrow=TRUE)
    )
  })
  return(S)
}

# Utility function to define the arguments needed to compute the cumulative distribution,
## with which to derive the survival function
args_surv <- function() {
    list(
      exp='list(t,rate=x[i,"rate"])',
      weibull='list(t,shape=x[i,"shape"],scale=x[i,"scale"])',
      weibullPH='list(t,shape=x[i,"shape"],scale=x[i,"scale"])',
      gamma='list(t,shape=x[i,"shape"],rate=x[i,"rate"])',
      gengamma='list(t,mu=x[i,"mu"],sigma=x[i,"sigma"],Q=x[i,"Q"])',
      gompertz='list(t,shape=x[i,"shape"],rate=x[i,"rate"])',
      lnorm='list(t,meanlog=x[i,"meanlog"],sdlog=x[i,"sdlog"])',
      llogis='list(t,shape=x[i,"shape"],scale=x[i,"scale"])',
      genf='list(t,mu=x[i,"mu"],sigma=x[i,"sigma"],Q=x[i,"Q"],P=x[i,"P])',
      rps='list(t,gamma=x[i,"gamma"],knots=m$knots,scale=scale,timescale=timescale,offset=offset,log=log)',
      survspline='list(t,gamma=x[i,"gamma"],knots=m$knots,scale=scale,timescale=timescale,offset=offset,log=log)'
  )
}
