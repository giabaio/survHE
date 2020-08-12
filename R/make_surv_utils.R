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
make_sim_mle <- function(fit,t,X,nsim,newdata) {
  # Simulates from the distribution of the model parameters
  # Extracts the model object from the survHE output
  m <- fit$models[[mod]]
  
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

make_sim_inla <- function(fit,t,X,nsim,newdata) {
  # Simulates from the distribution of the model parameters
  # Extracts the model object from the survHE output
  m <- fit$models[[mod]]
  
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

## Utility function to define the arguments needed to compute the cumulative distribution, 
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
