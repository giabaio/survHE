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
make_sim_mle <- function(fit,t,X,nsim,newdata,dist) {
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

make_sim_inla <- function(fit,t,X,nsim,newdata,dist) {
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
# Re-defines inla.contrib.sd --- in case it's not in the main INLA package anymore
inla.contrib.sd = function(model, nsamples=1000) {
  ## contributed by Gianluca Baio <gianluca@stats.ucl.ac.uk>
  
  ## Computes the sd for the random effects in an INLA model
  ## 1. Defines the precision (generates a matrix with bins and
  ## density on the precision scale)
  
  ## 2. Simulates replications from the posterior distributions of
  ## the quantities of interest
  
  ## Names of the variables associated with structured effects
  rand.effs <- names(model$marginals.hyperpar)
  for (i in 1:length(rand.effs)) {
    cmd <- paste("prec.marg.",i,"<-model$marginals.hyperpar$'",rand.effs[i],"'",sep="")
    eval(parse(text=cmd)) # marginal distribution of the precision, tau
    ## Simulation from the posterior marginal distribution for sigma = 1/sqrt(tau)
    cmd <- paste("sigma.", i,
                 "<- inla.rmarginal(nsamples,inla.marginal.transform(function(x) 1/sqrt(x), prec.marg.",
                 i,"))",sep="")
    eval(parse(text=cmd))
  }
  
  ## Outputs of the function
  mat <- matrix(NA, nsamples, length(rand.effs))
  for (i in 1:length(rand.effs)) {
    cmd <- paste("mat[,i] <- sigma.",i,sep="")
    eval(parse(text=cmd)) 
  }
  names2 <- gsub("Precision","sd",rand.effs)
  colnames(mat) <- names2
  
  tab <- matrix(NA,length(rand.effs),4)
  for (i in 1:length(rand.effs)) {
    tab[i,] <- c(mean(mat[,i]),sd(mat[,i]),quantile(mat[,i],.025),quantile(mat[,i],.975))
  }
  rownames(tab) <- names2
  colnames(tab) <- c("mean","sd","2.5%","97.5%")
  
  return (list(samples=mat, hyper=tab))
}


  # A function to rescale the parameters of a given model and then computes the survival curve
  rescale.inla <- function(m,linpred) {
    if (m$dlist$name=="weibull") {
      shape <- m$summary.hyperpar[1,1]
      # NB: As of Jan 11 2017, there's a little mistake in INLA and so need to minus the linpred HERE
      scale <- exp(-linpred)
      S <- lapply(1:length(scale), function(x) cbind(t,dweibull(t,shape,scale[x])/hweibull(t,shape,scale[x]))) 
    }
    if (m$dlist$name=="weibullPH") {
      shape <- m$summary.hyperpar[1,1]
      scale <- exp(linpred)
      S <- lapply(1:length(scale), function(x) cbind(t,dweibullPH(t,shape,scale[x])/hweibullPH(t,shape,scale[x]))) 
    }
    if (m$dlist$name=="exponential") {
      rate <- exp(linpred)
      S <- lapply(1:length(rate), function(x) cbind(t,dexp(t,rate[x])/hexp(t,rate[x])))  
    }
    if (m$dlist$name=="loglogistic") {
      shape <- m$summary.hyperpar[1,1]
      scale <- exp(linpred)
      S <- lapply(1:length(scale), function(x) cbind(t,dllogis(t,shape,scale[x])/hllogis(t,shape,scale[x]))) 
    }
    if (m$dlist$name=="lognormal") {
      mulog <- linpred
      sdlog <- inla.contrib.sd(m)$hyper[1,1]
      S <- lapply(1:length(mulog), function(x) cbind(t,dlnorm(t,mulog[x],sdlog)/hlnorm(t,mulog[x],sdlog))) 
    }
    return(S)
  }
  
  # Now computes the survival curves for the relevant case
  if (nsim==1) {
    S <- list()
    linpred <- apply(m$summary.fixed[,1]*t(X),2,sum)
    S[[1]] <- rescale.inla(m,linpred)
    sim <- NULL
  } else {
    jpost <- suppressWarnings(INLA::inla.posterior.sample(n=nsim,m))
    pos <- pmatch(rownames(m$summary.fixed),rownames(jpost[[1]]$latent))
    sim1 <- matrix(unlist(lapply(jpost,function(x) x$latent[pos,])),ncol=length(pos),byrow=T)
    colnames(sim1) <- m$names.fixed
    if (m$nhyper>0) {
      sim2 <- matrix(unlist(lapply(jpost,function(x) x$hyperpar)),ncol=m$nhyper,byrow=T)
      sim <- cbind(sim2,sim1)
      colnames(sim) <- c(paste0("hyperpar",1:m$nhyper),m$names.fixed)
    } else {
      sim <- sim1
    }
    linpred <- matrix(unlist(lapply(1:nsim,function(i) apply(sim1[i,]*t(X),2,sum))),nrow=nsim,byrow=T)
    if(m$dlist$name=="weibull") {
      shape <- sim[,1]
      scale <- exp(-linpred)
      S <- lapply(1:nsim,function(i) {
        lapply(1:dim(scale)[2],function(j) {
          cbind(t,dweibull(t,shape[i],scale[i,j])/hweibull(t,shape[i],scale[i,j]))
        })
      })
    }
    if(m$dlist$name=="weibullPH") {
      shape <- sim[,1]
      scale <- exp(linpred)
      S <- lapply(1:nsim,function(i) {
        lapply(1:dim(scale)[2],function(j) {
          cbind(t,dweibullPH(t,shape[i],scale[i,j])/hweibullPH(t,shape[i],scale[i,j]))
        })
      })
    }
    if(m$dlist$name=="exponential") {
      rate <- exp(linpred)
      S <- lapply(1:nsim,function(i) {
        lapply(1:dim(rate)[2],function(j) {
          cbind(t,dexp(t,rate[i,j])/hexp(t,rate[i,j]))
        })
      })
    }
    if (m$dlist$name=="loglogistic") {
      shape <- sim[,1]
      scale <- exp(linpred)
      S <- lapply(1:nsim,function(i) {
        lapply(1:dim(scale)[2],function(j) {
          cbind(t,dllogis(t,scale=scale[i,j],shape=shape[i])/hllogis(t,scale=scale[i,j],shape=shape[i]))
        })
      })
    }
    if (m$dlist$name=="lognormal") {
      mulog <- linpred
      sdlog <- sim[,1]
      S <- lapply(1:nsim,function(i) {
        lapply(1:dim(mulog)[2],function(j) {
          cbind(t,dlnorm(t,mulog[i,j],sdlog[i])/hlnorm(t,mulog[i,j],sdlog[i]))
        })
      })
    }
  }
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
