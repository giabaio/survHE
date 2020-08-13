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
  distr <- exArgs$name
  nsim <- exArgs$nsim
  
  # HERE DO THE LINPRED AND GET ALL THE PARAMS
  if (nsim==1 | nsim==nrow(beta)) {
    sim <- lapply(1:ncol(beta),function(x) do.call(paste0("rescale_hmc_",)))
  }
  if (nsim<nrow(beta)) {
    # HERE SAMPLE nsim FROM nrow(beta) AND DO THE SAME AS ABOVE
  }
  if (nsim>nrow(beta)) {
    stop()
  }
  return(sim)
}

rescale_hmc_exp <- function(m,X){
  # Rescales the original simulations to the list sim to be used by 'make.surv'
  # Exponential distribution
  sim <- exp(rstan::extract(m)$beta %*% t(X)); 
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

function(fit,mod,dist) {
  m <- fit$models[[mod]]
  beta <- rstan::extract(m)$beta
  if(dist%in%c("gam","gga","gef")) {
    covs <- fit$misc$data.stan[[mod]]$X_obs
  } else {
    covs <- fit$misc$data.stan[[mod]]$X
  }
  coefs=matrix(coefs[,apply(covmat,2,function(x) 1-all(x==0))==1],nrow=nrow(beta))
}

#   if(fit$method=="hmc") {
#     beta <- rstan::extract(m)$beta
#     coefs <- beta
#     if(fit$models[[mod]]@model_name%in%c("Gamma","GenGamma","GenF")) {
#       covmat <- fit$misc$data.stan[[mod]]$X_obs
#     } else {
#       covmat <- fit$misc$data.stan[[mod]]$X
#     }
#     coefs=matrix(coefs[,apply(covmat,2,function(x) 1-all(x==0))==1],nrow=nrow(beta))
#     # if (is.null(fit$misc$vars$factors) & is.null(fit$misc$vars$covs)) {
#     #   coefs <- matrix(beta[,1],nrow=nrow(beta),byrow=T)
#     # }
#     if(ncol(coefs)>0) {
#       if(dist!="RP") {
#         colnames(coefs) <- colnames(model.matrix(fit$misc$formula,fit$misc$data))
#       } else {
#         colnames(coefs) <- colnames(model.matrix(fit$misc$formula,fit$misc$data))[-1]
#       } 
#     }
#     basis <- function (knots, x) {
#       nx <- length(x)
#       if (!is.matrix(knots)) 
#         knots <- matrix(rep(knots, nx), byrow = TRUE, ncol = length(knots))
#       nk <- ncol(knots)
#       b <- matrix(nrow = length(x), ncol = nk)
#       if (nk > 0) {
#         b[, 1] <- 1
#         b[, 2] <- x
#       }
#       if (nk > 2) {
#         lam <- (knots[, nk] - knots)/(knots[, nk] - knots[, 1])
#         for (j in 1:(nk - 2)) {
#           b[, j + 2] <- pmax(x - knots[, j + 1], 0)^3 - lam[,j + 1] * pmax(x - knots[, 1], 0)^3 - 
#             (1 - lam[,j + 1]) * pmax(x - knots[, nk], 0)^3
#         }
#       }
#       b
#     }
#     
#     if (nsim==1) { # Computes the survival curve for the average value of all the parameters
#       S <- list()
#       sim <- NULL
#       coefs <- apply(coefs,2,mean)
#       if(dist=="Exponential") {
#         linpred <- exp(coefs%*%t(X))
#         s <- lapply(1:ncol(linpred),function(i) cbind(t,1-pexp(t,linpred[1,i])))
#       }
#       if (dist=="WeibullAF") {
#         shape <- mean(as.numeric(rstan::extract(m)$alpha))
#         linpred <- exp(coefs%*%t(X))
#         s <- lapply(1:ncol(linpred),function(j) cbind(t,1-pweibull(t,shape,linpred[1,j])))
#       }
#       if (dist=="WeibullPH") {
#         shape <- mean(as.numeric(rstan::extract(m)$alpha))
#         linpred <- exp(coefs%*%t(X))
#         s <- lapply(1:ncol(linpred),function(i) cbind(t,1-pweibullPH(t,shape,linpred[1,i])))
#       }
#       if (dist=="Gompertz") {
#         shape <- mean(as.numeric(rstan::extract(m)$alpha))
#         linpred <- exp(coefs%*%t(X))
#         s <- lapply(1:ncol(linpred),function(i) cbind(t,1-pgompertz(t,shape,linpred[1,i])))
#       }
#       if (dist=="Gamma") {
#         shape <- mean(as.numeric(rstan::extract(m)$alpha))
#         linpred <- exp(coefs%*%t(X))
#         s <- lapply(1:ncol(linpred),function(i) cbind(t,1-pgamma(t,shape,linpred[1,i])))
#       }
#       if (dist=="GenGamma") {
#         q <- mean(as.numeric(rstan::extract(m)$Q))
#         scale <- mean(as.numeric(rstan::extract(m)$sigma))
#         linpred <- (coefs%*%t(X))
#         s <- lapply(1:ncol(linpred),function(i) cbind(t,1-pgengamma(t,linpred[1,i],scale,q)))
#       }
#       if (dist=="GenF") {
#         Q <- mean(as.numeric(rstan::extract(m)$Q))
#         P <- mean(as.numeric(rstan::extract(m)$P))
#         sigma <- mean(as.numeric(rstan::extract(m)$sigma))
#         linpred <- (coefs%*%t(X))
#         s <- lapply(1:ncol(linpred),function(i) cbind(t,1-pgenf(t,linpred[1,i],sigma,Q,P)))
#       }
#       if (dist=="logNormal") {
#         sigma <- mean(as.numeric(rstan::extract(m)$alpha))
#         linpred <- (coefs%*%t(X))
#         s <- lapply(1:ncol(linpred),function(i) cbind(t,1-plnorm(t,linpred[1,i],sigma)))
#       }
#       if (dist=="logLogistic") {
#         sigma <- mean(as.numeric(rstan::extract(m)$alpha))
#         linpred <- exp(coefs%*%t(X))
#         s <- lapply(1:ncol(linpred),function(i) cbind(t,1-pllogis(t,scale=linpred[1,i],shape=sigma)))
#       }
#       if (dist=="RP") {
#         # Computes the knots wrt to the times selected for the analysis
#         # If there's a time=0, then add a little constant
#         t[t==0] <- min(0.00001,min(t[t>0]))
#         B <- basis(fit$misc$data.stan[[mod]]$knots,log(t))
#         gamma <- apply(rstan::extract(m)$gamma,2,mean)
#         coefs <- c(0,coefs)
#         if(nrow(X)==1) {
#           s <- cbind(t,1-flexsurv::psurvspline(q=t,gamma=gamma,beta=coefs,X=X,knots=fit$misc$data.stan[[mod]]$knots))
#         } else {
#           s <- lapply(1:ncol(X),function(i) 
#             cbind(t,1-flexsurv::psurvspline(q=t,gamma=gamma,beta=coefs,X=X[i,],knots=fit$misc$data.stan[[mod]]$knots)))
#         }
#       }
#       S[[1]] <- s
#     } else {
#       if (nsim>length(beta)) {nrow=length(beta)}
#       if(dist=="Exponential") {
#         linpred <- exp(coefs%*%t(X))
#         S <- lapply(1:nsim,function(i) {
#           lapply(1:ncol(linpred),function(j) {
#             cbind(t,1-pexp(t,linpred[i,j]))  
#           })
#         }) 
#         sim <- coefs[1:nsim,]
#       }
#       if (dist=="WeibullAF") {
#         shape <- as.numeric(rstan::extract(m)$alpha)
#         linpred <- exp(coefs%*%t(X))
#         S <- lapply(1:nsim,function(i) {
#           lapply(1:ncol(linpred),function(j) {
#             cbind(t,1-pweibull(t,shape[i],linpred[i,j]))  
#           })
#         }) 
#         sim <- cbind(coefs,shape)[1:nsim,]
#       }
#       if (dist=="WeibullPH") {
#         shape <- as.numeric(rstan::extract(m)$alpha)
#         linpred <- exp(coefs%*%t(X))
#         S <- lapply(1:nsim,function(i) {
#           lapply(1:ncol(linpred),function(j) {
#             cbind(t,1-pweibullPH(t,shape[i],linpred[i,j]))  
#           })
#         }) 
#         sim <- cbind(coefs,shape)[1:nsim,]
#       }
#       if (dist=="Gompertz") {
#         shape <- as.numeric(rstan::extract(m)$alpha)
#         linpred <- exp(coefs%*%t(X))
#         S <- lapply(1:nsim,function(i) {
#           lapply(1:ncol(linpred),function(j) {
#             cbind(t,1-pgompertz(t,shape[i],linpred[i,j]))  
#           })
#         }) 
#         sim <- cbind(coefs,shape)[1:nsim,]
#       }
#       if (dist=="Gamma") {
#         shape <- as.numeric(rstan::extract(m)$alpha)
#         linpred <- exp(coefs%*%t(X))
#         S <- lapply(1:nsim,function(i) {
#           lapply(1:ncol(linpred),function(j) {
#             cbind(t,1-pgamma(t,shape[i],linpred[i,j]))  
#           })
#         }) 
#         sim <- cbind(coefs,shape)[1:nsim,]
#       }
#       if (dist=="GenGamma") {
#         Q <- as.numeric(rstan::extract(m)$Q)
#         shape <- as.numeric(rstan::extract(m)$sigma)
#         linpred <- (coefs%*%t(X))
#         S <- lapply(1:nsim,function(i) {
#           lapply(1:ncol(linpred),function(j) {
#             cbind(t,1-pgengamma(t,linpred[i,j],shape[i],Q[i]))  
#           })
#         }) 
#         sim <- cbind(coefs,shape,Q)[1:nsim,]
#       }
#       if (dist=="GenF") {
#         Q <- as.numeric(rstan::extract(m)$Q)
#         P <- as.numeric(rstan::extract(m)$P)
#         sigma <- as.numeric(rstan::extract(m)$sigma)
#         linpred <- (coefs%*%t(X))
#         S <- lapply(1:nsim,function(i) {
#           lapply(1:ncol(linpred),function(j) {
#             cbind(t,1-pgenf(t,linpred[i,j],sigma[i],Q[i],P[i]))  
#           })
#         }) 
#         sim <- cbind(coefs,sigma,Q,P)[1:nsim,]
#       }
#       if (dist=="logNormal") {
#         sigma <- as.numeric(rstan::extract(m)$alpha)
#         linpred <- (coefs%*%t(X))
#         S <- lapply(1:nsim,function(i) {
#           lapply(1:ncol(linpred),function(j) {
#             cbind(t,1-plnorm(t,linpred[i,j],sigma[i]))  
#           })
#         }) 
#         sim <- cbind(coefs,sigma)[1:nsim,]
#       }
#       if (dist=="logLogistic") {
#         sigma=as.numeric(rstan::extract(m)$alpha)
#         linpred <- exp(coefs%*%t(X))
#         S <- lapply(1:nsim,function(i) {
#           lapply(1:ncol(linpred),function(j) {
#             cbind(t,1-pllogis(t,scale=linpred[i,j],shape=sigma[i]))  
#           })
#         }) 
#         sim <- cbind(coefs,sigma)[1:nsim,]
#       }
#       if (dist=="RP") {
#         # Computes the knots wrt to the times selected for the analysis
#         t[t==0] <- min(0.00001,min(t[t>0]))
#         B <- basis(fit$misc$data.stan[[mod]]$knots,log(t))
#         gamma <- rstan::extract(m)$gamma
#         coefs <- cbind(rep(0,nrow(coefs)),coefs)
#         if(nrow(X)==1) {
#           S <- lapply(1:nsim,function(i) {
#             lapply(1,function(j) {
#               cbind(t,1-flexsurv::psurvspline(q=t,gamma=gamma[i,],beta=coefs[i,],X=X,knots=fit$misc$data.stan[[mod]]$knots))
#             })
#           })
#         } else {
#           S <- lapply(1:nsim,function(i) {
#             lapply(1:ncol(X),function(j) {
#               cbind(t,1-flexsurv::psurvspline(q=t,gamma=gamma[i,],beta=coefs[i,],X=X[j,],knots=fit$misc$data.stan[[mod]]$knots))
#             })
#           })
#         }
#         sim <- cbind(coefs[,-1],gamma)[1:nsim,]
#       }
#     }
#   }

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
