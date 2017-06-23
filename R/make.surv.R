make.surv <- function(fit,mod=1,t=NULL,newdata=NULL,nsim=1,...) {
  ## Creates the survival curves for the fitted model(s)
  # fit = the result of the call to the fit.models function, containing the model fitting (and other relevant information)
  # mod = the index of the model. Default value is 1, but the user can choose which model fit to visualise, 
  #     if the call to fit.models has a vector argument for distr (so many models are fitted & stored in the same object)
  # t = the time framework to be used for the estimation of the survival curve
  # newdata = a list (of lists), specifiying the values of the covariates at which the computation is performed. For example
  #           'list(list(arm=0),list(arm=1))' will create two survival curves, one obtained by setting the covariate 'arm'
  #           to the value 0 and the other by setting it to the value 1. In line with 'flexsurv' notation, the user needs
  #           to either specify the value for *all* the covariates or for none (in which case, 'newdata=NULL', which is the
  #           default). If some value is specified and at least one of the covariates is continuous, then a single survival
  #           curve will be computed in correspondence of the average values of all the covariates (including the factors, 
  #           which in this case are expanded into indicators). The order of the variables in the list *must* be the same
  #           as in the formula used for the model
  # nsim = the number of simulations from the distribution of the survival curves. Default at nsim=1, in which case
  #          uses the point estimate for the relevant distributional parameters and computes the resulting survival curve
  # ... = additional options
  
  exArgs <- list(...)
  if(is.null(t)) {
    t <- sort(unique(fit$misc$km$time))
  }
  
  m <- fit$models[[mod]]                # Extracts the model object from the survHE output
  if (fit$method=="hmc") {dist <- m@model_name} else {dist <- fit$models[[mod]]$dlist$name}   # Extracts the name of the distribution fitted
  
  n.elements <- ifelse(is.null(newdata),0,length(newdata))
  n.provided <- unlist(lapply(newdata,function(x) length(x)))
  
  # If no newdata are provided then see what to do!
  data <- fit$misc$data
  test <- attributes(terms(fit$misc$formula))$term.labels
  ncovs <- length(test)
  formula.temp <- as.formula(gsub("inla.surv","Surv",deparse(fit$misc$formula)))
  Xraw <- model.frame(formula.temp,data=data)
  is.fac <- sapply(Xraw, is.factor)[-1]
  w <- (which(sapply(Xraw,is.factor)==1))-1
  X <- matrix(colMeans(model.matrix(formula.temp,data=data)), nrow = 1)
  if(fit$method=="inla") {
    colnames(X) <- rownames(m$summary.fixed)
  } else {
    colnames(X) <- colnames(model.matrix(formula.temp,data=data)) #c("Intercept",test)
  }
  
  # newdata is not given (ie = NULL); this implies n.provided=NULL
  if (n.elements==0) {
    # If all the covariates are factors and mode_factor = False, then get survival curves for all the combinations
    if(all(is.fac) & length(is.fac)>0 ) {
      X <- unique(model.matrix(formula.temp,data=data))
      nam <- as.matrix(unique(X))
      for (i in 2:ncol(nam)) {
        nam[, i] <- paste(colnames(nam)[i],nam[, i], sep = "=")
      }
      rownames(X) <- apply(nam, 1, paste, collapse = ",")
    }
  }
  # newdata is a list with many values for the individual profiles
  if (n.elements>=1) {
    if (!all(n.provided==ncovs)) {
      stop("You need to provide data for *all* the covariates specified in the model, in the list 'newdata'")
    } else {
      X <- matrix(rep(X,n.elements),nrow=n.elements,byrow=T)
      if(fit$method=="inla") {
        colnames(X) <- rownames(m$summary.fixed)
      } else {
        colnames(X) <- colnames(model.matrix(formula.temp,data=data))
      }
      # Just like flexsurv, if you provide values for the covariates, you have to do so for *all*!
      names <- unique(unlist(lapply(newdata,function(x) names(x))))
      positions <- lapply(1:length(names),function(i) which(grepl(names[i],colnames(X))))
      temp <- matrix(unlist(newdata),nrow=length(newdata),byrow=T)
      colnames(temp) <- names
      # Could change the value in X with the value in temp[-w] for the continuous variables
      contin <- (1:length(names))[-w]
      # Do this only if there're some continuous covariates!
      if (length(contin)>0) {
        for (i in 1:length(contin)) {
          for (j in 1:n.elements) {
            X[j,positions[[contin[i]]]] <- temp[j,contin[i]]
          }
        }
      }
      # And then change the value in X with the factor expansion for the categorical variables, if there are any
      if (length(w)>0) {
        for (i in 1:length(w)) {
          for (j in 1:n.elements) {
            check <- eval(parse(text=paste0("levels(as.factor(data$",names[w[i]],"))")))
            if (class(check)=="character") {
              # check will be 0 or 1 depending on which level of the factor is selected in newdata
              check <- as.numeric(grepl(temp[j,w[i]],check))
            } else {
              check <- eval(parse(text=paste0("temp[j,w[i]]==as.numeric(levels(as.factor(data$",names[w[i]],")))")))
            }
            X[j,positions[[w[i]]]] <- check[-1]
          }
        }
      }
      nam <- as.matrix(unique(X))
      for (i in 2:ncol(nam)) {
        nam[, i] <- paste(colnames(nam)[i],nam[, i], sep = "=")
      }
      rownames(X) <- apply(nam, 1, paste, collapse = ",")
    }
  }
  
  # Now does the simulations
  if(fit$method=="mle") {
    dist <- ifelse(dist=="weibull.quiet","weibull",dist)
    S <- sim <-list()
    if(nsim==1) {
      S <- lapply(1:n.elements,function(i) summary(m,t=t,newdata=newdata[[i]]))
      ###S[[1]] <- summary(m,t=t)
      sim <- NULL
    } else {
      if (is.null(newdata)){
        sim <- lapply(1:dim(X)[1],function(i) flexsurv::normboot.flexsurvreg(m,B=nsim,X=matrix(X[i,-1],nrow=1)))
      } else {
        sim <- lapply(1:n.elements,function(i) flexsurv::normboot.flexsurvreg(m,B=nsim,newdata=newdata[[i]]))
      }
      txt1 <- paste("x[",1:dim(sim[[1]])[2],"]",sep="",collapse=",")
      if (dist=="survspline") {
        if(exists("scale",where=exArgs)) {scale=exArgs$scale} else {scale="hazard"}
        if(exists("timescale",where=exArgs)) {timescale=exArgs$timescale} else {timescale="log"}
        if(exists("offset",where=exArgs)) {offset=exArgs$offset} else {offset=0}
        if(exists("log",where=exArgs)) {log=exArgs$log} else {log=FALSE}
        tmp <- lapply(1:length(sim), function(i) {
          matrix(unlist(
            lapply(1:nsim,function(j) {
              1-flexsurv::psurvspline(t,gamma=sim[[i]][j,],knots=m$knots,scale=scale,timescale=timescale,offset=offset,log=log)
            })
          ),nrow=nsim,byrow=T)
        })
      } else {
        tmp <- lapply(1:length(sim), function(i) {
          eval(parse(text=paste0("t(apply(sim[[",i,"]],1,function(x) d",dist,"(t,",txt1,")/h",dist,"(t,",txt1,")))")))     
        }) 
      }
      S <- list(list())
      S <- lapply(1:nsim,function(i) {
        lapply(1:length(sim),function(j) {
          cbind(t,tmp[[j]][i,])
        })
      })
    }
  } 
  
  # If the original model(s) have been fitted using INLA, then use the (summaries of the) posterior distributions to compute the survival curves
  if(fit$method=="inla") {
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
	    # NB: As of Jan 11 2017, there's a little mistake in INLA and so need to minus the linpred HERE
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
        sdlog <- INLA::inla.contrib.sd(m)$hyper[1,1]
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
	      # NB: As of Jan 11 2017, there's a little mistake in INLA and so need to minus the linpred HERE
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
  
  if(fit$method=="hmc") {
    beta <- rstan::extract(m)$beta
    coefs <- beta
    if(fit$models[[mod]]@model_name%in%c("Gamma","GenGamma","GenF")) {
      covmat <- fit$misc$data.stan[[mod]]$X_obs
    } else {
      covmat <- fit$misc$data.stan[[mod]]$X
    }
    coefs=matrix(coefs[,apply(covmat,2,function(x) 1-all(x==0))==1],nrow=nrow(beta))
    # if (is.null(fit$misc$vars$factors) & is.null(fit$misc$vars$covs)) {
    #   coefs <- matrix(beta[,1],nrow=nrow(beta),byrow=T)
    # }
    if(ncol(coefs)>0) {
      if(dist!="RP") {
        colnames(coefs) <- colnames(model.matrix(fit$misc$formula,fit$misc$data))
      } else {
        colnames(coefs) <- colnames(model.matrix(fit$misc$formula,fit$misc$data))[-1]
      } 
    }
    basis <- function (knots, x) {
      nx <- length(x)
      if (!is.matrix(knots)) 
        knots <- matrix(rep(knots, nx), byrow = TRUE, ncol = length(knots))
      nk <- ncol(knots)
      b <- matrix(nrow = length(x), ncol = nk)
      if (nk > 0) {
        b[, 1] <- 1
        b[, 2] <- x
      }
      if (nk > 2) {
        lam <- (knots[, nk] - knots)/(knots[, nk] - knots[, 1])
        for (j in 1:(nk - 2)) {
          b[, j + 2] <- pmax(x - knots[, j + 1], 0)^3 - lam[,j + 1] * pmax(x - knots[, 1], 0)^3 - 
            (1 - lam[,j + 1]) * pmax(x - knots[, nk], 0)^3
        }
      }
      b
    }
    
    if (nsim==1) { # Computes the survival curve for the average value of all the parameters
      S <- list()
      sim <- NULL
      coefs <- apply(coefs,2,mean)
      if(dist=="Exponential") {
        linpred <- exp(coefs%*%t(X))
        s <- lapply(1:ncol(linpred),function(i) cbind(t,1-pexp(t,linpred[1,i])))
      }
      if (dist=="WeibullAF") {
        shape <- mean(as.numeric(rstan::extract(m)$alpha))
        linpred <- exp(coefs%*%t(X))
        s <- lapply(1:ncol(linpred),function(j) cbind(t,1-pweibull(t,shape,linpred[1,j])))
      }
      if (dist=="WeibullPH") {
        shape <- mean(as.numeric(rstan::extract(m)$alpha))
        linpred <- exp(coefs%*%t(X))
        s <- lapply(1:ncol(linpred),function(i) cbind(t,1-pweibullPH(t,shape,linpred[1,i])))
      }
      if (dist=="Gompertz") {
        shape <- mean(as.numeric(rstan::extract(m)$alpha))
        linpred <- exp(coefs%*%t(X))
        s <- lapply(1:ncol(linpred),function(i) cbind(t,1-pgompertz(t,shape,linpred[1,i])))
      }
      if (dist=="Gamma") {
        shape <- mean(as.numeric(rstan::extract(m)$alpha))
        linpred <- exp(coefs%*%t(X))
        s <- lapply(1:ncol(linpred),function(i) cbind(t,1-pgamma(t,shape,linpred[1,i])))
      }
      if (dist=="GenGamma") {
        q <- mean(as.numeric(rstan::extract(m)$Q))
        scale <- mean(as.numeric(rstan::extract(m)$sigma))
        linpred <- (coefs%*%t(X))
        s <- lapply(1:ncol(linpred),function(i) cbind(t,1-pgengamma(t,linpred[1,i],scale,q)))
      }
      if (dist=="GenF") {
        Q <- mean(as.numeric(rstan::extract(m)$Q))
        P <- mean(as.numeric(rstan::extract(m)$P))
        sigma <- mean(as.numeric(rstan::extract(m)$sigma))
        linpred <- (coefs%*%t(X))
        s <- lapply(1:ncol(linpred),function(i) cbind(t,1-pgenf(t,linpred[1,i],sigma,Q,P)))
      }
      if (dist=="logNormal") {
        sigma <- mean(as.numeric(rstan::extract(m)$alpha))
        linpred <- (coefs%*%t(X))
        s <- lapply(1:ncol(linpred),function(i) cbind(t,1-plnorm(t,linpred[1,i],sigma)))
      }
      if (dist=="logLogistic") {
        sigma <- mean(as.numeric(rstan::extract(m)$alpha))
        linpred <- exp(coefs%*%t(X))
        s <- lapply(1:ncol(linpred),function(i) cbind(t,1-pllogis(t,scale=linpred[1,i],shape=sigma)))
      }
      if (dist=="RP") {
        # Computes the knots wrt to the times selected for the analysis
        # If there's a time=0, then add a little constant
        t[t==0] <- min(0.00001,min(t[t>0]))
        B <- basis(fit$misc$data.stan[[mod]]$knots,log(t))
        gamma <- apply(rstan::extract(m)$gamma,2,mean)
        coefs <- c(0,coefs)
        if(nrow(X)==1) {
          s <- cbind(t,1-flexsurv::psurvspline(q=t,gamma=gamma,beta=coefs,X=X,knots=fit$misc$data.stan[[mod]]$knots))
        } else {
          s <- lapply(1:ncol(X),function(i) 
            cbind(t,1-flexsurv::psurvspline(q=t,gamma=gamma,beta=coefs,X=X[i,],knots=fit$misc$data.stan[[mod]]$knots)))
        }
      }
      S[[1]] <- s
    } else {
      if (nsim>length(beta)) {nrow=length(beta)}
      if(dist=="Exponential") {
        linpred <- exp(coefs%*%t(X))
        S <- lapply(1:nsim,function(i) {
          lapply(1:ncol(linpred),function(j) {
            cbind(t,1-pexp(t,linpred[i,j]))  
          })
        }) 
        sim <- coefs[1:nsim,]
      }
      if (dist=="WeibullAF") {
        shape <- as.numeric(rstan::extract(m)$alpha)
        linpred <- exp(coefs%*%t(X))
        S <- lapply(1:nsim,function(i) {
          lapply(1:ncol(linpred),function(j) {
            cbind(t,1-pweibull(t,shape[i],linpred[i,j]))  
          })
        }) 
        sim <- cbind(coefs,shape)[1:nsim,]
      }
      if (dist=="WeibullPH") {
        shape <- as.numeric(rstan::extract(m)$alpha)
        linpred <- exp(coefs%*%t(X))
        S <- lapply(1:nsim,function(i) {
          lapply(1:ncol(linpred),function(j) {
            cbind(t,1-pweibullPH(t,shape[i],linpred[i,j]))  
          })
        }) 
        sim <- cbind(coefs,shape)[1:nsim,]
      }
      if (dist=="Gompertz") {
        shape <- as.numeric(rstan::extract(m)$alpha)
        linpred <- exp(coefs%*%t(X))
        S <- lapply(1:nsim,function(i) {
          lapply(1:ncol(linpred),function(j) {
            cbind(t,1-pgompertz(t,shape[i],linpred[i,j]))  
          })
        }) 
        sim <- cbind(coefs,shape)[1:nsim,]
      }
      if (dist=="Gamma") {
        shape <- as.numeric(rstan::extract(m)$alpha)
        linpred <- exp(coefs%*%t(X))
        S <- lapply(1:nsim,function(i) {
          lapply(1:ncol(linpred),function(j) {
            cbind(t,1-pgamma(t,shape[i],linpred[i,j]))  
          })
        }) 
        sim <- cbind(coefs,shape)[1:nsim,]
      }
      if (dist=="GenGamma") {
        Q <- as.numeric(rstan::extract(m)$Q)
        shape <- as.numeric(rstan::extract(m)$sigma)
        linpred <- (coefs%*%t(X))
        S <- lapply(1:nsim,function(i) {
          lapply(1:ncol(linpred),function(j) {
            cbind(t,1-pgengamma(t,linpred[i,j],shape[i],Q[i]))  
          })
        }) 
        sim <- cbind(coefs,shape,Q)[1:nsim,]
      }
      if (dist=="GenF") {
        Q <- as.numeric(rstan::extract(m)$Q)
        P <- as.numeric(rstan::extract(m)$P)
        sigma <- as.numeric(rstan::extract(m)$sigma)
        linpred <- (coefs%*%t(X))
        S <- lapply(1:nsim,function(i) {
          lapply(1:ncol(linpred),function(j) {
            cbind(t,1-pgenf(t,linpred[i,j],sigma[i],Q[i],P[i]))  
          })
        }) 
        sim <- cbind(coefs,sigma,Q,P)[1:nsim,]
      }
      if (dist=="logNormal") {
        sigma <- as.numeric(rstan::extract(m)$alpha)
        linpred <- (coefs%*%t(X))
        S <- lapply(1:nsim,function(i) {
          lapply(1:ncol(linpred),function(j) {
            cbind(t,1-plnorm(t,linpred[i,j],sigma[i]))  
          })
        }) 
        sim <- cbind(coefs,sigma)[1:nsim,]
      }
      if (dist=="logLogistic") {
        sigma=as.numeric(rstan::extract(m)$alpha)
        linpred <- exp(coefs%*%t(X))
        S <- lapply(1:nsim,function(i) {
          lapply(1:ncol(linpred),function(j) {
            cbind(t,1-pllogis(t,linpred[i,j],sigma[i]))  
          })
        }) 
        sim <- cbind(coefs,sigma)[1:nsim,]
      }
      if (dist=="RP") {
        # Computes the knots wrt to the times selected for the analysis
        t[t==0] <- min(0.00001,min(t[t>0]))
        B <- basis(fit$misc$data.stan[[mod]]$knots,log(t))
        gamma <- rstan::extract(m)$gamma
        coefs <- cbind(rep(0,nrow(coefs)),coefs)
        if(nrow(X)==1) {
          S <- lapply(1:nsim,function(i) {
            lapply(1,function(j) {
              cbind(t,1-flexsurv::psurvspline(q=t,gamma=gamma[i,],beta=coefs[i,],X=X,knots=fit$misc$data.stan[[mod]]$knots))
            })
          })
        } else {
          S <- lapply(1:nsim,function(i) {
            lapply(1:ncol(X),function(j) {
              cbind(t,1-flexsurv::psurvspline(q=t,gamma=gamma[i,],beta=coefs[i,],X=X[j,],knots=fit$misc$data.stan[[mod]]$knots))
            })
          })
        }
        sim <- cbind(coefs[,-1],gamma)[1:nsim,]
      }
    }
  }
  
  n.elements <- length(S[[1]]) 
  if (fit$method=="mle") {
    mat <- lapply(1:n.elements,function(j) matrix(unlist(lapply(1:nsim,function(i) S[[i]][[j]][,2])),nrow=nsim,byrow=T))
  }
  mat <- lapply(1:n.elements,function(j) matrix(unlist(lapply(1:nsim,function(i) S[[i]][[j]][,2])),nrow=nsim,byrow=T))
  
  ### des.mat = X; rownames(des.mat) = names(fit$misc$km$strata)
  
  # Now defines the output of the function
  # S = a list --- for each simulated value of the parameters, a list with the survival curves associated with the configuration of the covariates
  # sim = simulated values for the main parameters (eg scale, shape, rate, mean, sd) for each configuration of the covariates
  # nsim = the number of simulations saved
  # mat =  a list --- for each configuration of covariates a matrix with nsims rows and ntimes columns with the survival curves (to be read row-wise)
  # des.mat = a design matrix with the combination of the covariates used (each represents an element in the lists S and mat)
  
  list(S=S,sim=sim,nsim=nsim,mat=mat,des.mat=X,times=t)
}
