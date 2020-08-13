#' Engine for Probabilistic Sensitivity Analysis on the survival curves
#' 
#' Creates the survival curves for the fitted model(s)
#' 
#' 
#' @param fit the result of the call to the \code{fit.models} function,
#' containing the model fitting (and other relevant information)
#' @param mod the index of the model. Default value is 1, but the user can
#' choose which model fit to visualise, if the call to fit.models has a vector
#' argument for distr (so many models are fitted & stored in the same object)
#' @param t the time vector to be used for the estimation of the survival curve
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
#' @param nsim The number of simulations from the distribution of the survival
#' curves. Default at \code{nsim=1}, in which case uses the point estimate for
#' the relevant distributional parameters and computes the resulting survival
#' curve
#' @param ...  Additional options
#' @author Gianluca Baio
#' @seealso Something will go here
#' @references Something will go here
#' @keywords Survival models Bootstrap Probabilistic sensitivity analysis
#' @examples
#' 
#' # Loads an example dataset from 'flexsurv'
#' data(bc)
#' 
#' # Fits the same model using the 3 inference methods
#' mle = fit.models(formula=Surv(recyrs,censrec)~group,data=bc,
#'     distr="exp",method="mle")
#' p.mle = make.surv(mle)
#' 
#' @export make.surv
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
  
  # Defines list with optional parameters
  exArgs <- list(...)
  
  # Extracts the model object and the data from the survHE output
  m <- fit$models[[mod]]
  data <- fit$misc$data
  # Create a vector of times, if the user hasn't provided one, based on the observed data
  if(is.null(t)) {
    t <- sort(unique(fit$misc$km$time))
  }
  # Makes sure the distribution name(s) vector is in a useable format
  dist <- fit$misc$model_name
  # This extra input is needed for HMC only
  if(method=="hmc"){name <- fit$models[[mod]]@model_name} else {name <- NULL}
  
  # Now creates the profile of covariates for which to compute the survival curves
  X <- make_profile_surv(formula,data,newdata)
  
  # Draws a sample of nsim simulations from the distribution of the model parameters
  sim <- do.call(paste0("make_sim_",fit$method),
                 args=list(m=m,t=t,X=X,nsim=nsim,newdata=newdata,dist=dist)
         )
  # Computes the survival curves
  S <- do.call(compute_surv_curve,
               args=list(sim=sim,exArgs=exArgs,nsim=nsim,dist=dist,t=t) 
        )
  
  # Formats the output
  list(S=S,sim=sim,nsim=nsim,mat=mat,des.mat=X,times=t)
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
#   
#   n.elements <- length(S[[1]]) 
#   if (fit$method=="mle") {
#     mat <- lapply(1:n.elements,function(j) matrix(unlist(lapply(1:nsim,function(i) S[[i]][[j]][,2])),nrow=nsim,byrow=T))
#   }
#   mat <- lapply(1:n.elements,function(j) matrix(unlist(lapply(1:nsim,function(i) S[[i]][[j]][,2])),nrow=nsim,byrow=T))
#   
#   ### des.mat = X; rownames(des.mat) = names(fit$misc$km$strata)
#   
#   # Now defines the output of the function
#   # S = a list --- for each simulated value of the parameters, a list with the survival curves associated with the configuration of the covariates
#   # sim = simulated values for the main parameters (eg scale, shape, rate, mean, sd) for each configuration of the covariates
#   # nsim = the number of simulations saved
#   # mat =  a list --- for each configuration of covariates a matrix with nsims rows and ntimes columns with the survival curves (to be read row-wise)
#   # des.mat = a design matrix with the combination of the covariates used (each represents an element in the lists S and mat)
#   
#   list(S=S,sim=sim,nsim=nsim,mat=mat,des.mat=X,times=t)
# }
