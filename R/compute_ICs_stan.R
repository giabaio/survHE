#' Helper function to compute the information criteria statistics
#' when using hmc as the inferential engine. 'rstan' does not do
#' DIC automatically and AIC/BIC are also not standard for Bayesian
#' models, so can compute them post-hoc by manipulating the 
#' likelihood functions.
#' 
#' @param model The 'rstan' object with the model fit
#' @param distr3 The 'rstan' object with the model fit
#' @return \item{list}{A list containing the modified name of the 
#' distribution, the acronym (3-letters abbreviation), or the
#' labels (humane-readable name)}.
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso fit.models
#' @references Baio (2020). survHE
#' @keywords Parametric survival models Bayesian inference via Hamiltonian
#' Monte Carlo Bayesian inference via Integrated Nested Laplace Approximation
compute_ICs_stan <- function(model,distr3,data.stan) {
  # Computes the log-likelihood 
  beta <- rstan::extract(model)$beta
  # To safeguard against very asymmetric densities use the median (instead of the mean)
  beta.hat <- apply(beta,2,median)
  if (distr3 %in% c("exp", "wei", "wph", "gom", "lno", "llo","rps")) {
    linpred <- beta%*%t(data.stan$X)
    linpred.hat <- beta.hat%*%t(data.stan$X)
  } else {
    linpred <- list(
      lo=beta%*%t(data.stan$X_obs),
      lc=beta%*%t(data.stan$X_cens)
    )
    linpred.hat <- list(
      lo.bar=beta.hat%*%t(data.stan$X_obs),
      lc.bar=beta.hat%*%t(data.stan$X_cens)
    )
    #linpred <- linpred.hat <- NULL
  }
  
  # Now computes the densities using the helper functions
  out=do.call(what=paste0("lik_",distr3),args=list(distr3,linpred,linpred.hat,model,data.stan))
  # Extracts relevant variables from the list (See if there's a better way to do it!)
  logf=out$logf
  logf.hat=out$logf.hat
  npars=out$npars
  f=out$f
  f.bar=out$f.bar
  s=out$s
  s.bar=out$s.bar

  # Now computes the log-likelihood and then deviance and DIC, AIC, BIC
   if (distr3 %in% c("gam","gga","gef")) {
    loglik <- compute.loglik(f,s)
    D.theta <- -2*loglik 
    loglik.bar <- compute.loglik(f.bar,s.bar)
    D.bar <- -2*loglik.bar
    data.stan$n <- data.stan$n_obs+data.stan$n_cens
  } else {
    loglik <- apply(logf,1,sum)
    loglik.bar <- apply(logf.hat,1,sum)
  }
  D.theta <- -2*loglik
  D.bar <- -2*loglik.bar
  pD <- mean(D.theta) - D.bar
  pV <- 0.5*var(D.theta)
  dic <- mean(D.theta)+pD
  dic2 <- mean(D.theta)+pV
  # Approximates AIC & BIC using the mean deviance and the number of nominal parameters
  aic <- D.bar+2*npars                   #mean(D.theta)+2*pD
  bic <- D.bar+npars*log(data.stan$n)    #mean(D.theta)+pD*log(data.stan$n)
  
  list(aic=aic,bic=bic,dic=dic,dic2=dic2)
}

### Little function to compute the log-likelihood (for the obs vs censored cases)
compute.loglik <- function(f,s) {
  loglik <- (apply(log(f),1,sum) + apply(log(s),1,sum))
  return(loglik)
}

### These are utility functions to compute the log-density for the models fitted by hmc
lik_exp <- function(x,linpred,linpred.hat,model,data.stan) {
  logf <- matrix(
    unlist(lapply(1:nrow(linpred), function(i) {
      data.stan$d * log(hexp(data.stan$t, exp(linpred[i,]))) +
        log(1 - pexp(data.stan$t, exp(linpred[i,])))
    })),
    nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(
    data.stan$d * log(hexp(data.stan$t, exp(linpred.hat))) +
      log(1 - pexp(data.stan$t, exp(linpred.hat))),
    nrow = 1)
  # Number of parameters (for AIC): rate + covariates
  npars <- 1 + sum(1 - apply(data.stan$X, 2, function(x) all(x == 0)))
  list(logf=logf,logf.hat=logf.hat,npars=npars,f=NULL,f.bar=NULL,s=NULL,s.bar=NULL)
}

lik_wei <- function(x,linpred,linpred.hat,model,data.stan) {
  shape <- alpha <- as.numeric(rstan::extract(model)$alpha)
  shape.hat <- median(shape)
  logf <- matrix(
    unlist(lapply(1:nrow(linpred), function(i) {
      data.stan$d * log(hweibull(data.stan$t, shape[i], exp(linpred[i, ]))) + 
        log(1 - pweibull(data.stan$t, shape[i], exp(linpred[i, ])))
    })), 
    nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(
    data.stan$d * log(hweibull(data.stan$t, shape.hat, exp(linpred.hat))) + 
      log(1 - pweibull(data.stan$t, shape.hat, exp(linpred.hat))), 
    nrow = 1)
  # Number of parameters (for AIC): shape, scale + covariates
  npars <- 2 + sum(1 - apply(data.stan$X, 2, function(x) all(x == 0)))
  list(logf=logf,logf.hat=logf.hat,npars=npars,f=NULL,f.bar=NULL,s=NULL,s.bar=NULL)
}

lik_wph <- function(x,linpred,linpred.hat,model,data.stan) {
  shape <- alpha <- as.numeric(rstan::extract(model)$alpha)
  shape.hat = median(shape)
  logf <- matrix(
    unlist(lapply(1:nrow(linpred), function(i) {
      data.stan$d * log(hweibullPH(data.stan$t, shape[i], exp(linpred[i, ]))) +
        log(1 - pweibullPH(data.stan$t, shape[i], exp(linpred[i, ])))
    })), nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(
    data.stan$d * log(hweibullPH(data.stan$t, shape.hat, exp(linpred.hat))) +
      log(1 - pweibullPH(data.stan$t, shape.hat, exp(linpred.hat))), 
    nrow = 1)
  # Number of parameters (for AIC): shape, scale + covariates
  npars <- 2 + sum(1 - apply(data.stan$X, 2, function(x) all(x == 0)))
  list(logf=logf,logf.hat=logf.hat,npars=npars,f=NULL,f.bar=NULL,s=NULL,s.bar=NULL)
}

lik_gom <- function(x,linpred,linpred.hat,model,data.stan) {
  shape <- alpha <- as.numeric(rstan::extract(model)$alpha)
  shape.hat = median(shape)
  logf <- matrix(
    unlist(lapply(1:nrow(linpred), function(i) {
      data.stan$d * log(hgompertz(data.stan$t, shape = shape[i], rate = exp(linpred[i, ]))) +
        log(1 - pgompertz(shape = data.stan$t, shape[i], rate = exp(linpred[i, ])))
    })), 
    nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(
    data.stan$d * log(hgompertz(data.stan$t, shape.hat, exp(linpred.hat))) +
      log(1 - pgompertz(data.stan$t, shape.hat, exp(linpred.hat))),
    nrow = 1)
  # Number of parameters (for AIC): shape, rate + covariates
  npars <- 2 + sum(1 - apply(data.stan$X, 2, function(x) all(x == 0)))
  list(logf=logf,logf.hat=logf.hat,npars=npars,f=NULL,f.bar=NULL,s=NULL,s.bar=NULL)
}

lik_gam <- function(x,linpred,linpred.hat,model,data.stan) {
  shape <- alpha <- as.numeric(rstan::extract(model)$alpha)
  shape.bar <- median(shape)
  lo <- exp(linpred$lo)
  lc <- exp(linpred$lc)
  lo.bar <- exp(linpred.hat$lo.bar)
  lc.bar <- exp(linpred.hat$lc.bar)
  f = matrix(
    unlist(lapply(1:nrow(lo), function(i) 
      dgamma(data.stan$t, shape[i], lo[i, ]))),
    nrow = nrow(lo), byrow = T)
  f.bar = matrix(
    unlist(lapply(1:nrow(lo.bar), function(i)
      dgamma(data.stan$t, shape.bar, lo.bar[i, ]))), 
    nrow = 1, byrow = T)
  s = matrix(
    unlist(lapply(1:nrow(lc), function(i) 
      1 - pgamma(data.stan$d, shape[i], lc[i, ]))),
    nrow = nrow(lc), byrow = T)
  s.bar = matrix(
    unlist(lapply(1:nrow(lc.bar), function(i) 
      1 - pgamma(data.stan$d, shape.bar, lc.bar[i, ]))),
    nrow = 1, byrow = T)
  # Number of parameters (for AIC): shape, rate + covariates
  npars <- 2 + sum(1 - apply(data.stan$X_obs, 2, function(x) all(x == 0)))
  list(logf=NULL,logf.hat=NULL,npars=npars,f=f,f.bar=f.bar,s=s,s.bar=s)
}

lik_gga <- function(x,linpred,linpred.hat,model,data.stan) {
  q = as.numeric(rstan::extract(model)$Q)
  q.bar = median(q)
  scale = as.numeric(rstan::extract(model)$sigma)
  scale.bar = median(scale)
  lo <- exp(linpred$lo)
  lc <- exp(linpred$lc)
  lo.bar <- exp(linpred.hat$lo.bar)
  lc.bar <- exp(linpred.hat$lc.bar)
  f = matrix(
    unlist(lapply(1:nrow(lo), function(i)
      dgengamma(data.stan$t, lo[i, ], scale[i], q[i]))), 
    nrow = nrow(lo), byrow = T)
  f.bar = matrix(
    unlist(lapply(1:nrow(lo.bar), function(i)
      dgengamma(data.stan$t, lo.bar[i, ], scale.bar, q.bar))), 
    nrow = 1, byrow = T)
  s = matrix(
    unlist(lapply(1:nrow(lc), function(i)
      1 - pgengamma(data.stan$d, lc[i, ], scale[i], q[i]))), 
    nrow = nrow(lc), byrow = T)
  s.bar = matrix(
    unlist(lapply(1:nrow(lc.bar), function(i) 
      1 - pgengamma(data.stan$d, lc.bar[i, ], scale.bar, q.bar))), 
    nrow = 1, byrow = T)
  # Number of parameters (for AIC): mu, sigma, Q + covariates
  npars <- 3 + sum(1 - apply(data.stan$X_obs, 2, function(x) all(x == 0)))
  list(logf=NULL,logf.hat=NULL,npars=npars,f=f,f.bar=f.bar,s=s,s.bar=s)
}

lik_gef <- function(x,linpred,linpred.hat,model,data.stan) {
  Q = as.numeric(rstan::extract(model)$Q)
  Q.bar = median(Q)
  P = as.numeric(rstan::extract(model)$P)
  P.bar = median(P)
  sigma = as.numeric(rstan::extract(model)$sigma)
  sigma.bar = mean(sigma)
  lo <- exp(linpred$lo)
  lc <- exp(linpred$lc)
  lo.bar <- exp(linpred.hat$lo.bar)
  lc.bar <- exp(linpred.hat$lc.bar)
  f = matrix(
    unlist(lapply(1:nrow(lo), function(i) 
      dgenf(data.stan$t, lo[i, ], sigma[i], Q[i], P[i]))), 
    nrow = nrow(lo), byrow = T)
  f.bar = matrix(
    unlist(lapply(1:nrow(lo.bar), function(i) 
      dgenf(data.stan$t, lo.bar[i, ], sigma.bar, Q.bar, P.bar))), 
    nrow = 1, byrow = T)
  s = matrix(
    unlist(lapply(1:nrow(lc), function(i)
      1 - pgenf(data.stan$d, lc[i, ], sigma[i], Q[i], P[i]))),
    nrow = nrow(lc), byrow = T)
  s.bar = matrix(
    unlist(lapply(1:nrow(lc.bar), function(i) 
      1 - pgenf(data.stan$d, lc.bar[i, ], sigma.bar, Q.bar, P.bar))),
    nrow = 1, byrow = T)
  # Number of parameters (for AIC): mu, sigma, Q, P + covariates
  npars <- 4 + sum(
    1 - apply(data.stan$X_obs, 2, function(x) all(x == 0)))
  list(logf=NULL,logf.hat=NULL,npars=npars,f=f,f.bar=f.bar,s=s,s.bar=s)
}

lik_lno <- function(x,linpred,linpred.hat,model,data.stan) {
  sigma = as.numeric(rstan::extract(model)$alpha)
  sigma.hat = median(sigma)
  logf <- matrix(
    unlist(lapply(1:nrow(linpred), function(i) {
      data.stan$d * log(hlnorm(data.stan$t, (linpred[i, ]), sigma[i])) +
        log(1 - plnorm(data.stan$t, (linpred[i, ]), sigma[i]))
    })), nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(
    data.stan$d * log(hlnorm(data.stan$t, (linpred.hat), sigma.hat)) +
      log(1 - plnorm(data.stan$t, (linpred.hat), sigma.hat)),
    nrow = 1)
  # Number of parameters (for AIC): meanlog, sdlog + covariates
  npars <- 2 + sum(
    1 - apply(data.stan$X, 2, function(x) all(x == 0)))
  list(logf=logf,logf.hat=logf.hat,npars=npars,f=NULL,f.bar=NULL,s=NULL,s.bar=NULL)
}

lik_llo <- function(x,linpred,linpred.hat,model,data.stan) {
  sigma = as.numeric(rstan::extract(model)$alpha)
  sigma.hat = median(sigma)
  logf <- matrix(
    unlist(lapply(1:nrow(linpred), function(i) {
      data.stan$d * log(hllogis(data.stan$t, sigma[i], exp(linpred[i, ]))) +
        log(1 - pllogis(data.stan$t, sigma[i], exp(linpred[i, ])))
    })), nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(
    data.stan$d * log(hllogis(data.stan$t, sigma.hat, exp(linpred.hat))) +
      log(1 - pllogis(data.stan$t, sigma.hat, exp(linpred.hat))),
    nrow = 1)
  # Number of parameters (for AIC): shape, scale + covariates
  npars <- 2 + sum(
    1 - apply(data.stan$X, 2, function(x) all(x == 0)))
  list(logf=logf,logf.hat=logf.hat,npars=npars,f=NULL,f.bar=NULL,s=NULL,s.bar=NULL)
}

lik_rps <- function(x,linpred,linpred.hat,model,data.stan) {
  gamma <- rstan::extract(model)$gamma
  gamma.hat <- apply(gamma, 2, median)
  # Needs to reformat linpred.hat to a vector for consistency
  linpred.hat <- as.numeric(linpred.hat)
  logf <- data.stan$d * (-log(data.stan$t) + log(gamma %*% t(data.stan$DB)) + gamma %*% t(data.stan$B) + linpred) -
    exp(gamma %*% t(data.stan$B) + linpred)
  logf.hat <- t(
    data.stan$d * (-log(data.stan$t) + log(data.stan$DB %*% gamma.hat) + data.stan$B %*% gamma.hat + linpred.hat) -
      exp(data.stan$B %*% gamma.hat + linpred.hat))
  # Number of parameters (for AIC): gamma + covariates
  npars <- length(gamma.hat) + sum(
    apply(data.stan$X, 2, function(x) 1 - all(x == 0)))
  list(logf=logf,logf.hat=logf.hat,npars=npars,f=NULL,f.bar=NULL,s=NULL,s.bar=NULL)
}
