#' Fit Poly-Weibull model for survival analysis of mixture hazards
#' 
#' Runs the survival analysis using a Poly-Weibull model
#' 
#' On object in the class \code{survHE} containing the following elements
#' 
#' @param formula a list of formulae (one for each components of the mixture.
#' Can specify one single formula (in which case, the model is a simple Weibull
#' regression). For example, a valid call is using
#' \code{formula=list(Surv(time,event)~1,Surv(time,event)~arm)}
#' @param data A data frame containing the data to be used for the analysis.
#' This must contain data for the 'event' variable. In case there is no
#' censoring, then \code{event} is a column of 1s.
#' @param \dots Additional options (for INLA or HMC).
#' 
#' **HMC** specific options \code{chains} = number of chains to run in the HMC
#' (default = 2) \code{iter} = total number of iterations (default = 2000)
#' \code{warmup} = number of warmup iterations (default = iter/2) \code{thin} =
#' number of thinning (default = 1) \code{control} = a list specifying
#' Stan-related options, eg \code{control=list(adapt_delta=0.85)} (default =
#' NULL) \code{seed} = the random seed (to make things replicable) \code{pars}
#' = a vector of parameters (string, default = NA) \code{include} = a logical
#' indicator (if FALSE, then the pars are not saved; default = TRUE)
#' \code{priors} = a list (of lists) specifying the values for the parameters
#' of the prior distributions in the models \code{save.stan} = a logical
#' indicator (default = FALSE). If TRUE, then saves the data list for Stan and
#' the model file(s)
#' @return \item{models}{ A list containing the fitted models. These contain
#' the output from the original inference engine (\code{flexsurv}, \code{INLA}
#' or \code{rstan}). Can be processed using the methods specific to the
#' original packages, or via \code{survHE}-specific methods (such as
#' \code{plot}, \code{print}) or other specialised functions (eg to extrapolate
#' the survival curves, etc). } \item{model.fitting}{ A list containing the
#' output of the model-fit statistics (AIC, BIC, DIC). The AIC and BIC are
#' estimated for all methods, while the DIC is only estimated when using
#' Bayesian inference. } \item{method}{ A string indicating the method used to
#' fit the model, ie \code{'mle'}, \code{'inla'} or \code{'hmc'}.  }
#' \item{misc}{ A list containing the time needed to run the model(s) (in
#' seconds), the formula used, the results of the Kaplan-Meier analysis (which
#' is automatically performed using \code{npsurv}) and the original data frame.
#' }
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso Something will go here
#' @references Something will go here
#' @keywords Parametric survival models Bayesian inference via Hamiltonian
#' Monte Carlo Poly-Weibull model
#' @examples
#' 
#' ###
#' 
#' @export poly.weibull
poly.weibull <- function(formula=NULL,data,...) {
  # Fits the PolyWeibull model of Demiris et al (2015), SMMR 24(2), 287-301 to the data
  #
  # formula = a list of formulae to be used for each of the M components of the model
  # data = a data.frame containing the data
  # ... = optional arguments
  # chains = number of chains to run in the HMC (default = 2)
  # iter = total number of iterations (default = 2000)
  # warmup = number of warmup iterations (default = iter/2)
  # thin = number of thinning (default = 1)
  # control = a list specifying Stan-related options, eg control=list(adapt_delta=0.85) (default = NULL)
  # seed = the random seed (to make things replicable)
  # pars = a vector of parameters (string, default = NA)
  # include = a logical indicator (if FALSE, then the pars are not saved; default = TRUE)
  # cores = the number of CPU (cores) used to run the sampling procedure using rstan (default = 1)
  # priors = a list (of lists) specifying the values for the parameters of the prior distributions in the models
  # save.stan = a logical indicator (default = FALSE). If TRUE, then saves the data list for Stan and the model file(s)
  
  
  # This requires rstan, so if it's not installed asks for it!
  if(!isTRUE(requireNamespace("rstan",quietly=TRUE))) {
    stop("You need to install the R package 'rstan'. Please run in your R terminal:\n install.packages('rstan')")
  }
  
  # Optional arguments
  exArgs <- list(...)
  
  # Sets up defaults
  
  # Check whether the user has specified a list of formulae to be used in each compoenent of the model
  if (!is.null(formula)) {
    # a. The user has specified a formula, but has only given 1 element, then expand it
    #    so that formula becomes a list of formulae and sets up 2 components by default
    if (class(formula)=="formula") {
      formula <- list(formula)
    }
    M <- length(formula) 
    if(length(formula)<2){
      stop("Please speficy at least 2 components for the Poly-Weibull model by creating\n  a list of at least two formulae, eg: 'list(Surv(time,event)~1,Surv(time,event)~treatment)'")
    }
  }
  
  ## Stan options (the defaults are set in line with Stan's original)
  nlist <- NULL
  if(exists("chains",where=exArgs)) {chains <- exArgs$chains} else {chains <- 2} # DO WE WANT 4???
  if(exists("iter",where=exArgs)) {iter <- exArgs$iter} else {iter <- 2000}
  if(exists("warmup",where=exArgs)) {warmup <- exArgs$warmup} else {warmup <- floor(iter/2)}
  if(exists("thin",where=exArgs)) {thin <- exArgs$thin} else {thin <- 1}
  if(exists("control",where=exArgs)) {control <- exArgs$control} else {control <- NULL}
  if(exists("seed",where=exArgs)) {seed <- exArgs$seed} else {seed <- sample.int(.Machine$integer.max, 1)}
  if(exists("pars",where=exArgs)) {pars <- exArgs$pars} else {
    pars <- c("loglambda","lambda","lp__")
  }
  if(exists("include",where=exArgs)) {include <- exArgs$include} else {include <- FALSE}
  if(exists("cores",where=exArgs)) {cores <- exArgs$cores} else {cores <- 1}
  
  # Reconstructs the vars list based on the formula
  vars <- list()
  for (i in 1:M) {
    test <- attributes(terms(formula[[i]]))$term.labels
    ncovs <- length(test)
    time <- all.vars(formula[[i]],data)[1]
    event <- all.vars(formula[[i]],data)[2]
    if (ncovs>0) {
      Xraw <- model.frame(formula[[i]],data)
      w <- (which(sapply(Xraw,is.factor)==1))-1
      if (length(w)>=1) {
        factors <- gsub("as.factor[( )]","",test[w]) 
        factors <- gsub("[( )]","",factors)
        covs <- test[-w]
        if (length(covs)==0) {
          covs <- NULL
        }
      } else {
        factors <- NULL
        covs <- test
      }
    } else {
      covs <- factors <- NULL
    }
    # If there are covariates, creates a matrix and sets the dimension
    if(!is.null(covs)) {
      X <- data[,pmatch(covs,colnames(data))]
      K <- ifelse(!is.null(dim(X)[2]),dim(X)[2],1)
    }
    # If there are categorical covariates (in the vector 'factors'), makes sure they have the right form
    if(!is.null(factors)) {
      cols <- pmatch(factors,colnames(data))
      H <- length(cols)
      D <- numeric()
      for (i in 1:H) {
        data[,cols[i]] <- as.factor(data[,cols[i]])
        nlevs <- length(levels(data[,cols[i]]))
        D[i] <- nlevs
      } 
    } else {
      D <- 0
    }
    vars[[i]] <- list(time=time,event=event,factors=factors,covs=covs,nlevs=D)
  }
  
  # Loads the pre-compiled models
  dso <- stanmodels
  
  time2run <- numeric()
  # Selects the precompiled polyweibull model (CHECK IF THE ORDER IN availables.hmc CHANGES!!)
  dso <- dso[[6]] 
  
  data.stan <- list(t=data[,vars[[1]]$time], d=data[,vars[[1]]$event])
  data.stan$n <- length(data.stan$t)
  data.stan$M <- M
  X <- lapply(1:data.stan$M,function(i) model.matrix(formula[[i]],data))
  # max number of covariates in all the model formulae
  data.stan$H <- max(unlist(lapply(1:data.stan$M,function(i) ncol(X[[i]]))))
  # NB: Stan doesn't allow matrices with 1 column, so if there's only one covariate (eg intercept only), needs a little trick
  if (data.stan$H==1) {data.stan$H <- data.stan$H+1}
  X <- lapply(1:data.stan$M,function(i) {
    if(ncol(X[[i]])<data.stan$H) {
      X[[i]] <- cbind(X[[i]],matrix(0,nrow=nrow(X[[i]]),ncol=(data.stan$H-ncol(X[[i]]))))
    } else {
      X[[i]] <- X[[i]]
    }
  })
  data.stan$X <- array(NA,c(data.stan$M,data.stan$n,data.stan$H))
  for (m in 1:data.stan$M) {
    data.stan$X[m,,] <- X[[m]]
  }
  data.stan$mu_beta <- matrix(0,nrow=data.stan$H,ncol=data.stan$M)
  data.stan$sigma_beta <- matrix(10,data.stan$H,data.stan$M)
  
  # These are modified if the user gives values in the call to poly.weibull
  if(exists("priors",where=exArgs)) {
    priors <- exArgs$priors
    # Linear predictor coefficients
    if(!is.null(priors$mu_beta)) {
      data.stan$mu_beta <- priors$mu_beta
    }
    if(!is.null(priors$sigma_beta)) {
      data.stan$sigma_beta <- priors$sigma_beta
    }
  }
  
  # Now runs Stan to sample from the posterior distributions
  tic <- proc.time()
  out <- rstan::sampling(dso,data.stan,chains=chains,iter=iter,warmup=warmup,thin=thin,seed=seed,control=control,
                         pars=pars,include=include,cores=cores)
  toc <- proc.time()-tic
  time2run <- toc[3]
  list(out=out,data.stan=data.stan,time2run=time2run)
  
  if(exists("save.stan",where=exArgs)) {save.stan <- exArgs$save.stan} else {save.stan <- FALSE}
  time_survHE <- time2run
  time_stan <- sum(rstan::get_elapsed_time(out))
  time2run <- min(time_survHE,time_stan)
  names(time2run) <- "PolyWeibull"
  
  # Computes the log-likelihood 
  beta <- rstan::extract(out)$beta
  alpha <- rstan::extract(out)$alpha
  # NB: To make more robust estimate of AIC, BIC and DIC uses the median here (instead of the mean)
  #     This is likely to be a better approximation to the MLE when the underlying distributions are highly skewed!
  beta.hat <- apply(beta,c(2,3),median)
  alpha.hat <- apply(alpha,2,median)
  linpred <- lapply(1:data.stan$M,function(m) {
    beta[,m,]%*%t(data.stan$X[m,,])
  })
  linpred.hat <- lapply(1:data.stan$M,function(m) {
    beta.hat[m,]%*%t(data.stan$X[m,,])
  })
  
  h <- log_s <- array(NA,c(nrow(alpha),data.stan$n,data.stan$M))
  h_bar <- log_s_bar <- matrix(NA,data.stan$n,data.stan$M)
  for (m in 1:data.stan$M) {
    h_bar[,m] <- alpha.hat[m]*exp(linpred.hat[[m]])*data.stan$t^(alpha.hat[m]-1)
    log_s_bar[,m] <- exp(linpred.hat[[m]])*data.stan$t^(alpha.hat[m])
    for (i in 1:nrow(linpred[[m]])) {
      h[i,,m] <- alpha[i,m]*exp(linpred[[m]][i,])*data.stan$t^(alpha[i,m]-1)
      log_s[i,,m] <- exp(linpred[[m]][i,])*data.stan$t^(alpha[i,m])
    }
  }
  d_log_sum_h <- matrix(NA,nrow(alpha),data.stan$n)
  for (i in 1:nrow(alpha)) {
    d_log_sum_h[i,] <- data.stan$d * log(rowSums(h[i,,]))
  }
  loglik.bar <- sum(data.stan$d*log(rowSums(h_bar))-rowSums(log_s_bar))
  loglik <- rowSums(d_log_sum_h) - rowSums(log_s,2)
  D.theta <- -2*loglik
  D.bar <- -2*loglik.bar
  pD <- mean(D.theta) - D.bar
  pV <- 0.5*var(D.theta)  # Uses Gelman's definition of pD!
  dic <- mean(D.theta)+pD
  # Approximates AIC & BIC using the mean deviance and the number of nominal parameters
  npars <- data.stan$M + sum(unlist(lapply(1:data.stan$M,function(m) {sum(1-apply(data.stan$X[m,,],2,function(x) all(x==0)))})))
  aic <- D.bar+2*npars
  bic <- D.bar+npars*log(data.stan$n)
  
  # Now defines the outputs of the function
  model.fitting <- list(aic=aic,bic=bic,dic=dic)
  km <- list(time=data.stan$t)
  misc <- list(time2run=time2run,formula=formula,data=data,km=km)
  misc$vars <- vars; misc$data.stan=data.stan
  # If save.stan is set to TRUE, then saves the Stan model file(s) & data
  if(save.stan==TRUE) {
    model_code <- attr(mod$out@stanmodel,"model_code")
    con <- "PolyWeibull.stan"
    writeLines(model_code,con=con)
    txt <- paste0("Model code saved to the file: ",con,"\n")
    cat(txt)
  }
  mod <- list(out)
  names(mod) <- "PolyWeibull"
  # Finally prepares the output object
  res <- list(models=mod,model.fitting=model.fitting,method="hmc",misc=misc)
  # And sets its class attribute to "survHE"
  class(res) <- "survHE"
  return(res)
}
