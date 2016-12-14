## SET OF UTILITY FUNCTIONS TO INCLUDE SURVIVAL ANALYSIS RESULTS INTO A HEALTH ECONOMIC MODEL
## Gianluca Baio + Will Browne, 14 Dec 2016
##
# cat(paste0("\n# ",ls(),"\n"))
# bugs2survHE
# digitise
# fit.models
# make.ipd
# make.surv
# model.fit.plot
# plot.survHE
# print.survHE
# psa.plot
# test.linear.assumptions
# write.surv
# make.transition.probs


fit.models <- function(formula=NULL,data,distr=NULL,method="mle",...) {
  ## Main function - runs the survival analysis with several useful options
  ## formula = a formula specifying the model to be used, in the form
  ##           Surv(time,event)~treatment[+covariates] for flexsurv
  ##           inla.surv(time,event)~treatment[+covariates] for INLA
  ## data = a data frame containing the data to be used
  ## distr = a (vector of) string(s) containing the name(s) of the model(s) to be fitted
  ## method = a string specifying the inferential method ("mle", "inla" or "hmc) 
  ## 
  ## ... additional options (mainly to do with INLA & MCMC)
  ##
  ## **INLA** specific options 
  ## dz = defines the step length for the grid search over the hyperparameters space (default = 0.1)
  ## diff.logdens = defines the difference in the log-density for the hyperparameters to stop integration (default = 5)
  ## control.fixed = defines the default for the priors, unless specified by the user. Default values are
  ##                 prior mean = 0 for *all* fixed effects
  ##                 prior var = 1000 for *all* fixed effects
  ##                 prior mean = 0 for the intercept
  ##                 prior prec -> 0 for the intercept 
  ## control.family = a list of options. If distr is a vector, then can be provided as a named
  ##                  list of options, for example something like this: 
  ##                  control.family=list(weibull=list(param=c(.1,.1)),lognormal=list(initial=2))
  ##                  the names of the elements of the list need to be the same as those given
  ##                  in the vector distr
  ##### NB: COULD SPECIFY DIFFERENT control.family OPTIONS FOR THE INLA METHOD WITH MANY DISTRIBUTIONS
  ##
  ## max_splines 
  ##
  ## HMC (via Stan)
  ## chains = number of chains to run in the HMC (default = 2)
  ## iter = total number of iterations (default = 2000)
  ## warmup = number of warmup iterations (default = iter/2)
  ## thin = number of thinning (default = 1)
  ## control = a list specifying Stan-related options, eg control=list(adapt_delta=0.85) (default = NULL)
  ## seed = the random seed (to make things replicable)
  ## pars = a vector of parameters (string, default = NA)
  ## include = a logical indicator (if FALSE, then the pars are not saved; default = TRUE)
  ## k = number of knots (when using distr="rps") --- default = 0
  ## cores = the number of CPU (cores) used to run the sampling procedure using rstan (default = 1)
  ## priors = a list (of lists) specifying the values for the parameters of the prior distributions in the models
  ## save.stan = a logical indicator (default = FALSE). If TRUE, then saves the data list for Stan and the model file(s)
  
  # Lists all the additional inputs
  exArgs <- list(...)
  # Avoids the 'no visible binding for global variable' error, when compiling
  model <- NULL
  
  # Needs to specify either the formula or the list of variables!
  if(is.null(formula)) {
    stop("You need to specify a model 'formula', e.g. 'formula=Surv(time,event)~treat'")
  }
  # ensures method is lower case
  method <- tolower(method)
  # ensures method is one of "mle","inla", "mcmc" or "splines"
  if(!method %in% c("hmc","inla","mle")) {
    stop("Methods available for use are 'mle', 'hmc' or 'inla'")
  }
  # INLA can only do a limited set of models (for now) so if user has selected
  # one that is not available, then falls back on MLE analysis
  availables.mle <- c("genf", "genf.orig", "gengamma", "gengamma.orig", "exp", 
                      "weibull", "weibullPH", "lnorm", "gamma", "gompertz", 
                      "llogis", "exponential", "lognormal","rps")
  availables.inla <- c("exponential","weibull","lognormal","loglogistic")
  availables.hmc <- c("exponential","gamma","genf","gengamma","gompertz","polyweibull","rps",
                      "weibull","weibullPH","loglogistic","lognormal")
  ### NB: The Poly-Weibull model does work, but it needs some special function

  # Standardises labels for model names
  labs <- distr
  labs[pmatch("weibull",labs)] <- "Weibull (AFT)"
  labs[pmatch("weibullPH",labs)] <- "Weibull (PH)"
  labs[pmatch("exp",labs)] <- "Exponential"
  labs[pmatch("exponential",labs)] <- "Exponential"
  labs[pmatch("gamma",labs)] <- "Gamma"
  labs[pmatch("lnorm",labs)] <- "log-Normal"
  labs[pmatch("lognormal",labs)] <- "log-Normal"
  labs[pmatch("llogis",labs)] <- "log-Logistic"
  labs[pmatch("loglogistic",labs)] <- "log-Logistic"
  labs[pmatch("gengamma",labs)] <- "Gen. Gamma"
  labs[pmatch("genf",labs)] <- "Gen. F"
  labs[pmatch("gompz",labs)] <- "Gompertz"
  labs[pmatch("dloglogis",labs)] <- "log-Logistic"
  labs[pmatch("rps",labs)] <- "Royston-Parmar"
  
  if(method=="inla") {
    # Checks that the distribution name(s) are consistent with INLA
    user.distr <- distr
    distr[pmatch("llogis",user.distr)] <- "loglogistic"
    distr[pmatch("exp",user.distr)] <- "exponential"
    distr[pmatch("lnorm",user.distr)] <- "lognormal"
    # But if there still are some that are just not available then falls back on MLE
    if (any(is.na(pmatch(distr,availables.inla)))) {
      method <- "mle"
      cat("NB: INLA can only fit Exponential, Weibull, log-Logistic or log-Normal parametric survival models. \nFalling back on MLE analysis")
    }
  }
 
  if (method=="hmc") {
	 # Fixes the way in which the distribution is written up so it's consistent with flexsurv
	 user.distr <- distr
   distr[pmatch("exp",user.distr)] <- "exponential"
   distr[pmatch("lnorm",user.distr)] <- "lognormal"
   distr[pmatch("llogis",user.distr)] <- "loglogistic"
  }

  # Reconstructs the vars list based on the formula
  test <- attributes(terms(formula))$term.labels
  ncovs <- length(test)
  formula.temp <- as.formula(gsub("inla.surv","Surv",deparse(formula)))
  time <- all.vars(formula.temp,data)[1]
  event <- all.vars(formula.temp,data)[2]
  if (ncovs>0) {
    Xraw <- model.frame(formula.temp,data)
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
  vars <- list(time=time,event=event,factors=factors,covs=covs,nlevs=D)
  
  # Need to create a formula for the KM that is in the right format (flexsurv-like)
  chk <- is.na(pmatch("Surv",attributes(terms(formula))$variables[[2]][1])) # If TRUE, then it's in inla.surv() terms
  # If FALSE, then formula is in Surv() terms
  
  # If it's in inla.surv() needs to create a Surv() formula for the KM + check that the method isn't mle
  if(chk) {
    tmp <- deparse(formula)
    km.formula <- as.formula(gsub("inla.surv","Surv",tmp))
    # If the method was originally INLA but has been modified or if the method is mle but the formula is in the wrong
    # format, then need to change the formula to be consistent with the Surv() notation
    if (method=="mle") {
      formula <- km.formula
    }
  } else {
    km.formula <- formula
    # If the method is inla but the formula is in Surv() terms, then change the formula
    if (method=="inla") {
      tmp <- deparse(formula)
      formula <- as.formula(gsub("Surv","INLA::inla.surv",tmp))
    }
  }
  # Computes the Kaplan Meier curve using the package "rms"
  ObjSurvfit=rms::npsurv(        # Uses the function "npsurv" from the package "rms"
    formula=km.formula,          # to fit the model specified in the "formula" object
    data=data                    # to the dataset named "data"
  )
  
  # If method = MLE, then fits the model(s) using flexsurvreg
  if (method=="mle") {
    if(!isTRUE(requireNamespace("flexsurv",quietly=TRUE))) {
      stop("You need to install the R package 'flexsurv'. Please run in your R terminal:\n install.packages('flexsurv')")
    }
    # Checks that the distribution name(s) are consistent with flexsurv
    # The only problem here is if the user has specified a log-Logistic in INLA terminology
    user.distr <- distr
    distr[pmatch("loglogistic",user.distr)] <- "llogis"
    # Then run the model(s)
    runMLE <- function(distr) {
      tic <- proc.time()
      if(distr=="rps") {
          # If user selects RPS model, then could also provide some optional arguments - uses flexsurv defaults
          if(exists("bhazard",where=exArgs)) {bhazard=exArgs$bhazard} else {bhazard=NULL}
          if(exists("weights",where=exArgs)) {weights=exArgs$weights} else {weights=NULL}
          if(exists("subset",where=exArgs)) {subset=exArgs$subset} else {subset=NULL}
          if(exists("knots",where=exArgs)) {knots=exArgs$knots} else {knots=NULL}
          if(exists("k",where=exArgs)) {k=exArgs$k} else {k=0}
          if(exists("bknots",where=exArgs)) {bknots=exArgs$bknots} else {bknots=NULL}
          if(exists("scale",where=exArgs)) {scale=exArgs$scale} else {scale="hazard"}
          if(exists("timescale",where=exArgs)) {timescale=exArgs$scale} else {timescale="log"}
          model = flexsurv::flexsurvspline(formula=formula,data=data,k=k,knots=knots,bknots=bknots,scale=scale,timescale=timescale)
      } else {
          model <- flexsurv::flexsurvreg(formula=formula,data=data,dist=distr)
      }
      toc <- proc.time()-tic
      time2run=toc[3]
      list(model=model,time2run=time2run)
    }
    output <- lapply(distr,function(x) runMLE(x))
    mod <- lapply(output, function(x) x$model)
    time2run <- unlist(lapply(output, function(x) x$time2run)); names(time2run) <- labs
    aic <- unlist(lapply(mod,function(x) x$AIC))
    bic <- unlist(lapply(mod,function(x) -2*x$loglik+x$npars*log(x$N)))
    dic <- rep(NA,length(distr))
  }
  
  # If method = INLA, then fits model(s) using inla
  if (method=="inla") {
    # If INLA is not installed, then asks for it
    if (!isTRUE(requireNamespace("INLA", quietly = TRUE))) {
      stop("You need to install the packages 'INLA'. Please run in your R terminal:\n install.packages('INLA', repos='https://www.math.ntnu.no/inla/R/stable')")
    }
    # If INLA is installed, then attach the Namespace (so that all the relevant functions are available)
    if (isTRUE(requireNamespace("INLA", quietly = TRUE))) {
      if (!is.element("INLA", (.packages()))) {
        attachNamespace("INLA")
      }
      # Set up optional parameters to default values if the user hasn't done it themselves
      # 1. defines the step length for the grid search over the hyperparameters space
      if(exists("dz",where=exArgs)) {dz <- exArgs$dz} else {dz <- 0.1}
      # 2. defines the difference in the log-density for the hyperparameters to stop integration
      if(exists("diff.logdens",where=exArgs)) {diff.logdens <- exArgs$diff.logdens} else {diff.logdens <- 5}
      # 3. defines the default for the priors, unless specified by the user
      if(exists("control.fixed",where=exArgs)) {
        control.fixed <- exArgs$control.fixed
      } else {
        control.fixed <- INLA::inla.set.control.fixed.default()
        # prior mean = 0 for *all* fixed effects
        # prior var = 1000 for *all* fixed effects
        # prior mean = 0 for the intercept
        # prior prec -> 0 for the intercept 
        ## This makes the priors consistent with the defaults in HMC
        ## The available models all have sd=5 in HMC, which translates to a precision of 1/25!
        control.fixed$prec=control.fixed$prec.intercept=1/(5^2); 
      }
      if(exists("control.family",where=exArgs)) {
        control.family <- replicate(length(distr),list(INLA::inla.set.control.family.default()))
        names(control.family) <- distr
        cf <- exArgs$control.family
        news <- pmatch(names(cf),names(control.family))
        for (i in 1:length(news)) {
          control.family[[news[i]]] <- cf[[i]]
          names(control.family[[news[i]]]) <- names(cf[[i]])
        }
      } else {
        # If not specified, uses the default values in INLA (depending on the model selected)
        control.family <- replicate(length(distr),list(INLA::inla.set.control.family.default()))
        # And sets some more sensible/general values for some of the controls
        #       if(!is.na(pmatch("weibull",distr))) {
        #         control.family[[pmatch("weibull",distr)]]$param <- c(.1,.1)
        #       }
        #        for (i in 1:length(control.family)) {
        #          control.family[[i]]$initial=.1
        #        }
      }
      
      # 4. Finally runs INLA
      mod <- lapply(1:length(distr), function(x) {
        INLA::inla(formula,family=distr[x],data=data,control.compute=list(config=TRUE,dic=TRUE),
                   control.inla=list(int.strategy="grid",dz=dz,diff.logdens=diff.logdens),
                   control.fixed=control.fixed,control.family=control.family[[x]]
        )
      })
      # Now re-writes the formula in general terms (without linking to INLA::inla.surv)
      formula <- as.formula(gsub("INLA::inla.surv","Surv",deparse(formula)))
      time2run <- unlist(lapply(mod,function(x) x$cpu.used["Total"])); names(time2run) <- labs
      # NB Internally, DIC = model$dic$mean.deviance+model$dic$p.eff
      dic <- unlist(lapply(mod,function(x) x$dic$dic))
      # NB But to estimate the AIC & BIC is probably best to use model$dic$deviance.mean!
      aic <- unlist(lapply(1:length(mod), function(i) 2*mod[[i]]$dic$p.eff+mod[[i]]$dic$deviance.mean)); names(aic) <- NULL
      bic <- unlist(lapply(1:length(mod),function(i) 
        mod[[i]]$dic$deviance.mean+mod[[i]]$dic$p.eff*log(mod[[i]]$size.linear.predictor$n))); names(bic) <- NULL
      for (i in 1:length(distr)) {mod[[i]]$dlist$name <- distr[i]}
    }
  }
    
  if(method=="hmc") {
    if(!isTRUE(requireNamespace("rstan",quietly=TRUE))) {
      stop("You need to install the R package 'rstan'. Please run in your R terminal:\n install.packages('rstan')")
    }

    # Now runs the model
    ## Stan options (the defaults are set in line with Stan's original)
    nlist <- NULL
    if(exists("chains",where=exArgs)) {chains <- exArgs$chains} else {chains <- 2} # DO WE WANT 4???
    if(exists("iter",where=exArgs)) {iter <- exArgs$iter} else {iter <- 2000}
    if(exists("warmup",where=exArgs)) {warmup <- exArgs$warmup} else {warmup <- floor(iter/2)}
    if(exists("thin",where=exArgs)) {thin <- exArgs$thin} else {thin <- 1}
    if(exists("control",where=exArgs)) {
        check <- unlist(lapply(1:length(control),function(i) class(control[[i]])))
        nlists <- length(check)
        if (nlists==length(distr)) {
            control <- ifelse(nlist==1,list(exArgs$control),exArgs$control)
        } 
        if (nlists<length(check)) {
            control <- list(exArgs$control,replicate((nlists-1),list(NULL)))
        }
        if (nlists>length(check)) {
            msg <- paste0("Please provide at most ",length(distr)," lists in the 'control' argument")
            stop(msg)
        }
        control <- ifelse(class(exArgs$control[[1]])!="list",list(exArgs$control),exArgs$control)
    } else {
        control <- replicate(length(distr),list(NULL))
    }
    if(exists("seed",where=exArgs)) {seed <- exArgs$seed} else {seed <- sample.int(.Machine$integer.max, 1)}
    if(exists("pars",where=exArgs)) {pars <- exArgs$pars} else {
      pars <- c("lambda_cens","lambda_obs","cens","d","lp__","loglambda_cens","loglambda_obs","mu","logP","linpred")
    }
    if(exists("include",where=exArgs)) {include <- exArgs$include} else {include <- FALSE}
    if(exists("k",where=exArgs)) {k=exArgs$k} else {k=0}
    if(exists("cores",where=exArgs)) {cores=exArgs$cores} else {cores=1}

    non.on.log.scale <- c("genf","gengamma","lognormal")
    
    # Loads the pre-compiled models
    dso <- stanmodels
    
    ###############################################################################################
    ### THIS IS JUST A TEMPORARY LINE (UNTIL THE CORRECT MODELS ARE PRE-COMPILED!)
    ####dso = readRDS("~/Dropbox/UCL/Mapi/Projects/Survival/Stan_code/DSOs.rds")
    ###############################################################################################
    touse = time2run = numeric()
    for (i in 1:length(distr)) {
      touse[i] <- match(distr[i],availables.hmc)
    }
    dso <- lapply(touse,function(x) dso[[x]])

    mod <- lapply(1:length(distr), function(x) {
      # First makes the data list
      if (distr[x] %in% c("gamma","gengamma","genf")) {
        # If model is Gamma, GenGamma or GenF, then use the "obs vs" censored format
        data.stan = list(t=data[data[,vars$event]==1,vars$time],d=data[data[,vars$event]==0,vars$time])
        data.stan$n_obs <- length(data.stan$t); data.stan$n_cens <- length(data.stan$d)
        data.stan$X_obs = matrix(model.matrix(formula,data)[data[,vars$event]==1,],nrow=data.stan$n_obs,byrow=F)
        data.stan$X_cens = matrix(model.matrix(formula,data)[data[,vars$event]==0,],nrow=data.stan$n_cens,byrow=F)
        data.stan$H=ncol(data.stan$X_obs)
        # NB: Stan doesn't allow vectors of size 1, so if there's only one covariate (eg intercept only), needs a little trick
        if (data.stan$H==1) {
          data.stan$X_obs = cbind(data.stan$X_obs,rep(0,data.stan$n_obs))
          data.stan$X_cens = cbind(data.stan$X_cens,rep(0,data.stan$n_cens))
          data.stan$H = ncol(data.stan$X_obs)
        }
      } 
      if (distr[x] %in% c("exponential","gompertz","weibull","weibullPH","loglogistic","lognormal")) {
        # If it's one of the others (except polyweibull), use the "h,S" format
        data.stan = list(t=data[,vars$time], d=data[,vars$event]); data.stan$n = length(data.stan$t); 
        data.stan$X = model.matrix(formula,data)
        data.stan$H = ncol(data.stan$X)
        # NB: Stan doesn't allow vectors of size 1, so if there's only one covariate (eg intercept only), needs a little trick
        if (data.stan$H==1) {
          data.stan$X = cbind(data.stan$X,rep(0,data.stan$n))
          data.stan$H = ncol(data.stan$X)
        }
      }
      if (distr[x]=="rps"){
        # If it's Royston-Parmar splines, then gets the correct data 
        knots=quantile(log(data[data[,vars$event]==1,vars$time]), seq(0, 1, length = k+2))
        # Uses flexsurv to compute the basis and derivatives of the basis
        ######################################
        basis = function (knots, x) {
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
        dbasis = function (knots, x) {
          nx <- length(x)
          if (!is.matrix(knots)) 
            knots <- matrix(rep(knots, nx), byrow = TRUE, ncol = length(knots))
          nk <- ncol(knots)
          b <- matrix(nrow = length(x), ncol = nk)
          if (nk > 0) {
            b[, 1] <- 0
            b[, 2] <- 1
          }
          if (nk > 2) {
            lam <- (knots[, nk] - knots)/(knots[, nk] - knots[, 1])
            for (j in 3:nk) {
              b[, j] <- 3 * pmax(x - knots[, j - 1], 0)^2 - 3 * 
                lam[, j - 1] * pmax(x - knots[, 1], 0)^2 - 3 * 
                (1 - lam[, j - 1]) * pmax(x - knots[, nk], 0)^2
            }
          }
          b
        }
        ######################################
        B = basis(knots,log(data[,vars$time]))
        DB = dbasis(knots,log(data[,vars$time]))
        # Now checks to see whether the user wants to specify covariates and removes the intercept from the formula (for identifiability)
        mm = model.matrix(formula,data)[,-1]
        # a. if the formula is ~ 1, then adds two fictional covariates of all 0s
        if (length(mm)<1) {
            mm = matrix(rep(0,nrow(data)),nrow=nrow(data),ncol=2)
        }
        # b. in case there's only one covariate, then adds another fake covariate of all 0s
        if (is.null(dim(mm))) {
         mm = cbind(mm,rep(0,length(mm)))
        }
        data.stan=list(t=data[,vars$time], d=data[,vars$event], n=nrow(data),M=k,X=mm,H=ncol(mm),B=B,DB=DB,
                       mu_gamma=rep(0,k+2),sigma_gamma=rep(5,k+2),knots=knots) 
      }
      # ###########################################################################################################################
      # ### Poly-Weibull is in theory possible and pre-compiled, but it poses problems if the formula is a list
      # if (distr[x]=="polyweibull") {
      #   data.stan = list(t=data[,vars$time], d=data[,vars$event]); data.stan$n = length(data.stan$t); 
      #   data.stan$M = length(formula)
      #   X = lapply(1:data.stan$M,function(i) model.matrix(formula[[i]],data))
      #   data.stan$H = max(unlist(lapply(1:data.stan$M,function(i) ncol(X[[i]]))))
      #   X = lapply(1:data.stan$M,function(i) {
      #     if(ncol(X[[i]]<data.stan$H)) {
      #       X[[i]] = cbind(X[[i]],matrix(0,nrow=nrow(X[[i]]),ncol=(data.stan$H-ncol(X[[i]]))))
      #     }
      #   })
      #   data.stan$X=array(NA,c(data.stan$M,data.stan$n,data.stan$H))
      #   for (m in 1:data.stan$M) {
      #     data.stan$X[m,,] = X[[m]]
      #   }
      # }
      # # Linear predictor coefficients
      # if (distr[x]=="polyweibull") {
      #   data.stan$mu_beta=matrix(0,nrow=data.stan$H,ncol=data.stan$M) 
      #   data.stan$sigma_beta=matrix(5,nrow=data.stan$H,ncol=data.stan$M)
      # } else {
      data.stan$mu_beta=rep(0,data.stan$H)
      if (distr[x]%in%non.on.log.scale) {
        data.stan$sigma_beta = rep(100,data.stan$H)
      } else {
        data.stan$sigma_beta = rep(5,data.stan$H)
      }
      # }
      # Ancillary parameters
      if (distr[x]=="gamma") {data.stan$a_alpha=data.stan$b_alpha = 0.1}
      if (distr[x]=="genf") {
        data.stan$a_sigma=data.stan$b_sigma=0.1
        data.stan$mu_P=0; data.stan$sigma_P=0.5
        data.stan$mu_Q=0; data.stan$sigma_Q=2.5
      }
      if (distr[x]=="gengamma") {
        data.stan$a_sigma=data.stan$b_sigma=0.1
        data.stan$mu_Q=0; data.stan$sigma_Q=100
      }
      if (distr[x] %in% c("gompertz","loglogistic","weibull","weibullPH")) {data.stan$a_alpha=data.stan$b_alpha=0.1}
      if (distr[x]=="lognormal") {data.stan$a_alpha=0; data.stan$b_alpha=5}

      # These are modified if the user gives values in the call to fit.models
      if(exists("priors",where=exArgs)) {
        priors=exArgs$priors
        # If the user has not given values for all the distrs, then fill priors with empty lists
        if(length(priors)<length(distr)) {
          for (i in (length(priors)+1):length(distr)) {
            priors[[i]] = list()
          }
        }
        # Linear predictor coefficients
        if(!is.null(priors[[x]]$mu_beta)) {
          data.stan$mu_beta=priors[[x]]$mu_beta
        }
        if(!is.null(priors[[x]]$sigma_beta)) {
          data.stan$sigma_beta = priors[[x]]$sigma_beta
        }
        if(!is.null(priors[[x]]$mu_gamma) & distr[x]=="rps") {
          data.stan$mu_gamma = priors[[x]]$mu_gamma
        }
        if(!is.null(priors[[x]]$sigma_gamma) & distr[x]=="rps") {
            data.stan$sigma_gamma = priors[[x]]$sigma_gamma
        }
        # Ancillary parameters
        if(!is.null(priors[[x]]$a_sigma)) {a_sigma=priors[[x]]$a_sigma}
        if(!is.null(priors[[x]]$b_sigma)) {b_sigma=priors[[x]]$b_sigma}
        if(!is.null(priors[[x]]$mu_P)) {mu_P=priors[[x]]$mu_P}
        if(!is.null(priors[[x]]$sigma_P)) {sigma_P=priors[[x]]$sigma_P}
        if(!is.null(priors[[x]]$mu_Q)) {mu_Q=priors[[x]]$mu_Q}
        if(!is.null(priors[[x]]$sigma_Q)) {sigma_Q=priors[[x]]$sigma_Q}
        if(!is.null(priors[[x]]$a_alpha)) {a_alpha=priors[[x]]$a_alpha}
        if(!is.null(priors[[x]]$b_alpha)) {b_alpha=priors[[x]]$b_alpha}
      }
      
      # Now runs Stan to sample from the posterior distributions
      tic = proc.time()
      out=rstan::sampling(dso[[x]],data.stan,chains=chains,iter=iter,warmup=warmup,thin=thin,seed=seed,control=control[[x]],
                      pars=pars,include=include,cores=cores)
      toc = proc.time()-tic
      time2run[x]=toc[3]
      list(out=out,data.stan=data.stan,time2run=time2run)
    })
    if(exists("save.stan",where=exArgs)) {save.stan <- exArgs$save.stan} else {save.stan=FALSE}
    #   save.stan=ifelse(any(distr=="rps"),TRUE,FALSE)
    # }
    time_survHE = unlist(lapply(1:length(mod),function(x) mod[[x]]$time2run))
    time_stan <- unlist(lapply(1:length(mod),function(x) sum(rstan::get_elapsed_time(mod[[x]]$out))))
    time2run = pmin(time_survHE,time_stan)
    names(time2run) <- labs

	# Computes the log-likelihood 
    dic <- aic <- bic <- dic2 <- numeric()
	  for (i in 1:length(distr)) {
	    # Extracts the simulations for the relevant parameters
	    beta = rstan::extract(mod[[i]]$out)$beta
	    # To safeguard against very asymmetric densities use the median (instead of the mean)
	    beta.hat = apply(beta,2,median)
	    data.stan=mod[[i]]$data.stan
	    
	    if (distr[i] %in% c("exponential","weibull","weibullPH","gompertz","lognormal","loglogistic")) {
	      linpred = beta%*%t(data.stan$X)
	      linpred.hat = beta.hat%*%t(data.stan$X)
	    }

	    if(distr[i]=="exponential") {
	      logf = matrix(unlist(lapply(1:nrow(linpred),function(i) {
	        data.stan$d*log(hexp(data.stan$t,exp(linpred[i,]))) + log(1-pexp(data.stan$t,exp(linpred[i,])))
	      })),nrow=nrow(linpred),byrow=T)
	      logf.hat = matrix(data.stan$d*log(hexp(data.stan$t,exp(linpred.hat))) + log(1-pexp(data.stan$t,exp(linpred.hat))),nrow=1)
	      # Number of parameters (for AIC): rate + covariates
	      npars = 1+sum(1-apply(data.stan$X,2,function(x) all(x==0)))
	    }

	    if (distr[i]=="weibull") {
	        shape=as.numeric(rstan::extract(mod[[i]]$out)$alpha)
	        shape.hat=median(shape)
	        logf = matrix(unlist(lapply(1:nrow(linpred),function(i) {
	          data.stan$d*log(hweibull(data.stan$t,shape[i],exp(linpred[i,]))) + log(1-pweibull(data.stan$t,shape[i],exp(linpred[i,])))
	        })),nrow=nrow(linpred),byrow=T)
	        logf.hat = matrix(data.stan$d*log(hweibull(data.stan$t,shape.hat,exp(linpred.hat))) + log(1-pweibull(data.stan$t,shape.hat,exp(linpred.hat))),nrow=1)
	        # Number of parameters (for AIC): shape, scale + covariates
	        npars = 2+sum(1-apply(data.stan$X,2,function(x) all(x==0)))
	    }
	    
	    if (distr[i]=="weibullPH") {
	      shape=as.numeric(rstan::extract(mod[[i]]$out)$alpha)
	      shape.hat=median(shape)
	      logf = matrix(unlist(lapply(1:nrow(linpred),function(i) {
	        data.stan$d*log(hweibullPH(data.stan$t,shape[i],exp(linpred[i,]))) + 
	          log(1-pweibullPH(data.stan$t,shape[i],exp(linpred[i,])))
	      })),nrow=nrow(linpred),byrow=T)
	      logf.hat = matrix(data.stan$d*log(hweibullPH(data.stan$t,shape.hat,exp(linpred.hat)))+
	                          log(1-pweibullPH(data.stan$t,shape.hat,exp(linpred.hat))),nrow=1)
	      # Number of parameters (for AIC): shape, scale + covariates
	      npars = 2+sum(1-apply(data.stan$X,2,function(x) all(x==0)))
	    }

	    if (distr[i]=="gompertz") {
	      shape=as.numeric(rstan::extract(mod[[i]]$out)$alpha)
	      shape.hat=median(shape)
	      logf = matrix(unlist(lapply(1:nrow(linpred),function(i) {
	        data.stan$d*log(hgompertz(data.stan$t,shape[i],exp(linpred[i,]))) + 
	          log(1-pgompertz(data.stan$t,shape[i],exp(linpred[i,])))
	      })),nrow=nrow(linpred),byrow=T)
	      logf.hat = matrix(data.stan$d*log(hgompertz(data.stan$t,shape.hat,exp(linpred.hat)))+
	                          log(1-pgompertz(data.stan$t,shape.hat,exp(linpred.hat))),nrow=1)
	      # Number of parameters (for AIC): shape, rate + covariates
	      npars = 2+sum(1-apply(data.stan$X,2,function(x) all(x==0)))
	    }
	    
	    if (distr[i]=="gamma") {
	        shape=as.numeric(rstan::extract(mod[[i]]$out)$alpha)
	        shape.bar=median(shape)
	        lo=exp(beta%*%t(data.stan$X_obs))
	        lc=exp(beta%*%t(data.stan$X_cens))
	        lo.bar=exp(beta.hat%*%t(data.stan$X_obs))
	        lc.bar=exp(beta.hat%*%t(data.stan$X_cens))
	        f=matrix(unlist(lapply(1:nrow(lo),function(i) dgamma(data.stan$t,shape[i],lo[i,]))),nrow=nrow(lo),byrow=T)
	        f.bar=matrix(unlist(lapply(1:nrow(lo.bar),function(i) dgamma(data.stan$t,shape.bar,lo.bar[i,]))),nrow=1,byrow=T)
	        s=matrix(unlist(lapply(1:nrow(lc),function(i) 1-pgamma(data.stan$d,shape[i],lc[i,]))),nrow=nrow(lc),byrow=T)
	        s.bar=matrix(unlist(lapply(1:nrow(lc.bar),function(i) 1-pgamma(data.stan$d,shape.bar,lc.bar[i,]))),nrow=1,byrow=T)
	        # Number of parameters (for AIC): shape, rate + covariates
	        npars = 2+sum(1-apply(data.stan$X_obs,2,function(x) all(x==0)))
	    }
	    
	    if (distr[i]=="gengamma") {
	        q=as.numeric(rstan::extract(mod[[i]]$out)$Q)
	        q.bar=median(q)
	        scale=as.numeric(rstan::extract(mod[[i]]$out)$sigma)
	        scale.bar=median(scale)
	        lo=(beta%*%t(data.stan$X_obs))
	        lc=(beta%*%t(data.stan$X_cens))
	        lo.bar=(beta.hat%*%t(data.stan$X_obs))
	        lc.bar=(beta.hat%*%t(data.stan$X_cens))
	        f=matrix(unlist(lapply(1:nrow(lo),function(i) dgengamma(data.stan$t,lo[i,],scale[i],q[i]))),nrow=nrow(lo),byrow=T)
	        f.bar=matrix(unlist(lapply(1:nrow(lo.bar),function(i) dgengamma(data.stan$t,lo.bar[i,],scale.bar,q.bar))),nrow=1,byrow=T)
	        s=matrix(unlist(lapply(1:nrow(lc),function(i) 1-pgengamma(data.stan$d,lc[i,],scale[i],q[i]))),nrow=nrow(lc),byrow=T)
	        s.bar=matrix(unlist(lapply(1:nrow(lc.bar),function(i) 1-pgengamma(data.stan$d,lc.bar[i,],scale.bar,q.bar))),nrow=1,byrow=T)
	        # Number of parameters (for AIC): mu, sigma, Q + covariates
	        npars = 3+sum(1-apply(data.stan$X_obs,2,function(x) all(x==0)))
	    }
	    
	    if (distr[i]=="genf") {
	        Q=as.numeric(rstan::extract(mod[[i]]$out)$Q)
	        Q.bar=median(Q)
	        P=as.numeric(rstan::extract(mod[[i]]$out)$P)
	        P.bar=median(P)
	        sigma=as.numeric(rstan::extract(mod[[i]]$out)$sigma)
	        sigma.bar=mean(sigma)
	        lo=(beta%*%t(data.stan$X_obs))
	        lc=(beta%*%t(data.stan$X_cens))
	        lo.bar=(beta.hat%*%t(data.stan$X_obs))
	        lc.bar=(beta.hat%*%t(data.stan$X_cens))
	        f=matrix(unlist(lapply(1:nrow(lo),function(i) dgenf(data.stan$t,lo[i,],sigma[i],Q[i],P[i]))),nrow=nrow(lo),byrow=T)
	        f.bar=matrix(unlist(lapply(1:nrow(lo.bar),function(i) dgenf(data.stan$t,lo.bar[i,],sigma.bar,Q.bar,P.bar))),nrow=1,byrow=T)
	        s=matrix(unlist(lapply(1:nrow(lc),function(i) 1-pgenf(data.stan$d,lc[i,],sigma[i],Q[i],P[i]))),nrow=nrow(lc),byrow=T)
	        s.bar=matrix(unlist(lapply(1:nrow(lc.bar),function(i) 1-pgenf(data.stan$d,lc.bar[i,],sigma.bar,Q.bar,P.bar))),nrow=1,byrow=T)
	        # Number of parameters (for AIC): mu, sigma, Q, P + covariates
	        npars = 4+sum(1-apply(data.stan$X_obs,2,function(x) all(x==0)))
	    }
	    
	    if (distr[i]=="lognormal") {
	      sigma=as.numeric(rstan::extract(mod[[i]]$out)$alpha)
	      sigma.hat=median(sigma)
	      logf = matrix(unlist(lapply(1:nrow(linpred),function(i) {
	        data.stan$d*log(hlnorm(data.stan$t,(linpred[i,]),sigma[i])) + 
	          log(1-plnorm(data.stan$t,(linpred[i,]),sigma[i]))
	      })),nrow=nrow(linpred),byrow=T)
	      logf.hat = matrix(data.stan$d*log(hlnorm(data.stan$t,(linpred.hat),sigma.hat))+
	                          log(1-plnorm(data.stan$t,(linpred.hat),sigma.hat)),nrow=1)
	      # Number of parameters (for AIC): meanlog, sdlog + covariates
	      npars = 2+sum(1-apply(data.stan$X,2,function(x) all(x==0)))
	    }

	    if (distr[i]=="loglogistic") {
	      sigma=as.numeric(rstan::extract(mod[[i]]$out)$alpha)
	      sigma.hat=median(sigma)
	      logf = matrix(unlist(lapply(1:nrow(linpred),function(i) {
	        data.stan$d*log(hllogis(data.stan$t,sigma[i],exp(linpred[i,]))) + 
	          log(1-pllogis(data.stan$t,sigma[i],exp(linpred[i,])))
	      })),nrow=nrow(linpred),byrow=T)
	      logf.hat = matrix(data.stan$d*log(hllogis(data.stan$t,sigma.hat,exp(linpred.hat)))+
	                          log(1-pllogis(data.stan$t,sigma.hat,exp(linpred.hat))),nrow=1)
	      # Number of parameters (for AIC): shape, scale + covariates
	      npars = 2+sum(1-apply(data.stan$X,2,function(x) all(x==0)))
	    }	      

	    if (distr[i]=="rps") {
	        gamma = rstan::extract(mod[[i]]$out)$gamma
	        gamma.hat = apply(gamma,2,median)
	        logf = data.stan$d*(-log(data.stan$t)+log(gamma%*%t(data.stan$DB)) + gamma%*%t(data.stan$B)+ beta%*%t(data.stan$X)) -
	                     exp(gamma%*%t(data.stan$B)+ beta%*%t(data.stan$X))
	        logf.hat = t(data.stan$d*(-log(data.stan$t)+log(data.stan$DB%*%gamma.hat)+data.stan$B%*%gamma.hat + data.stan$X%*%beta.hat) - 
	            exp(data.stan$B%*%gamma.hat + data.stan$X%*%beta.hat))
	        # Number of parameters (for AIC): gamma + covariates
	        npars = length(gamma.hat)+sum(apply(data.stan$X,2,function(x) 1-all(x==0)))
	    }	  
	    
	    # Now computes the log-likelihood and then deviance and DIC, AIC, BIC
	    # Little function to compute the log-likelihood (for the obs vs censored cases)
	    compute.loglik <- function(f,s) {
	      loglik = (apply(log(f),1,sum)+apply(log(s),1,sum))
	      return(loglik)
	    }
	    if (distr[i] %in% c("gamma","gengamma","genf")) {
	      loglik = compute.loglik(f,s); D.theta=-2*loglik 
	      loglik.bar = compute.loglik(f.bar,s.bar); D.bar=-2*loglik.bar
	      data.stan$n = data.stan$n_obs+data.stan$n_cens
	    } else {
	      loglik = apply(logf,1,sum)
	      loglik.bar = apply(logf.hat,1,sum)
	    }
	    D.theta=-2*loglik
	    D.bar=-2*loglik.bar
	    pD = mean(D.theta) - D.bar
	    pV = .5*var(D.theta)
	    dic[i] = mean(D.theta)+pD
	    dic2[i] = mean(D.theta) + pV
	    # Approximates AIC & BIC using the mean deviance and the number of nominal parameters
	    aic[i] = D.bar+2*npars                   #mean(D.theta)+2*pD
	    bic[i] = D.bar+npars*log(data.stan$n)    #mean(D.theta)+pD*log(data.stan$n)
	  }
  }
  
  # Now defines the outputs of the function
  model.fitting <- list(aic=aic,bic=bic,dic=dic)
  misc <- list(time2run=time2run,formula=formula,km=ObjSurvfit,data=data)
  if(method=="hmc") {
    misc$vars <- vars
    misc$data.stan=data.stan
	model.fitting$dic2=dic2
    # If save.stan is set to TRUE, then saves the Stan model file(s) & data
    if(save.stan==TRUE) {
      write_model <- lapply(1:length(distr),function(i) {
        model_code <- attr(mod[[i]]$out@stanmodel,"model_code")
        con <- paste0(distr[i],".stan")
        writeLines(model_code, con=con)
        txt <- paste0("Model code saved to the file: ",con,"\n")
        cat(txt)
      })
    }
    mod = lapply(1:length(mod),function(i) mod[[i]]$out)
  }
  # Names the models list
  names(mod) = names(misc$time2run)
  # Finally prepares the output object
  res <- list(models=mod,model.fitting=model.fitting,method=method,misc=misc)
  # And sets its class attribute to "survHE"
  class(res) <- "survHE"
  return(res)
}


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
            check = eval(parse(text=paste0("levels(as.factor(data$",names[w[i]],"))")))
            if (class(check)=="character") {
              # check will be 0 or 1 depending on which level of the factor is selected in newdata
              check = as.numeric(grepl(temp[j,w[i]],check))
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
              tmp = lapply(1:length(sim), function(i) {
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
        scale <- exp(linpred)^(1/-shape)
        S <- lapply(1:length(scale), function(x) cbind(t,dweibull(t,shape,scale[x])/hweibull(t,shape,scale[x]))) 
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
        scale <- exp(linpred)^(1/-shape)
        S <- lapply(1:nsim,function(i) {
          lapply(1:dim(scale)[2],function(j) {
            cbind(t,dweibull(t,shape[i],scale[i,j])/hweibull(t,shape[i],scale[i,j]))
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
            cbind(t,dlogis(log(t),scale[i,j],1/shape[i])/hllogis(log(t),scale[i,j],shape[i]))
          })
        })
      }
      if (m$dlist$name=="lognormal") {
        mulog <- linpred
        sdlog <- sim[,1]
        S <- lapply(1:nsim,function(i) {
          lapply(1:dim(scale)[2],function(j) {
            cbind(t,dlnorm(t,mulog[i,j],sdlog[i])/hlnorm(t,mulog[i,j],sdlog[i]))
          })
        })
      }
    }
  }
  
  if(fit$method=="hmc") {
      beta <- rstan::extract(m)$beta
      coefs = beta
      if(fit$models[[mod]]@model_name%in%c("Gamma","GenGamma","GenF")) {
        covmat = fit$misc$data.stan$X_obs
      } else {
        covmat = fit$misc$data.stan$X
      }
      coefs=matrix(coefs[,apply(covmat,2,function(x) 1-all(x==0))==1],nrow=nrow(beta))
      # if (is.null(fit$misc$vars$factors) & is.null(fit$misc$vars$covs)) {
      #   coefs = matrix(beta[,1],nrow=nrow(beta),byrow=T)
      # }
      if(ncol(coefs)>0) {
        if(dist!="RP") {
          colnames(coefs) = colnames(model.matrix(fit$misc$formula,fit$misc$data))
        } else {
          colnames(coefs) = colnames(model.matrix(fit$misc$formula,fit$misc$data))[-1]
        } 
      }
      basis = function (knots, x) {
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
              linpred=exp(coefs%*%t(X))
              s <- lapply(1:ncol(linpred),function(i) cbind(t,1-pexp(t,linpred[1,i])))
          }
          if (dist=="WeibullAF") {
              shape=mean(as.numeric(rstan::extract(m)$alpha))
              linpred=exp(coefs%*%t(X))
              s=lapply(1:ncol(linpred),function(j) cbind(t,1-pweibull(t,shape,linpred[1,j])))
          }
          if (dist=="WeibullPH") {
              shape=mean(as.numeric(rstan::extract(m)$alpha))
              linpred=exp(coefs%*%t(X))
              s=lapply(1:ncol(linpred),function(i) cbind(t,1-pweibullPH(t,shape,linpred[1,i])))
          }
          if (dist=="Gompertz") {
              shape=mean(as.numeric(rstan::extract(m)$alpha))
              linpred=exp(coefs%*%t(X))
              s=lapply(1:ncol(linpred),function(i) cbind(t,1-pgompertz(t,shape,linpred[1,i])))
          }
          if (dist=="Gamma") {
              shape=mean(as.numeric(rstan::extract(m)$alpha))
              linpred=exp(coefs%*%t(X))
              s=lapply(1:ncol(linpred),function(i) cbind(t,1-pgamma(t,shape,linpred[1,i])))
          }
          if (dist=="GenGamma") {
              q=mean(as.numeric(rstan::extract(m)$Q))
              scale=mean(as.numeric(rstan::extract(m)$sigma))
              linpred=(coefs%*%t(X))
              s=lapply(1:ncol(linpred),function(i) cbind(t,1-pgengamma(t,linpred[1,i],scale,q)))
          }
          if (dist=="GenF") {
              Q=mean(as.numeric(rstan::extract(m)$Q))
              P=mean(as.numeric(rstan::extract(m)$P))
              sigma=mean(as.numeric(rstan::extract(m)$sigma))
              linpred=(coefs%*%t(X))
              s=lapply(1:ncol(linpred),function(i) cbind(t,1-pgenf(t,linpred[1,i],sigma,Q,P)))
          }
          if (dist=="logNormal") {
              sigma=mean(as.numeric(rstan::extract(m)$alpha))
              linpred=(coefs%*%t(X))
              s=lapply(1:ncol(linpred),function(i) cbind(t,1-plnorm(t,linpred[1,i],sigma)))
          }
          if (dist=="logLogistic") {
              sigma=mean(as.numeric(rstan::extract(m)$alpha))
              linpred=exp(coefs%*%t(X))
              s=lapply(1:ncol(linpred),function(i) cbind(t,1-pllogis(t,scale=linpred[1,i],shape=sigma)))
          }
          if (dist=="RP") {
            # Computes the knots wrt to the times selected for the analysis
            # If there's a time=0, then add a little constant
            t[t==0] = min(0.00001,min(t[t>0]))
            B = basis(fit$misc$data.stan$knots,log(t))
            gamma=apply(rstan::extract(m)$gamma,2,mean)
            coefs=c(0,coefs)
            if(nrow(X)==1) {
              s=cbind(t,1-flexsurv::psurvspline(q=t,gamma=gamma,beta=coefs,X=X,knots=fit$misc$data.stan$knots))
            } else {
              s=lapply(1:ncol(X),function(i) 
                cbind(t,1-flexsurv::psurvspline(q=t,gamma=gamma,beta=coefs,X=X[i,],knots=fit$misc$data.stan$knots)))
            }
          }
          S[[1]] <- s
      } else {
          if (nsim>length(beta)) {nrow=length(beta)}
          if(dist=="Exponential") {
              linpred=exp(coefs%*%t(X))
              S=lapply(1:nsim,function(i) {
                  lapply(1:ncol(linpred),function(j) {
                      cbind(t,1-pexp(t,linpred[i,j]))  
                  })
              }) 
              sim = coefs[1:nsim,]
          }
          if (dist=="WeibullAF") {
              shape=as.numeric(rstan::extract(m)$alpha)
              linpred=exp(coefs%*%t(X))
              S=lapply(1:nsim,function(i) {
                  lapply(1:ncol(linpred),function(j) {
                      cbind(t,1-pweibull(t,shape[i],linpred[i,j]))  
                  })
              }) 
              sim = cbind(coefs,shape)[1:nsim,]
          }
          if (dist=="WeibullPH") {
              shape=as.numeric(rstan::extract(m)$alpha)
              linpred=exp(coefs%*%t(X))
              S=lapply(1:nsim,function(i) {
                  lapply(1:ncol(linpred),function(j) {
                      cbind(t,1-pweibullPH(t,shape[i],linpred[i,j]))  
                  })
              }) 
              sim = cbind(coefs,shape)[1:nsim,]
          }
          if (dist=="Gompertz") {
              shape=as.numeric(rstan::extract(m)$alpha)
              linpred=exp(coefs%*%t(X))
              S=lapply(1:nsim,function(i) {
                  lapply(1:ncol(linpred),function(j) {
                      cbind(t,1-pgompertz(t,shape[i],linpred[i,j]))  
                  })
              }) 
              sim = cbind(coefs,shape)[1:nsim,]
          }
          if (dist=="Gamma") {
              shape=as.numeric(rstan::extract(m)$alpha)
              linpred=exp(coefs%*%t(X))
              S=lapply(1:nsim,function(i) {
                  lapply(1:ncol(linpred),function(j) {
                      cbind(t,1-pgamma(t,shape[i],linpred[i,j]))  
                  })
              }) 
              sim = cbind(coefs,shape)[1:nsim,]
          }
          if (dist=="GenGamma") {
              Q=as.numeric(rstan::extract(m)$Q)
              shape=as.numeric(rstan::extract(m)$sigma)
              linpred=(coefs%*%t(X))
              S=lapply(1:nsim,function(i) {
                  lapply(1:ncol(linpred),function(j) {
                      cbind(t,1-pgengamma(t,linpred[i,j],shape[i],Q[i]))  
                  })
              }) 
              sim = cbind(coefs,shape,Q)[1:nsim,]
          }
          if (dist=="GenF") {
              Q=as.numeric(rstan::extract(m)$Q)
              P=as.numeric(rstan::extract(m)$P)
              sigma=as.numeric(rstan::extract(m)$sigma)
              linpred=(coefs%*%t(X))
              S=lapply(1:nsim,function(i) {
                  lapply(1:ncol(linpred),function(j) {
                      cbind(t,1-pgenf(t,linpred[i,j],sigma[i],Q[i],P[i]))  
                  })
              }) 
              sim = cbind(coefs,sigma,Q,P)[1:nsim,]
          }
          if (dist=="logNormal") {
              sigma=as.numeric(rstan::extract(m)$alpha)
              linpred=(coefs%*%t(X))
              S=lapply(1:nsim,function(i) {
                  lapply(1:ncol(linpred),function(j) {
                      cbind(t,1-plnorm(t,linpred[i,j],sigma[i]))  
                  })
              }) 
              sim = cbind(coefs,sigma)[1:nsim,]
          }
          if (dist=="logLogistic") {
              sigma=as.numeric(rstan::extract(m)$alpha)
              linpred=exp(coefs%*%t(X))
              S=lapply(1:nsim,function(i) {
                  lapply(1:ncol(linpred),function(j) {
                      cbind(t,1-pllogis(t,linpred[i,j],sigma[i]))  
                  })
              }) 
              sim = cbind(coefs,sigma)[1:nsim,]
          }
          if (dist=="RP") {
            # Computes the knots wrt to the times selected for the analysis
            t[t==0] = min(0.00001,min(t[t>0]))
            B = basis(fit$misc$data.stan$knots,log(t))
            gamma=rstan::extract(m)$gamma
            coefs=cbind(rep(0,nrow(coefs)),coefs)
            if(nrow(X)==1) {
              S=lapply(1:nsim,function(i) {
                lapply(1,function(j) {
                  cbind(t,1-flexsurv::psurvspline(q=t,gamma=gamma[i,],beta=coefs[i,],X=X,knots=fit$misc$data.stan$knots))
                })
              })
            } else {
              S = lapply(1:nsim,function(i) {
                lapply(1:ncol(X),function(j) {
                  cbind(t,1-flexsurv::psurvspline(q=t,gamma=gamma[i,],beta=coefs[i,],X=X[j,],knots=fit$misc$data.stan$knots))
                })
              })
            }
            sim = cbind(coefs[,-1],gamma)[1:nsim,]
          }
      }
  }

  n.elements <- length(S[[1]]) 
  if (fit$method=="mle") {
    mat = lapply(1:n.elements,function(j) matrix(unlist(lapply(1:nsim,function(i) S[[i]][[j]][,2])),nrow=nsim,byrow=T))
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



print.survHE <- function(x,mod=1,...) {
  # Creates a print method for the objects in the class survHE
  # x is the survHE object (the output of the call to fit.models)
  # mod is the index of the model. Default value is 1, but the user can choose which model fit to visualise, 
  #     if the call to fit.models has a vector argument for distr (so many models are fitted & stored in the same object)
  # ... optional arguments
  # digits = number of significant digits to be shown in the summary table (default = 6)
  # nsim = number of simulations from the joint posterior for INLA (default = 100)
  # original = a flag to say whether the *original* table from either INLA or MCMC should be printed
  
  exArgs <- list(...)
  ##  if(exists("original",where=exArgs)) {original=exArgs$original} else {original=FALSE}
  
  # Available models
  availables.mle <- c("genf", "genf.orig", "gengamma", "gengamma.orig", "exp", 
                      "weibull", "weibull.quiet","weibullPH", "lnorm", "gamma", "gompertz", 
                      "llogis", "exponential", "lognormal","survspline")
  availables.inla <- c("exponential","weibull","lognormal","loglogistic")
  availables.hmc <- c("Exponential","Gamma","GenF","GenGamma","Gompertz","PolyWeibull","RP",
                      "WeibullAF","WeibullPH","logLogistic","logNormal")
  # If the distribution specified is not-standard (eg user-defined in MLE, or using random effects or non-standard
  # distributions in Stan), then sets original=TRUE and gives the original version of the print table.
  if (exists("original",where=exArgs)) {original=exArgs$original} else {
    if (x$method=="mle") {
      if (x$models[[mod]]$dlist$name %in% availables.mle) {original=FALSE} else {original=TRUE}
    }
    if (x$method=="inla") {
      if (x$models[[mod]]$dlist$name %in% availables.inla) {original=FALSE} else {original=TRUE}
    }
    if (x$method=="hmc") {
      if (x$models[[mod]]@model_name %in% availables.hmc) {original=FALSE} else {original=TRUE}
    }
  }
  # Can select the number of digits to be printed in the output table
  if(!exists("digits",where=exArgs)){digits=6} else {digits=exArgs$digits}
  
  if(x$method =="mle") {
    res <- x$models[[mod]]$res[,c(1,4,2,3)]
    if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
  }
  if(x$method=="inla" & original==FALSE) {
    # Rescales the parameters to make the estimates comparable with flexsurvreg
    if(!exists("nsim",where=exArgs)){nsim <- 100} else {nsim=exArgs$nsim} 
    
    # This is a rescaling function for the built-in models (that INLA can do by default)
    rescale.print.inla <- function(x,mod,nsim) {
      # Simulates from the joint posterior of *all* parameters & hyperparameters
      jpost <- suppressWarnings(INLA::inla.posterior.sample(n=nsim,x$models[[mod]]))
      # This finds the position of the hyperparameters in the simulations from the joint posterior
      pos <- pmatch(rownames(x$models[[mod]]$summary.fixed),rownames(jpost[[1]]$latent))
      
      if(x$models[[mod]]$dlist=="weibull") {
        shape <- unlist(lapply(jpost,function(x) x$hyperpar)); names(shape) <- NULL
        scale <- exp(unlist(lapply(jpost,function(x) x$latent[pos[1],])))^(1/-shape)
        effects <- matrix(NA,nrow=(length(pos)-1),ncol=nsim)
        if(length(attributes(terms(x$misc$formula))$term.labels)>0) {
          for (j in 2:length(pos)) {
            effects[(j-1),] <- log(exp(unlist(lapply(jpost,function(x) x$latent[pos[j],])))^(1/-shape))
          }
          rownames(effects) <- x$models[[mod]]$names.fixed[-1]
        }
        tab <- rbind(shape,scale,effects)
      }
      if(x$models[[mod]]$dlist=="exponential") {
        rate <- exp(unlist(lapply(jpost,function(x) x$latent[pos[1],])))
        effects <- matrix(NA,nrow=(length(pos)-1),ncol=nsim)
        if(length(attributes(terms(x$misc$formula))$term.labels)>0) {
          for (j in 2:length(pos)) {
            effects[(j-1),] <- unlist(lapply(jpost,function(x) x$latent[pos[j],]))
          }
          rownames(effects) <- x$models[[mod]]$names.fixed[-1]
        }
        tab <- rbind(rate,effects)
      }
      if(x$models[[mod]]$dlist=="lognormal") {
        prec <- unlist(lapply(jpost,function(x) x$hyperpar)); names(prec) <- NULL
        sdlog <- 1/sqrt(prec)
        meanlog <- unlist(lapply(jpost,function(x) x$latent[pos[1],]))
        effects <- matrix(NA,nrow=(length(pos)-1),ncol=nsim)
        if(length(attributes(terms(x$misc$formula))$term.labels)>0) {
          for (j in 2:length(pos)) {
            effects[(j-1),] <- unlist(lapply(jpost,function(x) x$latent[pos[j],]))
          }
          rownames(effects) <- x$models[[mod]]$names.fixed[-1]
        }
        tab <- rbind(meanlog,sdlog,effects)
      }
      if(x$models[[mod]]$dlist=="loglogistic") {
        shape <- unlist(lapply(jpost,function(x) x$hyperpar)); names(shape) <- NULL
        scale <- exp(unlist(lapply(jpost,function(x) x$latent[pos[1],])))
        effects <- matrix(NA,nrow=(length(pos)-1),ncol=nsim)
        if(length(attributes(terms(x$misc$formula))$term.labels)>0) {
          for (j in 2:length(pos)) {
            effects[(j-1),] <- unlist(lapply(jpost,function(x) x$latent[pos[j],]))
          }
          rownames(effects) <- x$models[[mod]]$names.fixed[-1]
        }
        tab <- rbind(shape,scale,effects)
      }
      return(tab)
    }
    # The user could specify a rescale.print function for their own specific model and that would be used instead
    if(exists("rescale.print",where=exArgs)) {
      func <- exArgs$rescale.print
      if(exists("inputs",where=exArgs)) {
        inputs=exArgs$inputs
      } else {
        inputs=list()
      }
    } else {
      func <- rescale.print.inla
      inputs <- list(x,mod,nsim)
    }
    tab <- do.call(what=func,args=inputs)
    
    res <- t(apply(tab,1,function(x) c(mean(x),sd(x),quantile(x,.025),quantile(x,.975))))
    if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
  }
  
  if(x$method=="hmc") {
    quiet <- function(x) { 
      sink(tempfile()) 
      on.exit(sink()) 
      invisible(force(x)) 
    } 
    # If the model is intercept only or only one covariate, then gets rid of unnecessary beta's
    quiet(print(x$models[[mod]]))
    table <- cbind(x$models[[mod]]@.MISC$summary$msd,x$models[[mod]]@.MISC$summary$quan[,c("2.5%","97.5%")])
    take.out = which(rownames(table)=="lp__")
    betas = grep("beta",rownames(table))
    if(x$models[[mod]]@model_name%in%c("Gamma","GenGamma","GenF")) {
      covmat = x$misc$data.stan$X_obs
    } else {
      covmat = x$misc$data.stan$X
    }
    take.out = c(take.out,betas[apply(covmat,2,function(x) all(x==0))])
    # if (is.null(x$misc$vars$factors) & is.null(x$misc$vars$covs)) {
    #   take.out = c(take.out,which(rownames(table)=="beta[2]"))
    # }
    table=table[-take.out,]
  
    if (original==FALSE) {
      if (x$models[[mod]]@model_name=="Exponential") {
        rate <- matrix(table[grep("rate",rownames(table)),],ncol=4)
        rownames(rate) <- "rate"
        if(length(attributes(terms(x$misc$formula))$term.labels)>0) {
          effects=matrix(table[-which(rownames(table) %in% c("rate")),][-1,],ncol=4,byrow=F)
          rownames(effects) = colnames(model.matrix(x$misc$formula,x$misc$data))[-1]
        } else {
          effects = matrix(NA,nrow=0,ncol=4)
        }
        res = rbind(rate,effects)
        if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
      }
      
      if (x$models[[mod]]@model_name=="Gamma") {
        rate <- matrix(table[grep("rate",rownames(table)),],ncol=4)
        rownames(rate) <- "rate"
        shape <- matrix(table[grep("alpha",rownames(table)),],ncol=4)
        rownames(shape) <- "shape"
        if(length(attributes(terms(x$misc$formula))$term.labels)>0) {
          effects=matrix(table[-which(rownames(table) %in% c("rate","alpha")),][-1,],ncol=4,byrow=F)
          rownames(effects) = colnames(model.matrix(x$misc$formula,x$misc$data))[-1]
        } else {
          effects = matrix(NA,nrow=0,ncol=4)
        }
        res = rbind(shape,rate,effects)
        if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
      }
      
      if (x$models[[mod]]@model_name=="WeibullAF" | x$models[[mod]]@model_name=="WeibullPH") {
        scale <- matrix(table[grep("scale",rownames(table)),],ncol=4)
        rownames(scale) <- "scale"
        shape <- matrix(table[grep("alpha",rownames(table)),],ncol=4)
        rownames(shape) <- "shape"
        if(length(attributes(terms(x$misc$formula))$term.labels)>0) {
          effects=matrix(table[-which(rownames(table) %in% c("scale","alpha")),][-1,],ncol=4,byrow=F)
          rownames(effects) = colnames(model.matrix(x$misc$formula,x$misc$data))[-1]
        } else {
          effects = matrix(NA,nrow=0,ncol=4)
        }
        res = rbind(shape,scale,effects)
        if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
      }
      
      if (x$models[[mod]]@model_name=="Gompertz") {
        rate <- matrix(table[grep("rate",rownames(table)),],ncol=4)
        rownames(rate) <- "rate"
        shape <- matrix(table[grep("alpha",rownames(table)),],ncol=4)
        rownames(shape) <- "shape"
        if(length(attributes(terms(x$misc$formula))$term.labels)>0) {
          effects=matrix(table[-which(rownames(table) %in% c("rate","alpha")),][-1,],ncol=4,byrow=F)
          rownames(effects) = colnames(model.matrix(x$misc$formula,x$misc$data))[-1]
        } else {
          effects = matrix(NA,nrow=0,ncol=4)
        }
        res = rbind(shape,rate,effects)
        if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
      }
      
      if (x$models[[mod]]@model_name=="logNormal") {
        meanlog <- matrix(table[grep("meanlog",rownames(table)),],ncol=4)
        rownames(meanlog) <- "meanlog"
        sigma <- matrix(table[grep("alpha",rownames(table)),],ncol=4)
        rownames(sigma) <- "sdlog"
        if(length(attributes(terms(x$misc$formula))$term.labels)>0) {
          effects=matrix(table[-which(rownames(table) %in% c("meanlog","alpha")),][-1,],ncol=4,byrow=F)
          rownames(effects) = colnames(model.matrix(x$misc$formula,x$misc$data))[-1]
        } else {
          effects = matrix(NA,nrow=0,ncol=4)
        }
        res = rbind(meanlog,sigma,effects)
        if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
      }
      
      if (x$models[[mod]]@model_name=="logLogistic") {
        rate <- matrix(table[grep("rate",rownames(table)),],ncol=4)
        rownames(rate) <- "rate"
        shape <- matrix(table[grep("alpha",rownames(table)),],ncol=4)
        rownames(shape) <- "shape"
        if(length(attributes(terms(x$misc$formula))$term.labels)>0) {
          effects=matrix(table[-which(rownames(table) %in% c("rate","alpha")),][-1,],ncol=4,byrow=F)
          rownames(effects) = colnames(model.matrix(x$misc$formula,x$misc$data))[-1]
        } else {
          effects = matrix(NA,nrow=0,ncol=4)
        }
        res = rbind(shape,rate,effects)
        if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
      }
      
      if (x$models[[mod]]@model_name=="GenF") {
        mu <- matrix(table[grep("beta",rownames(table)),],ncol=4,nrow=1)
        rownames(mu) <- "mu"
        sigma <- matrix(table[grep("sigma",rownames(table)),],ncol=4)
        rownames(sigma) <- "sigma"
        Q <- matrix(table[grep("Q",rownames(table)),],ncol=4)
        rownames(Q) <- "Q"
        P <- matrix(table[match("P",rownames(table)),],ncol=4)
        rownames(P) <- "P"
        if(length(attributes(terms(x$misc$formula))$term.labels)>0) {
          effects=matrix(table[-which(rownames(table) %in% c("beta[1]","sigma","Q","P")),],ncol=4,byrow=F)
          rownames(effects) = colnames(model.matrix(x$misc$formula,x$misc$data))[-1]
        } else {
          effects = matrix(NA,nrow=0,ncol=4)
        }
        res <- rbind(mu,sigma,Q,P,effects)
        if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
      }
      
      if (x$models[[mod]]@model_name=="GenGamma") {
        mu <- matrix(table[grep("beta",rownames(table)),][1,],ncol=4,nrow=1)
        rownames(mu) <- "mu"
        sigma <- matrix(table[grep("sigma",rownames(table)),],ncol=4)
        rownames(sigma) <- "sigma"
        Q <- matrix(table[grep("Q",rownames(table)),],ncol=4)
        rownames(Q) <- "Q"
        if(length(attributes(terms(x$misc$formula))$term.labels)>0) {
          effects=matrix(table[-which(rownames(table) %in% c("beta[1]","Q","sigma")),],ncol=4,byrow=F)
          rownames(effects) = colnames(model.matrix(x$misc$formula,x$misc$data))[-1]
        } else {
          effects = matrix(NA,nrow=0,ncol=4)
        }
        res <- rbind(mu,sigma,Q,effects)
        if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
      }
      if (x$models[[mod]]@model_name=="RP") {
        # # First remove "fake covariates" (used to trick Stan into having a formula with only 1 or 0 covariates for Xbeta)
        # betas = grep("beta",rownames(table))
        # take.out = betas[apply(x$misc$data.stan$X,2,function(x) all(x==0))]
        # if(length(take.out)>0) {table = table[-take.out,]}
        # Now formats the gammas
        gamma <- matrix(table[grep("gamma",rownames(table)),],ncol=4)
        rownames(gamma) <- paste0("gamma",0:(nrow(gamma)-1))
        if(length(attributes(terms(x$misc$formula))$term.labels)>0) {
          effects=matrix(table[-grep("gamma",rownames(table)),],ncol=4,byrow=F)
          rownames(effects) = colnames(model.matrix(x$misc$formula,x$misc$data))[-1]
        } else {
          effects = matrix(NA,nrow=0,ncol=4)
        }
        res <- rbind(gamma,effects)
        if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
      }
      if (x$models[[mod]]@model_name=="PolyWeibull") {
        alpha = matrix(table[grep("alpha",rownames(table)),],ncol=4)
        rownames(alpha) <- paste0("shape_",1:x$misc$data.stan$M)
        to.rm=matrix(unlist(lapply(1:length(x$misc$formula),function(m) apply(x$misc$data.stan$X[m,,],2,function(x) all(x==0)))),
                     nrow=length(x$misc$formula),byrow=T)
		nmatch = length(which(to.rm==T))
		idx = matrix(unlist(lapply(1:nmatch,function(i) {
			paste0(which(to.rm==TRUE,arr.ind=T)[i,],collapse=",")
		})))
		if (!is.null(nrow(idx))) {
			take.out = match(paste0("beta[",idx,"]"),rownames(table))
		}
        if(all(!is.na(take.out))) {table=table[-take.out,]}
        effects=table[-grep("alpha",rownames(table)),]
        rownames(effects) = unlist(lapply(1:x$misc$data.stan$M,function(m) {
          paste0(colnames(model.matrix(x$misc$formula[[m]],x$misc$data)),"_",m)
        }))
        res = rbind(alpha,effects)
        if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
      }
    }
  }

  # Finally creates the table
  # Original formatting of the tables from INLA & Stan
  if(original==TRUE) {
    if (x$method=="mle") {
      print(x$models[[mod]])
    }
    if (x$method=="inla") {
      print(summary(x$models[[mod]]))
    }
    if (x$method=="hmc") {
	  if (x$models[[mod]]@model_name=="PolyWeibull") {
	      take.out = betas[unlist(lapply(1:length(x$misc$formula),function(m) apply(x$misc$data.stan$X[m,,],2,function(x) all(x==0))))]
	  } else {
	      take.out = betas[unlist(lapply(1:length(x$misc$formula),function(m) apply(covmat,2,function(x) all(x==0))))]
	  }
      take.out = c(take.out,grep("lp__",rownames(rstan::summary(x$models[[mod]])$summary)))
      tab=rstan::summary(x$models[[mod]],probs=c(.025,.975))$summary[-take.out,]
      n_kept <- x$models[[mod]]@sim$n_save - x$models[[mod]]@sim$warmup2
      cat("Inference for Stan model: ", x$models[[mod]]@model_name, ".\n", sep = "")
      cat(x$models[[mod]]@sim$chains, " chains, each with iter=", x$models[[mod]]@sim$iter, 
          "; warmup=", x$models[[mod]]@sim$warmup, "; thin=", x$models[[mod]]@sim$thin, "; \n", 
          "post-warmup draws per chain=", n_kept[1], ", ", "total post-warmup draws=", 
          sum(n_kept), ".\n\n", sep = "")
      print(tab,digits=digits)
      sampler <- attr(x$models[[mod]]@sim$samples[[1]], "args")$sampler_t
      cat("\nSamples were drawn using ", sampler, " at ", x$models[[mod]]@date, 
          ".\n", "For each parameter, n_eff is a crude measure of effective sample size,\n", 
          "and Rhat is the potential scale reduction factor on split chains (at \n", 
          "convergence, Rhat=1).\n", sep = "")
    }
  } else {
    # FORMATS THE TABLE
    # Now recodes the model name to a standardised string
    if (x$method=="hmc") {
      label <- x$models[[mod]]@model_name
      if (label=="RP") {label="Royston & Parmar splines"}
      label.met="Stan (Bayesian inference via \nHamiltonian Monte Carlo)"
    } else {
      if(x$models[[mod]]$dlist$name=="exp" | x$models[[mod]]$dlist$name=="exponential") {label="Exponential"}
      if(x$models[[mod]]$dlist$name=="gamma") {label="Gamma"}
      if(x$models[[mod]]$dlist$name=="lognormal" | x$models[[mod]]$dlist$name=="lnorm") {label="log-Normal"}
      if(x$models[[mod]]$dlist$name=="llogis" | x$models[[mod]]$dlist$name=="loglogistic") {label="log-Logistic"}
      if(x$models[[mod]]$dlist$name=="gengamma") {label="Generalised Gamma"}
      if(x$models[[mod]]$dlist$name=="weibull" | x$models[[mod]]$dlist$name=="weibull.quiet" | x$models[[mod]]$dlist$name=="weibullPH") {label="Weibull"}
      if(x$models[[mod]]$dlist$name=="genf") {label="Generalised F"}
      if(x$models[[mod]]$dlist$name=="gompertz") {label="Gompertz"}
      if(x$models[[mod]]$dlist$name=="survspline") {label="Royston & Parmar splines"}
    }
    if(x$method=="mle") {label.met="Flexsurvreg \n(Maximum Likelihood Estimate)"}
    if(x$method=="inla") {label.met="INLA (Bayesian inference via \nIntegrated Nested Laplace Approximation)"}

    cat("\n")
    cat(paste0("Model fit for the ",label," model, obtained using ",label.met,". Running time: ",
               format(x$misc$time2run[[mod]],digits=5,nsmall=3)," seconds"))
    cat("\n\n")
    print(res,quote=F,digits=digits,justify="center")
    cat("\n")
    cat("Model fitting summaries\n")
    cat(paste0("Akaike Information Criterion (AIC)....: ",format(x$model.fitting$aic[[mod]],digits=6,nsmall=3)))
    cat("\n")
    cat(paste0("Bayesian Information Criterion (BIC)..: ",format(x$model.fitting$bic[[mod]],digits=6,nsmall=3)))
    if(x$method=="inla" | x$method=="hmc") {
      cat("\n")
      cat(paste0("Deviance Information Criterion (DIC)..: ",format(x$model.fitting$dic[[mod]],digits=6,nsmall=3)))
    }
    cat("\n\n")
  }
}


psa.plot <- function(psa,...) {
  # Plots the survival curves for all the PSA simulations
  # psa = the result of the call to the function make.surv
  # ... = additional arguments
  # xlab = label for the x-axis
  # ylab = label for the y-axis
  # col = vector of colours with which to plot the curves
  # alpha = parameter to determine the transparency (default = 0.1)
  # main = a string to write the title
  # labs = logical (default = TRUE): should text to identify which profile has been plotted printed on the graph?
  # xpos = the point on the x-axis in which to write the legend with the profiles (default = 65% of the x-axis)
  # ypos = the point on the y-axis in which to write the legend with the profiles (default = 100% of the y-axis)
  # cex.txt = the factor by which to write the text (default = .75)
  # offset = how much space between text for each line? (default = .35)
  # nsmall = number of decimal places (default = 2)
  # digits = number of digits used for the numerical values in the labels
  
  n.elements <- length(psa$S[[1]]) 
  times <- psa$S[[1]][[1]][,1]
  exArgs <- list(...)
  if(!exists("xlab",where=exArgs)) {xlab="Time"} else {xlab=exArgs$xlab}
  if(!exists("ylab",where=exArgs)) {ylab="Survival"} else {ylab=exArgs$ylab}
  if(!exists("col",where=exArgs)) {col=sample(colors(),n.elements)} else {col=exArgs$col}
  if(!exists("alpha",where=exArgs)) {alpha=0.1} else {alpha=exArgs$alpha}
  if(!exists("main",where=exArgs)) {main=""} else {main=exArgs$main}
  if(!exists("labs",where=exArgs)) {labs=TRUE} else {labs=exArgs$labs}
  if(!exists("xpos",where=exArgs)) {xpos=max(times)*0.65} else {xpos=exArgs$xpos}
  if(!exists("ypos",where=exArgs)) {ypos=1} else {ypos=exArgs$ypos}
  if(!exists("cex.txt",where=exArgs)) {cex.txt=.75} else {cex.txt=exArgs$cex.txt}
  if(!exists("offset",where=exArgs)) {off=seq(1,nrow(psa$des.mat))*.35} else {off=seq(1,nrow(psa$des.mat))*exArgs$offset}
  if(!exists("nsmall",where=exArgs)) {nsmall=2} else {nsmall=exArgs$nsmall}
  if(!exists("digits",where=exArgs)) {digits=5} else {digits=exArgs$digits}

  # If there's only the average value for the survival curve, simpler plot
  if (psa$nsim==1) {
    alpha <- 1
    plot(psa$S[[1]][[1]][,1:2],t="l",xlab=xlab,ylab=ylab,col=adjustcolor(col[1],alpha),ylim=c(0,1),xlim=range(pretty(times)),
         main=main,axes=F)
    if (n.elements>1) {
      pts2 <- lapply(2:n.elements,function(i) points(psa$S[[1]][[i]],t="l",col=adjustcolor(col[i],alpha)))
    }
  }
  
  # If there are nsim simulations from the survival curves, then more complex plot
  if (psa$nsim>1) {
    tmp <- lapply(1:n.elements,function(j) matrix(unlist(lapply(1:psa$nsim,function(i) psa$S[[i]][[j]][,2])),nrow=psa$nsim,byrow=T))
    q025 <- lapply(1:n.elements, function(j) apply(tmp[[j]],2,function(x) quantile(x,.025)))
    q500 <- lapply(1:n.elements, function(j) apply(tmp[[j]],2,function(x) quantile(x,.5))) 
    q975 <- lapply(1:n.elements, function(j) apply(tmp[[j]],2,function(x) quantile(x,.975))) 
    plot(psa$S[[1]][[1]][,1],q500[[1]],col=adjustcolor(col[1],1),t="l",xlab=xlab,ylab=ylab,ylim=c(0,1),xlim=range(pretty(times)),lwd=2,main=main,axes=F)
    polygon(c(psa$S[[1]][[1]][,1],rev(psa$S[[1]][[1]][,1])),c(q975[[1]],rev(q025[[1]])),col=adjustcolor(col[1],alpha),border=NA)
    if (n.elements==1) {
    }
    if (n.elements>1) {
      lapply(2:n.elements, function(i) {
        pts1 <- points(psa$S[[1]][[i]][,1],q500[[i]],col=adjustcolor(col[i],1),t="l",lwd=2) 
        pts2 <- polygon(c(psa$S[[1]][[i]][,1],rev(psa$S[[1]][[i]][,1])),c(q975[[i]],rev(q025[[i]])),col=adjustcolor(col[i],alpha),border=NA)
      })
    }
  }
  axis(1);axis(2)
  if (labs==TRUE) {
    txt1=lapply(1:ncol(psa$des.mat),function(i) {
      text(xpos,ypos-(i-1)/40,paste0(colnames(psa$des.mat)[i]," : "),cex=cex.txt,pos=2,col="black")
    })
    txt2=lapply(1:ncol(psa$des.mat),function(i) {
      lapply(1:nrow(psa$des.mat),function(j) {
        text(xpos+off[j],ypos-(i-1)/40,format(psa$des.mat[j,i],nsmall=nsmall,digits=digits),cex=cex.txt,pos=2,col=col[j])
      })
    })
  }
}



plot.survHE <- function(...) {
    ## Plots the KM + the results of the model fitted by fit.models()
    ## Uses different commands, depending on which method has been used to fit the models
    #
    # x = the result of the call to the fit.model function. Can be x,y,z,... (each survHE objects)
    #
    # mod = a numeric vector --- selects the models to plot (so mod=c(1,3) only selects the 1st and 3rd arguments)
    # xlab
    # ylab
    # lab.trt
    # cex.trt
    # n.risk
    # xlim
    # colors
    # labs
    # add.km = TRUE (whether to also add the Kaplan Meier estimates of the data)
    # newdata = a list (of lists), specifiying the values of the covariates at which the computation is performed. For example
    #           'list(list(arm=0),list(arm=1))' will create two survival curves, one obtained by setting the covariate 'arm'
    #           to the value 0 and the other by setting it to the value 1. In line with 'flexsurv' notation, the user needs
    #           to either specify the value for *all* the covariates or for none (in which case, 'newdata=NULL', which is the
    #           default). If some value is specified and at least one of the covariates is continuous, then a single survival
    #           curve will be computed in correspondence of the average values of all the covariates (including the factors, 
    #           which in this case are expanded into indicators). 
    
    exArgs <- list(...) 		# Lists all the additional inputs
    nexArgs <- length(exArgs)
    classes <- unlist(lapply(1:nexArgs,function(i) class(exArgs[[i]])))
    w=which(classes=="survHE")
    original.method=unlist(lapply(w,function(i) exArgs[[i]]$method))
    if(length(w)==0) {
        stop("You need to input at least one 'survHE' object to run this function!")
    }
    if(length(w)==1) {
      totmodels <- unlist(lapply(w,function(i) length(exArgs[[i]]$models)))
      mods <- exArgs[[w]]$models
      method <- rep(exArgs[[w]]$method,totmodels) 
      aic <- unlist(exArgs[[w]]$model.fitting$aic)
      bic <- unlist(exArgs[[w]]$model.fitting$bic)
      dic <- unlist(exArgs[[w]]$model.fitting$dic)
      if(totmodels>1){
        if (!is.null(exArgs$mod)) {which.model <- exArgs$mod} else {which.model <- 1:length(mods)}
        mods <- lapply(which.model,function(i) mods[[i]])
        method <- method[which.model]
        aic <- aic[which.model]
        bic <- bic[which.model]
        dic <- dic[which.model]
      } 
    }
    if (length(w)>1) {
      mods <- unlist(lapply(w,function(i) exArgs[[i]]$models),recursive = F)
      totmodels <- unlist(lapply(w,function(i) length(exArgs[[i]]$models)))
      method <- unlist(lapply(w,function(i) rep(exArgs[[i]]$method,totmodels[i])))
      aic <- unlist(lapply(w,function(i) exArgs[[i]]$model.fitting$aic))
      bic <- unlist(lapply(w,function(i) exArgs[[i]]$model.fitting$bic))
      dic <- unlist(lapply(w,function(i) exArgs[[i]]$model.fitting$dic))
      if (!is.null(exArgs$mod)) {which.model <- exArgs$mod} else {which.model <- 1:length(mods)}
      mods <- lapply(which.model,function(i) mods[[i]])
      method <- method[which.model]
      aic <- aic[which.model]
      bic <- bic[which.model]
      dic <- dic[which.model]
    }
    model.fitting <- list(aic=aic,bic=bic,dic=dic)
    x <- list(); x$models <- mods; class(x) <- "survHE"
    x$model.fitting <- model.fitting
    ## Needs to include in the misc object the element vars (which is used for HMC models)
    if (any(method=="hmc")) {
      x$misc=exArgs[[min(which(original.method=="hmc"))]]$misc
    } else {
      # If none of the survHE objects are HMC, then just use the first
      x$misc=exArgs[[1]]$misc
    }

    nmodels <- length(x$models)  # Number of models fitted by fit.models()
    
    # Checks that extra options are specified
    if (is.null(exArgs$t)) {t=sort(unique(x$misc$km$time))} else {t=exArgs$t}
    if (is.null(exArgs$xlab)) {xl="time"} else {xl=exArgs$xlab}
    if (is.null(exArgs$ylab)) {yl="Survival"} else {yl=exArgs$ylab}
    if (is.null(exArgs$lab.trt)) {lab.trt=names(x$misc$km$strata)} else {lab.trt=names(x$km$strata)<-exArgs$lab.trt}
    if (is.null(exArgs$cex.trt)) {cex.trt=.8} else {cex.trt=exArgs$cex.trt}
    if (is.null(exArgs$n.risk)) {nrisk=FALSE} else {nrisk=exArgs$n.risk}
    if (is.null(exArgs$main)) {main=""} else {main=exArgs$main}
    if (is.null(exArgs$newdata)) {newdata = NULL} else {newdata=exArgs$newdata}
    if (is.null(exArgs$cex.lab)) {cex.lab = .8} else {cex.lab=exArgs$cex.lab}
    
    if (is.null(exArgs$xlim) & is.null(exArgs$t)) {
        xlm=range(pretty(x$misc$km$time))
    } 
    if (is.null(exArgs$xlim) & !is.null(exArgs$t)) {
        xlm=range(pretty(t))
    }
    if (!is.null(exArgs$xlim) & is.null(exArgs$t)) {
        xlm <- exArgs$xlim
    }
    if (!is.null(exArgs$xlim) & !is.null(exArgs$t)) {
        xlm <- exArgs$xlim
    }
    
    if (is.null(exArgs$colors)) {
        if (nmodels>1) {colors=(2:(nmodels+1))} else {colors=2}
    } else {colors=exArgs$colors}
    if(is.null(exArgs$axes)){axes=TRUE} else {axes=exArgs$axes}
    if (is.null(exArgs$labs)) {
        labs <- unlist(lapply(1:length(x$models),function(i) {
            if(class(x$models[[i]])=="stanfit") {tolower(x$models[[i]]@model_name)} else {x$models[[i]]$dlist$name}
        }))
        labs[labs %in% c("weibull.quiet","weibull","weibullaf","weibullph")] = "Weibull"
        labs[labs %in% c("exp","exponential")] = "Exponential"
        labs[labs %in% "gamma"] <- "Gamma"
        labs[labs %in% c("lnorm","lognormal")] = "log-Normal"
        labs[labs %in% c("llogis","loglogistic","loglogis")] = "log-Logistic"
        labs[labs %in% "gengamma"] <- "Gen. Gamma"
        labs[labs %in% "genf"] <- "Gen. F"
        labs[labs %in% "gompertz"] <- "Gompertz"
        labs[labs %in% c("survspline","rp")] <- "Royston & Parmar splines"
    } else {labs=exArgs$labs}
    labs <- c("Kaplan Meier",labs)
    if (is.null(exArgs$add.km)) {add.km=TRUE} else {add.km=exArgs$add.km}
    
    # Now plots the KM curve using "rms" if add.km is set to TRUE
    if (add.km==TRUE & is.null(newdata)) {
        rms::survplot(x$misc$km,                                     # Specialised plot from "rms" 
                      xlab=xl,ylab=yl,		                           # x- and y- labels
                      label.curves=list(labels=lab.trt,cex=cex.trt), # specifies curve labels
                      n.risk=nrisk,   	                             # tells R to show number at risk 
                      lwd=2,xlim=xlm  	                             # defines the size of the lines (2 pts)
        )
      col <- c("black",colors)
      title(main)
    } else {
        labs <- labs[-1]
        if(class(colors)!="character") {colors <- colors-1}
        plot(0,0,col="white",xlab=xl,ylab=yl,axes=F,xlim=xlm,ylim=c(0,1),main=main)
        if(axes==TRUE) {axis(1);axis(2)}
        col <- colors
    }
    res <- lapply(1:nmodels,function(i) {
        x$method=method[i]
        make.surv(x,nsim=1,t=t,mod=i,newdata=newdata)
    })
    
    if (!is.null(newdata)) {
      # Needs to distinguish between mle and non-mle because of how make.surv saves the S list
      options(digits=5,nsmall=2)
      pts <- list()
      for (i in 1:nmodels) {
        if (method[i]=="mle") {
          pts[[i]] = lapply(1:length(newdata),function(j) {
            tmp=matrix(unlist(res[[i]]$S[[j]]),ncol=4)
            cbind(tmp[,1],tmp[,2])
          })
        } else {
          pts[[i]] = lapply(1:length(newdata),function(j) {
            res[[i]]$S[[1]][[j]]
          })
        }
      }
      colors = 1:nmodels
      leg.txt = character()
      for (i in 1:nmodels) {
        for (j in 1:length(newdata)) {
          points(pts[[i]][[j]],t="l",col=colors[i],lty=j)
          leg.txt[j] = paste0(names(newdata[[j]]),"=",prettyNum(newdata[[j]],format="fg"),collapse=", ")
        }
      }
      legend("topright",legend=leg.txt,bty="n",lty=1:length(newdata),cex=cex.lab)
    }
    if(is.null(newdata)) {
      # With no newdata this works!
      for (i in 1:nmodels) {
        pts <- lapply(res[[i]]$S[[1]],function(m) cbind(m[,1],m[,2]))
        lapply(1:length(pts), function(x) points(pts[[x]],t="l",col=colors[i],lty=x))
      }
      legend(x="topright",legend=labs,lwd=2,bty="n",col=col,cex=cex.lab)
    }
}


# model.checking <- function(x,mod=1,...) {
#   # x = a survHE object with the results of the call to the fit.models function
#   # mod = the model to analyse
#   # ... = additional options
#   
#   exArgs = list(...)
#   if (x$method=="mle") {
#     stop("Model checking is available only for Bayesian models")
#   }
#   if (x$method=="inla") {
#     
#   }
#   if (x$method=="hmc") {
#     rstan::traceplot(x$models[[mod]])
#     rstan::stan_ac(x$models[[mod]])
#   }
# }


model.fit.plot <- function(...,type="aic") {
    ## Plots a summary of the model fit for all the models 
    ## Can also combine several survHE objects each containing the fit for one model
    
  exArgs <- list(...) 		# Lists all the additional inputs
  nexArgs <- length(exArgs)
  classes <- unlist(lapply(1:nexArgs,function(i) class(exArgs[[i]])))
  w=which(classes=="survHE")
  original.method=unlist(lapply(w,function(i) exArgs[[i]]$method))
  if(length(w)==0) {
    stop("You need to input at least one 'survHE' object to run this function!")
  }
  if(length(w)==1) {
    totmodels <- unlist(lapply(w,function(i) length(exArgs[[i]]$models)))
    mods <- exArgs[[w]]$models
    method <- rep(exArgs[[w]]$method,totmodels) 
    aic <- unlist(exArgs[[w]]$model.fitting$aic)
    bic <- unlist(exArgs[[w]]$model.fitting$bic)
    dic <- unlist(exArgs[[w]]$model.fitting$dic)
    if(totmodels>1){
      if (!is.null(exArgs$mod)) {which.model <- exArgs$mod} else {which.model <- 1:length(mods)}
      mods <- lapply(which.model,function(i) mods[[i]])
      method <- method[which.model]
      aic <- aic[which.model]
      bic <- bic[which.model]
      dic <- dic[which.model]
    } 
  }
  if (length(w)>1) {
    mods <- unlist(lapply(w,function(i) exArgs[[i]]$models),recursive = F)
    totmodels <- unlist(lapply(w,function(i) length(exArgs[[i]]$models)))
    method <- unlist(lapply(w,function(i) rep(exArgs[[i]]$method,totmodels[i])))
    aic <- unlist(lapply(w,function(i) exArgs[[i]]$model.fitting$aic))
    bic <- unlist(lapply(w,function(i) exArgs[[i]]$model.fitting$bic))
    dic <- unlist(lapply(w,function(i) exArgs[[i]]$model.fitting$dic))
    if (!is.null(exArgs$mod)) {which.model <- exArgs$mod} else {which.model <- 1:length(mods)}
    mods <- lapply(which.model,function(i) mods[[i]])
    method <- method[which.model]
    aic <- aic[which.model]
    bic <- bic[which.model]
    dic <- dic[which.model]
  }
  model.fitting <- list(aic=aic,bic=bic,dic=dic)
  fit <- list(); fit$models <- mods; class(fit) <- "survHE"
  fit$model.fitting <- model.fitting
  ## Needs to include in the misc object the element vars (which is used for HMC models)
  if (any(method=="hmc")) {
    fit$misc=exArgs[[min(which(original.method=="hmc"))]]$misc
  } else {
    # If none of the survHE objects are HMC, then just use the first
    fit$misc=exArgs[[1]]$misc
  }
  
    if (is.null(exArgs$models)) {
        models <- unlist(lapply(1:length(fit$models),function(i) {
            if(class(fit$models[[i]])=="stanfit") {tolower(fit$models[[i]]@model_name)} else {fit$models[[i]]$dlist$name}
        }))
        models[models %in% c("weibull.quiet","weibull","weibullaf","weibullph")] = "Weibull"
        models[models %in% c("exp","exponential")] = "Exponential"
        models[models %in% "gamma"] = "Gamma"
        models[models %in% c("lnorm","lognormal")] = "log-Normal"
        models[models %in% c("llogis","loglogistic","loglogis")] = "log-Logistic"
        models[models %in% "gengamma"] = "Gen. Gamma"
        models[models %in% "genf"] = "Gen. F"
    } else {
        models=exArgs$models 
    }
    
    # Defines the data to be plotted
    if (type=="aic" | type=="AIC" | type=="a" | type=="A") {
        mf <- data.frame(model=models,AIC=fit$model.fitting$aic)
        lab.type <- "AIC"
    } else if (type=="bic" | type=="BIC" | type=="b" | type=="B") {
        mf <- data.frame(model=models,BIC=fit$model.fitting$bic)
        lab.type <- "BIC"
    } else if (type=="dic" | type=="DIC" | type=="d" | type=="D") {
        mf <- data.frame(model=models,DIC=fit$model.fitting$dic)
        lab.type <- "DIC"
    }
    
    # Finally do the plot
    if (is.null(exArgs$xlim)) {xlm=range(pretty(mf[,2]))} else {xlm=exArgs$xlim}
    if (is.null(exArgs$digits)) {digits=7} else {digits=exArgs$digits}
    if (is.null(exArgs$nsmall)) {nsmall=3} else {nsmall=exArgs$nsmall}
    if (is.null(exArgs$main)) {main=paste0("Model comparison based on ",lab.type)} else {main=exArgs$main}
    if (is.null(exArgs$mar)) {mar=c(4,6,3,1.3)} else {mar=exArgs$mar}
    if (is.null(exArgs$cex.names)) {cex.names=0.8} else {cex.names=exArgs$cex.names}
    par(mar=mar)                                           # Bottom,left,top & right margins
    b <- barplot(                                          # Function to draw a barplot (see BMS NICE submission)
        mf[,2],  	                                         # Makes a barplot using the values of the AIC or BIC
        names.arg=mf$model,	                               # Names of the models (can be formatted differently)
        xlab=lab.type,                                     # Label for the x-axis
        xlim=xlm,
        xpd=F,                                             # Graphical parameter to clip at the lowest end of the range
        horiz=T,                                           # Plots the graph horizontally (better readability)
        las=1,                                             # Rotates the labels on the y-axis (better readability)
        cex.names=cex.names,                               # Rescales the labels on the y-axis to 80% of normal size
        main=main
    )
    # And then adds the actual value of the AIC/BIC for each of the models
    text(mf[,2],		                                       # Position of the text on the x-axis
         b,                                                # Position of the text on the y-axis
         format(mf[,2],digits=digits,nsmall=nsmall),       # Formats the values of the AICs/BICs/DICs, using 3 dp
         pos=4,                                            # Puts the text to the right of the bars
         cex=.8                                            # Rescales the labels on the y-axis to 80% of normal size
    )
}


test.linear.assumptions <- function(fit,mod=1,coxph=TRUE,label = FALSE,...){
  ## THIS IS INTERESTING, BUT NEEDS TO COMPLETE WITH THE OTHER DISTRIBUTIONS!!!
  exArgs <- list(...)
  
  m <- fit$models[[mod]]
  dist <- ifelse(fit$method=="hmc",tolower(m@model_name),m$dlist$name)
  split_vector <- c(1)
  for(i in 2:length(fit$misc$km$time)){
    if(fit$misc$km$time[i]<fit$misc$km$time[i-1]){
      split_vector <- c(split_vector,i-1,i)
    }
  }
  split_vector <- c(split_vector,length(fit$misc$km$time))
  split_mat <- matrix(split_vector,length(split_vector)/2,2,byrow = T)
  times <- lapply(1:dim(split_mat)[1],function(x) fit$misc$km$time[split_mat[x,][1]:split_mat[x,][2]])
  survs <- lapply(1:dim(split_mat)[1],function(x) fit$misc$km$surv[split_mat[x,][1]:split_mat[x,][2]])
  if (dist %in% c("Exponential","exp","exponential")){
    plot(0,0,col="white",xlab='time',ylab='log(S(t))',axes=F,xlim=range(pretty(fit$misc$km$time)))
    axis(1)
    axis(2)
    pts <- lapply(1:dim(split_mat)[1],function(m) cbind(times[[m]],log(survs[[m]])))
    lapply(1:length(pts), function(x) points(pts[[x]],t="l",lty=x))
    if (label){legend('topright','Exponential distributional assumption',bty='n')}
    
    #text(max(pts[[1]][,1]),max(pts[[1]][,2]), 'Exponential linear assumption', cex=0.6, col="red")
    
  }
  if (dist %in% c("weibull","weibullPH","weibull.quiet","weibullaf","weibullph")){
    plot(0,0,col="white",xlab='log(time)',ylab='log(-log(S(t))) = log cumulative hazard',
         axes=F,xlim=range(pretty(log(fit$misc$km$time))), 
         ylim =range(pretty(log(-log(survs[[1]])))))
    axis(1)
    axis(2)
    pts <- lapply(1:dim(split_mat)[1],function(m) cbind(log(times[[m]]),log(-log(survs[[m]]))))
    lapply(1:length(pts), function(x) points(pts[[x]],t="l",lty=x))
    if (label){legend('topright','Weibull distributional assumption',bty='n')}
  }
  if (dist %in% c("llogis","loglogistic")){
    plot(0,0,col="white",xlab='time',ylab='log(S(t)/(1-S(t)))',axes=F,xlim=range(pretty(log(fit$misc$km$time))), ylim =range(pretty(log(survs[[1]]/(1-survs[[1]])))))
    axis(1)
    axis(2)
    pts <- lapply(1:dim(split_mat)[1],function(m) cbind(log(times[[m]]),log(survs[[m]]/(1-survs[[m]]))))
    lapply(1:length(pts), function(x) points(pts[[x]],t="l",lty=x))
    if (label){legend('topright','log-Logistic distributional assumption',bty='n')}
  }
  if (dist %in% c("lognormal","lnorm")){
    ### add log normal 
    plot(0,0,col="white",xlab='time',ylab='log(S(t))',axes=F,xlim=range(pretty(log(fit$misc$km$time)))   )#, ylim =range(qnorm(1-survs[[1]])))
    axis(1)
    axis(2)
    pts <- lapply(1:dim(split_mat)[1],function(m) cbind(times[[m]],qnorm(1-survs[[m]])))
    lapply(1:length(pts), function(x) points(pts[[2]],t="l",lty=x))
    if (label){legend('topright','lognormal distributional assumption',bty='n')}
    
  }
  if (dist == "gompertz"){
    estimate.h <- function(s,t){
      denom <- t-c(t[-1],max(t)+1)
      print(denom)
      numerator <- log(s) - log(c(s[-1],0))
      print(numerator)
      return(-numerator/denom)
    }
    plot(0,0,col="white",xlab='log(time)',ylab='h(t)',axes=F,xlim=range(pretty(fit$misc$km$time)), ylim =range(pretty(estimate.h(survs[[1]],times[[1]]))))
    axis(1)
    axis(2)
    ### NEED TO CHECK --- WHAT IS V2???
    pts <- lapply(1:dim(split_mat)[1],function(m) data.table::data.table(cbind(times[[m]],estimate.h(survs[[m]],times[[m]])))[V2!=0,])
    lapply(1:length(pts), function(x) points(pts[[x]],t="l",lty=x))
    if (label){legend('topright','Gompertz distributional assumption',bty='n')}
    
  }
}


#### THIS IS IMPORTANT --- NEED TO DECIDE WHAT AND HOW WE WANT TO EXPORT THE PSA RESULTS!!!
write.surv <- function(object,file,sheet=NULL,what="surv") {
  # Writes the survival summary to an excel file (helpful to then call the values in the Markov model)
  # object = a summary.flexsurvreg object containing the survival curves (with times, estimates and interval limits)
  # file = a string with the full path to the file name to be saved
  # sheet = a string with the name of the sheet to be created
  # what = the object to be exported. Possible values are:
  #        'surv' = a matrix with nsim rows and ntimes columns with the survival curve (one such matix for each configuration of the covariates)
  #        'sim' = a matrix with nsim rows and simulations for the survival parameters (scale, shape, rate, etc)
  
  # If xlsx is not installed, then request installation
  if(!isTRUE(requireNamespace("xlsx",quietly=TRUE))) {
    stop("You need to install the R package 'xlsx'. Please run in your R terminal:\n install.packages('xlsx')")
  }
  # But if it is installed, check if it's loaded and if not make its Namespace available so that all its functions are available
  if (isTRUE(requireNamespace("xlsx",quietly=TRUE))) {
    if (!is.element("xlsx", (.packages()))) {
      attachNamespace("xlsx")
    }
    # Extracts the relevant component of the make.surv output
    if(what=="surv") {
      export <- object$mat
      export.lab <- paste0(object$nsim," simulation(s) for the survival curve:")
      nobjs <- length(export)
      profile.lab <- unlist(lapply(1:nobjs,function(i) paste0("[[",i,"]] = ",rownames(object$des.mat)[i],"\n")))
      dims <- dim(export[[1]])
      # Finds the total number of rows necessary to write the simulations to the output file
      tot.rows <- dims[1]*nobjs + nobjs
      cn <- as.character(object$S[[1]][[1]][,1])
      for (i in 1:nobjs) {
        colnames(export[[i]])=paste0("t_",cn)
      }
    }
    
    # # Gets the extension of the file --- then decides whether to do write.csv or xlsx (do we want more formats??)
    # exts <- tools::file_ext(file)
    # if(exts=="csv") {
    #   out=export[[1]]
    #   if (nobjs>1) {
    #     for (i in 2:nobjs) {
    #       out = rbind(out,rep(NA,ncol(out)),export[[i]])
    #     }
    #   }
    #   write.csv(out,file=file)
    # }
    
    if(is.null(sheet)) {sheet="Sheet 1"}
    
    # If it already exists, we need to append the data to a different sheet
    if (file.exists(file)) {
      wb <- xlsx::loadWorkbook(file)
      # If worksheet already exists needs to replace it & overwrite it
      if (sheet %in% names(xlsx::getSheets(wb))) {xlsx::removeSheet(wb,sheetName=sheet)}
      sheet <- xlsx::createSheet(wb,sheet)
      sr <- seq(from=1,by=(dims[1]+2),to=tot.rows)
      ex <- lapply(1:nobjs,function(i) xlsx::addDataFrame(export[[i]],sheet=sheet,startRow=sr[i],startColumn=1,row.names=T,col.names=T))
      # sc <- seq(from=1,by=5,length.out=nobjs)
      # for (i in 1:nobjs) {
      #   xlsx::addDataFrame(export[[i]],sheet=sheet,startRow=1,startColumn=sc[i],row.names=F)
      # }
      xlsx::saveWorkbook(wb,file)
    }
    
    # But if file does not exist, then create it
    if (!file.exists(file)) {
      exts <- tools::file_ext(file)
      ## Should put some restriction as to what file extensions we want?
      wb <- xlsx::createWorkbook(type=exts)
      sheet <- xlsx::createSheet(wb,sheet)
      sr <- seq(from=1,by=(dims[1]+2),to=tot.rows)
      ex <- lapply(1:nobjs,function(i) xlsx::addDataFrame(export[[i]],sheet=sheet,startRow=sr[i],startColumn=1,row.names=T,col.names=T))
      # sc <- seq(from=1,by=5,length.out=nobjs)
      # for (i in 1:nobjs) {
      #   xlsx::addDataFrame(object[[i]],sheet=sheet,startRow=1,startColumn=sc[i],row.names=F)
      # }
      xlsx::saveWorkbook(wb,file)
    }
    
    msg <- paste0("Written to file: ",file)
    cat(export.lab,"\n",profile.lab,"\n",msg)
  }
}

digitise <- function(surv_inp,nrisk_inp,km_output="KMdata.txt",ipd_output="IPDdata.txt") {
  # Post-process the data obtained by DigitizeIT to obtain the KM data and the individual level data
  # surv_inp = a txt file obtained by DigitizeIT and containing the input survival times from graph reading
  # nrisk_inp = a txt file obtained by DigitizeIT and containing the reported number at risk
  # km_output = the name of the file to which the KM data will be written
  # ipd_output = the name of the file to which the individual level data data will be written
  # Adapted from Patricia Guyot (2012)
  
  # Defines the working directory (same as the one where the DigitizeIT data are)
  working.dir <- dirname(surv_inp)
  ####  setwd(working.dir); working.dir <- paste0(getwd(),"/")
  tot.events<-"NA"  #tot.events = total no. of events reported. If not reported, then tot.events="NA"
  arm.id<-1         #arm indicator
  
  #Read in survival times read by digizeit
  digizeit <- read.table(surv_inp,header=TRUE,row.names=NULL)
  t.S<-digizeit[,2]     # times recorded from DigitizeIT
  S<-digizeit[,3]       # survival from DigitizeIT
  
  #Read in published numbers at risk, n.risk, at time, t.risk, lower and upper indexes for time interval
  pub.risk<-read.table(nrisk_inp,header=TRUE,row.names=NULL)
  ## Needs to get rid of possible time intervals with no digitised observations
  pub.risk <- pub.risk[pub.risk[,4]>0,]
  ## Needs to recode the first ever occurrence to 1??
  if (!(pub.risk[1,3]==1)) {pub.risk[1,3] <- 1}
  
  # Defines the variables needed for the algorithm
  t.risk<-pub.risk[,2]
  lower<-pub.risk[,3]
  upper<-pub.risk[,4]
  n.risk<-pub.risk[,5]
  n.int<-length(n.risk)
  n.t<- upper[n.int]
  
  #Initialise vectors
  arm <- rep(arm.id,n.risk[1])
  n.censor <- rep(0,(n.int-1))
  n.hat <- rep(n.risk[1]+1,n.t)
  cen <- d <- rep(0,n.t)
  KM.hat <- rep(1,n.t)
  last.i <- rep(1,n.int)
  sumdL <- 0
  
  # Executes Patricia's algorithm to determine censoring
  if (n.int > 1){
    #Time intervals 1,...,(n.int-1)
    for (i in 1:(n.int-1)){
      #First approximation of no. censored on interval i
      n.censor[i]<- round(n.risk[i]*S[lower[i+1]]/S[lower[i]]- n.risk[i+1])
      #Adjust tot. no. censored until n.hat = n.risk at start of interval (i+1)
      while((n.hat[lower[i+1]]>n.risk[i+1])||((n.hat[lower[i+1]]<n.risk[i+1])&&(n.censor[i]>0))){
        if (n.censor[i]<=0){
          cen[lower[i]:upper[i]]<-0
          n.censor[i]<-0
        }
        if (n.censor[i]>0){
          cen.t<-rep(0,n.censor[i])
          for (j in 1:n.censor[i]){
            cen.t[j]<- t.S[lower[i]] +
              j*(t.S[lower[(i+1)]]-t.S[lower[i]])/(n.censor[i]+1)
          }
          #Distribute censored observations evenly over time. Find no. censored on each time interval.
          cen[lower[i]:upper[i]]<-hist(cen.t,breaks=t.S[lower[i]:lower[(i+1)]],plot=F)$counts
        }
        #Find no. events and no. at risk on each interval to agree with K-M estimates read from curves
        n.hat[lower[i]]<-n.risk[i]
        last<-last.i[i]
        for (k in lower[i]:upper[i]){
          if (i==1 & k==lower[i]){
            d[k]<-0
            KM.hat[k]<-1
          }
          else {
            d[k]<-round(n.hat[k]*(1-(S[k]/KM.hat[last])))
            KM.hat[k]<-KM.hat[last]*(1-(d[k]/n.hat[k]))
          }
          n.hat[k+1]<-n.hat[k]-d[k]-cen[k]
          if (d[k] != 0) last<-k
        }
        n.censor[i]<- n.censor[i]+(n.hat[lower[i+1]]-n.risk[i+1])
      }
      if (n.hat[lower[i+1]]<n.risk[i+1]) n.risk[i+1]<-n.hat[lower[i+1]]
      last.i[(i+1)]<-last
    }
  }
  #Time interval n.int.
  if (n.int>1){
    #Assume same censor rate as average over previous time intervals.
    n.censor[n.int]<- min(round(sum(n.censor[1:(n.int-1)])*(t.S[upper[n.int]]-
                                                              t.S[lower[n.int]])/(t.S[upper[(n.int-1)]]-t.S[lower[1]])), n.risk[n.int])
  }
  if (n.int==1){n.censor[n.int]<-0}
  if (n.censor[n.int] <= 0){
    cen[lower[n.int]:(upper[n.int]-1)]<-0
    n.censor[n.int]<-0
  }
  if (n.censor[n.int]>0){
    cen.t<-rep(0,n.censor[n.int])
    for (j in 1:n.censor[n.int]){
      cen.t[j]<- t.S[lower[n.int]] +
        j*(t.S[upper[n.int]]-t.S[lower[n.int]])/(n.censor[n.int]+1)
    }
    cen[lower[n.int]:(upper[n.int]-1)]<-hist(cen.t,breaks=t.S[lower[n.int]:upper[n.int]],plot=F)$counts
  }
  #Find no. events and no. at risk on each interval to agree with K-M estimates read from curves
  n.hat[lower[n.int]]<-n.risk[n.int]
  last<-last.i[n.int]
  for (k in lower[n.int]:upper[n.int]){
    if(KM.hat[last] !=0){
      d[k]<-round(n.hat[k]*(1-(S[k]/KM.hat[last])))} else {d[k]<-0}
    KM.hat[k]<-KM.hat[last]*(1-(d[k]/n.hat[k]))
    n.hat[k+1]<-n.hat[k]-d[k]-cen[k]
    #No. at risk cannot be negative
    if (n.hat[k+1] < 0) {
      n.hat[k+1]<-0
      cen[k]<-n.hat[k] - d[k]
    }
    if (d[k] != 0) last<-k
  }
  #If total no. of events reported, adjust no. censored so that total no. of events agrees.
  if (tot.events != "NA"){
    if (n.int>1){
      sumdL<-sum(d[1:upper[(n.int-1)]])
      #If total no. events already too big, then set events and censoring = 0 on all further time intervals
      if (sumdL >= tot.events){
        d[lower[n.int]:upper[n.int]]<- rep(0,(upper[n.int]-lower[n.int]+1))
        cen[lower[n.int]:(upper[n.int]-1)]<- rep(0,(upper[n.int]-lower[n.int]))
        n.hat[(lower[n.int]+1):(upper[n.int]+1)]<- rep(n.risk[n.int],(upper[n.int]+1-lower[n.int]))
      }
    }
    #Otherwise adjust no. censored to give correct total no. events
    if ((sumdL < tot.events)|| (n.int==1)){
      sumd<-sum(d[1:upper[n.int]])
      while ((sumd > tot.events)||((sumd< tot.events)&&(n.censor[n.int]>0))){
        n.censor[n.int]<- n.censor[n.int] + (sumd - tot.events)
        if (n.censor[n.int]<=0){
          cen[lower[n.int]:(upper[n.int]-1)]<-0
          n.censor[n.int]<-0
        }
        if (n.censor[n.int]>0){
          cen.t<-rep(0,n.censor[n.int])
          for (j in 1:n.censor[n.int]){
            cen.t[j]<- t.S[lower[n.int]] +
              j*(t.S[upper[n.int]]-t.S[lower[n.int]])/(n.censor[n.int]+1)
          }
          cen[lower[n.int]:(upper[n.int]-1)]<-hist(cen.t,breaks=t.S[lower[n.int]:upper[n.int]],plot=F)$counts
        }
        n.hat[lower[n.int]]<-n.risk[n.int]
        last<-last.i[n.int]
        for (k in lower[n.int]:upper[n.int]){
          d[k]<-round(n.hat[k]*(1-(S[k]/KM.hat[last])))
          KM.hat[k]<-KM.hat[last]*(1-(d[k]/n.hat[k]))
          if (k != upper[n.int]){
            n.hat[k+1]<-n.hat[k]-d[k]-cen[k]
            #No. at risk cannot be negative
            if (n.hat[k+1] < 0) {
              n.hat[k+1]<-0
              cen[k]<-n.hat[k] - d[k]
            }
          }
          if (d[k] != 0) last<-k
        }
        sumd<- sum(d[1:upper[n.int]])
      }
    }
  }
  
  # Now writes the results to the output files
  KMdata <- data.frame(time=t.S,n.risk=n.hat[1:n.t],n.event=d,n.censored=cen)
  write.table(KMdata,km_output,sep="\t",row.names=FALSE,col.names=TRUE)
  
  # And forms IPD data
  #Initialise vectors
  t.IPD <- rep(t.S[n.t],n.risk[1])
  event.IPD <- rep(0,n.risk[1])
  #Write event time and event indicator (=1) for each event, as separate row in t.IPD and event.IPD
  k=1
  for (j in 1:n.t){
    if(d[j]!=0){
      t.IPD[k:(k+d[j]-1)]<- rep(t.S[j],d[j])
      event.IPD[k:(k+d[j]-1)]<- rep(1,d[j])
      k<-k+d[j]
    }
  }
  #Write censor time and event indicator (=0) for each censor, as separate row in t.IPD and event.IPD
  for (j in 1:(n.t-1)){
    if(cen[j]!=0){
      t.IPD[k:(k+cen[j]-1)]<- rep(((t.S[j]+t.S[j+1])/2),cen[j])
      event.IPD[k:(k+cen[j]-1)]<- rep(0,cen[j])
      k<-k+cen[j]
    }
  }
  #Output IPD
  IPD <- data.frame(time=t.IPD,event=event.IPD,arm)
  write.table(IPD,ipd_output,sep="\t",row.names=FALSE,col.names=TRUE)
  
  if (dirname(km_output)==".") {
    cat("\n")
    cat(paste0("Kaplan Meier data written to file: ",working.dir,km_output))    
  } else {
    cat("\n")
    cat(paste0("Kaplan Meier data written to file: ",km_output))    
  }
  if (dirname(ipd_output)==".") {
    cat("\n")
    cat(paste0("IPD data written to file: ",working.dir,ipd_output))    
  } else {
    cat("\n")
    cat(paste0("IPD data written to file: ",ipd_output))
  }
}



make.ipd <- function(ipd_files,ctr=1,var.labs=c("time","event","arm")) {
  ## Piles in the simulated IPD resulting from running digitise for more than one treatment arm  
  ## ipd_files = a list including the names of the IPD files created as output of digitise
  ## ctr = the index of the file associated with the control arm (default, the first file).
  ##       This will be coded as 0
  ## var.labs = a vector of labels for the column of the resulting data matrix. NB these
  ##            should match the arguments to the formula specified for fit.models. The
  ##            user can specify values. These should be 3 elements (TIME, EVENT, ARM)
  
  # Identifies the number of arms (= number of IPD files)
  n_arms <- length(ipd_files)
  index <- 1:n_arms
  active <- index[-ctr]
  data <- read.table(ipd_files[[ctr]],header=TRUE,row.names=NULL)
  data[,"arm"] <- 0 # sets the value of "arm" to 0, for the control group
  arm.ind <- 1
  for (i in active) {
    tmp <- read.table(ipd_files[[index[i]]],header=TRUE,row.names=NULL)
    tmp[,"arm"] <- arm.ind
    data <- rbind(data,tmp)
    arm.ind <- arm.ind+1
  }
  colnames(data) <- var.labs
  return(data)
}



make.transition.probs <- function(x,...) {
  # Computes the transition probabilities (to feed a discrete-time Markov model), based on the output of make.surv
  # Uses the formula p(t)=1-S(t+k)/S(t) where k is the MM cycle length and t is a generic time
  # x = an object obtained as output of the call to make.surv
  # ... = additional arguments. 
  #       labs = a string vector of names for the elements of the list (strata for the survival analysis)

  exArgs=list(...)
  
  tp = lapply(1:length(x$mat),function(i) {
    matrix(unlist(lapply(1:nrow(x$mat[[1]]),function(j) {
      1-(x$mat[[i]][j,2:ncol(x$mat[[1]])]/x$mat[[i]][j,1:ncol(x$mat[[1]])-1])
    })),byrow=T,nrow=nrow(x$mat[[1]]))
  })
  
  if (exists("labs",where=exArgs)) {
    if(length(exArgs$labs)==length(x$mat)) {
      names(tp) = exArgs$labs
    }
  } else {
    names(tp) = rownames(x$des.mat)
  }
  
  # Defines the column labels (to specify what the times refer to)
  col.labs = character()
  for (i in 1:length(x$times)-1) {
    col.labs[i] = paste0(x$times[i],"-",x$times[i+1])
  }
  for (i in 1:length(tp)) {
    colnames(tp[[i]]) = col.labs
  }

  # Output
  # tpa: list of matrices with nsim rows and length(time) columns with the simulations for the transition probabilities between
  #      consecutive times.
  return(tp)
}


poly.weibull = function(formula=NULL,data,...) {
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
  exArgs = list(...)
  
  # Sets up defaults
  
  # Check whether the user has specified a list of formulae to be used in each compoenent of the model
  if (!is.null(formula)) {
    # a. The user has specified a formula, but has only given 1 element, then expand it
    #    so that formula becomes a list of formulae and sets up 2 components by default
    if (class(formula)=="formula") {
      formula = list(formula)
    }
    M = length(formula) 
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
  if(exists("control",where=exArgs)) {control=exArgs$control} else {control=NULL}
  if(exists("seed",where=exArgs)) {seed <- exArgs$seed} else {seed <- sample.int(.Machine$integer.max, 1)}
  if(exists("pars",where=exArgs)) {pars <- exArgs$pars} else {
    pars <- c("loglambda","lambda","lp__")
  }
  if(exists("include",where=exArgs)) {include <- exArgs$include} else {include <- FALSE}
  if(exists("cores",where=exArgs)) {cores=exArgs$cores} else {cores=1}

  # Reconstructs the vars list based on the formula
  vars = list()
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
    vars[[i]] = list(time=time,event=event,factors=factors,covs=covs,nlevs=D)
  }
  
  # Loads the pre-compiled models
  dso <- stanmodels
  
  ###############################################################################################
  ### THIS IS JUST A TEMPORARY LINE (UNTIL THE CORRECT MODELS ARE PRE-COMPILED!)
  ####dso = readRDS("~/Dropbox/UCL/Mapi/Projects/Survival/Stan_code/DSOs.rds")
  ###############################################################################################
  time2run = numeric()
  # Selects the precompiled polyweibull model (CHECK IF THE ORDER IN availables.hmc CHANGES!!)
  dso <- dso[[6]] 
  
  data.stan = list(t=data[,vars[[1]]$time], d=data[,vars[[1]]$event]); data.stan$n = length(data.stan$t); 
  data.stan$M = M;
  X = lapply(1:data.stan$M,function(i) model.matrix(formula[[i]],data))
  # max number of covariates in all the model formulae
  data.stan$H = max(unlist(lapply(1:data.stan$M,function(i) ncol(X[[i]]))))
  # NB: Stan doesn't allow matrices with 1 column, so if there's only one covariate (eg intercept only), needs a little trick
  if (data.stan$H==1) {data.stan$H=data.stan$H+1}
  X = lapply(1:data.stan$M,function(i) {
    if(ncol(X[[i]])<data.stan$H) {
      X[[i]] = cbind(X[[i]],matrix(0,nrow=nrow(X[[i]]),ncol=(data.stan$H-ncol(X[[i]]))))
    } else {
      X[[i]] = X[[i]]
    }
  })
  data.stan$X=array(NA,c(data.stan$M,data.stan$n,data.stan$H))
  for (m in 1:data.stan$M) {
    data.stan$X[m,,] = X[[m]]
  }
  data.stan$mu_beta=matrix(0,nrow=data.stan$H,ncol=data.stan$M); data.stan$sigma_beta=matrix(10,data.stan$H,data.stan$M)
  
  # These are modified if the user gives values in the call to poly.weibull
  if(exists("priors",where=exArgs)) {
    priors=exArgs$priors
    # Linear predictor coefficients
    if(!is.null(priors$mu_beta)) {
      data.stan$mu_beta=priors$mu_beta
    }
    if(!is.null(priors$sigma_beta)) {
      data.stan$sigma_beta = priors$sigma_beta
    }
  }

  # Now runs Stan to sample from the posterior distributions
  tic = proc.time()
  out=rstan::sampling(dso,data.stan,chains=chains,iter=iter,warmup=warmup,thin=thin,seed=seed,control=control,
                      pars=pars,include=include,cores=cores)
  toc = proc.time()-tic
  time2run=toc[3]
  list(out=out,data.stan=data.stan,time2run=time2run)
  
  if(exists("save.stan",where=exArgs)) {save.stan <- exArgs$save.stan} else {save.stan=FALSE}
  time_survHE = time2run
  time_stan <- sum(rstan::get_elapsed_time(out))
  time2run = min(time_survHE,time_stan)
  names(time2run) <- "PolyWeibull"
    
  # Computes the log-likelihood 
  beta = rstan::extract(out)$beta
  alpha = rstan::extract(out)$alpha
  # NB: To make more robust estimate of AIC, BIC and DIC uses the median here (instead of the mean)
  #     This is likely to be a better approximation to the MLE when the underlying distributions are highly skewed!
  beta.hat = apply(beta,c(2,3),median)
  alpha.hat = apply(alpha,2,median)
  linpred = lapply(1:data.stan$M,function(m) {
    beta[,m,]%*%t(data.stan$X[m,,])
  })
  linpred.hat = lapply(1:data.stan$M,function(m) {
    beta.hat[m,]%*%t(data.stan$X[m,,])
  })
  
  h=log_s=array(NA,c(nrow(alpha),data.stan$n,data.stan$M))
  h_bar = log_s_bar = matrix(NA,data.stan$n,data.stan$M)
  for (m in 1:data.stan$M) {
    h_bar[,m]=alpha.hat[m]*exp(linpred.hat[[m]])*data.stan$t^(alpha.hat[m]-1)
    log_s_bar[,m] = exp(linpred.hat[[m]])*data.stan$t^(alpha.hat[m])
    for (i in 1:nrow(linpred[[m]])) {
      h[i,,m] = alpha[i,m]*exp(linpred[[m]][i,])*data.stan$t^(alpha[i,m]-1)
      log_s[i,,m] = exp(linpred[[m]][i,])*data.stan$t^(alpha[i,m])
    }
  }
  d_log_sum_h = matrix(NA,nrow(alpha),data.stan$n)
  for (i in 1:nrow(alpha)) {
    d_log_sum_h[i,] = data.stan$d * log(rowSums(h[i,,]))
  }
  loglik.bar = sum(data.stan$d*log(rowSums(h_bar))-rowSums(log_s_bar))
  loglik = rowSums(d_log_sum_h) - rowSums(log_s,2)
  D.theta=-2*loglik
  D.bar=-2*loglik.bar
  pD = mean(D.theta) - D.bar
  pV = .5*var(D.theta)  # Uses Gelman's definition of pD!
  dic = mean(D.theta)+pD
  # Approximates AIC & BIC using the mean deviance and the number of nominal parameters
  npars = data.stan$M + sum(unlist(lapply(1:data.stan$M,function(m) {sum(1-apply(data.stan$X[m,,],2,function(x) all(x==0)))})))
  aic = D.bar+2*npars
  bic = D.bar+npars*log(data.stan$n)
  
  # Now defines the outputs of the function
  model.fitting <- list(aic=aic,bic=bic,dic=dic)
  km = list(time=data.stan$t)
  misc <- list(time2run=time2run,formula=formula,data=data,km=km)
  misc$vars <- vars; misc$data.stan=data.stan
  # If save.stan is set to TRUE, then saves the Stan model file(s) & data
  if(save.stan==TRUE) {
	model_code = attr(mod$out@stanmodel,"model_code")
    con = "PolyWeibull.stan"
    writeLines(model_code,con=con)
	txt = paste0("Model code saved to the file: ",con,"\n")
    cat(txt)
  }
  mod = list(out)
  names(mod)="PolyWeibull"
  # Finally prepares the output object
  res <- list(models=mod,model.fitting=model.fitting,method="hmc",misc=misc)
  # And sets its class attribute to "survHE"
  class(res) <- "survHE"
  return(res)
}



##
summary.survHE = function(object,mod=1,t=NULL,nsim=1000,...) {
  # Computes the estimated mean survival as the area under the survival curve
  # This is obtained using the trapezoidal method by taking the average of the "left" and "right" y-values.
  # object: is the output from a fit.models call
  # mod: the model to be analysed (default = 1)
  # t: the vector of times to be used in the computation. Default = NULL, which means the observed times will be used.
  #     NB: the vector of times should be: i) long enough so that S(t) goes to 0; and ii) dense enough so that
  #         the approximation to the AUC is sufficiently precise
  # nsim: number of simulations from the survival curve distributions to be used (to compute interval estimates)
  # stats: a logical value. If TRUE, also shows a table 
  # ...: optional arguments
  # newdata = a list (of lists), specifiying the values of the covariates at which the computation is performed. For example
  #           'list(list(arm=0),list(arm=1))' will create two survival curves, one obtained by setting the covariate 'arm'
  #           to the value 0 and the other by setting it to the value 1. In line with 'flexsurv' notation, the user needs
  #           to either specify the value for *all* the covariates or for none (in which case, 'newdata=NULL', which is the
  #           default). If some value is specified and at least one of the covariates is continuous, then a single survival
  #           curve will be computed in correspondence of the average values of all the covariates (including the factors, 
  #           which in this case are expanded into indicators). The order of the variables in the list *must* be the same
  #           as in the formula used for the model
  # labs: a vector of strings giving the names of the "profile" of covariates for which the mean survival times are computed
  #
  # NB: NEED TO FIX THIS FOR THE POLY-WEIBULL
  #
  # Defines the utility function to compute the stats table
  make.stats = function(x, dim = 2) {
    bugs.stats = function(x) {
      c(mean(x), sd(x), quantile(x, 0.025), median(x), quantile(x, 0.975))
    }
    if (is.null(dim(x)) == TRUE) {
      tab <- bugs.stats(x)
      names(tab) <- c("mean", "sd", "2.5%", "median", "97.5%")
    }
    if (is.null(dim(x)) == FALSE) {
      tab <- t(apply(x, dim, function(x) bugs.stats(x)))
      colnames(tab) <- c("mean", "sd", "2.5%", "median", "97.5%")
    }
    return(tab)
  }
  
  exArgs = list(...)
  if (!exists("newdata",where=exArgs)) {newdata=NULL} else {newdata=exArgs$newdata}
  if (!exists("labs",where=exArgs)) {labs=NULL} else {labs=exArgs$labs}
  if (is.null(t)) {t=object$misc$km$time} else {t=t}
  
  psa = make.surv(object,mod=mod,t=t,nsim=nsim,newdata=newdata)
  rlabs = rownames(psa$des.mat)
  if (!is.null(rlabs)) {
    rlabs=gsub("^1,","",rlabs)
  } else {
    rlabs = ""
  }
  if(!is.null(labs) & length(labs)==length(rlabs)) {rlabs=labs}
  
  mean.surv = matrix(unlist(lapply(1:psa$nsim,function(i) {
    lapply(1:length(psa$S[[1]]),function(j) {
      xvar = psa$S[[i]][[j]][,1]
      yvar = psa$S[[i]][[j]][,2]
      sum(diff(xvar) * (head(yvar,-1)+tail(yvar,-1)), na.rm=T)/2
    })
  })),nrow=psa$nsim,byrow=T)
  if (ncol(mean.surv)==length(names(object$misc$km$strata))) {
    colnames(mean.surv) = names(object$misc$km$strata)
  }
  
  tab = NULL
  if(psa$nsim>1) {
    tab = make.stats(mean.surv)
    rownames(tab) = rlabs
    # if(!is.null(names(x$misc$km$strata))) {
    #   if (ncol(mean.surv)==length(names(x$misc$km$strata))) {
    #     rownames(tab) = names(x$misc$km$strata)
    #   } else {
    #     rownames(tab) = rlabs
    #   }
    # } else {
    #   rownames(tab) = rlabs
    # }
    cat("\nEstimated average survival time distribution* \n")
    print(tab)
    cat(paste0("\n*Computed over the range: [",paste(format(range(t),digits=4,nsmall=3),collapse="-"),"] using ",psa$nsim," simulations.\nNB: Check that the survival curves tend to 0 over this range!\n"))
  }
  invisible(list(mean.surv=mean.surv,tab=tab))
}
