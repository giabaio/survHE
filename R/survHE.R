## SET OF UTILITY FUNCTIONS TO INCLUDE SURVIVAL ANALYSIS RESULTS INTO A HEALTH ECONOMIC MODEL
## Gianluca Baio + Will Browne + Peter Konings (10 Jan 2017)
##
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
  availables.inla <- c("exponential","weibull","weibullPH","lognormal","loglogistic")
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
      message("NB: INLA can only fit Exponential, Weibull, log-Logistic or log-Normal parametric survival models. Falling back on MLE analysis")
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
  ObjSurvfit <- rms::npsurv(        # Uses the function "npsurv" from the package "rms"
    formula = km.formula,          # to fit the model specified in the "formula" object
    data = data                    # to the dataset named "data"
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
          if(exists("bhazard",where=exArgs)) {bhazard <- exArgs$bhazard} else {bhazard <-NULL}
          if(exists("weights",where=exArgs)) {weights <- exArgs$weights} else {weights <- NULL}
          if(exists("subset",where=exArgs)) {subset <- exArgs$subset} else {subset <- NULL}
          if(exists("knots",where=exArgs)) {knots <- exArgs$knots} else {knots <- NULL}
          if(exists("k",where=exArgs)) {k <- exArgs$k} else {k <- 0}
          if(exists("bknots",where=exArgs)) {bknots <- exArgs$bknots} else {bknots <- NULL}
          if(exists("scale",where=exArgs)) {scale <- exArgs$scale} else {scale <- "hazard"}
          if(exists("timescale",where=exArgs)) {timescale <- exArgs$scale} else {timescale <- "log"}
          model <- flexsurv::flexsurvspline(formula=formula,data=data,k=k,knots=knots,bknots=bknots,scale=scale,timescale=timescale)
      } else {
          model <- flexsurv::flexsurvreg(formula=formula,data=data,dist=distr)
      }
      toc <- proc.time()-tic
      time2run <- toc[3]
      list(model=model,time2run=time2run)
    }
    output <- lapply(distr,function(x) runMLE(x))
    mod <- lapply(output, function(x) x$model)
    time2run <- unlist(lapply(output, function(x) x$time2run))
    names(time2run) <- labs
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
        control.fixed$prec <- control.fixed$prec.intercept <- 1/(5^2)
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
	# As of 9 Jan 2017, INLA is creating new distribution names for survival models
	# so needs to update the name
	if(distr[x] %in% c("exponential","lognormal")) {distr[x]=paste0(distr[x],"surv")}
        if(distr[x]=="weibullPH") {
		distr[x]="weibullsurv"
		control.family[[x]]$variant=0
	}
	if(distr[x]=="weibull") {
		distr[x]="weibullsurv"
		control.family[[x]]$variant=1
	}
	####
	## Workaround to load the libraries (needed in LINUX????)
	## INLA:::inla.dynload.workaround()
	####
        INLA::inla(formula,family=distr[x],data=data,control.compute=list(config=TRUE,dic=TRUE),
                   control.inla=list(int.strategy="grid",dz=dz,diff.logdens=diff.logdens),
                   control.fixed=control.fixed,control.family=control.family[[x]]
        )
      })
      # Now re-writes the formula in general terms (without linking to INLA::inla.surv)
      formula <- as.formula(gsub("INLA::inla.surv","Surv",deparse(formula)))
      time2run <- unlist(lapply(mod,function(x) x$cpu.used["Total"]))
      names(time2run) <- labs
      # NB Internally, DIC = model$dic$mean.deviance+model$dic$p.eff
      dic <- unlist(lapply(mod,function(x) x$dic$dic))
      # NB But to estimate the AIC & BIC is probably best to use model$dic$deviance.mean!
      aic <- unlist(lapply(1:length(mod), function(i) 2*mod[[i]]$dic$p.eff+mod[[i]]$dic$deviance.mean))
      names(aic) <- NULL
      bic <- unlist(lapply(1:length(mod),function(i) 
        mod[[i]]$dic$deviance.mean+mod[[i]]$dic$p.eff*log(mod[[i]]$size.linear.predictor$n)))
      names(bic) <- NULL
      for (i in 1:length(distr)) {
        mod[[i]]$dlist$name <- distr[i]}
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
    if(exists("k",where=exArgs)) {k <- exArgs$k} else {k <- 0}
    if(exists("cores",where=exArgs)) {cores <- exArgs$cores} else {cores <- 1}

    non.on.log.scale <- c("genf","gengamma","lognormal")
    
    # Loads the pre-compiled models
    dso <- stanmodels
    
    touse <- time2run <- numeric()
    for (i in 1:length(distr)) {
      touse[i] <- match(distr[i],availables.hmc)
    }
    dso <- lapply(touse,function(x) dso[[x]])

    mod <- lapply(1:length(distr), function(x) {
      # First makes the data list
      if (distr[x] %in% c("gamma","gengamma","genf")) {
        # If model is Gamma, GenGamma or GenF, then use the "obs vs" censored format
        data.stan <- list(t=data[data[,vars$event]==1,vars$time],d=data[data[,vars$event]==0,vars$time])
        data.stan$n_obs <- length(data.stan$t)
        data.stan$n_cens <- length(data.stan$d)
        data.stan$X_obs <- matrix(model.matrix(formula,data)[data[,vars$event]==1,],nrow=data.stan$n_obs,byrow=F)
        data.stan$X_cens <- matrix(model.matrix(formula,data)[data[,vars$event]==0,],nrow=data.stan$n_cens,byrow=F)
        data.stan$H=ncol(data.stan$X_obs)
        # NB: Stan doesn't allow vectors of size 1, so if there's only one covariate (eg intercept only), needs a little trick
        if (data.stan$H==1) {
          data.stan$X_obs <- cbind(data.stan$X_obs,rep(0,data.stan$n_obs))
          data.stan$X_cens <- cbind(data.stan$X_cens,rep(0,data.stan$n_cens))
          data.stan$H <- ncol(data.stan$X_obs)
        }
      } 
      if (distr[x] %in% c("exponential","gompertz","weibull","weibullPH","loglogistic","lognormal")) {
        # If it's one of the others (except polyweibull), use the "h,S" format
        data.stan <- list(t=data[,vars$time], d=data[,vars$event])
        data.stan$n <- length(data.stan$t) 
        data.stan$X <- model.matrix(formula,data)
        data.stan$H <- ncol(data.stan$X)
        # NB: Stan doesn't allow vectors of size 1, so if there's only one covariate (eg intercept only), needs a little trick
        if (data.stan$H==1) {
          data.stan$X <- cbind(data.stan$X,rep(0,data.stan$n))
          data.stan$H <- ncol(data.stan$X)
        }
      }
      if (distr[x]=="rps"){
        # If it's Royston-Parmar splines, then gets the correct data 
        knots <- quantile(log(data[data[,vars$event]==1,vars$time]), seq(0, 1, length = k+2))
        # Uses flexsurv to compute the basis and derivatives of the basis
        ######################################
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
        dbasis <- function (knots, x) {
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
        B <- basis(knots,log(data[,vars$time]))
        DB <- dbasis(knots,log(data[,vars$time]))
        # Now checks to see whether the user wants to specify covariates and removes the intercept from the formula (for identifiability)
        mm <- model.matrix(formula,data)[,-1]
        # a. if the formula is ~ 1, then adds two fictional covariates of all 0s
        if (length(mm)<1) {
            mm <- matrix(rep(0,nrow(data)),nrow=nrow(data),ncol=2)
        }
        # b. in case there's only one covariate, then adds another fake covariate of all 0s
        if (is.null(dim(mm))) {
         mm <- cbind(mm,rep(0,length(mm)))
        }
        data.stan=list(t=data[,vars$time], d=data[,vars$event], n=nrow(data),M=k,X=mm,H=ncol(mm),B=B,DB=DB,
                       mu_gamma=rep(0,k+2),sigma_gamma=rep(5,k+2),knots=knots) 
      }
      # ###########################################################################################################################
      # ### Poly-Weibull is in theory possible and pre-compiled, but it poses problems if the formula is a list
      # if (distr[x]=="polyweibull") {
      #   data.stan <- list(t=data[,vars$time], d=data[,vars$event]); data.stan$n <- length(data.stan$t); 
      #   data.stan$M <- length(formula)
      #   X <- lapply(1:data.stan$M,function(i) model.matrix(formula[[i]],data))
      #   data.stan$H <- max(unlist(lapply(1:data.stan$M,function(i) ncol(X[[i]]))))
      #   X <- lapply(1:data.stan$M,function(i) {
      #     if(ncol(X[[i]]<data.stan$H)) {
      #       X[[i]] <- cbind(X[[i]],matrix(0,nrow=nrow(X[[i]]),ncol=(data.stan$H-ncol(X[[i]]))))
      #     }
      #   })
      #   data.stan$X=array(NA,c(data.stan$M,data.stan$n,data.stan$H))
      #   for (m in 1:data.stan$M) {
      #     data.stan$X[m,,] <- X[[m]]
      #   }
      # }
      # # Linear predictor coefficients
      # if (distr[x]=="polyweibull") {
      #   data.stan$mu_beta=matrix(0,nrow=data.stan$H,ncol=data.stan$M) 
      #   data.stan$sigma_beta=matrix(5,nrow=data.stan$H,ncol=data.stan$M)
      # } else {
      data.stan$mu_beta=rep(0,data.stan$H)
      if (distr[x]%in%non.on.log.scale) {
        data.stan$sigma_beta <- rep(100,data.stan$H)
      } else {
        data.stan$sigma_beta <- rep(5,data.stan$H)
      }
      # }
      # Ancillary parameters
      if (distr[x]=="gamma") {data.stan$a_alpha=data.stan$b_alpha <- 0.1}
      if (distr[x]=="genf") {
        data.stan$a_sigma=data.stan$b_sigma=0.1
        data.stan$mu_P=0
        data.stan$sigma_P=0.5
        data.stan$mu_Q=0
        data.stan$sigma_Q=2.5
      }
      if (distr[x]=="gengamma") {
        data.stan$a_sigma=data.stan$b_sigma=0.1
        data.stan$mu_Q=0
        data.stan$sigma_Q=100
      }
      if (distr[x] %in% c("gompertz","loglogistic","weibull","weibullPH")) {data.stan$a_alpha=data.stan$b_alpha=0.1}
      if (distr[x]=="lognormal") {
        data.stan$a_alpha=0
        data.stan$b_alpha=5}

      # These are modified if the user gives values in the call to fit.models
      if(exists("priors",where=exArgs)) {
        priors=exArgs$priors
        # If the user has not given values for all the distrs, then fill priors with empty lists
        if(length(priors)<length(distr)) {
          for (i in (length(priors)+1):length(distr)) {
            priors[[i]] <- list()
          }
        }
        # Linear predictor coefficients
        if(!is.null(priors[[x]]$mu_beta)) {
          data.stan$mu_beta=priors[[x]]$mu_beta
        }
        if(!is.null(priors[[x]]$sigma_beta)) {
          data.stan$sigma_beta <- priors[[x]]$sigma_beta
        }
        if(!is.null(priors[[x]]$mu_gamma) & distr[x]=="rps") {
          data.stan$mu_gamma <- priors[[x]]$mu_gamma
        }
        if(!is.null(priors[[x]]$sigma_gamma) & distr[x]=="rps") {
            data.stan$sigma_gamma <- priors[[x]]$sigma_gamma
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
      tic <- proc.time()
      out=rstan::sampling(dso[[x]],data.stan,chains=chains,iter=iter,warmup=warmup,thin=thin,seed=seed,control=control[[x]],
                      pars=pars,include=include,cores=cores)
      toc <- proc.time()-tic
      time2run[x]=toc[3]
      list(out=out,data.stan=data.stan,time2run=time2run)
    })
    if(exists("save.stan",where=exArgs)) {save.stan <- exArgs$save.stan} else {save.stan=FALSE}
    #   save.stan=ifelse(any(distr=="rps"),TRUE,FALSE)
    # }
    time_survHE <- unlist(lapply(1:length(mod),function(x) mod[[x]]$time2run))
    time_stan <- unlist(lapply(1:length(mod),function(x) sum(rstan::get_elapsed_time(mod[[x]]$out))))
    time2run <- pmin(time_survHE,time_stan)
    names(time2run) <- labs

	# Computes the log-likelihood 
    dic <- aic <- bic <- dic2 <- numeric()
	  for (i in 1:length(distr)) {
	    # Extracts the simulations for the relevant parameters
	    beta <- rstan::extract(mod[[i]]$out)$beta
	    # To safeguard against very asymmetric densities use the median (instead of the mean)
	    beta.hat <- apply(beta,2,median)
	    data.stan=mod[[i]]$data.stan
	    
	    if (distr[i] %in% c("exponential","weibull","weibullPH","gompertz","lognormal","loglogistic")) {
	      linpred <- beta%*%t(data.stan$X)
	      linpred.hat <- beta.hat%*%t(data.stan$X)
	    }

	    if(distr[i]=="exponential") {
	      logf <- matrix(unlist(lapply(1:nrow(linpred),function(i) {
	        data.stan$d*log(hexp(data.stan$t,exp(linpred[i,]))) + log(1-pexp(data.stan$t,exp(linpred[i,])))
	      })),nrow=nrow(linpred),byrow=T)
	      logf.hat <- matrix(data.stan$d*log(hexp(data.stan$t,exp(linpred.hat))) + log(1-pexp(data.stan$t,exp(linpred.hat))),nrow=1)
	      # Number of parameters (for AIC): rate + covariates
	      npars <- 1+sum(1-apply(data.stan$X,2,function(x) all(x==0)))
	    }

	    if (distr[i]=="weibull") {
	        shape <- as.numeric(rstan::extract(mod[[i]]$out)$alpha)
	        shape.hat <- median(shape)
	        logf <- matrix(unlist(lapply(1:nrow(linpred),function(i) {
	          data.stan$d*log(hweibull(data.stan$t,shape[i],exp(linpred[i,]))) + log(1-pweibull(data.stan$t,shape[i],exp(linpred[i,])))
	        })),nrow=nrow(linpred),byrow=T)
	        logf.hat <- matrix(data.stan$d*log(hweibull(data.stan$t,shape.hat,exp(linpred.hat))) + log(1-pweibull(data.stan$t,shape.hat,exp(linpred.hat))),nrow=1)
	        # Number of parameters (for AIC): shape, scale + covariates
	        npars <- 2+sum(1-apply(data.stan$X,2,function(x) all(x==0)))
	    }
	    
	    if (distr[i]=="weibullPH") {
	      shape <- as.numeric(rstan::extract(mod[[i]]$out)$alpha)
	      shape.hat=median(shape)
	      logf <- matrix(unlist(lapply(1:nrow(linpred),function(i) {
	        data.stan$d*log(hweibullPH(data.stan$t,shape[i],exp(linpred[i,]))) + 
	          log(1-pweibullPH(data.stan$t,shape[i],exp(linpred[i,])))
	      })),nrow=nrow(linpred),byrow=T)
	      logf.hat <- matrix(data.stan$d*log(hweibullPH(data.stan$t,shape.hat,exp(linpred.hat)))+
	                          log(1-pweibullPH(data.stan$t,shape.hat,exp(linpred.hat))),nrow=1)
	      # Number of parameters (for AIC): shape, scale + covariates
	      npars <- 2+sum(1-apply(data.stan$X,2,function(x) all(x==0)))
	    }

	    if (distr[i]=="gompertz") {
	      shape <- as.numeric(rstan::extract(mod[[i]]$out)$alpha)
	      shape.hat=median(shape)
	      logf <- matrix(unlist(lapply(1:nrow(linpred),function(i) {
	        data.stan$d*log(hgompertz(data.stan$t,shape[i],exp(linpred[i,]))) + 
	          log(1-pgompertz(data.stan$t,shape[i],exp(linpred[i,])))
	      })),nrow=nrow(linpred),byrow=T)
	      logf.hat <- matrix(data.stan$d*log(hgompertz(data.stan$t,shape.hat,exp(linpred.hat)))+
	                          log(1-pgompertz(data.stan$t,shape.hat,exp(linpred.hat))),nrow=1)
	      # Number of parameters (for AIC): shape, rate + covariates
	      npars <- 2+sum(1-apply(data.stan$X,2,function(x) all(x==0)))
	    }
	    
	    if (distr[i]=="gamma") {
	        shape <- as.numeric(rstan::extract(mod[[i]]$out)$alpha)
	        shape.bar <- median(shape)
	        lo <- exp(beta%*%t(data.stan$X_obs))
	        lc <- exp(beta%*%t(data.stan$X_cens))
	        lo.bar <- exp(beta.hat%*%t(data.stan$X_obs))
	        lc.bar <- exp(beta.hat%*%t(data.stan$X_cens))
	        f=matrix(unlist(lapply(1:nrow(lo),function(i) dgamma(data.stan$t,shape[i],lo[i,]))),nrow=nrow(lo),byrow=T)
	        f.bar=matrix(unlist(lapply(1:nrow(lo.bar),function(i) dgamma(data.stan$t,shape.bar,lo.bar[i,]))),nrow=1,byrow=T)
	        s=matrix(unlist(lapply(1:nrow(lc),function(i) 1-pgamma(data.stan$d,shape[i],lc[i,]))),nrow=nrow(lc),byrow=T)
	        s.bar=matrix(unlist(lapply(1:nrow(lc.bar),function(i) 1-pgamma(data.stan$d,shape.bar,lc.bar[i,]))),nrow=1,byrow=T)
	        # Number of parameters (for AIC): shape, rate + covariates
	        npars <- 2+sum(1-apply(data.stan$X_obs,2,function(x) all(x==0)))
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
	        npars <- 3+sum(1-apply(data.stan$X_obs,2,function(x) all(x==0)))
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
	        npars <- 4+sum(1-apply(data.stan$X_obs,2,function(x) all(x==0)))
	    }
	    
	    if (distr[i]=="lognormal") {
	      sigma=as.numeric(rstan::extract(mod[[i]]$out)$alpha)
	      sigma.hat=median(sigma)
	      logf <- matrix(unlist(lapply(1:nrow(linpred),function(i) {
	        data.stan$d*log(hlnorm(data.stan$t,(linpred[i,]),sigma[i])) + 
	          log(1-plnorm(data.stan$t,(linpred[i,]),sigma[i]))
	      })),nrow=nrow(linpred),byrow=T)
	      logf.hat <- matrix(data.stan$d*log(hlnorm(data.stan$t,(linpred.hat),sigma.hat))+
	                          log(1-plnorm(data.stan$t,(linpred.hat),sigma.hat)),nrow=1)
	      # Number of parameters (for AIC): meanlog, sdlog + covariates
	      npars <- 2+sum(1-apply(data.stan$X,2,function(x) all(x==0)))
	    }

	    if (distr[i]=="loglogistic") {
	      sigma=as.numeric(rstan::extract(mod[[i]]$out)$alpha)
	      sigma.hat=median(sigma)
	      logf <- matrix(unlist(lapply(1:nrow(linpred),function(i) {
	        data.stan$d*log(hllogis(data.stan$t,sigma[i],exp(linpred[i,]))) + 
	          log(1-pllogis(data.stan$t,sigma[i],exp(linpred[i,])))
	      })),nrow=nrow(linpred),byrow=T)
	      logf.hat <- matrix(data.stan$d*log(hllogis(data.stan$t,sigma.hat,exp(linpred.hat)))+
	                          log(1-pllogis(data.stan$t,sigma.hat,exp(linpred.hat))),nrow=1)
	      # Number of parameters (for AIC): shape, scale + covariates
	      npars <- 2+sum(1-apply(data.stan$X,2,function(x) all(x==0)))
	    }	      

	    if (distr[i]=="rps") {
	        gamma <- rstan::extract(mod[[i]]$out)$gamma
	        gamma.hat <- apply(gamma,2,median)
	        logf <- data.stan$d*(-log(data.stan$t)+log(gamma%*%t(data.stan$DB)) + gamma%*%t(data.stan$B)+ beta%*%t(data.stan$X)) -
	                     exp(gamma%*%t(data.stan$B)+ beta%*%t(data.stan$X))
	        logf.hat <- t(data.stan$d*(-log(data.stan$t)+log(data.stan$DB%*%gamma.hat)+data.stan$B%*%gamma.hat + data.stan$X%*%beta.hat) - 
	            exp(data.stan$B%*%gamma.hat + data.stan$X%*%beta.hat))
	        # Number of parameters (for AIC): gamma + covariates
	        npars <- length(gamma.hat)+sum(apply(data.stan$X,2,function(x) 1-all(x==0)))
	    }	  
	    
	    # Now computes the log-likelihood and then deviance and DIC, AIC, BIC
	    # Little function to compute the log-likelihood (for the obs vs censored cases)
	    compute.loglik <- function(f,s) {
	      loglik <- (apply(log(f),1,sum)+apply(log(s),1,sum))
	      return(loglik)
	    }
	    if (distr[i] %in% c("gamma","gengamma","genf")) {
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
	    dic[i] <- mean(D.theta)+pD
	    dic2[i] <- mean(D.theta)+pV
	    # Approximates AIC & BIC using the mean deviance and the number of nominal parameters
	    aic[i] <- D.bar+2*npars                   #mean(D.theta)+2*pD
	    bic[i] <- D.bar+npars*log(data.stan$n)    #mean(D.theta)+pD*log(data.stan$n)
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
    mod <- lapply(1:length(mod),function(i) mod[[i]]$out)
  }
  # Names the models list
  names(mod) <- names(misc$time2run)
  # Finally prepares the output object
  res <- list(models=mod,model.fitting=model.fitting,method=method,misc=misc)
  # And sets its class attribute to "survHE"
  class(res) <- "survHE"
  return(res)
}
