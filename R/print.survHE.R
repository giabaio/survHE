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
  availables.inla <- c("exponential","weibull","weibullPH","lognormal","loglogistic")
  availables.hmc <- c("Exponential","Gamma","GenF","GenGamma","Gompertz","PolyWeibull","RP",
                      "WeibullAF","WeibullPH","logLogistic","logNormal")
  # If the distribution specified is not-standard (eg user-defined in MLE, or using random effects or non-standard
  # distributions in Stan), then sets original=TRUE and gives the original version of the print table.
  if (exists("original",where=exArgs)) {original=exArgs$original} else {
    if (x$method=="mle") {
      if (x$models[[mod]]$dlist$name %in% availables.mle) {original=FALSE} else {original=TRUE}
    }
    if (x$method=="inla") {
      # Needs to attache the namespace as it uses a hidden function 
      # (that is actually only required when doing knitr and printing the tables)
      if (!is.element("INLA", (.packages()))) {
        suppressPackageStartupMessages(attachNamespace("INLA"))
      }
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
	shape <- unlist(lapply(jpost,function(x) x$hyperpar))
	names(shape) <- NULL
	## NB: As of Jan 11 2017, there's a mistake in INLA and so needs to minus the argument of the exp here!
	scale <- exp(-unlist(lapply(jpost,function(x) x$latent[pos[1],])))
	effects <- matrix(NA,nrow=(length(pos)-1),ncol=nsim)
	if(length(attributes(terms(x$misc$formula))$term.labels)>0) {
	  for (j in 2:length(pos)) {
	    effects[(j-1),] <- log(exp(-unlist(lapply(jpost,function(x) x$latent[pos[j],]))))
	  }
	  rownames(effects) <- x$models[[mod]]$names.fixed[-1]
	}
	tab <- rbind(shape,scale,effects)
      }
      if(x$models[[mod]]$dlist=="weibullPH") {
	shape <- unlist(lapply(jpost,function(x) x$hyperpar))
	names(shape) <- NULL
	scale <- exp(unlist(lapply(jpost,function(x) x$latent[pos[1],])))
	effects <- matrix(NA,nrow=(length(pos)-1),ncol=nsim)
	if(length(attributes(terms(x$misc$formula))$term.labels)>0) {
	  for (j in 2:length(pos)) {
	    effects[(j-1),] <- log(exp(unlist(lapply(jpost,function(x) x$latent[pos[j],]))))
	  }
	  rownames(effects) <- x$models[[mod]]$names.fixed[-1]
	}
	tab <- rbind(shape,scale,effects)
      }
      ### OLD CODE --- worked on older versions of INLA!
#      if(x$models[[mod]]$dlist=="weibull") {
#        shape <- unlist(lapply(jpost,function(x) x$hyperpar))
#        names(shape) <- NULL
#        scale <- exp(unlist(lapply(jpost,function(x) x$latent[pos[1],])))^(1/-shape)
#        effects <- matrix(NA,nrow=(length(pos)-1),ncol=nsim)
#        if(length(attributes(terms(x$misc$formula))$term.labels)>0) {
#          for (j in 2:length(pos)) {
#            effects[(j-1),] <- log(exp(unlist(lapply(jpost,function(x) x$latent[pos[j],])))^(1/-shape))
#          }
#          rownames(effects) <- x$models[[mod]]$names.fixed[-1]
#        }
#        tab <- rbind(shape,scale,effects)
#      }
      ###
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
        prec <- unlist(lapply(jpost,function(x) x$hyperpar))
        names(prec) <- NULL
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
        shape <- unlist(lapply(jpost,function(x) x$hyperpar))
        names(shape) <- NULL
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
    take.out <- which(rownames(table)=="lp__")
    betas <- grep("beta",rownames(table))
    if(x$models[[mod]]@model_name%in%c("Gamma","GenGamma","GenF")) {
      covmat <- x$misc$data.stan[[mod]]$X_obs
    } else {
      covmat <- x$misc$data.stan[[mod]]$X
    }
    take.out <- c(take.out,betas[apply(covmat,2,function(x) all(x==0))])
    # if (is.null(x$misc$vars$factors) & is.null(x$misc$vars$covs)) {
    #   take.out = c(take.out,which(rownames(table)=="beta[2]"))
    # }
    table <- table[-take.out,]
    
    if (original==FALSE) {
      if (x$models[[mod]]@model_name=="Exponential") {
        rate <- matrix(table[grep("rate",rownames(table)),],ncol=4)
        rownames(rate) <- "rate"
        if(length(attributes(terms(x$misc$formula))$term.labels)>0) {
          effects <- matrix(table[-which(rownames(table) %in% c("rate")),][-1,],ncol=4,byrow=F)
          rownames(effects) <- colnames(model.matrix(x$misc$formula,x$misc$data))[-1]
        } else {
          effects <- matrix(NA,nrow=0,ncol=4)
        }
        res <- rbind(rate,effects)
        if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
      }
      
      if (x$models[[mod]]@model_name=="Gamma") {
        rate <- matrix(table[grep("rate",rownames(table)),],ncol=4)
        rownames(rate) <- "rate"
        shape <- matrix(table[grep("alpha",rownames(table)),],ncol=4)
        rownames(shape) <- "shape"
        if(length(attributes(terms(x$misc$formula))$term.labels)>0) {
          effects <- matrix(table[-which(rownames(table) %in% c("rate","alpha")),][-1,],ncol=4,byrow=F)
          rownames(effects) <- colnames(model.matrix(x$misc$formula,x$misc$data))[-1]
        } else {
          effects <- matrix(NA,nrow=0,ncol=4)
        }
        res <- rbind(shape,rate,effects)
        if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
      }
      
      if (x$models[[mod]]@model_name=="WeibullAF" | x$models[[mod]]@model_name=="WeibullPH") {
        scale <- matrix(table[grep("scale",rownames(table)),],ncol=4)
        rownames(scale) <- "scale"
        shape <- matrix(table[grep("alpha",rownames(table)),],ncol=4)
        rownames(shape) <- "shape"
        if(length(attributes(terms(x$misc$formula))$term.labels)>0) {
          effects <- matrix(table[-which(rownames(table) %in% c("scale","alpha")),][-1,],ncol=4,byrow=F)
          rownames(effects) <- colnames(model.matrix(x$misc$formula,x$misc$data))[-1]
        } else {
          effects <- matrix(NA,nrow=0,ncol=4)
        }
        res <- rbind(shape,scale,effects)
        if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
      }
      
      if (x$models[[mod]]@model_name=="Gompertz") {
        rate <- matrix(table[grep("rate",rownames(table)),],ncol=4)
        rownames(rate) <- "rate"
        shape <- matrix(table[grep("alpha",rownames(table)),],ncol=4)
        rownames(shape) <- "shape"
        if(length(attributes(terms(x$misc$formula))$term.labels)>0) {
          effects <- matrix(table[-which(rownames(table) %in% c("rate","alpha")),][-1,],ncol=4,byrow=F)
          rownames(effects) <- colnames(model.matrix(x$misc$formula,x$misc$data))[-1]
        } else {
          effects <- matrix(NA,nrow=0,ncol=4)
        }
        res <- rbind(shape,rate,effects)
        if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
      }
      
      if (x$models[[mod]]@model_name=="logNormal") {
        meanlog <- matrix(table[grep("meanlog",rownames(table)),],ncol=4)
        rownames(meanlog) <- "meanlog"
        sigma <- matrix(table[grep("alpha",rownames(table)),],ncol=4)
        rownames(sigma) <- "sdlog"
        if(length(attributes(terms(x$misc$formula))$term.labels)>0) {
          effects=matrix(table[-which(rownames(table) %in% c("meanlog","alpha")),][-1,],ncol=4,byrow=F)
          rownames(effects) <- colnames(model.matrix(x$misc$formula,x$misc$data))[-1]
        } else {
          effects <- matrix(NA,nrow=0,ncol=4)
        }
        res <- rbind(meanlog,sigma,effects)
        if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
      }
      
      if (x$models[[mod]]@model_name=="logLogistic") {
        rate <- matrix(table[grep("rate",rownames(table)),],ncol=4)
        rownames(rate) <- "rate"
        shape <- matrix(table[grep("alpha",rownames(table)),],ncol=4)
        rownames(shape) <- "shape"
        if(length(attributes(terms(x$misc$formula))$term.labels)>0) {
          effects=matrix(table[-which(rownames(table) %in% c("rate","alpha")),][-1,],ncol=4,byrow=F)
          rownames(effects) <- colnames(model.matrix(x$misc$formula,x$misc$data))[-1]
        } else {
          effects <- matrix(NA,nrow=0,ncol=4)
        }
        res <- rbind(shape,rate,effects)
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
          rownames(effects) <- colnames(model.matrix(x$misc$formula,x$misc$data))[-1]
        } else {
          effects <- matrix(NA,nrow=0,ncol=4)
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
          rownames(effects) <- colnames(model.matrix(x$misc$formula,x$misc$data))[-1]
        } else {
          effects <- matrix(NA,nrow=0,ncol=4)
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
          effects <- matrix(table[-grep("gamma",rownames(table)),],ncol=4,byrow=F)
          rownames(effects) <- colnames(model.matrix(x$misc$formula,x$misc$data))[-1]
        } else {
          effects <- matrix(NA,nrow=0,ncol=4)
        }
        res <- rbind(gamma,effects)
        if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
      }
      if (x$models[[mod]]@model_name=="PolyWeibull") {
        alpha <- matrix(table[grep("alpha",rownames(table)),],ncol=4)
        rownames(alpha) <- paste0("shape_",1:x$misc$data.stan$M)
        to.rm=matrix(unlist(lapply(1:length(x$misc$formula),function(m) apply(x$misc$data.stan$X[m,,],2,function(x) all(x==0)))),
                     nrow=length(x$misc$formula),byrow=T)
        nmatch <- length(which(to.rm==T))
        idx <- matrix(unlist(lapply(1:nmatch,function(i) {
          paste0(which(to.rm==TRUE,arr.ind=T)[i,],collapse=",")
        })))
        if (!is.null(nrow(idx))) {
          take.out <- match(paste0("beta[",idx,"]"),rownames(table))
        }
        if(all(!is.na(take.out))) {table=table[-take.out,]}
        effects=table[-grep("alpha",rownames(table)),]
        rownames(effects) <- unlist(lapply(1:x$misc$data.stan$M,function(m) {
          paste0(colnames(model.matrix(x$misc$formula[[m]],x$misc$data)),"_",m)
        }))
        res <- rbind(alpha,effects)
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
        take.out <- betas[unlist(lapply(1:length(x$misc$formula),function(m) apply(x$misc$data.stan$X[m,,],2,function(x) all(x==0))))]
      } else {
        take.out <- betas[unlist(lapply(1:length(x$misc$formula),function(m) apply(covmat,2,function(x) all(x==0))))]
      }
      take.out <- c(take.out,grep("lp__",rownames(rstan::summary(x$models[[mod]])$summary)))
      tab <- rstan::summary(x$models[[mod]],probs=c(.025,.975))$summary[-take.out,]
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
      if (label=="RP") {label <- "Royston & Parmar splines"}
      label.met <- "Stan (Bayesian inference via \nHamiltonian Monte Carlo)"
    } else {
      if(x$models[[mod]]$dlist$name=="exp" | x$models[[mod]]$dlist$name=="exponential") {label <- "Exponential"}
      if(x$models[[mod]]$dlist$name=="gamma") {label <- "Gamma"}
      if(x$models[[mod]]$dlist$name=="lognormal" | x$models[[mod]]$dlist$name=="lnorm") {label <- "log-Normal"}
      if(x$models[[mod]]$dlist$name=="llogis" | x$models[[mod]]$dlist$name=="loglogistic") {label <-"log-Logistic"}
      if(x$models[[mod]]$dlist$name=="gengamma") {label <- "Generalised Gamma"}
      if(x$models[[mod]]$dlist$name=="weibull" | x$models[[mod]]$dlist$name=="weibull.quiet" | x$models[[mod]]$dlist$name=="weibullPH") {label <- "Weibull"}
      if(x$models[[mod]]$dlist$name=="genf") {label <- "Generalised F"}
      if(x$models[[mod]]$dlist$name=="gompertz") {label <- "Gompertz"}
      if(x$models[[mod]]$dlist$name=="survspline") {label <- "Royston & Parmar splines"}
    }
    if(x$method=="mle") {label.met <- "Flexsurvreg \n(Maximum Likelihood Estimate)"}
    if(x$method=="inla") {label.met <- "INLA (Bayesian inference via \nIntegrated Nested Laplace Approximation)"}
    
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
