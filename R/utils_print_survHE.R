#' Helper function to get the relevant stats to print the summary table
#' 
#' @param x The 'survHE' object with the fitted model
#' @param mod A number identifying which of the models is to be used
#' @return \item{res}{The resulting stats}
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords MLE
#' @noRd 
get_stats_mle <- function(x,mod) {
  # Can use directly 'flexsurv' output to make the results table
  res=x$models[[mod]]$res[,c(1,4,2,3),drop=FALSE]
  colnames(res)=c("mean","se","L95%","U95%")
  return(res)
}


#' Helper function to get the relevant stats to print the summary table
#' 
#' @param x The 'survHE' object with the fitted model
#' @param mod A number identifying which of the models is to be used
#' @return \item{res}{The resulting stats}
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords INLA
#' @noRd 
get_stats_inla <- function(x,mod) {
  # Calls the helper functions to make the results table
  res=do.call(paste0("rescale_stats_inla_",x$misc$model_name[mod]),
              args=list(x,mod))
  return(res)
}


#' Helper function to get the relevant stats to print the summary table
#' 
#' @param x The 'survHE' object with the fitted model
#' @param mod A number identifying which of the models is to be used
#' @return \item{res}{The resulting stats}
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords HMC
#' @noRd 
get_stats_hmc <- function(x,mod) {
  ######quiet(print(x$models[[mod]]))
  # Gets the original summary stats from the 'rstan' run
  table = rstan::summary(x$models[[mod]])$summary[,c("mean","sd","2.5%","97.5%")]
  ###table <- cbind(x$models[[mod]]@.MISC$summary$msd,x$models[[mod]]@.MISC$summary$quan[,c("2.5%","97.5%")])
  # Removes the node 'lp___'
  table=table[-grep("lp__",rownames(table)),]
  # If the model is intercept only, removes the unnecessary covariates created to suit 'stan' format
  if("X_obs" %in% names(x$misc$data.stan[[1]])) {
    if(any(apply(x$misc$data.stan[[1]]$X_obs,2,function(x) all(x==0)))) {
      table=table[-grep("beta\\[2\\]",rownames(table)),]
    }
  } else {
    if(any(apply(x$misc$data.stan[[1]]$X,2,function(x) all(x==0)))) {
      table=table[-grep("beta\\[2\\]",rownames(table)),]
    }
  }
  # Now calls the helper functions to make the results table
  res=do.call(paste0("rescale_stats_hmc_",x$misc$model_name[mod]),
              args=list(table=table,x=x))
  return(res)
}

#' Helper function to rescale the stats for the Exponential model
#' 
#' @param table The table with the relevant values for the model 
#' parameters
#' @param x The original 'survHE' object
#' @return \item{res}{The resulting stats}
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords HMC Exponential
#' @noRd 
rescale_stats_hmc_exp <- function(table,x) {
  rate <- matrix(table[grep("rate",rownames(table)),],ncol=4)
  rownames(rate) <- "rate"
  effects=add_effects_hmc(table,x)
  res <- rbind(rate,effects)
  if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
  return(res)
}

#' Helper function to rescale the stats for the Weibull AFT model
#' 
#' @param table The table with the relevant values for the model 
#' parameters
#' @param x The original 'survHE' object
#' @return \item{res}{The resulting stats}
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords HMC WeibullAFT
#' @noRd 
rescale_stats_hmc_wei <- function(table,x) {
  scale <- matrix(table[grep("scale",rownames(table)),],ncol=4)
  rownames(scale) <- "scale"
  shape <- matrix(table[grep("alpha",rownames(table)),],ncol=4)
  rownames(shape) <- "shape"
  effects=add_effects_hmc(table,x)
  res <- rbind(shape,scale,effects)
  if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
  return(res)
}

#' Helper function to rescale the stats for the Weibull PH model
#' 
#' @param table The table with the relevant values for the model 
#' parameters
#' @param x The original 'survHE' object
#' @return \item{res}{The resulting stats}
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords HMC WeibullPH
#' @noRd 
rescale_stats_hmc_wph <- function(table,x) {
  scale <- matrix(table[grep("scale",rownames(table)),],ncol=4)
  rownames(scale) <- "scale"
  shape <- matrix(table[grep("alpha",rownames(table)),],ncol=4)
  rownames(shape) <- "shape"
  effects=add_effects_hmc(table,x)
  res <- rbind(shape,scale,effects)
  if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
  return(res)
}

#' Helper function to rescale the stats for the Gompertz model
#' 
#' @param table The table with the relevant values for the model 
#' parameters
#' @param x The original 'survHE' object
#' @return \item{res}{The resulting stats}
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords HMC Gompertz
#' @noRd 
rescale_stats_hmc_gom <- function(table,x) {
  rate <- matrix(table[grep("rate",rownames(table)),],ncol=4)
  rownames(rate) <- "rate"
  shape <- matrix(table[grep("alpha",rownames(table)),],ncol=4)
  rownames(shape) <- "shape"
  effects=add_effects_hmc(table,x)
  res <- rbind(shape,rate,effects)
  if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
  return(res)
}

#' Helper function to rescale the stats for the logNormal model
#' 
#' @param table The table with the relevant values for the model 
#' parameters
#' @param x The original 'survHE' object
#' @return \item{res}{The resulting stats}
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords HMC logNormal
#' @noRd 
rescale_stats_hmc_lno <- function(table,x) {
  meanlog <- matrix(table[grep("meanlog",rownames(table)),],ncol=4)
  rownames(meanlog) <- "meanlog"
  sdlog <- matrix(table[grep("alpha",rownames(table)),],ncol=4)
  rownames(sdlog) <- "sdlog"
  effects=add_effects_hmc(table,x)
  res <- rbind(meanlog,sdlog,effects)
  if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
  return(res)
}

#' Helper function to rescale the stats for the Gamma model
#' 
#' @param table The table with the relevant values for the model 
#' parameters
#' @param x The original 'survHE' object
#' @return \item{res}{The resulting stats}
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords HMC Gamma
#' @noRd 
rescale_stats_hmc_gam <- function(table,x) {
  rate <- matrix(table[grep("rate",rownames(table)),],ncol=4)
  rownames(rate) <- "rate"
  shape <- matrix(table[grep("alpha",rownames(table)),],ncol=4)
  rownames(shape) <- "shape"
  effects=add_effects_hmc(table)
  res <- rbind(shape,rate,effects,x)
  if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
  return(res)
}

#' Helper function to rescale the stats for the logLogistic model
#' 
#' @param table The table with the relevant values for the model 
#' parameters
#' @param x The original 'survHE' object
#' @return \item{res}{The resulting stats}
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords HMC logLogistic
#' @noRd 
rescale_stats_hmc_llo <- function(table,x) {
  rate <- matrix(table[grep("rate",rownames(table)),],ncol=4)
  rownames(rate) <- "scale"
  shape <- matrix(table[grep("alpha",rownames(table)),],ncol=4)
  rownames(shape) <- "shape"
  effects=add_effects_hmc(table,x)
  res <- rbind(shape,rate,effects)
  if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
  return(res)
}

#' Helper function to rescale the stats for the Gen F model
#' 
#' @param table The table with the relevant values for the model 
#' parameters
#' @param x The original 'survHE' object
#' @return \item{res}{The resulting stats}
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords HMC GenF
#' @noRd 
rescale_stats_hmc_gef <- function(table,x) {
  mu <- matrix(table[grep("beta",rownames(table)),],ncol=4,nrow=1)
  rownames(mu) <- "mu"
  sigma <- matrix(table[grep("sigma",rownames(table)),],ncol=4)
  rownames(sigma) <- "sigma"
  Q <- matrix(table[grep("Q",rownames(table)),],ncol=4)
  rownames(Q) <- "Q"
  P <- matrix(table[match("P",rownames(table)),],ncol=4)
  rownames(P) <- "P"
  effects=add_effects_hmc(table,x)
  res <- rbind(mu,sigma,Q,P,effects)
  if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
  return(res)
}

#' Helper function to rescale the stats for the Gen Gamma model
#' 
#' @param table The table with the relevant values for the model 
#' parameters
#' @param x The original 'survHE' object
#' @return \item{res}{The resulting stats}
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords HMC GenGamma
#' @noRd 
rescale_stats_hmc_gga <- function(table,x) {
  mu <- matrix(table[grep("beta",rownames(table)),,drop=FALSE][1,],ncol=4,nrow=1)
  rownames(mu) <- "mu"
  sigma <- matrix(table[grep("sigma",rownames(table)),],ncol=4)
  rownames(sigma) <- "sigma"
  Q <- matrix(table[grep("Q",rownames(table)),],ncol=4)
  rownames(Q) <- "Q"
  effects=add_effects_hmc(table,x)
  res <- rbind(mu,sigma,Q,effects)
  if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
  return(res)
}

#' Helper function to rescale the stats for the RPS model
#' 
#' @param table The table with the relevant values for the model 
#' parameters
#' @param x The original 'survHE' object
#' @return \item{res}{The resulting stats}
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords HMC Royston-Parmar splines
#' @noRd 
rescale_stats_hmc_rps <- function(table,x) {
  gamma <- matrix(table[grep("gamma",rownames(table)),],ncol=4)
  rownames(gamma) <- paste0("gamma",0:(nrow(gamma)-1))
  # If there covariates adds their effects
  if(length(grep("beta",rownames(table)))>0) {
    effects <- matrix(table[grep("beta",rownames(table)),],ncol=4)
    cn=colnames(model.matrix(x$misc$formula,x$misc$data))
    rownames(effects) <- cn[-grep("Intercept",cn),drop=FALSE]
  } else {
    effects <- matrix(NA,nrow=0,ncol=4)
  }
  res <- rbind(gamma,effects)
  if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
  return(res)
}

#' Helper function to rescale the stats for the Poly-Weibull model
#' 
#' @param table The table with the relevant values for the model 
#' parameters
#' @param x The original 'survHE' object
#' @return \item{res}{The resulting stats}
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords HMC Poly-Weibull
#' @noRd 
rescale_stats_hmc_pow <- function(table,x) {
  rownames(table)[grep("alpha",rownames(table))]=paste0("shape_",1:length(grep("alpha",rownames(table))))

  # Figures out which beta coefficients should be removed (because they are multiplied by a covariate that is constantly 0)
  to.rm=matrix(unlist(lapply(1:length(x$misc$formula),function(m) apply(x$misc$data.stan[[1]]$X[m,,],2,function(x) all(x==0)))),
               nrow=length(x$misc$formula),byrow=T)
  nmatch <- length(which(to.rm==T))
  if(nmatch>0){
    idx <- matrix(unlist(lapply(1:nmatch,function(i) {
      paste0(which(to.rm==TRUE,arr.ind=T)[i,],collapse=",")
    })))  
  } else {idx=NULL}
  if (!is.null(nrow(idx))) {
    take.out <- match(paste0("beta[",idx,"]"),rownames(table))
  } else {take.out=NULL}
  if(all(!is.null(take.out))) {table=table[-take.out,]}
  effects=table[-grep("shape",rownames(table)),]
  rownames(effects) <- unlist(lapply(1:x$misc$data.stan[[1]]$M,function(m) {
    paste0(colnames(model.matrix(x$misc$formula[[m]],x$misc$data)),"_",m)
  }))
  res <- rbind(table[grep("shape",rownames(table)),],effects)
  if (is.null(dim(res))) {names(res) <- c("mean","se","L95%","U95%")} else {colnames(res) <- c("mean","se","L95%","U95%")}
  return(res)
}

#' Helper function to rescale the stats for the Weibull AFT model
#' 
#' @param table The table with the relevant values for the model 
#' parameters
#' @return \item{res}{The resulting stats}
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords INLA WeibullAFT
#' @noRd 
rescale_stats_inla_wei <- function(x,mod,nsim=1000) {
  # The scale and effects are computed as a *non linear* function of the AFT effects and the shape
  # But for simplicity can approximate this using 'inla.rmarginal'
  shape_sim=INLA::inla.rmarginal(nsim,x$models[[mod]]$marginals.hyperpar[[1]])
  fixeff_sim=lapply(1:nrow(x$models[[mod]]$summary.fixed),function(i) {
    INLA::inla.rmarginal(nsim,x$models[[mod]]$marginals.fixed[[i]])
  })
  shape=shape_sim %>% make_stats %>% matrix(.,ncol=4)
  ## NB: INLA has a weird parameterisation and with Weibull AFT, the coefficients have the wrong sign
  if(attributes(terms(x$misc$formula))$intercept==1) {
    scale=exp(-fixeff_sim[[1]]+log(max(x$misc$km$time))) %>% make_stats%>% matrix(.,ncol=4)
  }
  rownames(scale) <- "scale"
  rownames(shape) <- "shape"
  res=rbind(shape,scale)
  # If there are covariates then add them too
  if(length(fixeff_sim)>1) {
    effects=lapply(2:nrow(x$models[[mod]]$summary.fixed),function(i) {
      -fixeff_sim[[i]] 
    })
    effects=matrix(unlist(lapply(effects,function(i) i %>% make_stats)),
                   nrow=length(fixeff_sim)-1,ncol=4,byrow=T)
    rownames(effects) <- x$models[[mod]]$names.fixed[-1]
    res=rbind(res,effects)
  }
  colnames(res)=c("mean","se","L95%","U95%")
  return(res)
}

#' Helper function to rescale the stats for the Weibull PH model
#' 
#' @param table The table with the relevant values for the model 
#' parameters
#' @return \item{res}{The resulting stats}
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords INLA WeibullPH
#' @noRd 
rescale_stats_inla_wph <- function(x,mod,nsim=1000) {
  # The scale and effects are computed as a *non linear* function of the AFT effects and the shape
  # But for simplicity can approximate this using 'inla.rmarginal'
  shape_sim=INLA::inla.rmarginal(nsim,x$models[[mod]]$marginals.hyperpar[[1]])
  fixeff_sim=lapply(1:nrow(x$models[[mod]]$summary.fixed),function(i) {
    INLA::inla.rmarginal(nsim,x$models[[mod]]$marginals.fixed[[i]])
  })
  shape=shape_sim %>% make_stats %>% matrix(.,ncol=4)
  if(attributes(terms(x$misc$formula))$intercept==1) {
    scale=exp(fixeff_sim[[1]]+log(max(x$misc$km$time)))^(-shape_sim) %>% make_stats%>% matrix(.,ncol=4)
  }
  rownames(scale) <- "scale"
  rownames(shape) <- "shape"
  res=rbind(shape,scale)
  # If there are covariates then add them too
  if(length(fixeff_sim)>1) {
    effects=lapply(2:nrow(x$models[[mod]]$summary.fixed),function(i) {
      fixeff_sim[[i]]
    })
    effects=matrix(unlist(lapply(effects,function(i) i %>% make_stats)),
                   nrow=length(fixeff_sim)-1,ncol=4,byrow=T)
    rownames(effects) <- x$models[[mod]]$names.fixed[-1]
    res=rbind(res,effects)
  }
  colnames(res)=c("mean","se","L95%","U95%")
  return(res)
}

#' Helper function to rescale the stats for the Exponential model
#' 
#' @param table The table with the relevant values for the model 
#' parameters
#' @return \item{res}{The resulting stats}
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords INLA Exponential
#' @noRd 
rescale_stats_inla_exp <- function(x,mod,nsim=1000) {
  fixeff_sim=lapply(1:nrow(x$models[[mod]]$summary.fixed),function(i) {
    INLA::inla.rmarginal(nsim,x$models[[mod]]$marginals.fixed[[i]])
  })
  if(attributes(terms(x$misc$formula))$intercept==1) {
    rate=exp(fixeff_sim[[1]]-log(max(x$misc$km$time))) %>% make_stats %>% matrix(.,ncol=4)
  }
  rownames(rate)="rate"
  res=rate
  if(length(fixeff_sim)>1) {
    effects=lapply(2:nrow(x$models[[mod]]$summary.fixed),function(i) {
      fixeff_sim[[i]]
    })
    effects=matrix(unlist(lapply(effects,function(i) i %>% make_stats)),
                   nrow=length(fixeff_sim)-1,ncol=4,byrow=T)
    rownames(effects) <- x$models[[mod]]$names.fixed[-1]
    res=rbind(res,effects)
  }
  colnames(res)=c("mean","se","L95%","U95%")
  return(res)
}

#' Helper function to rescale the stats for the logNormal model
#' 
#' @param table The table with the relevant values for the model 
#' parameters
#' @return \item{res}{The resulting stats}
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords INLA logNormal
#' @noRd 
rescale_stats_inla_lno <- function(x,mod,nsim=1000) {
  prec_sim=INLA::inla.rmarginal(nsim,x$models[[mod]]$marginals.hyperpar[[1]])
  fixeff_sim=lapply(1:nrow(x$models[[mod]]$summary.fixed),function(i) {
    INLA::inla.rmarginal(nsim,x$models[[mod]]$marginals.fixed[[i]])
  })
  if(attributes(terms(x$misc$formula))$intercept==1) {
    meanlog=(fixeff_sim[[1]]+log(max(x$misc$km$time))) %>% make_stats %>% matrix(.,ncol=4)
  }
  sdlog=sqrt(1/prec_sim) %>% make_stats %>% matrix(.,ncol=4)
  rownames(meanlog)="meanlog"
  rownames(sdlog)="sdlog"
  res=rbind(meanlog,sdlog)
  if(length(fixeff_sim)>1) {
    effects=lapply(2:nrow(x$models[[mod]]$summary.fixed),function(i) {
      fixeff_sim[[i]]
    })
    effects=matrix(unlist(lapply(effects,function(i) i %>% make_stats)),
                   nrow=length(fixeff_sim)-1,ncol=4,byrow=T)
    rownames(effects) <- x$models[[mod]]$names.fixed[-1]
    res=rbind(res,effects)
  }
  colnames(res)=c("mean","se","L95%","U95%")
  return(res)
}

#' Helper function to rescale the stats for the logLogistic model
#' 
#' @param table The table with the relevant values for the model 
#' parameters
#' @return \item{res}{The resulting stats}
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords INLA logLogistic
#' @noRd 
rescale_stats_inla_llo <- function(x,mod,nsim=1000) {
  # Uses 'variant=1' in INLA
  shape_sim=INLA::inla.rmarginal(nsim,x$models[[mod]]$marginals.hyperpar[[1]])
  fixeff_sim=lapply(1:nrow(x$models[[mod]]$summary.fixed),function(i) {
    INLA::inla.rmarginal(nsim,x$models[[mod]]$marginals.fixed[[i]])
  })
  shape=shape_sim %>% make_stats %>% matrix(.,ncol=4)
  if(attributes(terms(x$misc$formula))$intercept==1) {
    scale=exp(-fixeff_sim[[1]]+log(max(x$misc$km$time))) %>% make_stats %>% matrix(.,ncol=4)
  }
  rownames(shape)="shape"
  rownames(scale)="scale"
  res=rbind(shape,scale)
  if(length(fixeff_sim)>1) {
    effects=lapply(2:nrow(x$models[[mod]]$summary.fixed),function(i) {
      -fixeff_sim[[i]]
    })
    effects=matrix(unlist(lapply(effects,function(i) i %>% make_stats)),
                   nrow=length(fixeff_sim)-1,ncol=4,byrow=T)
    rownames(effects) <- x$models[[mod]]$names.fixed[-1]
    res=rbind(res,effects)
  }
  colnames(res)=c("mean","se","L95%","U95%")
  return(res)
}

#' Helper function to rescale the stats for the Gompertz model
#' 
#' @param table The table with the relevant values for the model 
#' parameters
#' @return \item{res}{The resulting stats}
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords INLA Gompertz
#' @noRd 
rescale_stats_inla_gom <- function(x,mod,nsim=1000) {
  shape_sim=INLA::inla.rmarginal(nsim,x$models[[mod]]$marginals.hyperpar[[1]])/max(x$misc$km$time)
  fixeff_sim=lapply(1:nrow(x$models[[mod]]$summary.fixed),function(i) {
    INLA::inla.rmarginal(nsim,x$models[[mod]]$marginals.fixed[[i]])
  })
  shape=shape_sim %>% make_stats %>% matrix(.,ncol=4)
  # Need to rescale only if there is an intercept
  if(attributes(terms(x$misc$formula))$intercept==1) {
   rate=exp(fixeff_sim[[1]]-log(max(x$misc$km$time))) %>% make_stats %>% matrix(.,ncol=4)
  }
  rownames(shape)="shape"
  rownames(rate)="rate"
  res=rbind(shape,rate)
  if(length(fixeff_sim)>1) {
    effects=lapply(2:nrow(x$models[[mod]]$summary.fixed),function(i) {
      fixeff_sim[[i]]
    })
    effects=matrix(unlist(lapply(effects,function(i) i %>% make_stats)),
                   nrow=length(fixeff_sim)-1,ncol=4,byrow=T)
    rownames(effects) <- x$models[[mod]]$names.fixed[-1]
    res=rbind(res,effects)
  }
  colnames(res)=c("mean","se","L95%","U95%")
  return(res)
}

#' Helper function to create summary stats
#' 
#' @param x A vector of simulations
#' @return \item{tab}{The resulting stats}
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords Summaries Print
#' @noRd 
make_stats <- function(x) {
  tab=c(mean(x), sd(x), quantile(x, 0.025), quantile(x,0.975))
  return(tab)
}

#' Helper function to for Stan, which needs to first print the output of 
#' the model before you can access the elements in the object 
#' '@.MISC$summary', so can use this function to print quietly...
#' 
#' @param x A 'survHE' object
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords HMC Stan
#' @noRd 
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

#' Helper function to checks whether covariates effects should be 
#' included in the 'res' table for HMC
#' 
#' @param table The table with the summary statistics
#' @param x The original 'survHE' object
#' @return \item{effects}{The effects}
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords HMC Stan
#' @noRd 
add_effects_hmc <- function(table,x) {
  # If there's more than one beta, then there are "effects" (otherwise it's only intercept)
  if(length(grep("beta",rownames(table)))>1) {
    effects <- matrix(table[grep("beta",rownames(table)),],ncol=4)
    rownames(effects) <- colnames(model.matrix(x$misc$formula,x$misc$data))
    # Now removes the line with the intercept (which is already rescaled to the shape/rate/mean parameter)
    effects=effects[-grep("Intercept",rownames(effects)),,drop=FALSE]
  } else {
    effects <- matrix(NA,nrow=0,ncol=4)
  }
  return(effects)
}

#' Helper function to create the original summary table
#' 
#' @param x The 'survHE' model
#' @param mod Which of the models to be used
#' @param digits The number of digits to print
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords MLE
#' @noRd 
original_table_mle <- function(x,mod,digits) {
  print(x$models[[mod]],digits=digits)
}

#' Helper function to create the original summary table
#' 
#' @param x The 'survHE' model
#' @param mod Which of the models to be used
#' @param digits The number of digits to print
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords INLA
#' @noRd 
original_table_inla <- function(x,mod,digits) {
  print(summary(x$models[[mod]]),digits=digits)
  cat("\n")
  cat("NB: notice that INLA models are fitted to data rescaled in [0-1] for computational stability.")
  cat("\nThe estimates are rescaled on the original scale, applying a suitable back-transformation.") 
  cat("\nThe numbers shown when 'original=TRUE' will be different than those shown in the 'survHE' format.")
}

#' Helper function to create the original summary table
#' 
#' @param x The 'survHE' model
#' @param mod Which of the models to be used
#' @param digits The number of digits to print
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords HMC
#' @noRd 
original_table_hmc <- function(x,mod,digits) {
  print(x$models[[mod]],digits=digits)
}


#' Helper function to format the summary table with the model parameters
#' 
#' @param x The 'survHE' model
#' @param mod Which of the models to be used
#' @param res The output table
#' @param digits The number of digits to print
#' @author Gianluca Baio
#' @seealso print.survHE
#' @references Baio (2020). survHE
#' @keywords Table formatting Print
#' @noRd 
format_table <- function(x,mod,res,digits){
  # First re-format some of the labels (eg model names)
  if(x$misc$model_name[mod]=="exp") {label <- "Exponential"}
  if(x$misc$model_name[mod]=="gam") {label <- "Gamma"}
  if(x$misc$model_name[mod]=="lno") {label <- "log-Normal"}
  if(x$misc$model_name[mod]=="llo") {label <-"log-Logistic"}
  if(x$misc$model_name[mod]=="gga") {label <- "Generalised Gamma"}
  if(x$misc$model_name[mod]=="wei") {label <- "Weibull AF"}
  if(x$misc$model_name[mod]=="wph") {label <- "Weibull PH"}
  if(x$misc$model_name[mod]=="gef") {label <- "Generalised F"}
  if(x$misc$model_name[mod]=="gom") {label <- "Gompertz"}
  if(x$misc$model_name[mod]=="rps") {label <- "Royston & Parmar splines"}
  if(x$misc$model_name[mod]=="pow") {label <- "Poly-Weibull"}
  
  # Creates label of the method used
  label.met <- ifelse(
    x$method=="mle","Flexsurvreg \n(Maximum Likelihood Estimate)",
    ifelse(x$method=="inla","INLA (Bayesian inference via \nIntegrated Nested Laplace Approximation)",
           "Stan (Bayesian inference via \nHamiltonian Monte Carlo)")
  )
  
  # Finally prints the formatted table
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