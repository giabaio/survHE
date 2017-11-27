#' Plot survival curves for the models fitted using \code{fit.models}
#' 
#' Plots the results of model fit.
#' 
#' 
#' @param ...  Must include at least one result object saved as the call to the
#' \code{fit.models} function.  Other possibilities are additional (mainly
#' graphical) options. These are: \code{xlab} = a string with the label for the
#' x-axis (default = "time") \code{ylab} = a string with the label for the
#' y-axis (default = "Survival") \code{lab.trt} = a (vector of) string(s)
#' indicating the labels associated with the strata defining the different
#' survival curves to plot. Default to the value used by the Kaplan Meier
#' estimate given in \code{fit.models} \code{cex.trt} = factor by which the
#' size of the font used to write the strata is resized (default = 0.8)
#' \code{n.risk} = logical. If TRUE (defaults) writes the number at risk at
#' different time points (as determined by the Kaplan Meier estimate)
#' \code{newdata} = a list (of lists) providing the values for the relevant
#' covariates If NULL, then will use the mean values for the covariates if at
#' least one is a continuous variable, or the combination of the categorical
#' covariates. \code{xlim} = a vector determining the limits for the x-axis
#' \code{colors} = a vector of characters defining the colours in which to plot
#' the different survival curves \code{labs} = a vector of characters defining
#' the names of the models fitted \code{add.km} = TRUE (whether to also add the
#' Kaplan Meier estimates of the data)
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso Something will go here
#' @references Something will go here
#' @keywords Parametric survival models
#' @examples
#' 
#' data(bc)
#' 
#' mle = fit.models(formula=Surv(recyrs,censrec)~group,data=bc,
#'     distr="exp",method="mle")
#' plot(mle)
#' 
#' @export plot.survHE
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
  w <- which(classes=="survHE")
  original.method <- unlist(lapply(w,function(i) exArgs[[i]]$method))
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
  x <- list()
  x$models <- mods
  nmodels <- length(x$models)  # Number of models fitted by fit.models()
  class(x) <- "survHE"
  x$model.fitting <- model.fitting
  ## Needs to include in the misc object the element vars (which is used for HMC models)
  if (any(method=="hmc")) {
   x$misc <- exArgs[[min(which(original.method=="hmc"))]]$misc
   x$misc$data.stan=x$misc$data.stan[[1]]
   if (exists("X",x$misc$data.stan)) {
     x$misc$data.stan$X_obs <- x$misc$data.stan$X
   } else {
     x$misc$data.stan$X <- x$misc$data.stan$X_obs
   }
   x$misc$data.stan <- rep(list(x$misc$data.stan),nmodels)
  } else {
   # If none of the survHE objects are HMC, then just use the first
   x$misc <- exArgs[[1]]$misc
  }
  
  # Checks that extra options are specified
  if (is.null(exArgs$t)) {t <- sort(unique(x$misc$km$time))} else {t <- exArgs$t}
  if (is.null(exArgs$xlab)) {xl <- "time"} else {xl <- exArgs$xlab}
  if (is.null(exArgs$ylab)) {yl <- "Survival"} else {yl <- exArgs$ylab}
  if (is.null(exArgs$lab.trt)) {lab.trt <- names(x$misc$km$strata)} else {lab.trt<- names(x$km$strata)<-exArgs$lab.trt}
  if (is.null(exArgs$cex.trt)) {cex.trt <- 0.8} else {cex.trt <- exArgs$cex.trt}
  if (is.null(exArgs$n.risk)) {nrisk <- FALSE} else {nrisk <- exArgs$n.risk}
  if (is.null(exArgs$main)) {main <- ""} else {main <- exArgs$main}
  if (is.null(exArgs$newdata)) {newdata <- NULL} else {newdata <- exArgs$newdata}
  if (is.null(exArgs$cex.lab)) {cex.lab <- 0.8} else {cex.lab <- exArgs$cex.lab}
  
  if (is.null(exArgs$xlim) & is.null(exArgs$t)) {
    xlm <- range(pretty(x$misc$km$time))
  } 
  if (is.null(exArgs$xlim) & !is.null(exArgs$t)) {
    xlm <- range(pretty(t))
  }
  if (!is.null(exArgs$xlim) & is.null(exArgs$t)) {
    xlm <- exArgs$xlim
  }
  if (!is.null(exArgs$xlim) & !is.null(exArgs$t)) {
    xlm <- exArgs$xlim
  }
  
  if (is.null(exArgs$colors)) {
    if (nmodels>1) {colors <- (2:(nmodels+1))} else {colors <- 2}
  } else {colors <- exArgs$colors}
  if(is.null(exArgs$axes)){axes <- TRUE} else {axes <- exArgs$axes}
  if (is.null(exArgs$labs)) {
    labs <- unlist(lapply(1:length(x$models),function(i) {
      if(class(x$models[[i]])=="stanfit") {tolower(x$models[[i]]@model_name)} else {x$models[[i]]$dlist$name}
    }))
    labs[labs %in% c("weibull.quiet","weibull","weibullaf","weibullph")] <- "Weibull"
    labs[labs %in% c("exp","exponential")] <- "Exponential"
    labs[labs %in% "gamma"] <- "Gamma"
    labs[labs %in% c("lnorm","lognormal")] <- "log-Normal"
    labs[labs %in% c("llogis","loglogistic","loglogis")] <- "log-Logistic"
    labs[labs %in% "gengamma"] <- "Gen. Gamma"
    labs[labs %in% "genf"] <- "Gen. F"
    labs[labs %in% "gompertz"] <- "Gompertz"
    labs[labs %in% c("survspline","rp")] <- "Royston & Parmar splines"
  } else {labs <- exArgs$labs}
  labs <- c("Kaplan Meier",labs)
  if (is.null(exArgs$add.km)) {add.km <- TRUE} else {add.km <- exArgs$add.km}
  
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
    if(axes==TRUE) {
      axis(1)
      axis(2)}
    col <- colors
  }
  res <- lapply(1:nmodels,function(i) {
    x$method <- method[i]
    make.surv(x,nsim=1,t=t,mod=i,newdata=newdata)
  })
  
  if (!is.null(newdata)) {
    # Needs to distinguish between mle and non-mle because of how make.surv saves the S list
    options(digits=5,nsmall=2)
    pts <- list()
    for (i in 1:nmodels) {
      if (method[i]=="mle") {
        pts[[i]] <- lapply(1:length(newdata),function(j) {
          tmp <- matrix(unlist(res[[i]]$S[[j]]),ncol=4)
          cbind(tmp[,1],tmp[,2])
        })
      } else {
        pts[[i]] <- lapply(1:length(newdata),function(j) {
          res[[i]]$S[[1]][[j]]
        })
      }
    }
    colors <- 1:nmodels
    leg.txt <- character()
    for (i in 1:nmodels) {
      for (j in 1:length(newdata)) {
        points(pts[[i]][[j]],t="l",col=colors[i],lty=j)
        leg.txt[j] <- paste0(names(newdata[[j]]),"=",prettyNum(newdata[[j]],format="fg"),collapse=", ")
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
