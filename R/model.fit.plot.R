#' Graphical representation of the measures of model fitting based on
#' Information Criteria
#' 
#' Plots a summary of the model fit for all the models fitted
#' 
#' Something will go here
#' 
#' @param ...  Optional inputs. Must include at least one \code{survHE} object.
#' @param type should the AIC, the BIC or the DIC plotted? (values = \code{"aic"},
#' \code{"bic"} or \code{"dic"})
#' @param scale If \code{scale='absolute'} (default), then plot the absolute value 
#' of the *IC. If \code{scale='relative'} then plot a rescaled version taking
#' the percentage increase in the *IC in comparison with the best-fitting model
#' @return Something will go here
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso Something will go here
#' @references Something will go here
#' @keywords Model fitting Parametric survival models
#' @examples
#' 
#' # Something will go here
#' 
#' @export model.fit.plot
model.fit.plot <- function(...,type="aic",scale="absolute") {
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
  fit <- list()
  fit$models <- mods
  class(fit) <- "survHE"
  fit$model.fitting <- model.fitting
  ## Needs to include in the misc object the element vars (which is used for HMC models)
  if (any(method=="hmc")) {
    fit$misc <- exArgs[[min(which(original.method=="hmc"))]]$misc
  } else {
    # If none of the survHE objects are HMC, then just use the first
    fit$misc <- exArgs[[1]]$misc
  }
  
  if (is.null(exArgs$models)) {
    models <- unlist(lapply(1:length(fit$models),function(i) {
      if(class(fit$models[[i]])=="stanfit") {tolower(fit$models[[i]]@model_name)} else {fit$models[[i]]$dlist$name}
    }))
    models[models %in% c("weibull.quiet","weibull","weibullaf","weibullph")] <- "Weibull"
    models[models %in% c("exp","exponential")] <- "Exponential"
    models[models %in% "gamma"] <- "Gamma"
    models[models %in% "gompertz"] <- "Gompertz"
    models[models %in% c("lnorm","lognormal")] <- "log-Normal"
    models[models %in% c("llogis","loglogistic","loglogis")] <- "log-Logistic"
    models[models %in% "gengamma"] <- "Gen. Gamma"
    models[models %in% "genf"] <- "Gen. F"
  } else {
    models <- exArgs$models 
  }
  
  # Defines the data to be plotted
  if (type=="aic" | type=="AIC" | type=="a" | type=="A") {
    mf <- data.frame(model=models,AIC=fit$model.fitting$aic,AIC=fit$model.fitting$aic)
    lab.type <- "AIC"
  } else if (type=="bic" | type=="BIC" | type=="b" | type=="B") {
    mf <- data.frame(model=models,BIC=fit$model.fitting$bic,BIC=fit$model.fitting$bic)
    lab.type <- "BIC"
  } else if (type=="dic" | type=="DIC" | type=="d" | type=="D") {
    mf <- data.frame(model=models,DIC=fit$model.fitting$dic,DIC=fit$model.fitting$dic)
    lab.type <- "DIC"
  }
  if(scale=="rel" | scale=="relative") {
    # Version of the plot with percentage increase in *IC in comparison to the best fitting model
    if (type=="aic" | type=="AIC" | type=="a" | type=="A") {
      mf <- data.frame(model=models,AIC=100*(fit$model.fitting$aic-min(fit$model.fitting$aic))/min(fit$model.fitting$aic),
                       AIC=fit$model.fitting$aic)
      lab.type <- "Percentage difference in AIC in comparison to the 'best-fitting' model"
      lab.type1="AIC"
    } else if (type=="bic" | type=="BIC" | type=="b" | type=="B") {
      mf <- data.frame(model=models,BIC=100*(fit$model.fitting$bic-min(fit$model.fitting$bic))/min(fit$model.fitting$bic),
                       BIC=fit$model.fitting$bic)
      lab.type <- "Percentage difference in BIC in comparison to the 'best-fitting' model"
      lab.type1="BIC"
    } else if (type=="dic" | type=="DIC" | type=="d" | type=="D") {
      mf <- data.frame(model=models,DIC=100*(fit$model.fitting$dic-min(fit$model.fitting$dic))/min(fit$model.fitting$dic),
                       DIC=fit$model.fitting$dic)
      lab.type <- "Percentage difference in DIC in comparison to the 'best-fitting' model"
      lab.type1="DIC"
    }
  }
  
  # Finally do the plot
  if (is.null(exArgs$xlim)) {xlm <- range(pretty(mf[,2]))} else {xlm <- exArgs$xlim}
  if (is.null(exArgs$digits)) {digits <- 7} else {digits <- exArgs$digits}
  if (is.null(exArgs$nsmall)) {nsmall <- 3} else {nsmall <- exArgs$nsmall}
  if (is.null(exArgs$main)) {main <- paste0("Model comparison based on ",lab.type1)} else {main <- exArgs$main}
  if (is.null(exArgs$mar)) {mar <- c(4,6,3,1.3)} else {mar <- exArgs$mar}
  if (is.null(exArgs$cex.names)) {cex.names <- 0.8} else {cex.names <- exArgs$cex.names}
  par(mar=mar)                                         # Bottom,left,top & right margins
  b <- barplot(                                        # Function to draw a barplot (see BMS NICE submission)
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
       format(mf[,3],digits=digits,nsmall=nsmall),       # Formats the values of the AICs/BICs/DICs, using 3 dp
       pos=4,                                            # Puts the text to the right of the bars
       cex=.8                                            # Rescales the labels on the y-axis to 80% of normal size
  )
}
