#' Graphical depiction of the probabilistic sensitivity analysis for the
#' survival curves
#' 
#' Plots the survival curves for all the PSA simulations
#' 
#' 
#' @param psa the result of the call to the function \code{make.surv}
#' @param ...  Optional graphical parameters, such as: \code{xlab} = label for
#' the x-axis \code{ylab} = label for the y-axis \code{col} = (vector) of
#' colors for the lines to be plotted \code{alpha} = the level of transparency
#' for the curves (default = 0.1)
#' @author Gianluca Baio
#' @keywords Survival models Bootstrap Probabilistic sensitivity analysis
#' @examples
#' 
#' data(bc)
#' 
#' # Fits the same model using the 3 inference methods
#' mle = fit.models(formula=Surv(recyrs,censrec)~group,data=bc,
#'     distr="exp",method="mle")
#' p.mle = make.surv(mle,nsim=100)
#' psa.plot(p.mle)
#' 
#' @export psa.plot
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
  # xlim = a vector of limits for the times
  # scale.time = a factor by which to rescale the time variable (eg 1/52 rescales weeks to years)
  # round.fact = the number of decimal places to include in the x-axis for the time variable (when rescaling!)
  
  n.elements <- length(psa$S[[1]]) 
  times <- psa$S[[1]][[1]][,1]
  exArgs <- list(...)
  if(!exists("xlab",where=exArgs)) {xlab <- "Time"} else {xlab <- exArgs$xlab}
  if(!exists("ylab",where=exArgs)) {ylab="Survival"} else {ylab <- exArgs$ylab}
  if(!exists("col",where=exArgs)) {col <- sample(colors(),n.elements)} else {col <- exArgs$col}
  if(!exists("alpha",where=exArgs)) {alpha <- 0.1} else {alpha <- exArgs$alpha}
  if(!exists("main",where=exArgs)) {main <- ""} else {main <- exArgs$main}
  if(!exists("labs",where=exArgs)) {labs <- TRUE} else {labs <- exArgs$labs}
  if(!exists("xpos",where=exArgs)) {xpos <- max(times)*0.65} else {xpos <- exArgs$xpos}
  if(!exists("ypos",where=exArgs)) {ypos <- 1} else {ypos <- exArgs$ypos}
  if(!exists("cex.txt",where=exArgs)) {cex.txt <- 0.75} else {cex.txt <- exArgs$cex.txt}
  if(!exists("offset",where=exArgs)) {off <- seq(1,nrow(psa$des.mat))*.35} else {off <- seq(1,nrow(psa$des.mat))*exArgs$offset}
  if(!exists("nsmall",where=exArgs)) {nsmall <- 2} else {nsmall <- exArgs$nsmall}
  if(!exists("digits",where=exArgs)) {digits <- 5} else {digits <- exArgs$digits}
  if(exists("xlim",where=exArgs)) {xlim=exArgs$xlim} else {xlim <- range(pretty(times))}
  
  # If there's only the average value for the survival curve, simpler plot
  if (psa$nsim==1) {
    alpha <- 1
    plot(psa$S[[1]][[1]][,1:2],t="l",xlab=xlab,ylab=ylab,col=adjustcolor(col[1],alpha),ylim=c(0,1),xlim=xlim,
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
    plot(psa$S[[1]][[1]][,1],q500[[1]],col=adjustcolor(col[1],1),t="l",xlab=xlab,ylab=ylab,ylim=c(0,1),xlim=xlim,lwd=2,main=main,axes=F)
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
  # If the time has been re-scaled, then use the new scale as labels on the x-axis. Otherwise, use original scale
  if(exists("scale.time",where=exArgs)) {
      if(exists("round.fact",where=exArgs)) {round.fact=exArgs$round.fact} else {round.fact=1}
      axis(1,at=pretty(times),labels=round(pretty(times)*exArgs$scale.time,round.fact))
  } else {
      axis(1)
  }
  axis(2)
  if (class(labs)=="logical") {
    if(labs){ 
       txt1 <- lapply(1:ncol(psa$des.mat),function(i) {
          text(xpos,ypos-(i-1)/40,paste0(colnames(psa$des.mat)[i]," :  "),cex=cex.txt,pos=2,col="black")
       })
       txt2 <- lapply(1:ncol(psa$des.mat),function(i) {
          lapply(1:nrow(psa$des.mat),function(j) {
             text(xpos+off[j],ypos-(i-1)/40,format(psa$des.mat[j,i],nsmall=nsmall,digits=digits),cex=cex.txt,pos=2,col=col[j])
          })
       })
    }
  } else {
     legend("topright",legend=labs,col=col,bty="n",lty=1,cex=cex.txt)
  }
}
