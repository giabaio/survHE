test.linear.assumptions <- function(fit,mod=1,coxph=TRUE,label=FALSE,...){
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
    pts <- lapply(1:dim(split_mat)[1],function(m) cbind(times[[m]],estimate.h(survs[[m]],times[[m]])))[V2!=0,]
    lapply(1:length(pts), function(x) points(pts[[x]],t="l",lty=x))
    if (label){legend('topright','Gompertz distributional assumption',bty='n')}
    
  }
}
