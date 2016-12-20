make.transition.probs <- function(x,...) {
  # Computes the transition probabilities (to feed a discrete-time Markov model), based on the output of make.surv
  # Uses the formula p(t)=1-S(t+k)/S(t) where k is the MM cycle length and t is a generic time
  # x = an object obtained as output of the call to make.surv
  # ... = additional arguments. 
  #       labs = a string vector of names for the elements of the list (strata for the survival analysis)
  
  exArgs <- list(...)
  
  tp <- lapply(1:length(x$mat),function(i) {
    matrix(unlist(lapply(1:nrow(x$mat[[1]]),function(j) {
      1-(x$mat[[i]][j,2:ncol(x$mat[[1]])]/x$mat[[i]][j,1:ncol(x$mat[[1]])-1])
    })),byrow=T,nrow=nrow(x$mat[[1]]))
  })
  
  if (exists("labs",where=exArgs)) {
    if(length(exArgs$labs)==length(x$mat)) {
      names(tp) <- exArgs$labs
    }
  } else {
    names(tp) <- rownames(x$des.mat)
  }
  
  # Defines the column labels (to specify what the times refer to)
  col.labs <- character()
  for (i in 1:length(x$times)-1) {
    col.labs[i] <- paste0(x$times[i],"-",x$times[i+1])
  }
  for (i in 1:length(tp)) {
    colnames(tp[[i]]) <- col.labs
  }
  
  # Output
  # tpa: list of matrices with nsim rows and length(time) columns with the simulations for the transition probabilities between
  #      consecutive times.
  return(tp)
}
