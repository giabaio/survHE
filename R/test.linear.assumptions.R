test.linear.assumptions <- function(fit, mod = 1, label = FALSE, ...) {
  ## THIS IS INTERESTING, BUT NEEDS TO COMPLETE WITH THE OTHER DISTRIBUTIONS!!!
  
  m <- fit$models[[mod]]
  dist <- tolower(ifelse(fit$method == "hmc", m@model_name, m$dlist$name))
  split_vector <- c(1)
  for (i in 2:length(fit$misc$km$time)) {
    if (fit$misc$km$time[i] < fit$misc$km$time[i - 1]) {
      split_vector <- c(split_vector, i - 1, i)
    }
  }
  split_vector <- c(split_vector, length(fit$misc$km$time))
  split_mat <- matrix(split_vector,
                      length(split_vector)/2,
                      2,
                      byrow = TRUE)
  times <- lapply(1:dim(split_mat)[1],
                  function(x)
                    fit$misc$km$time[split_mat[x, ][1]:split_mat[x, ][2]])
  survs <- lapply(1:dim(split_mat)[1],
                  function(x)
                    fit$misc$km$surv[split_mat[x, ][1]:split_mat[x, ][2]])
  
  if (dist %in% c("exp", "exponential")) {
    xlab <- "time"
    ylab <- "log(S(t))"
    xlim <- range(pretty(fit$misc$km$time))
    ylim <- NULL
    legend_text <- "Exponential distributional assumption"
    pts <- lapply(1:dim(split_mat)[1],
                  function(m)
                    cbind(times[[m]],
                          log(survs[[m]])))
  }
  
  if (dist %in% c("weibull", "weibull.quiet", "weibullaf", "weibullph")) {
    xlab <- "log(time)"
    ylab <- "log(-log(S(t))) = log cumulative hazard"
    xlim <- range(pretty(log(fit$misc$km$time)))
    ylim <- range(pretty(log(-log(survs[[1]]))))
    legend_text <- "Weibull distributional assumption"
    pts <- lapply(1:dim(split_mat)[1],
                  function(m)
                    cbind(log(times[[m]]),
                          log(-log(survs[[m]]))))
  }
  
  
  if (dist %in% c("llogis", "loglogistic")) {
    xlab <- "time"
    ylab <- "log(S(t)/(1-S(t)))"
    xlim <- range(pretty(log(fit$misc$km$time))) 
    ylim <- range(pretty(log(survs[[1]]/(1 - survs[[1]]))))
    legend_text <- "log-Logistic distributional assumption"
    pts <- lapply(1:dim(split_mat)[1],
                  function(m)
                    cbind(log(times[[m]]),
                          log(survs[[m]]/(1 - survs[[m]]))))
  }
  
  if (dist %in% c("lognormal", "lnorm")) {
    xlab <- "time"
    ylab <- "log(S(t))"
    axes <- FALSE
    xlim <- range(pretty(log(fit$misc$km$time)))
    legend_text <- "lognormal distributional assumption"
    pts <- lapply(1:dim(split_mat)[1],
                  function(m)
                    cbind(times[[m]],
                          qnorm(1 - survs[[m]])))
  }
  
  if (dist == "gompertz") {
    
    estimate.h <- function(s, t) {
      denom <- t - c(t[-1], max(t) + 1)
      numerator <- log(s) - log(c(s[-1], 0))
      return(-numerator/denom)
    }
    
    xlab <- "log(time)"
    ylab <- "h(t)"
    xlim = range(pretty(fit$misc$km$time))
    ylim = range(pretty(estimate.h(survs[[1]], times[[1]])))
    legend_text <- "Gompertz distributional assumption"
    ### NEED TO CHECK --- WHAT IS V2???
    pts <- lapply(1:dim(split_mat)[1],
                  function(m)
                    cbind(times[[m]],
                          estimate.h(survs[[m]],
                                     times[[m]])))[V2 != 0, ]
  }
  
  # actual plot
  plot(x = 0,
       y = 0,
       type = "n",
       axes = FALSE,
       xlab = xlab,
       ylab = ylab,
       xlim = xlim,
       ylim = ylim)
  
  axis(1)
  axis(2)
  
  lapply(1:length(pts),
         function(x)
           points(pts[[x]],
                  t = "l",
                  lty = x,
                  ...))
  
  if (isTRUE(label)) {
    legend("topright", legend_text, bty = "n")
  }
  
  return(invisible())
}
