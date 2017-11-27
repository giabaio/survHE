#' Tests the linear assumptions for the parametric model
#' 
#' Tests the linear assumptions for the parametric model
#' 
#' 
#' @param fit an object of class survHE
#' @param mod index or name of a model in fit. Defaults to 1.
#' @param label_plot if TRUE, labels assumptions. Defaults to FALSE.
#' @param \dots further arguments, passed on to points()
#' @return A diagnostic plot
#' @author William Browne
#' @keywords survival hplot
#' @export test.linear.assumptions
test.linear.assumptions <- function(fit, mod = 1, label_plot = FALSE, ...) {
  ## THIS IS INTERESTING, BUT NEEDS TO COMPLETE WITH THE OTHER DISTRIBUTIONS!!!
  
  stopifnot(length(mod) == 1)
  if(is.numeric(mod)) stopifnot(mod <= length(fit$models))
  if(is.character(mod)) stopifnot(mod %in% names(fit$models))
  
  m <- fit$models[[mod]]
  dist <- tolower(ifelse(fit$method == "hmc", m@model_name, m$dlist$name))
  
  stopifnot(dist %in% c("exp", "exponential", "weibull", "weibull.quiet", "weibullaf", "weibullph",
                        "llogis", "loglogistic", "lognormal", "lnorm", "gompertz"))
  
  model_strata <- rep(x = names(fit$misc$km$strata),
                      times = fit$misc$km$strata)

  times <- tapply(fit$misc$km$time,
                  model_strata,
                  I)
  survs <- tapply(fit$misc$km$surv,
                  model_strata,
                  I)

  
  if (dist %in% c("exp", "exponential")) {
    xlab <- "time"
    ylab <- "log(S(t))"
    legend_text <- "Exponential distributional assumption"
    x <- times
    y <- lapply(survs, log)
  }
  
  if (dist %in% c("weibull", "weibull.quiet", "weibullaf", "weibullph")) {
    xlab <- "log(time)"
    ylab <- "log(-log(S(t))) = log cumulative hazard"
    legend_text <- "Weibull distributional assumption"
    x <- lapply(times, log)
    y <- lapply(survs, function(x) log(-log(x)))
  }
  
  
  if (dist %in% c("llogis", "loglogistic")) {
    xlab <- "time"
    ylab <- "log(S(t)/(1-S(t)))"
    legend_text <- "log-Logistic distributional assumption"
    x <- lapply(times, log)
    y <- lapply(survs, function(x) log(x) / (1 - log(x)))
  }
  
  if (dist %in% c("lognormal", "lnorm")) {
    xlab <- "time"
    ylab <- "log(S(t))"
    axes <- FALSE
    legend_text <- "lognormal distributional assumption"
    x <- times
    y <- lapply(survs, function(x) qnorm(1 - x))
  }
  
  if (dist == "gompertz") {
    # placeholder
    warning("Gompertz models are not yet implemented in test.linear.assumptions()")
    x <- 0
    y <- 0
    # estimate.h <- function(s, t) {
    #   denom <- t - c(t[-1], max(t) + 1)
    #   numerator <- log(s) - log(c(s[-1], 0))
    #   return(-numerator/denom)
    # }
    
    xlab <- "log(time)"
    ylab <- "h(t)"
    legend_text <- "Gompertz distributional assumption"
    ### NEED TO CHECK --- WHAT IS V2???
    # pts <- lapply(1:dim(split_mat)[1],
    #               function(m)
    #                 cbind(times[[m]],
    #                       estimate.h(survs[[m]],
    #                                  times[[m]])))[V2 != 0, ]
  }
  
  # actual plot
  pts <- mapply(cbind, x, y)

  plot(x = 0,
       y = 0,
       type = "n",
       axes = FALSE,
       xlab = xlab,
       ylab = ylab,
       xlim = range(pretty(unlist(x))),
       ylim = range(pretty(unlist(y))))
  
  axis(1)
  axis(2)
  
  lapply(1:length(pts),
         function(x)
           points(pts[[x]],
                  t = "l",
                  lty = x,
                  ...))
  
  if (isTRUE(label_plot)) {
    legend("topright", legend_text, bty = "n")
  }
  
  return(invisible())
}
