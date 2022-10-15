#' Linearly transformed survival plots
#' 
#' Tests the linear assumptions for the parametric model
#' 
#' 
#' @param fit an object of class survHE
#' @param mod index or name of a model in fit. Defaults to 1.
#' @param label_plot if TRUE, labels assumptions. Defaults to \code{FALSE}.
#' @param \dots further arguments, passed on to plot.
#' @return Diagnostic plot
#' @author William Browne
#' @keywords survival hplot
#' @export test.linear.assumptions
#' @noRd 
#' @examples 
#' 
#' data(bc)
#' 
#' # exponential distribution
#'
#' mle <- fit.models(formula = Surv(recyrs, censrec) ~ group, data = bc,
#'                   distr = "exp", method = "mle")
#' test.linear.assumptions(mle)
#'                  
test.linear.assumptions <- function(fit, mod = 1, label_plot = FALSE, ...) {
  
  dots <- list(...)
  
  if (length(mod) != 1)
    stop("")
  
  if (!is.numeric(mod) && mod <= length(fit$models))
    stop("")
  
  if (!is.character(mod) && mod %in% names(fit$models))
    stop("")
  
  m <- fit$models[[mod]]
  
  get_dist <- function(fit)
    tolower(ifelse(fit$method == "hmc", m@model_name, m$dlist$name))
  
  dist <- get_dist(fit)
  
  if (!dist %in% c("exp", "exponential", "weibull", "weibull.quiet", "weibullaf", "weibullph",
                   "llogis", "loglogistic", "lognormal", "lnorm", "gompertz"))
    stop("Distribution not available.")
  
  n_strata <- length(fit$misc$km$strata)
  
  model_strata <- rep(x = names(fit$misc$km$strata),
                      times = fit$misc$km$strata)
  
  times <- split(fit$misc$km$time,
                 model_strata)
  
  survs <- split(fit$misc$km$surv,
                 model_strata)
  
  params <- list()
  
  if (dist %in% c("exp", "exponential")) {
    params <- list(
      FUN = "lines",
      xlab = "time",
      ylab = "log(S(t))",
      main = "Exponential distributional assumption",
      x = times,
      y = lapply(survs, log),
      lty = 1:n_strata,
      col = 1:n_strata,
      type = "l")
  }
  
  if (dist %in% c("weibull", "weibull.quiet", "weibullaf", "weibullph")) {
    params <- list(
      FUN = "lines",
      xlab = "log(time)",
      ylab = "log(-log(S(t))) = log cumulative hazard",
      main = "Weibull distributional assumption",
      x = lapply(times, log),
      y = lapply(survs, function(x) log(-log(x))))
  }
  
  
  if (dist %in% c("llogis", "loglogistic")) {
    params <- list(
      FUN = "lines",
      xlab = "time",
      ylab = "log(S(t)/(1-S(t)))",
      main = "log-Logistic distributional assumption",
      x = lapply(times, log),
      y = lapply(survs, function(x) log(x/(1 - x))))
  }
  
  if (dist %in% c("lognormal", "lnorm")) {
    params <- list(
      FUN = "lines",
      xlab = "time",
      ylab = "log(S(t))",
      axes = FALSE,
      main = "lognormal distributional assumption",
      x = times,
      y = lapply(survs, function(x) qnorm(1 - x)))
  }
  
  if (dist == "gompertz") {
    # placeholder
    warning("Gompertz models are not yet implemented in test.linear.assumptions()")
    
    params <- list(
      x = 0,
      y = 0,
      # estimate.h <- function(s, t) {
      #   denom <- t - c(t[-1], max(t) + 1)
      #   numerator <- log(s) - log(c(s[-1], 0))
      #   return(-numerator/denom)
      # }
      
      xlab = "log(time)",
      ylab = "h(t)",
      main = "Gompertz distributional assumption")
    ### NEED TO CHECK --- WHAT IS V2???
    # pts <- lapply(1:dim(split_mat)[1],
    #               function(m)
    #                 cbind(times[[m]],
    #                       estimate.h(survs[[m]],
    #                                  times[[m]])))[V2 != 0, ]
  }
  
  default_pars <- list(
    x = NULL,
    type = "n",
    axes = FALSE,
    xlab = params$xlab,
    ylab = params$ylab,
    main = params$main,
    xlim = range(pretty(unlist(params$x))),
    ylim = range(pretty(unlist(params$y))))
  
  setup_pars <- modifyList(
    default_pars, dots[names(default_pars)])
  
  # empty plot
  do.call(plot, setup_pars)
  
  axis(1); axis(2)
  
  # plot lines
  do.call(mapply, modifyList(params, dots))
  
  # if (isTRUE(label_plot)) {
  #   legend("topright", legend_text, bty = "n")
  # }
}
