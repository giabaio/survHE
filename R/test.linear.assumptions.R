
#' Linearly transformed survival plots
#' 
#' Tests the linear assumptions for the parametric model
#' 
#' @param fit an object of class survHE
#' @param mod index or name of a model in fit. Defaults to 1.
#' @param label_plot if \code{TRUE}, labels assumptions. Defaults to \code{FALSE}.
#' @param graph Type of plot: base or ggplot2.
#' @param \dots further arguments, passed on to plot.
#' @return Diagnostic plot
#' @author William Browne, Nathan Green
#' @keywords survival hplot
#' @export
#' @examples 
#' 
#' data(bc)
#' 
#' # exponential distribution
#' fit_exp <- fit.models(formula = Surv(recyrs, censrec) ~ group, data = bc,
#'                   distr = "exp", method = "mle")
#' test.linear.assumptions(fit_exp)
#' test.linear.assumptions(fit_exp, graph = "ggplot2")
#'                  
#' # weibull distribution
#' fit_wei <- fit.models(formula = Surv(recyrs, censrec) ~ group, data = bc,
#'                   distr = "weibull", method = "mle")
#' test.linear.assumptions(fit_wei)
#' test.linear.assumptions(fit_wei, graph = "ggplot2")
#'                  
#' # loglogistic distribution
#' fit_llog <- fit.models(formula = Surv(recyrs, censrec) ~ group, data = bc,
#'                   distr = "loglogistic", method = "mle")
#' test.linear.assumptions(fit_llog)
#' test.linear.assumptions(fit_llog, graph = "ggplot2")
#'                  
test.linear.assumptions <- function(fit, mod = 1, label_plot = FALSE,
                                    graph = "base", ...) {
  
  dots <- list(...)
  
  graph <- match.arg(graph, c("base", "ggplot2"))
  
  if (length(mod) != 1)
    stop("")
  
  if (!is.numeric(mod) && mod <= length(fit$models))
    stop("")
  
  if (!is.character(mod) && mod %in% names(fit$models))
    stop("")
  
  dist <- get_distribution(fit, mod)
  
  distn_names <- list(
    "exp" = c("exp", "exponential"),
    "weibull" = c("weibull", "weibull.quiet", "weibullaf", "weibullph"),
    "loglogistic" = c("llogis", "loglogistic"),
    "lognormal" = c("lognormal", "lnorm"),
    "gompertz" = "gompertz")
  
  if (!dist %in% unname(unlist(distn_names)))
    stop("Distribution not available.")
  
  fit_km <- fit$misc$km
  
  n_strata <- length(fit_km$strata)
  
  model_strata <- rep(x = names(fit_km$strata),
                      times = fit_km$strata)
  
  times <- split(fit_km$time, model_strata)
  survs <- split(fit_km$surv, model_strata)
  
  params <- list()
  
  if (dist %in% distn_names[["exp"]]) {
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
  
  if (dist %in% distn_names[["weibull"]]) {
    params <- list(
      FUN = "lines",
      xlab = "log(time)",
      ylab = "log(-log(S(t))) = log cumulative hazard",
      main = "Weibull distributional assumption",
      x = lapply(times, log),
      y = lapply(survs, function(x) log(-log(x))),
      lty = 1:n_strata,
      col = 1:n_strata,
      type = "l")
  }
  
  if (dist %in% distn_names[["loglogistic"]]) {
    params <- list(
      FUN = "lines",
      xlab = "time",
      ylab = "log(S(t)/(1-S(t)))",
      main = "log-Logistic distributional assumption",
      x = lapply(times, log),
      y = lapply(survs, function(x) log(x/(1 - x))),
      lty = 1:n_strata,
      col = 1:n_strata,
      type = "l")
  }
  
  if (dist %in% distn_names[["lognormal"]]) {
    params <- list(
      FUN = "lines",
      xlab = "time",
      ylab = expression(Phi^-1 ~ (1 - S(t))),
      main = "lognormal distributional assumption",
      x = times,
      y = lapply(survs, function(x) qnorm(1 - x)),
      lty = 1:n_strata,
      col = 1:n_strata,
      type = "l")
  }
  
  if (dist %in% distn_names[["gompertz"]]) {
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
  
  if (graph == "base") {
    
    # empty plot
    do.call(plot, setup_pars)
    
    axis(1); axis(2)
    
    # plot lines
    do.call(mapply, modifyList(params, dots))
    
    if (isTRUE(label_plot)) {
      legend("topright", names(survs), col = params$col,
             lty = params$lty, bty = "n")
    }
  }
  
  if (graph == "ggplot2") {
    
    ggdata <- 
      data.frame(time = unlist(params$x),
                 y = unlist(params$y)) |>
      tibble::rownames_to_column("Group") |> 
      mutate(Group = gsub("\\d+", "", Group))
    
    p <- 
      ggplot(ggdata, aes(x = time, y = y, group = Group, col = Group)) +
      geom_line()
    
    print(p)
  }
  
  invisible(params)
}


#'
get_distribution <- function(fit, mod) {
  m <- fit$models[[mod]]
  tolower(ifelse(fit$method == "hmc", m@model_name, m$dlist$name))
}

