
#' Plot to assess suitability of parametric model
#' 
#' Perform an exploratory investigation for linearity of
#' transformed survival models.
#' 
#' For the Weibull, twice taking logs of the survivor function
#' 
#' \deqn{log(-log S(t)) = log \lambda + \gamma log t}
#' 
#' A plot of \eqn{log(-log S(t))} against \eqn{log(t)} would give an approximately
#' straight line if the Weibull assumption is reasonable.
#' The plot could also be used to give a rough estimate of the parameters.
#' 
#' Similarly, for the log-logistic distribution
#' 
#' \deqn{logS(t)/(1 - S(t)) = \theta - \kappa log t}
#' 
#' For the log-normal distribution
#' 
#' \deqn{\Phi^{-1} (1 - S(t)) = (log t - \mu) / \sigma}
#' 
#' We can also check the assumption made with using the Cox regression model
#' of proportional hazards by inspecting the log-cumulative hazard plot.
#' 
#' \deqn{log H_i(t) = \beta x_i + log H_0(t)}
#' 
#' The transformed curves for different values of the explanatory variables
#' will be parallel if PH holds.
#' 
#' @param fit An object of class survHE.
#' @param mod Index or name of a model in fit. Defaults to 1.
#' @param add_legend If \code{TRUE}, labels assumptions. Defaults to \code{FALSE}.
#' @param graph Type of plot: base or ggplot2.
#' @param \dots Further arguments, passed on to plot.
#' @return Diagnostic plot
#' @author William Browne, Nathan Green
#' @references Collett (2015) Modelling Survival Data in Medical Research, CRC Press
#' @keywords survival hplot
#' @export
#' @examples 
#' 
#' data(bc)
#' form <- formula("Surv(recyrs, censrec) ~ group")
#' 
#' # exponential distribution
#' fit_exp <- fit.models(form, data = bc,
#'                   distr = "exp", method = "mle")
#' plot_transformed_km(fit_exp)
#' plot_transformed_km(fit_exp, graph = "ggplot2")
#'                  
#' # weibull distribution
#' fit_wei <- fit.models(form, data = bc,
#'                   distr = "weibull", method = "mle")
#' plot_transformed_km(fit_wei)
#' plot_transformed_km(fit_wei, graph = "ggplot2")
#'                  
#' # loglogistic distribution
#' fit_llog <- fit.models(form, data = bc,
#'                   distr = "loglogistic", method = "mle")
#' plot_transformed_km(fit_llog)
#' plot_transformed_km(fit_llog, graph = "ggplot2")
#'                  
#' # lognormal distribution
#' fit_lnorm <- fit.models(form, data = bc,
#'                   distr = "lognormal", method = "mle")
#' plot_transformed_km(fit_lnorm)
#' plot_transformed_km(fit_lnorm, graph = "ggplot2")
#'                  
plot_transformed_km <- function(fit, mod = 1, add_legend = FALSE,
                                graph = "base", ...) {
  
  dots <- list(...)
  
  graph <- match.arg(graph, c("base", "ggplot2"))
  
  if (length(mod) != 1)
    stop("Please provide at most one model index.")
  
  if (is.numeric(mod) && !mod <= length(fit$models))
    stop("More model names provided than available in list of model fits provided.")
  
  if (is.character(mod) && !mod %in% names(fit$models))
    stop("Model name not available in list of model fits provided.")
  
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
      ylab = "log(-log(S(t))) i.e. log cumulative hazard",
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
      xlab = "log(time)",
      ylab = expression(Phi^-1 ~ (1 - S(t))),
      main = "Log-normal distributional assumption",
      x = lapply(times, log),
      y = lapply(survs, function(x) qnorm(1 - x)),
      lty = 1:n_strata,
      col = 1:n_strata,
      type = "l")
  }
  
  if (dist %in% distn_names[["gompertz"]]) {
    stop("Gompertz not yet implemented.")
         
    #   estimate.h <- function(s, t) {
    #     denom <- t - c(t[-1], max(t) + 1)
    #     numerator <- log(s) - log(c(s[-1], 0))
    #     -numerator/denom
    #     }
    
    # params <- list(
    #   x = lapply(times, log),
    #   y = estimate.h(survs, times),
    #   xlab = "log(time)",
    #   ylab = "h(t)",
    #   main = "Gompertz distributional assumption",
    # lty = 1:n_strata,
    # col = 1:n_strata,
    # type = "l")
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
    
    if (isTRUE(add_legend)) {
      legend("topright", names(survs), col = params$col,
             lty = params$lty, bty = "n")
    }
  }
  
  if (graph == "ggplot2") {
    
    if (!add_legend) {
      pos.legend <- "none"
    }
    
    ggdata <- 
      data.frame(time = unlist(params$x),
                 y = unlist(params$y)) |>
      tibble::rownames_to_column("Group") |> 
      mutate(Group = gsub("\\d+", "", Group))
    
    p <- 
      ggplot(ggdata, aes(x = .data$time, y = .data$y,
                         group = .data$Group, col = .data$Group)) +
      geom_line() +
      do.call(labs,
              list(title = setup_pars$main,
                   x = setup_pars$xlab,
                   y = setup_pars$ylab)) +
      theme_bw() +
      theme(legend.position = pos.legend)
    
    print(p)
  }
  
  invisible(params)
}


#'
get_distribution <- function(fit, mod) {
  m <- fit$models[[mod]]
  tolower(ifelse(fit$method == "hmc", m@model_name, m$dlist$name))
}

