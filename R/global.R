#' This defines variables, mainly used with 'dplyr' and pipes that
#' throw a 'no visible binding for global variable' during 
#' 'R CMD check' (before CRAN submission)
#' 
#' @noRd 
utils::globalVariables(c(
  ".",
  "low",
  "upp",
  "time",
  "event",
  "(Intercept)",
  "S",
  "model_name",
  "strata",
  "object_name",
  "lower",
  "upper",
  "value",
  "lab",
  "aic",
  "bic",
  "dic",
  "mods_id",
  "xmin",
  "xmax",
  "ymin",
  "ymax",
  "where"
))