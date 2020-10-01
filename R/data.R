#' A fictional survival trial.
#'
#' A dataset containing fictional data from a trial, where
#' the main outcome is in terms of time-to-event and 
#' censoring indicator and with additional covariates.
#'
#' @format A data frame with 367 rows and 8 variables:
#' \describe{
#'   \item{ID_patient}{The individual level identifier}
#'   \item{time}{The observed time at which the event happens}
#'   \item{censored}{An indicator to describe whether the 
#'   event is fully observed or censored}
#'   \item{arm}{An indicator for the treatment arm, with
#'   0 = control and 1 = active treatment}
#'   \item{sex}{An indicator for the individual's sex, with
#'   0 = male and 1 = female}
#'   \item{age}{A numeric variable with the individual's age}
#'   \item{imd}{A categorical variable representing a measure
#'   of area-level social deprivation}
#'   \item{ethnic}{A categorical variable representing the 
#'   individual's ethnic group, as measured from a Census}
#' }
"data"