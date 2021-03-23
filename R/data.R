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


#' NICE TA174 dataset.
#'
#' A dataset containing the data used for NICE TA174, as
#' made publicly available as part of the supplementary
#' material for Williams et al (2017). Medical Decision
#' Making, 37;427-439.
#'
#' @format A tibble with 810 rows and 8 variables:
#' \describe{
#'   \item{patid}{A numeric patient identifier}
#'   \item{treat}{The treatment indicator. 1=rituximab 
#'   in combination with fludarabine andcyclophosphamide 
#'   (RFC); 0=fludarabine and cyclo-phosphamide alone (FC)}
#'   \item{prog}{An indicator to describe whether the 
#'   patient has experience a progression}
#'   \item{death}{An indicator to describe whether the 
#'   patient has experience death}
#'   \item{prog_t}{The observed time at progression, or
#'   the time at which the patient has been censored; 
#'   measured in months}
#'   \item{death_t}{The observed time at death, or
#'   the time at which the patient has been censored; 
#'   measured in months}
#'   \item{prog_ty}{The observed time at progression, or
#'   the time at which the patient has been censored; 
#'   measured in years}
#'   \item{death_ty}{The observed time at death, or
#'   the time at which the patient has been censored; 
#'   measured in years}
#'   }
"ta174"


#' NICE TA174 dataset in multi-state format.
#'
#' These are the same data contained in NICE TA174, as
#' made publicly available as part of the supplementary
#' material for Williams et al (2017). Medical Decision
#' Making, 37;427-439. However, the data have been 
#' restructured (by using the function \code{make_data_multi_state()})
#' to be used for multi-state analysis
#'
#' @format A tibble with 1868 rows and 16 variables:
#' \describe{
#'   \item{id}{A numeric patient identifier}
#'   \item{from}{An indicator of the starting state. 
#'   1=Pre-progression; 2=Progression; 3=Death}
#'   \item{to}{An indicator for the arriving state}
#'   \item{trans}{A code for the actual transition considered.
#'   1=Pre-progression -> Progression; 2=Pre-progression ->
#'   Death; 3=Progression -> Death}
#'   \item{Tstart}{The time of entry into the observation}
#'   \item{Tstop}{The time of exit from observation}
#'   \item{time}{The observed time until even (progression or
#'   death), or censoring occurs}
#'   \item{status}{The event indicator; takes value 1 if the
#'   underlying event (which varies depending on which 
#'   transition is being considered) happens and 0 otherwise}
#'   \item{treat}{The treatment indicator. 1=rituximab 
#'   in combination with fludarabine andcyclophosphamide 
#'   (RFC); 0=fludarabine and cyclo-phosphamide alone (FC)}
#'   \item{patid}{The original numeric patient identifier}
#'   \item{prog}{The original indicator to describe whether 
#'   the patient has experience a progression}
#'   \item{death}{The original indicator to describe whether 
#'   the patient has experience death}
#'   \item{prog_t}{The original observed time at progression, 
#'   or the time at which the patient has been censored; 
#'   measured in months}
#'   \item{death_t}{The original observed time at death, or
#'   the time at which the patient has been censored; 
#'   measured in months}
#'   \item{prog_ty}{The original observed time at progression, 
#'   or the time at which the patient has been censored; 
#'   measured in years}
#'   \item{death_ty}{The original observed time at death, or
#'   the time at which the patient has been censored; 
#'   measured in years}
#'   }
"msmdata"