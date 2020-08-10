#' Helper function to provide a list of models available in each method
#' 
#' @return \item{availables}{A list of models available in each method}.
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso fit.models
#' @references Baio (2020). survHE
#' @keywords Parametric survival models Bayesian inference via Hamiltonian
#' Monte Carlo Bayesian inference via Integrated Nested Laplace Approximation
load_availables <- function() {
  # INLA can only do a limited set of models (for now) so if user has selected
  # one that is not available, then falls back on MLE analysis
  availables=list(
    mle=c("genf" = "gef",
          "genf.orig" = "gof",
          "gengamma" = "gga",
          "gengamma.orig" = "ggo",
          "exp" = "exp",
          "weibull" = "wei",
          "weibullPH" = "wph",
          "lnorm" = "lno",
          "gamma" = "gam",
          "gompertz" = "gom",
          "llogis" = "llo",
          "lognormal" = "lno",
          "rps" = "rps"
    ),
    inla=c("exponential" = "exp",
           "weibull" = "wei",
           "weibullPH" = "wph",
           "lognormal" = "lno",
           "loglogistic" = "llo",
           "rps" = "rps"
    ),
    hmc=c("Exponential" = "exp",
          "Gamma" = "gam",
          "GenF" = "gef",
          "GenGamma" = "gga",
          "Gompertz" = "gom",
          "PolyWeibull" = "pow",
          "RP" = "rps",
          "WeibullAF" = "wei",
          "WeibullPH" = "wph",
          "logLogistic" = "llo",
          "logNormal" = "lno"
    )
  )
  return(availables)
}