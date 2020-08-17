#' Helper function to manipulate the strings of text defining the 
#' distributions selected by the user so they are consistent with the
#' various methods
#' 
#' @param x A string with the distribution name selected by the user.
#' @return \item{list}{A list containing the modified name of the 
#' distribution, the acronym (3-letters abbreviation), or the
#' labels (humane-readable name)}.
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso fit.models
#' @references Baio (2020). survHE
#' @keywords Parametric survival models Bayesian inference via Hamiltonian
#' Monte Carlo Bayesian inference via Integrated Nested Laplace Approximation
manipulate_distributions <- function(x){
  # selected model checks -----
  matchTable = list(
    "exp" = c("exponential", "exp"),
    "wei" = c("weibull", "weibullaft", "weiaft", "waft", "weibullaf", "weiaf", "waf", "wei"),
    "wph" = c("weibullph", "weiph", "wph"),
    "gam" = c("gamma", "gam", "gma"),
    "lno" = c("lognormal", "lnormal", "lnorm", "lognorm", "lno"),
    "llo" = c("loglogistic", "loglog", "llogistic", "llogis", "llo", "llogist"),
    "gga" = c("generalisedgamma", "generalizedgamma", "ggamma", "gengamma", "gga", "ggam"),
    "ggo" = c("gengamma.orig", "ggo"),
    "gef" = c("generalisedf", "generalizedf", "genf", "gef"),
    "gof" = c("genf.orig", "gof"),
    "gom" = c("gompertz", "gpz", "gomp", "gompz", "gom"),
    "rps" = c("roystonparmar", "roystonparmarsplines", "roystonparmarspline", "spline", "splines", "rps")
  )
  # Human readable label
  labelTable = c(
    "exp" = "Exponential",
    "wei" = "Weibull (AFT)",
    "wph" = "Weibull (PH)",
    "gam" = "Gamma",
    "lno" = "log-Normal", 
    "llo" = "log-Logistic",
    "gga" = "Gen. Gamma", "ggo" = "Gen. Gamma (orig parametrisation)",
    "gef" = "Gen. F", "gof" = "Gen. F (orig parametrisation)",
    "gom" = "Gompertz",
    "rps" = "Royston-Parmar",
    "pow" = "Poly-Weibull")
  # Labels used by R to define p..., r... and d... commands
  labelR = c(
    "exp" = "exp",
    "wei" = "weibull", 
    "wph" = "weibullPH", 
    "gam" = "gamma",
    "lno" = "lnorm", 
    "llo" = "llogis",
    "gga" = "gengamma", 
    "ggo" = "gengamma.orig",
    "gef" = "genf",
    "gof" = "genf.orig",
    "gom" = "gompertz",
    "rps" = "survspline"
  )
  
  distr = gsub("[ ]*[-]*", "", tolower(x))
  isDistrUnmatched = which(!sapply(
    1:length(distr),
    '%in%',
    unname(unlist(sapply(matchTable, match, distr)))))
  if (length(isDistrUnmatched) > 0) {
    stop(paste0("Distribution ", paste(distr[isDistrUnmatched], collapse = ", "), " could not be matched."))
  }
  
  distr3 <- numeric()
  for (i in 1:length(distr)) {
    distr3[i] <- names(which(unlist(lapply(matchTable,function(x) distr[i]%in%x))))  
  }
  labs <- unname(labelTable[distr3])
  distr <- unname(labelR[distr3])
  
  list(distr=distr,distr3=distr3,labs=labs)
}
