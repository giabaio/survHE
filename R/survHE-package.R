

#' Internal objects used by stan
#' 
#' Objects, C++ and S4 classes used by stan to fit parametric survival models.
#' 
#' 
#' @name internal_stan
#' @aliases internal_stan stanmodels stan_files model_Exponential model_Gamma
#' model_GenF model_GenGamma model_Gompertz model_logLogistic model_logNormal
#' model_PolyWeibull model_RP model_WeibullAF model_WeibullPH
#' Rcpp_model_Exponential-class Rcpp_model_Gamma-class Rcpp_model_GenF-class
#' Rcpp_model_GenGamma-class Rcpp_model_Gompertz-class
#' Rcpp_model_logLogistic-class Rcpp_model_logNormal-class
#' Rcpp_model_PolyWeibull-class Rcpp_model_RP-class Rcpp_model_WeibullAF-class
#' Rcpp_model_WeibullPH-class
#' @docType data
#' @keywords internal
NULL





#' Survival Analysis in Health Economic Evaluation
#' 
#' Contains a suite of functions to perform survival analysis with the aim of
#' aiding in health economic modelling (extrapolation, model checking and PSA)
#' 
#' \tabular{ll}{ Package: \tab survHE\cr Type: \tab Package\cr Version: \tab
#' 1.1.1\cr Date: \tab 2020-10-01\cr License: \tab GPL2 \cr LazyLoad: \tab
#' yes\cr } Contains a suite of functions to perform survival analysis with the
#' aim of aiding in health economic modelling (extrapolation, model checking
#' and PSA)
#' 
#' @name survHE-package
#' @aliases survHE-package survHE
#' @docType package
#' @author Gianluca Baio
#' 
#' Maintainer: Gianluca Baio
#' @template refs
#' @keywords Survival Modelling Health Economic Evaluation
#' @examples
#' 
#' # Loads some survival data
#' data(bc)
#' # Fits a parametric model
#' m <- fit.models(formula=Surv(recyrs,censrec)~group,data=bc,
#'     distr="exp",method="mle")
#' # Print output in tabular format
#' print(m)
#' # Visualise output in terms of survival curves
#' plot(m)
NULL



