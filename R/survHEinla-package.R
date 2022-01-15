#' Survival Analysis in Health Economic Evaluation using INLA
#' 
#' Contains a suite of functions to perform survival analysis with the aim of
#' aiding in health economic modelling (extrapolation, model checking and PSA)
#' 
#' \tabular{ll}{ Package: \tab survHEinla\cr Type: \tab Package\cr Version: \tab
#' 0.0.1\cr Date: \tab 2022-01-15\cr License: \tab GPL2 \cr LazyLoad: \tab
#' yes\cr } Module to expand the facilities of survHE to run Bayesian models using
#' INLA
#' 
#' @name survHEinla-package
#' @aliases survHEinla-package survHEinla
#' @docType package
#' @author Gianluca Baio
#' 
#' Maintainer: Gianluca Baio
#' @template refs
#' @keywords Survival Modelling Health Economic Evaluation using INLA
#' @examples
#' \dontrun{ 
#' # Loads some survival data
#' library(survHE)
#' data(bc)
#' # Fits a parametric model using INLA
#' m <- fit.models(formula=Surv(recyrs,censrec)~group,data=bc,
#'     distr="exp",method="inla")
#' # Print output in tabular format
#' print(m)
#' # Visualise output in terms of survival curves
#' plot(m)
#' }
NULL
