#' Survival Analysis in Health Economic Evaluation
#' 
#' Contains a suite of functions to perform survival analysis with the aim of
#' aiding in health economic modelling (extrapolation, model checking and PSA)
#' 
#' \tabular{ll}{ Package: \tab survHE\cr Type: \tab Package\cr Version: \tab
#' 2.0.1cr Date: \tab 2022-11-10\cr License: \tab GPL2 \cr LazyLoad: \tab
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
#' \dontrun{ 
#' # Loads some survival data
#' data(bc)
#' # Fits a parametric model
#' m <- fit.models(formula=Surv(recyrs,censrec)~group,data=bc,
#'     distr="exp",method="mle")
#' # Print output in tabular format
#' print(m)
#' # Visualise output in terms of survival curves
#' plot(m)
#' }
NULL



