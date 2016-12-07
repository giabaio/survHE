\name{poly.weibull}
\alias{poly.weibull}
\title{
Fit Poly-Weibull model for survival analysis of mixture hazards
}
\description{
Runs the survival analysis using a Poly-Weibull model
}
\usage{
poly.weibull(formula=NULL,data,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{formula}{
a list of formulae (one for each components of the mixture. Can specify one single formula
(in which case, the model is a simple Weibull regression). For example, a valid call is using
\code{formula=list(Surv(time,event)~1,Surv(time,event)~arm)}
}
\item{data}{
A data frame containing the data to be used for the analysis. This must contain data for
the 'event' variable. In case there is no censoring, then \code{event} is a column of 1s.
}
\item{\dots}{
Additional options (for INLA or HMC). 

**HMC** specific options
\code{chains} = number of chains to run in the HMC (default = 2)
\code{iter} = total number of iterations (default = 2000)
\code{warmup} = number of warmup iterations (default = iter/2)
\code{thin} = number of thinning (default = 1)
\code{control} = a list specifying Stan-related options, eg 
\code{control=list(adapt_delta=0.85)} (default = NULL)
\code{seed} = the random seed (to make things replicable)
\code{pars} = a vector of parameters (string, default = NA)
\code{include} = a logical indicator (if FALSE, then the pars are not saved; 
default = TRUE)
\code{priors} = a list (of lists) specifying the values for the parameters of the 
prior distributions in the models
\code{save.stan} = a logical indicator (default = FALSE). If TRUE, then saves the data 
list for Stan and the model file(s)}
}
\details{
On object in the class \code{survHE} containing the following elements
}
\value{
\item{models}{
A list containing the fitted models. These contain the output from the original 
inference engine (\code{flexsurv}, \code{INLA} or \code{rstan}). Can be processed
using the methods specific to the original packages, or via \code{survHE}-specific
methods (such as \code{plot}, \code{print}) or other specialised functions (eg to
extrapolate the survival curves, etc).
}
\item{model.fitting}{
A list containing the output of the model-fit statistics (AIC, BIC, DIC). The AIC and
BIC are estimated for all methods, while the DIC is only estimated when using Bayesian
inference.
}
\item{method}{
A string indicating the method used to fit the model, ie \code{'mle'}, \code{'inla'} or
\code{'hmc'}. 
}
\item{misc}{
A list containing the time needed to run the model(s) (in seconds), the formula used, the
results of the Kaplan-Meier analysis (which is automatically performed using \code{npsurv})
and the original data frame.
}
}
\references{
Something will go here
}
\author{
Gianluca Baio
}
\note{
Something will go here
}
\seealso{
Something will go here
}
\examples{
###
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Parametric survival models}
\keyword{Bayesian inference via Hamiltonian Monte Carlo}
\keyword{Poly-Weibull model}

