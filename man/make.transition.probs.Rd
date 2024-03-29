% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make.transition.probs.R
\name{make.transition.probs}
\alias{make.transition.probs}
\title{make.transition.probs}
\usage{
make.transition.probs(fit, labs = NULL, ...)
}
\arguments{
\item{fit}{an object obtained as output of the call to \code{fit.models}}

\item{labs}{a vector with labels to identify the 'profiles' ie the
combination of covariates that have been passed onto the model formula.
If 'NULL' (default), then figures it out from the 'survHE' object.}

\item{...}{additional arguments. Includes the standard inputs to the
call to \code{make.surv}, so \code{mod} (the index of the possibly many
models stored in the 'survHE' object), \code{t} (the vector of times
over which to compute the survival curves), \code{newdata} (a list that
defines the profile of covariates) and \code{nsim} (the number of
simulations to use - default is \code{nsim}=1)}
}
\value{
A tibble 'lambda' with an indicator for the treatment arm,
the times at which the probabilities have been computed and \code{nsim}
columns each with a simulation of the transition probabilities for
all the times specified by the user
}
\description{
Computes the transition probabilities (to be passed to a Markov model) from
the cumulative hazard curves obtained using \code{fit.models}, using the formula
p(t)=1-exp(H(t-k)/H(t)), where k is the Markov model cycle length (or the
difference across two consecutive times) and t is a generic time
}
\note{
Something will go here
}
\examples{
\dontrun{
# Something will go here
}

}
\references{
Something will go here
}
\seealso{
\code{\link{make.surv}}
}
\author{
Gianluca Baio
}
\keyword{Markov}
\keyword{Transition}
\keyword{models}
\keyword{probabilities}
