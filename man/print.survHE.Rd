\name{print.survHE}
\alias{print.survHE}
\title{
Print a summary of the survival model(s) fitted by \code{fit.models}
}
\description{
Prints the summary table for the model(s) fitted, with the estimate of the 
parameters
}
\usage{
\method{print}{survHE}(x,mod=1,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
the \code{survHE} object (the output of the call to \code{fit.models})
}
\item{mod}{
is the index of the model. Default value is 1, but the user can choose which model 
fit to visualise, if the call to fit.models has a vector argument for distr (so many 
models are fitted & stored in the same object)
}
\item{\dots}{
additional options, including:
\code{digits} = number of significant digits to be shown in the summary table 
(default = 6)
\code{nsim} = number of simulations from the joint posterior for \code{INLA} 
(default = 100)
\code{original} = a flag to say whether the *original* table from either \code{flexsurv}
or \code{INLA} or \code{rstan} should be printed
}
}
\references{
Something will go here
}
\author{
Gianluca Baio
}
\examples{
data(bc)

mle = fit.models(formula=Surv(recyrs,censrec)~group,data=bc,
    distr="exp",method="mle")
print(mle)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Parametric survival models}

