\name{mean.surv}
\alias{mean.surv}
\title{
Computes the mean survival time for a given model and data
}
\description{
Calculates the mean survival time as the area under the survival curve
}
\usage{
mean.surv(x,mod=1,t=NULL,nsim=1000,stats=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
a \code{survHE} object (resulting from the call to \code{fit.models}}
\item{mod}{
the model to be analysed (default = 1)
}
\item{t}{the vector of times to be used in the computation. Default = NULL, 
which means the observed times will be used. NB: the vector of times should be: 
i) long enough so that S(t) goes to 0; and ii) dense enough so that
the approximation to the AUC is sufficiently precise}
\item{nsim}{the number of simulations from the survival curve distributions to be used 
(to compute interval estimates)}
\item{stats}{a logical value. If TRUE, also print a table with summary statistics from the 
distribution of the mean survival times}
\item{\dots}{
Additional options}
}
\details{
A list comprising of the following elements
}
\value{
\item{mean.surv}{
A matrix with the simulated values for the mean survival times
}
\item{tab}{
A summary table
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
\keyword{Mean survival time}

