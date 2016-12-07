\name{make.transition.probs}
\alias{make.transition.probs }
\title{
make.transition.probs 
}
\description{
Computes the transition probabilities (to be passed to a Markov model) from the survival
curves obtained using \code{fit.models} and \code{make.surv}, using the formula 
p(t)=1-S(t+k)/S(t), where k is the Markov model cycle length and t is a generic time
}
\usage{
make.transition.probs(x,...)
}
\arguments{
\item{x}{
an object obtained as output of the call to \code{make.surv}
}
\item{...}{
additional arguments. Includes labs = a string vector of names for the elements of 
the list (strata for the survival analysis)
}
}
\details{
Something will go here
}
\value{
Something will go here
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
# Something will go here
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Transition probabilities}
\keyword{Markov models}

