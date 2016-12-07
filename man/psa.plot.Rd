\name{psa.plot}
\alias{psa.plot}
\title{
Graphical depiction of the probabilistic sensitivity analysis for the survival curves 
}
\description{
Plots the survival curves for all the PSA simulations
}
\usage{
psa.plot(psa,...)
}
\arguments{
\item{psa}{
the result of the call to the function \code{make.surv}
}
\item{...}{
Optional graphical parameters, such as:
\code{xlab} = label for the x-axis
\code{ylab} = label for the y-axis
\code{col} = (vector) of colors for the lines to be plotted
\code{alpha} = the level of transparency for the curves (default = 0.1)
}
}
\author{
Gianluca Baio
}
\examples{
data(bc)

# Fits the same model using the 3 inference methods
mle = fit.models(formula=Surv(recyrs,censrec)~group,data=bc,
    distr="exp",method="mle")
p.mle = make.surv(mle,nsim=100)
psa.plot(p.mle)
}
\keyword{Survival models}
\keyword{Bootstrap}
\keyword{Probabilistic sensitivity analysis}
