\name{test.linear.assumptions}
\alias{test.linear.assumptions}

\title{
Tests the linear assumptions for the parametric model
}
\description{
Tests the linear assumptions for the parametric model
}
\usage{
test.linear.assumptions(fit, mod = 1, label = FALSE, ...)
}

\arguments{
  \item{fit}{
an object of class survHE
}
  \item{mod}{
index or name of a model in fit. Defaults to 1.
}
  \item{label}{
if TRUE, labels assumptions. Defaults to FALSE.
}
  \item{\dots}{
further arguments, passed on to points()
}
}
\value{
A diagnostic plot
}
\author{
William Browne
}
\keyword{survival}
\keyword{hplot}
