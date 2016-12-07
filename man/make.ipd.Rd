\name{make.ipd}
\alias{make.ipd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create an individual level dataset from digitised data
}
\description{
Piles in the simulated IPD resulting from running digitise for more than one 
treatment arm 
}
\usage{
make.ipd(ipd_files,ctr=1,var.labs=c("time","event","arm"))
}
\arguments{
\item{ipd_files}{
a list including the names of the IPD files created as output of digitise
}
\item{ctr}{
the index of the file associated with the control arm (default, the first file).
This will be coded as 0
}
\item{var.labs}{
a vector of labels for the column of the resulting data matrix. NB these
should match the arguments to the formula specified for fit.models. The
user can specify values. These should be 4 elements (ID, TIME, EVENT, ARM)
}
}
\references{
Something will go here
}
\author{
Gianluca Baio
}
\seealso{
Something will go here
}
\examples{
# Something will go here
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Digitized Kaplan Meier curve}

