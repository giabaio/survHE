\name{digitise}
\alias{digitise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Format digitised data for use in survival analysis
}
\description{
Produces txt files with Kaplan Meier and individual level survival data from digitised 
Kaplan Meier curves obtained by DigitizeIT
}
\usage{
digitise(surv_inp, nrisk_inp, km_output = "KMdata.txt", 
ipd_output = "IPDdata.txt")
}
\arguments{
\item{surv_inp}{
a txt file obtained for example by DigitizeIT and containing the input survival times 
from graph reading
}
\item{nrisk_inp}{
a txt file obtained by DigitizeIT and containing the reported number at risk
}
\item{km_output}{
the name of the file to which the KM data will be written
}
\item{ipd_output}{
the name of the file to which the individual level data data will be written
}
}
\references{
Something will go here
}
\author{
Patricia Guyot and Gianluca Baio
}
\seealso{
}
\examples{
# Something will go here
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Digitized Kaplan Meier curve}

