\name{cmatB.AB}
\alias{cmatB.AB}
\docType{data}
\title{ Contrast matrix describing the experimental questions}
\description{
This matrix of numerical constants describes the experimental question,  say p. Each experimental question is
described by a single contrast vector (a single row in cmat) or by a contrast matrix
(several rows in cmat). The ordering of the columns corresponds to that of the associated
design matrix X. Thus, typically the first two elements in a row of cmat are reserved for
for the two dyes. E.g. the first line of the matrix cmatB.AB describes the main effect B.
}
\usage{data(cmatB.AB)}
\format{
}
\details{
}
\source{
}
\references{
Bretz, F and Landgrebe J and Brunner E (2003):"Design and analysis of
        two colour factorial microarray experiments", submitted.
}
\examples{
data(cmatB.AB)
}
\keyword{datasets}
