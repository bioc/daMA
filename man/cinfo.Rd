\name{cinfo}
\alias{cinfo}
\docType{data}
\title{Vector indexing the matrix cmat}
\description{
This vector is used to describe the structure of the rows of the contrast
matrix cmat. The number of entries in cinfo mirrors the number of experimental
questions. "1" indicates a contrast in vectorial form, integers n > 1 indicate \eqn{n}
contrasts given in matrix form.
}
\usage{Such a vector has to be created for every contrast martix. sum(cinfo) has to be 
equal to nrow(cmat). }

\format{}
\details{}
\source{}
\references{
Bretz, F and Landgrebe J and Brunner E (2003):"Design and analysis of
        two colour factorial microarray experiments", submitted.
}
\examples{
data(cinfo)
}
\keyword{datasets}
