\name{core}
\alias{core}
\title{Internal function of analyseMA}
\description{
This internal function of analyseMA computes the statistics and estimators
that are organised and given out by the main function analyseMA.
}
\usage{
core(vector, design, cmat, cinfo, tol)
}
\arguments{
  \item{vector}{a simple help variable for the apply call}
  \item{design}{ the design matrix of size \eqn{N \times (K+2)}, where K is the number of
experimental conditions. This is the design matrix X known from linear model theory and
its elements are typically 0, 1, or -1. A 0 means that the associated parameter does not
apply for the corresponding observation (i.e., row). The first two columns are reserved
for the two dyes and are usually filled up with 1 and -1, respectively. }
  \item{cmat}{a matrix describing the p experimental questions (contrasts) to be analysed in
	the experiment. The matrix can be composed of vectorial contrasts (a single row of the
	matrix) and of contrasts in matrix form (several rows of the matrix),
	e.g. an \eqn{A \times B} interaction effect in a \eqn{3 \times 2} design. All
	contrasts have to be combined into one matrix (using rbind for instance). }
  \item{cinfo}{ a vector of length p describing the grouping of the contrast matrix rows in vector or matrix 
	form. E.g. if the design matrix contains three contrasts in vector form, cinfo =
	rep(1,3), if it contains two vectorial contratst and one as matrix with three rows,
	cinfo=c(1,1,3). }
  \item{tol}{A value indicating the tolerance for contrast estimability check }
}
\details{
}
\value{
}
\references{Bretz, F and Landgrebe J and Brunner E (2003):"Design and analysis of
        two colour factorial microarray experiments", submitted.
	        \url{http://www.microarrays.med.uni-goettingen.de/}
		}
\author{Jobst Landgrebe (jlandgr1@gwdg.de) and Frank Bretz (bretz@bioinf.uni-hannover.de) }
\note{}

\seealso{}
\examples{}
\keyword{}
