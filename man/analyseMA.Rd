\name{analyseMA}
\alias{analyseMA}

\title{ANALYSIS OF FACTORIAL MICROARRAY EXPERIMENTS}
\description{
analyseMA is used for the ananlysis of factorial two-colour microarray experiments based
on the experimental design, a user-defined matrix containing the experimental question in
contrast form and a vector to discern vectorial contrasts from contrasts given in matrix
form. 
}
\usage{
analyseMA( data, design, id, cmat, cinfo, padj=c("none","bonferroni","fdr"), tol=1e-06 ) 
}

\arguments{
  \item{data}{ a matrix of size \eqn{G \times N} containing the normalized and/or standardized
data to be analyzed, where G is the number of spots under investigation and N is the
number of arrays used in the experiment. The matrix should contain one row for each spot.
The matrix should contain as many columns as arrays involved in the experiment, such that
each column contains the data for one single array. The matrix should not contain any ID variables, 
which are entered separately. Missing values should be entered as NA.}
  \item{design}{ the design matrix of size \eqn{N \times (K+2)}, where K is the number of
experimental conditions. This is the design matrix X known from linear model theory and
its elements are typically 0, 1, or -1. A 0 means that the associated parameter does not
apply for the corresponding observation (i.e., row). The first two columns are reserved
for the two dyes and are usually filled up with 1 and -1, respectively. }
  \item{id}{an ID vector of length G for the identification of the spots.}
  \item{cmat}{a matrix describing the p experimental questions (contrasts) to be analysed in
	the experiment. The matrix can be composed of vectorial contrasts (a single row of the
	matrix) and of contrasts in matrix form (several rows of the matrix),
	e.g. an \eqn{A \times B} interaction effect in a \eqn{3 \times 2} design. All
	contrasts have to be combined into one matrix (using rbind for instance). }
  \item{cinfo}{ a vector of length p describing the grouping of the contrast matrix rows in vector or matrix 
	form. E.g. if the design matrix contains three contrasts in vector form, cinfo =
	rep(1,3), if it contains two vectorial contratst and one as matrix with three rows,
	cinfo=c(1,1,3). }
  \item{padj}{a quoted string indicating the multiplicity adjustment that should be used.
"none" - no multiplicity adjustment, "bonferroni" - Bonferroni single step adjustment,
"fdr" - linear step-up procedure of Benjamini and Hochberg.}
  \item{tol}{A value indicating the tolerance for contrast estimability check }
}
\details{
  The analysis is perfomed separately for each spot. For each spot, arrays with NA values
are dropped. Then, for each experimental question (either contrast vector or contrast
matrix) a check on the estimabilty of the resulting linear function is done. If the linear
function of interest is estimable, t- or F-tests (whichever is appropriate) are computed
and the associated unadjusted $P-$values are computed. Multiplicity adjustment is done
over the number of spots only.
}
\value{
  a \eqn{G \times (4p+3)} matrix with the following row-wise components. 
  \item{(i)}{the first column contains the ID}
  \item{(ii)}{columns 2 though p+1 contain the estimates of the linear function (in case
of vectorial contrasts) or the dregrees of freedom for the quadratic form in the numerator
(in case of contrasts given in matrix form and that F-tests are used), depending on cinfo.}
  \item{(iii)}{columns p+2 through 2p+1 contain the test statistics (either t- or F-tests,
depending on cinfo)}
  \item{(iv)}{columns 2p+2 through 3p+1 contain the raw P-values, associated to the t- and
F-tests}
  \item{(v)}{column 3p+2 contains the mean square error}
  \item{(vi)}{column 3p+3 contains the residual degrees of freedom}
  \item{(vii)}{columns 3p+4 through 4p+3 contain the multiplicity adjusted P-values,
associated to the raw P-values, as long as a multiplicty adjustment method has been selected}
}
\references{Bretz, F and Landgrebe J and Brunner E (2003):"Design and analysis of
	two colour factorial microarray experiments", submitted.
	\url{http://www.microarrays.med.uni-goettingen.de/} 
}
\author{Jobst Landgrebe (jlandgr1@gwdg.de) and Frank Bretz (bretz@bioinf.uni-hannover.de) }
\note{}

\seealso{}
\examples{
	\dontrun{      result <-	analyseMA( data=data.3x2, design=designs.composite$BSBSBS, id=id.3x2,
		cmat=cmatB.AB, cinfo=c(1,3), padj=c("fdr"), tol=1e-06 ) # analyse a dataset with
									# 30012 spots and 18 arrays. The design
									# is 3x2 with 3 replicates, the
									# contrasts of interest are the main effect
									# B and the interaction effect AxB.
	}
}
\keyword{design}
\keyword{models}
