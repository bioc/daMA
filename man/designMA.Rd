\name{designMA}
\alias{designMA}

\title{DESIGN OF FACTORIAL MICROARRAY EXPERIMENTS}

\description{designMA computes efficient factorial microarray experimental designs
	for two-colour microarrays based on a list of user-defined design matrices, a
	matrix describing the experimental questions (contrasts), a vector to discern
	vectorial contrasts from contrasts given in matrix form
	and a design optimality criterion.}

\usage{
designMA(design.list, cmat, cinfo, type = c("d", "e", "t"), tol = 1e-06)
}

\arguments{
  \item{design.list}{a named list of design matrices. Each design matrix should have nrow = number of
	arrays and ncol= number of experimental conditions. With p columns, the first two columns describe
	the dye labeling (green and red), the remaining columns describe the experimental 
	conditions.
	} 
  \item{cmat}{a matrix describing the experimental questions (contrasts) to be analysed in
	the experiment. The matrix can be composed of vectorial contrasts (a single row of the
	matrix) and of contrasts in matrix form (several rows of the matrix),
	e.g. an \eqn{A \times B} interaction effect in a \eqn{3 \times 2} design. All
	contrasts have to be combined into one matrix (using rbind for instance).
	}

  \item{cinfo}{ a vector describing the grouping of the contrast matrix rows in vector or matrix 
	form. E.g. if the design matrix contains three contrasts in vector form, cinfo =
	rep(1,3), if it contains two vectorial contratst and one as matrix with three rows,
	cinfo=c(1,1,3). }	

  \item{type}{a quoted letter indicating the optimality criterion that shoul be used. "d"
	- determinant, "e" - eigenvalue, "t" - trace.}

  \item{tol}{A value indicating the tolerance for contrast estimability check.}
}

\details{The choice of the optimality criterion influences the design defined as best. We
	propose the trace criterion because of its straightforward interpretability. For a
	detailed description of optimality criteria cf. Pukelsheim, F. "Optimal Design of
	Experiments", New York 1993.
}
\value{ a list with the following components
  \item{alleff}{a matrix giving the absolute efficiency values (cols) for each contrast
	(rows). NA if contrast is not estimatable.
       }
  \item{alleffrel}{a matrix giving the relative efficiency values (cols) for each
	contrast (rows). The values are obtained by dividing the absolute values by the
	by the maximal efficiency value for a given contrast. NA if contrast is not estimatable.
       }
  \item{alleffave}{a vector giving the average efficiency for each design  over all
	contrasts.}
  \item{effdesign}{the name of the design with the highest alleffave value.}
  \item{df}{a vector with the degrees of freedom of the F-statistics obtained by the designs.}
}

\references{Bretz, F and Landgrebe J and Brunner E (2003):"Design and analysis of
	two colour factorial microarray experiments", submitted.
	\url{http://www.microarrays.med.uni-goettingen.de/} 
}

\author{Jobst Landgrebe (jlandgr1@gwdg.de) and Frank Bretz (bretz@bioinf.uni-hannover.de) }

\note{}

\seealso{}

\examples{ 	\dontrun{designs.basic}	# look at typical basic designs
		\dontrun{designs.composite } #look at comlpex composite designs
		\dontrun{t.eff.3x2.B.AB <-  designMA(designs.composite,
				   	    cmatB.AB,cinfoB.AB,type="t")# compute design efficiencies for
									# a \eqn{3 \times 2} factorial experiment
									# using 18 microarrays and asking for 
									# the main effect B and the interaction effect \eqn{A \times B}
		}
		\dontrun{t.eff.3x2.all <-  designMA(designs.composite,
					 	    cmat,cinfo,type="t")
		 }							#compute design efficiencies design for
									# a \eqn{3 \times 2} factorial
									# experiment using 18
                                                                        # microarrays and asking for 
                                                                        # the the simple B
                                                                        # effects, the main effects
                                                                        # A, B and the interaction
                                                                        # effect \eqn{A  \times B}
}

\keyword{design}
\keyword{models}
