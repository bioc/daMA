\name{designs.composite}
\alias{designs.composite}
\docType{data}
\title{Composite designs for two-colour factorial 3 x 2 microarray data}
\description{A list of matrices describing composite designs for two-colour factorial
microarray data of size 3 x 2 using 18 microarrays each. The design matrices are made up of
basic designs. Matrix rows represent microarrays, matrix columns represent parameters.
}
\usage{data(designs.composite)}
\format{
List of 10 matrices of size 18 x 9.
}
\details{
The matrix names reflect the basic designs they are made up from. The first two digits of
the names abbreviated the first basic design, the second two the second design etc.
The basic design abbreviations are: BS - swap over B, AL - A loop, XL - crossed loop, CL -
circle loop, RS - star swap, TL - triangular loop, CR - common reference.
BSBSBS is a tripled basic BS design, CLCLTL is a double circle loop design combined with a
triangular design and so on.
}
\source{
cDNA microarray lab of the University of Goettingen, Germany.
\url{http://www.microarrays.med.uni-goettingen.de}
}

\examples{
data(designs.composite)
}
\keyword{datasets}
