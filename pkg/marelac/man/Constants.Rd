\name{Constants}
\docType{data}
\alias{Constants}
\title{Useful Physical and Chemical Constants}
\description{
  Physical and chemical constants useful for marine and lacustrine
  sciences.
}
\usage{Constants}
\format{
  A list with, for each constant the value, the units, and a
  description.
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{

data.frame(cbind(acronym = names(Constants),
           matrix(ncol = 3, byrow = TRUE, data = unlist(Constants),
           dimnames = list(NULL, c("value", "units", "description")))))
}
\references{Mohr, P.J. and B.N. Taylor, 2005.
CODATA recommended values of the fundamental physical constants: 2002,
Review of Modern Physics 77, pp. 1 - 107}
\keyword{datasets}
