\name{Oceans}
\docType{data}
\alias{Oceans}
\title{Useful Characteristics of the Oceans}
\description{
  Surfaces and volumes of the world's oceans
}
\usage{Oceans}
\format{
  A list with, for each characteristic the value, the units, and a
  description.
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{

data.frame(cbind(acronym = names(Oceans),
           matrix(ncol = 3, byrow = TRUE, data = unlist(Oceans),
           dimnames = list(NULL, c("value", "units", "description")))))
}
\references{
  Sarmiento JL and Gruber N, 2006. Ocean Biogeochemical
  Dynamics. Princeton University Press, Princeton. p 85.
}
\seealso{
\code{\link{AtomicWeight}}, \code{\link{Bathymetry}},
\code{\link{Constants}}, \code{\link{earth_surf}}
}
\keyword{datasets}
