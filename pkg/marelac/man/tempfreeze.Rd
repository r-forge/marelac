\name{tempfreeze}
\alias{tempfreeze}
\title{Freezing Temperature of Sea Water}
\description{
  Estimates the freezing temperature of seawater, using the UNESCO 1983
  polynomial.
  
  Valid for salinity 4-40
}
\usage{
tempfreeze(S=35, P = 1.013253, hydroP = max(0, P-1.013253))
}
\arguments{
  \item{S }{salinity, -,}
  \item{P }{pressure, bar,}
  \item{hydroP }{hydrostatic pressure, pressure relative to the local
    atmospheric pressure, bar.
  }
}
\value{
  Temperature, \eqn{^\circ}{degrees }C
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
tempfreeze(S = 40,hydroP = 50)
}
\note{
  \code{hydroP} is hydrostatic pressure, 0 bar at sea surface.

  Tempfreeze for S = 40, \code{hydroP} = 50 is -2.588567 dgC.
}
\references{
  Fofonoff, P. and Millard, R.C. Jr. Unesco 1983. Algorithms for
  computation of fundamental properties of seawater. Unesco
  Tech. Pap. in Mar. Sci., No. 44, 53 pp.
}
\keyword{utilities}

