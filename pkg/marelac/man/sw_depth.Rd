\name{sw_depth}
\alias{sw_depth}
\title{Water Depth}
\description{
  Computes the water depth for water of salinity 35, and temperature 0
  dg, based on latitude and hydrostatic pressure, using the UNESCO 1983
  polynomial.
}
\usage{
sw_depth(p = P-1.013253, P=1.013253, lat = 0)
}
\arguments{
  \item{p }{gauge or applied pressure, pressure relative to the local atmospheric pressure, bar}
  \item{P }{true pressure, bar}
  \item{lat }{latitude (-90 - +90), -,}
}
\value{
  Water depth in m.
}
\note{

  \code{sw_depth} for  p = 1000, lat = 30 is 9712.653 m.
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
sw_depth(p = 1000, lat = 30)
}
\references{
  Fofonoff, P. and Millard, R.C. Jr. Unesco 1983. Algorithms for
  computation of fundamental properties of seawater. Unesco
  Tech. Pap. in Mar. Sci., No. 44, 53 pp.
}
\keyword{utilities}

