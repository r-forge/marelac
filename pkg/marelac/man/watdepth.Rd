\name{watdepth}
\alias{watdepth}
\title{water depth}
\description{Computes the water depth for water of salinity 35, and temperature 0 dg,
based on latitude and hydrostatic pressure, using the UNESCO 1983 polynomial}
\usage{watdepth(P=1.013253,hydroP=max(0,P-1.013253), lat = 0)}
\arguments{
  \item{P }{pressure, bar}
  \item{hydroP }{hydrostatic pressure, pressure relative to the local atmospheric pressure, bar}
  \item{lat }{latitude (-90 - +90), -}
}
\value{
  water depth, in m
}
\note{hydroP is hydrostatic pressure, 0 bar at sea surface.

depth for  hydroP=1000, lat=30 is 9712.653 m}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{watdepth(hydroP=1000,lat=30)}
\references{Fofonoff, P. and Millard, R.C. Jr. Unesco 1983. Algorithms for computation of fundamental properties of
 seawater. Unesco Tech. Pap. in Mar. Sci., No. 44, 53 pp.}
\keyword{utilities}

