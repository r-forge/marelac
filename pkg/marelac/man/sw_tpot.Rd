\name{sw_tpot}
\alias{sw_tpot}
\title{Potential Temperature of Sea Water}
\description{
  Estimates the potential temperature of seawater, using the UNESCO 1983
  polynomial.

  It is the temperature an element of seawater would have if raised
  adiabatically with no change of salinity, to atmospheric pressure.
}
\usage{sw_tpot(S = 35, t = 25, p, pref=0)}
\arguments{
  \item{t }{temperature, \eqn{^\circ}{degrees }C,}
  \item{S }{salinity, -,}
  \item{p }{gauge or applied pressure, pressure relative to the local atmospheric pressure, bar}
  \item{pref }{reference hydrostatic pressure, bar.}
}
\value{
  Temperature, \eqn{^\circ}{degrees }C.
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
sw_tpot(S = 40, t = 40, p = 1000)
}
\references{
  Fofonoff, P. and Millard, R.C. Jr. Unesco 1983. Algorithms for
  computation of fundamental properties of seawater. Unesco
  Tech. Pap. in Mar. Sci., No. 44, 53 pp.
}
\note{
\code{sw_tpot} for S = 40, t = 40, p = 1000 is 36.89073 dgC
}
\keyword{utilities}

