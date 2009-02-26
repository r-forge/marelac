\name{convert_RtoS}
\alias{convert_RtoS}
\title{Conductivity-Salinity Conversion}
\description{
  Estimates salinity from conductivity, temperature and pressure.

  The equation is valid over ranges of temperature -2 to 35
  \eqn{^\circ}{degrees }C, pressure of 0 - 1000 bar and salinity 2-42 in the world ocean.
}
\usage{
convert_RtoS(R = 1, t = 25, p = max(0, P-1.013253), P = 1.013253)
}
\arguments{
  \item{R }{Conductivity ratio, -, the conductivity at (S, t, P) divided
    by the conductivity at S = 35, t = 25, p = 0}
  \item{t }{Temperature, \eqn{^\circ}{degrees }C}
  \item{p }{gauge (or applied) pressure, the pressure relative to the local
    atmospheric pressure, bar.}
  \item{P }{True pressure, bar}
}
\value{
  The salinity.
}
\note{
  The conductivity ratio for Salinity = 40.0000, t = 40, p = 1000 is
  1.888091.
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
convert_RtoS(R = 1.888091, t = 40, p = 1000)
## Salinity = 40.0000, t = 40, p = 1000, cond = 1.888091
}
\references{
  Fofonoff, P. and Millard, R.C. Jr. Unesco 1983. Algorithms for
  computation of fundamental properties of seawater. Unesco
  Tech. Pap. in Mar. Sci., No. 44, 53 pp.
}
\keyword{utilities}

