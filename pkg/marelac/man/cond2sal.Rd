\name{cond2sal}
\alias{cond2sal}
\title{Conductivity-Salinity Conversion}
\description{
  Estimates salinity from conductivity, temperature and pressure.

  The equation is valid over ranges of temperature -2 to 35
  \eqn{^\circ}{degrees }C, and salinity 2-42 in the world ocean.
}
\usage{
cond2sal(R = 1, T = 25, P = 1.013253, hydroP = max(0, P-1.013253))
}
\arguments{
  \item{R }{Conductivity ratio, -, the conductivity at T, S, P divided
    by the conductivity at S = 35, T = 25, P = 1,}
  \item{T }{Temperature, \eqn{^\circ}{degrees }C,}
  \item{P }{Pressure, bar,}
  \item{hydroP }{hydrostatic pressure, pressure relative to the local
    atmospheric pressure, bar.}
}
\value{
  The salinity.
}
\note{
  Pressure here is true pressure, 1 bar (at sea surface), in contrast to
  hydrostatic pressure in dbar of original formula.

  The conductivity ratio for Salinity = 40.0000, T = 40, P = 1001 is
  1.888091.
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
cond2sal(R = 1.888091, T = 40, hydroP = 1000)
## Salinity = 40.0000, T = 40, P = 1000, cond = 1.888091
}
\references{
  Fofonoff, P. and Millard, R.C. Jr. Unesco 1983. Algorithms for
  computation of fundamental properties of seawater. Unesco
  Tech. Pap. in Mar. Sci., No. 44, 53 pp.
}
\keyword{utilities}

