\name{convert_StoR}
\alias{convert_StoR}
\title{Salinity-Conductivity Conversion}
\description{
  Estimates the conductivity ratio from salinity, temperature and
  pressure.

  The equation is valid over ranges of temperature -2 to 35
  \eqn{^\circ}{degrees }C, pressure of 0 - 1000 bar and salinity 2-42 in the world ocean.
}
\usage{
convert_StoR(S = 35, t = 25, P = 1.013253, p = max(0, P-1.013253))
}
\arguments{
  \item{S }{(practical) Salinity, -,}
  \item{t }{Temperature, \eqn{^\circ}{degrees }C,}
  \item{P }{True Pressure, bar,}
  \item{p }{gauge (or applied) pressure, the pressure relative to the local
    atmospheric pressure, bar.}
}
\value{
  The conductivity ratio (-), this is the conductivity at T, S, P
  divided by the conductivity at S = 35, T = 25, \code{p} = 0.
}
\note{
  Pressure here is true pressure, 1 bar (at sea surface), in contrast to
  hydrostatic pressure in dbar of original formula.

  The conductivity ratio for Salinity = 40.0000, T = 40, \code{p} =
  1000 is 1.888091.
}

\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
convert_StoR(S = 40, t = 40, p = 1000)
}
\references{
  Fofonoff, P. and Millard, R.C. Jr. Unesco 1983. Algorithms for
  computation of fundamental properties of seawater. Unesco
  Tech. Pap. in Mar. Sci., No. 44, 53 pp.
}
\keyword{utilities}

