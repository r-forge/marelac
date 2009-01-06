\name{temppot}
\alias{temppot}
\title{Potential temperature of sea water}
\description{Estimates the potential temperature of seawater, using the UNESCO 1983 polynomial.

It is the temperature an element of seawater would have if raised adiabatically with no change
of salinity, to atmospheric pressure}
\usage{temppot(T=25, S=35, hydroP, hydroPref = 0)}
\arguments{
  \item{T }{Temperature, \eqn{^\circ}{degrees }C}
  \item{S }{Salinity, -}
  \item{hydroP }{(hydrostatic) Pressure, bar}
  \item{hydroPref }{Reference hydrostatic pressure, bar}
}
\value{
  Temperature, \eqn{^\circ}{degrees }C
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{temppot(S=40,T=40,hydroP=1000)}
\references{Fofonoff, P. and Millard, R.C. Jr. Unesco 1983. Algorithms for computation of fundamental properties of
 seawater. Unesco Tech. Pap. in Mar. Sci., No. 44, 53 pp.}
\note{

temppot for S=40, T=40, hydroP=1000 is 36.89073 dgC}
\keyword{utilities}

