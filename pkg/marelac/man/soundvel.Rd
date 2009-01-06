\name{soundvel}
\alias{soundvel}
\title{velocity of the sound in seawater}
\description{Computes the velocity of the sound in seawater, using the UNESCO 1983 polynomial.

Valid for salinity from 0 to 40, temperature from 0 to 40 dgC, pressure from 1 to 1000 bars.}
\usage{soundvel(T=25, S=35, P=1.013253, hydroP=max(0,P-1.013253))}
\arguments{
  \item{S }{Salinity, -}
  \item{T }{Temperature, \eqn{^\circ}{degrees }C}
  \item{P }{Pressure, bar}
  \item{hydroP }{hydrostatic pressure, pressure relative to the local atmospheric pressure, bar}
}
\value{
  sound velocity, in m / sec
}
\note{hydro P is hydrostatic pressure, 0 bar at sea surface.

sound velocity for S=40, T=40, hydroP=1000 is 1731.995}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{soundvel(T=40,S=40,hydroP=1000)}
\references{Fofonoff, P. and Millard, R.C. Jr. Unesco 1983. Algorithms for computation of fundamental properties of
 seawater. Unesco Tech. Pap. in Mar. Sci., No. 44, 53 pp.}
\keyword{utilities}

