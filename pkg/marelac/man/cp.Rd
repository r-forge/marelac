\name{cp}
\alias{cp}
\title{Heat capacity of sea water}
\description{Estimates the heat capacity of seawater, using the UNESCO 1983 polynomial.

Valid for S=0 to 40, T = 0 to 35 dg C}
\usage{cp(T=25, S=35, P=1.013253, hydroP=max(0,P-1.013253))}
\arguments{
  \item{T }{Temperature, \eqn{^\circ}{degrees }C}
  \item{S }{Salinity, -}
  \item{P }{Pressure, bar}
  \item{hydroP }{hydrostatic pressure, pressure relative to the local atmospheric pressure, bar}
}
\value{
  Heat capacity, in \eqn{J kg^{-1} dgC^{-1}}
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{cp(S=40,T=40,hydroP=1000)}
\references{Fofonoff, P. and Millard, R.C. Jr. Unesco 1983. Algorithms for computation of fundamental properties of
 seawater. Unesco Tech. Pap. in Mar. Sci., No. 44, 53 pp.}
\note{hydroP is hydrostatic pressure, 0 bar at sea surface.

cp for S=40, T=40, P=1000 is 3849.5 J/(kg dgC)}
\keyword{utilities}

