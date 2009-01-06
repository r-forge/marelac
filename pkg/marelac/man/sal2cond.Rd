\name{sal2cond}
\alias{sal2cond}
\title{Salinity-Conductivity conversion}
\description{Estimates the conductivity ratio from salinity, temperature and pressure

The equation is valid over ranges of temperature -2 to 35 \eqn{^\circ}{degrees }C, and
salinity 2-42 in the world ocean.
}
\usage{sal2cond(T=25, S=35, P=1.013253, hydroP=max(0,P-1.013253))}
\arguments{
  \item{S }{Salinity, -}
  \item{T }{Temperature, \eqn{^\circ}{degrees }C}
  \item{P }{Pressure, bar}
  \item{hydroP }{hydrostatic pressure, pressure relative to the local atmospheric pressure, bar}
}
\value{
  The conductivity ratio, -
  
this is the conductivity at T, S, P divided by the conductivity at S=35, T=25, hydroP=0.

}
\note{pressure here is true pressure, 1 bar (at sea surface), in contrast to hydrostatic pressure in dbar of original formula

The conductivity ratio for Salinity = 40.0000, T=40, hydroP= 1000 is 1.888091}

\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{sal2cond(S=40, T=40, hydroP=1000)
}
\references{Fofonoff, P. and Millard, R.C. Jr. Unesco 1983. Algorithms for computation of fundamental properties of
 seawater. Unesco Tech. Pap. in Mar. Sci., No. 44, 53 pp.}
\keyword{utilities}

