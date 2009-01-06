\name{adtgrad}
\alias{adtgrad}
\title{adiabatic temperature gradient in seawater}
\description{Computes the adiabatic temperature gradient in seawater, using the UNESCO 1983 polynomial.

Also known as the adiabatic lapse rate, the change of temperature per unit pressure for an adiabatic change
of pressure of an element of seawater.
It is assumed that no heat or salt is exchanged with the surroundings.}
\usage{adtgrad(T=25, S=35, P=1.013253, hydroP=max(0,P-1.013253))}
\arguments{
  \item{S }{Salinity, -}
  \item{T }{Temperature, \eqn{^\circ}{degrees }C}
  \item{P }{Pressure, bar}
  \item{hydroP }{hydrostatic pressure, pressure relative to the local atmospheric pressure, bar}
}
\value{
  adiabatic temperature gradient, in dg K / bar
}
\note{hydroP is hydrostatic pressure, 0 bar at sea surface.

Note: in the original formula, the units of adtgrad are dg K/dbar (here: dg K/bar)

adtgrad for S=40, T=40, hydroP=1000 is 3.255976E-3}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{adtgrad(T=40,S=40,hydroP=1000)}
\references{Fofonoff, P. and Millard, R.C. Jr. Unesco 1983. Algorithms for computation of fundamental properties of
 seawater. Unesco Tech. Pap. in Mar. Sci., No. 44, 53 pp.}
\keyword{utilities}

