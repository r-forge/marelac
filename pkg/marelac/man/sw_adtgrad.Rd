\name{sw_adtgrad}
\alias{sw_adtgrad}
\title{Adiabatic Temperature Gradient in Seawater}

\description{
  Computes the adiabatic temperature gradient in seawater,
  using the UNESCO 1983 polynomial.\cr

  Also known as the adiabatic lapse rate, the change of temperature per
  unit pressure for an adiabatic change of pressure of an element of
  seawater. It is assumed that no heat or salt is exchanged with the
  surroundings.
}
\usage{
sw_adtgrad(S = 35, t = 25, p = P-1.013253, P=1.013253 )
}
\arguments{
  \item{S }{Salinity, -}
  \item{t }{Temperature, \eqn{^\circ}{degrees }C}
  \item{p }{gauge or applied pressure, pressure relative to the local atmospheric pressure, bar}
  \item{P }{true pressure, bar}
}
\value{
  adiabatic temperature gradient, in dg K / bar
}
\note{

  Note: in the original formula, the units of \code{sw_adtgrad} are dg
  K/dbar (here: dg K/bar).

  \code{sw_adtgrad} for S = 40, t = 40, p = 1000 is 3.255976E-3

}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl> }
\examples{
sw_adtgrad(t = 40, S = 40, p = 1000)
}
\references{
  Fofonoff, P. and Millard, R.C. Jr. Unesco 1983. Algorithms for computation of fundamental properties of
  seawater. Unesco Tech. Pap. in Mar. Sci., No. 44, 53 pp.
}
\keyword{utilities}

