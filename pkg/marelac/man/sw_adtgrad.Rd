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
  \item{S }{Practical salinity (-),}
  \item{t }{Temperature, \eqn{^\circ}{degrees }C}
  \item{p }{gauge or applied pressure, pressure relative to the local
    atmospheric pressure, bar}
  \item{P }{true pressure, bar}
}
\value{
  adiabatic temperature gradient, in dg K / bar
}
\note{

  Note: in the original formula, the units of \code{sw\_adtgrad} are dg
  K/dbar (here: dg K/bar).

  \code{sw\_adtgrad} for S = 40, t = 40, p = 1000 is 3.255976e-3

}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl> }
\examples{
sw_adtgrad(t = 40, S = 40, p = 1000)
}
\references{
  Fofonoff NP and Millard RC Jr, 1983. Algorithms for computation of
  fundamental properties of seawater. UNESCO technical papers in marine
  science, 44, 53 pp.

  \url{http://unesdoc.unesco.org/images/0005/000598/059832EB.pdf}
}
\seealso{
\code{\link{sw_alpha}}, \code{\link{sw_beta}},
\code{\link{sw_comp}}, \code{\link{sw_conserv}}, \code{\link{sw_cp}},
\code{\link{sw_dens}},
\code{\link{sw_depth}}, \code{\link{sw_enthalpy}}, \code{\link{sw_entropy}},
\code{\link{sw_gibbs}}, \code{\link{sw_kappa}},
\code{\link{sw_kappa_t}}, \code{\link{sw_sfac}},
\code{\link{sw_svel}}, \code{\link{sw_tfreeze}}, \code{\link{sw_tpot}}

}
\keyword{utilities}

