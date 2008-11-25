\name{rhoH2O}
\alias{rhoH2O}
\title{Density of Water}
\description{
  Density of water in \eqn{kg m^{-3}}{kg m^{-3}}
}
\usage{
rhoH2O(S=35, T=25, P=0, method=c("Millero", "Chen", "marine", "limnological"))
}
\arguments{
  \item{S}{salinity in PSU,}
  \item{T}{temperature in degree centigrade,}
  \item{P}{Hydrostatic pressure in bar (surface =0),}
  \item{method}{"Millero" for marine systems (default), "Chen" for the limnological
                range (i.e. fresh water systems). Abbreviation (e.g. "M" or "m")
                is possible.}
}

\value{
  Density of water in \eqn{kg m^{-3}}{kg m^{-3}}.
}
\details{
  To avoid confusion between the arguments (S, T, P) it is advisable to use
  named arguments in general (e.g. \code{rho(T=4, method="Chen")}.
  The marine formula is imported from package \pkg{seacarb}.)
}
\references{ 
Chen, Ch.-T. and F.J. Millero (1986). Precise thermodynamic properties of natural 
waters covering only the limnological range. Limnol. Oceanogr. 31 No. 3, 657 - 662.

Fofonoff N.P. and R.C. Millard (1983). Algorithms for computation of fundamental
properties of seawater. UNESCO technical papers in marine science, 44, 53 pp.

Millero, F. J. and Poisson, A. (1981). International one-atmosphere equation of state
of seawater. Deep-Sea Research 28A, 625-629.
}

\seealso{\code{\link[seacarb]{rho}} in package seacarb.}
\examples{
rhoH2O(T=4)
T <- 0:30
plot(T, rhoH2O(T=T))
S <- 0:40
plot(S, rhoH2O(S=S, T=4, method="Millero"))
lines(S, rhoH2O(S=S, T=4, method="Chen"))

}
\keyword{ misc }
