\name{rhoH2O}
\alias{rhoH2O}
\title{Density of Water}
\description{
  Density of water in \eqn{kg m^{-3}}{kg m^{-3}}
}
\usage{
rhoH2O(S = 35, t = 25, p = max(0, P-1.013253), P = 1.013253,
       method = c("Millero", "Chen", "marine", "limnological"))
}
\arguments{
  \item{S }{salinity in PSU,}
  \item{t }{temperature in \eqn{^\circ}{degrees }C,}
  \item{P }{true pressure in bar,}
  \item{p }{gauge or applied pressure, pressure relative to the
            local atmospheric pressure, bar.
  }
  \item{method}{"Millero" for marine systems (default), "Chen" for the
    limnological range (i.e. fresh water systems). Abbreviation
    (e.g. "M" or "m") is possible.
  }
}

\value{
  Density of water in \eqn{kg m^{-3}}{kg m^-3}.
}
\details{
  To avoid confusion between the arguments (t, S, P) it is advisable to
  use named arguments in general (e.g. \code{rho(t = 4, method =
  "Chen")}.  The marine formula is imported from package \pkg{seacarb}.)
}
\note{
  Pressure used here is 1 bar (true pressure), in contrast to
  hydrostatic pressure (0 bar at surface) in original formula.
}

\references{
  
  Chen Ch.-T. and Millero FJ, 1986. Precise thermodynamic properties
  of natural waters covering only the limnological
  range. Limnol. Oceanogr. 31 No. 3, 657 - 662.

  Fofonoff NP and Millard RC Jr, 1983. Algorithms for computation of
  fundamental properties of seawater. UNESCO technical papers in marine
  science, 44, 53 pp.

  \url{http://unesdoc.unesco.org/images/0005/000598/059832EB.pdf}

  Millero FJ and Poisson A, 1981. International one-atmosphere
  equation of state of seawater. Deep-Sea Research 28A, 625-629.
}

\seealso{\code{\link[seacarb]{rho}} in package seacarb.}
\examples{
rhoH2O(t = 4)
T <- 0:30
plot(T, rhoH2O(t = T))
S <- 0:40
plot(S, rhoH2O(S = S, t = 4, method = "Millero"))
lines(S, rhoH2O(S = S, t = 4, method = "Chen"))
}
\keyword{ misc }
