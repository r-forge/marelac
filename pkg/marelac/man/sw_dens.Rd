\name{sw_dens}
\alias{sw_dens}
\title{Density of Sea Water}
\description{
  Density of sea water in \eqn{kg m^{-3}}{kg m^{-3}}
}
\usage{
sw_dens(S = 35, t = 25, p = max(0,P-1.013253), P=1.013253, UNESCO=FALSE)
}
\arguments{
  \item{S }{Salinity (dimensionless),}
  \item{t }{Temperature, \eqn{^\circ}{degrees }C,}
  \item{p }{gauge or applied pressure, pressure relative to the
            local atmospheric pressure, bar}
  \item{P }{true pressure, bar}
  \item{UNESCO }{When TRUE, uses the UNESCO (1983) polynomial,
                else using the gibbs function according to Feistel, 2008 }
}

\value{
  Density of water in \eqn{kg m^{-3}}{kg m^-3}.
}
\details{
  To avoid confusion between the arguments (S, t, p) it is advisable to
  use named arguments in general (e.g. \code{rho(t = 4)}.
  The UNESCO formula is imported from package \pkg{seacarb}.)
}
\note{
  Pressure used here is 1 bar (true pressure), in contrast to
  hydrostatic pressure (0 bar at surface) in original formula.

  The coefficients from McDougall et al., 2009 were used.
  For temperature, they differ slightly from Feistel 2003 and Feistel 2008,
  which is why, for temparatures different from 0, there is a slight offset
  from the estimates as from table 22 or 21 from Feistel (2008).

}

\references{

  Fofonoff, P. and Millard, R.C. Jr. Unesco 1983. Algorithms for
  computation of fundamental properties of seawater. Unesco
  Tech. Pap. in Mar. Sci., No. 44, 53 pp.


  Feistel, R., 2008: A Gibbs function for seawater thermodynamics for -6 to 80
  dgC and salinity up to 120 g/kg. Deep-Sea Res. I, 55, 1639-1671.

  McDougall, T.J., R. Feistel, F. J. Millero, D. R. Jackett, D. G. Wright,
  B. A. King, G. M. Marion, C.-T. A. Chen and P. Spitzer, 2009: Calculation of
  the Thermophysical Properties of Seawater, Global Ship-based Repeat
  Hydrography Manual, IOCCP Report No. 14, ICPO Publication Series no. 134.

}

\seealso{\code{\link[seacarb]{rho}} in package seacarb.}
\examples{
# table 22 Feistel 2008
sw_dens(0,0,0)                #0.999843086e3
sw_dens(0,79.85,0)            #0.97188383e3   - deviates
sw_dens(0,0,998.98675)        #0.104527796e4

# table 21 Feistel 2008
sw_dens(35.16504,0,0)         #0.10281072e4
sw_dens(100,79.85,0)          #0.102985888e4
sw_dens(35.16504,0,998.98675) #0.10709264e4

sw_dens(35.7,25.5,102.3)      #1027.95249315662

S <- 0:40
plot(S, sw_dens(S = S, t = 4, UNESCO=TRUE))

lines(S, sw_dens(S = S, t = 4, UNESCO=FALSE),col="red")
}
\keyword{ misc }
