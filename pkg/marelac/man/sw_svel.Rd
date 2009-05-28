\name{sw_svel}
\alias{sw_svel}
\title{Velocity of the Sound in Seawater}
\description{
  Computes the velocity of the sound in seawater,
  using the UNESCO 1983 polynomial or based on the Gibbs function.

  Valid for salinity from 0 to 40, temperature from 0 to 40 dgC, pressure
  from 1 to 1000 bars.
}
\usage{
sw_svel(S = 35, t = 25, p = P-1.013253, P = 1.013253, UNESCO = FALSE)
}
\arguments{
  \item{S }{Salinity (dimensionless),}
  \item{t }{Temperature, \eqn{^\circ}{degrees }C,}
  \item{p }{gauge or applied pressure, pressure relative to the local
    atmospheric pressure, bar
  }
  \item{P }{true pressure, bar}
  \item{UNESCO }{When TRUE, uses the UNESCO (1983) polynomial,
    else using the gibbs function according to Feistel, 2008
  }
}
\value{
  Sound velocity, in m / sec.
}
\note{
  Sound velocity for S = 40, t = 40, p = 1000 is 1731.995 using UNESCO
  polynomial.
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{

sw_svel(t = 40, S = 40, p = 1000, UNESCO = TRUE)  # 1731.995
sw_svel(S = 35.7, t = 25.5, p = 102.3)            # 1552.93372863425

}
\references{
  Fofonoff, P. and Millard, R.C. Jr. Unesco 1983. Algorithms for
  computation of fundamental properties of seawater. Unesco
  Tech. Pap. in Mar. Sci., No. 44, 53 pp.


  Feistel, R., 2008: A Gibbs function for seawater thermodynamics for -6
  to 80 dgC and salinity up to 120 g/kg. Deep-Sea Res. I, 55, 1639-1671.

  McDougall, T.J., R. Feistel, F. J. Millero, D. R. Jackett, D. G. Wright,
  B. A. King, G. M. Marion, C.-T. A. Chen and P. Spitzer, 2009: Calculation
  of the Thermophysical Properties of Seawater, Global Ship-based Repeat
  Hydrography Manual, IOCCP Report No. 14, ICPO Publication Series no. 134.


}
\keyword{utilities}

