\name{sw_cp}
\alias{sw_cp}
\title{Heat Capacity of Sea Water}
\description{
  Estimates the heat capacity of seawater.

  Valid for S = 0 to 40, T = 0 to 35 dg C
}
\usage{
sw_cp(S = 35, t = 25, p = P-1.013253, P=1.013253, UNESCO=FALSE)
}
\arguments{
  \item{S }{Salinity (dimensionless),}
  \item{t }{Temperature, \eqn{^\circ}{degrees }C,}
  \item{p }{gauge or applied pressure, pressure relative to the local atmospheric pressure, bar}
  \item{P }{true pressure, bar}
  \item{UNESCO }{When TRUE, uses the UNESCO (1983) polynomial,
                else using the gibbs function according to Feistel, 2008 }
}
\value{
  Heat capacity, in \eqn{J kg^{-1} dgC^{-1}}
}
\author{Karline Soetaert < k.soetaert@nioo.knaw.nl > }
\examples{
sw_cp(S = 40, t = 40, p = 1000, UNESCO=TRUE)  # 3849.5

sw_cp(35.7,25.5,102.3)#3974.42541259729

}
\references{
  Fofonoff, P. and Millard, R.C. Jr. Unesco 1983. Algorithms for
  computation of fundamental properties of seawater. Unesco
  Tech. Pap. in Mar. Sci., No. 44, 53 pp.


  Feistel, R., 2008: A Gibbs function for seawater thermodynamics for -6 to 80 dgC and salinity up
to 120 g /kg. Deep-Sea Res. I, 55, 1639--1671.

  McDougall, T.J., R. Feistel, F. J. Millero, D. R. Jackett, D. G. Wright, B. A. King,
  G. M. Marion, C.--T. A. Chen
and P. Spitzer, 2009: Calculation of the Thermophysical Properties of Seawater,
Global Ship-based Repeat
Hydrography Manual, IOCCP Report No. 14, ICPO Publication Series no. 134.

}
\note{
  \code{p} is applied pressure, 0 bar at sea surface.

  when using UNESCO polynomial, \code{cp} for S = 40, T = 40, P = 1000 is 3849.5 J/(kg dg C).
}
\keyword{utilities}

