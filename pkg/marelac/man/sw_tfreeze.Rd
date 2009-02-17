\name{sw_tfreeze}
\alias{sw_tfreeze}
\title{Freezing Temperature of Sea Water}
\description{
  Estimates the freezing temperature of seawater, using the UNESCO 1983
  polynomial.
  
  Valid for salinity 4-40
}
\usage{
sw_tfreeze(S=35, p = P-1.013253, P=1.013253 )
}
\arguments{
  \item{S }{salinity, -,}
  \item{p }{gauge or applied pressure, pressure relative to the local atmospheric pressure, bar}
  \item{P }{true pressure, bar}
}
\value{
  Temperature, \eqn{^\circ}{degrees }C
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
sw_tfreeze(S = 40,p = 50)
}
\note{

  freezing temperature for S = 40, \code{p} = 50 is -2.588567 dgC.
}
\references{
  Fofonoff, P. and Millard, R.C. Jr. Unesco 1983. Algorithms for
  computation of fundamental properties of seawater. Unesco
  Tech. Pap. in Mar. Sci., No. 44, 53 pp.
}
\keyword{utilities}

