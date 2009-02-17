\name{gravity}
\alias{gravity}
\title{Gravity on Earth}
\description{
  Computes the gravity, based on latitude.
}
\usage{
gravity(lat = 0, UNESCO =FALSE)
}
\arguments{
  \item{lat }{latitude (-90 - +90).}
  \item{UNESCO }{When TRUE, uses the UNESCO (1983) polynomial, else according to Moritz, 2000 }
}
\value{
  Gravity, in \eqn{\rm m\,sec^{-2}}{m/s2}.
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
gravity(lat = 30)}
\references{
The UNESCO polynomial:

  Fofonoff, P. and Millard, R.C. Jr. Unesco 1983. Algorithms for
  computation of fundamental properties of seawater. Unesco
  Tech. Pap. in Mar. Sci., No. 44, 53 pp.

Moritz, H., 2000: Geodetic reference system 1980. Journal of Geodesy, 74, 128-133.
}
\keyword{utilities}

