\name{gravity}
\alias{gravity}
\title{Gravity on Earth}
\description{
  Computes the gravity, based on latitude, using the UNESCO 1983
  polynomial.
}
\usage{
gravity(lat = 0)
}
\arguments{
  \item{lat }{latitude (-90 - +90).}
}
\value{
  Gravity, in \eqn{\rm m\,sec^{-2}}{m/s2}.
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
gravity(lat = 30)}
\references{
  Fofonoff, P. and Millard, R.C. Jr. Unesco 1983. Algorithms for
  computation of fundamental properties of seawater. Unesco
  Tech. Pap. in Mar. Sci., No. 44, 53 pp.
}
\keyword{utilities}

