\name{vapor}
\alias{vapor}
\title{Saturation Water Vapor Pressure}
\description{
  The partial pressure of water in saturated air (pH20/P), as in Weiss
  and Price (1980), where P is the total atmospheric pressure, (1
  atmosphere), and pH2O is the partial pressure of the water vapor.
}
\usage{vapor(S = 35, T = 25)}
\arguments{
  \item{S }{Salinity (-),}
  \item{T }{Temperature, \eqn{^\circ}{degrees }C.}
}
\value{
  The saturation vapor pressure (-).
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
plot(0:30, vapor(T = 0:30), xlab = "Temperature, dgC", ylab = "pH2O/P")
}
\references{
  Sarmiento J.L. and N. Gruber, 2006. Ocean Biogeochemical
  Dynamics. Princeton University Press, Princeton. p 74

  Weiss, R. and B.A. Price, 1980. Nitrous oxide solubility in wate and
  sewater. Mar. Chem. 8, 347-359.
}
\keyword{utilities}

