\name{sw_entropy}
\alias{sw_entropy}
\title{Specific Entropy of Seawater}
\description{
  Computes the seawater specific entropy, J/(kg*K)
}
\usage{
sw_entropy(S = 35, t = 25, p = P-1.013253, P=1.013253)
}
\arguments{
  \item{S }{Salinity (dimensionless),}
  \item{t }{Temperature, \eqn{^\circ}{degrees }C,}
  \item{p }{gauge or applied pressure, pressure relative to the
            local atmospheric pressure, bar}
  \item{P }{true pressure, bar}
}
\value{
  Specific entropy, J/(kg*K).
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
sw_entropy(35.7,25.5,102.3) #352.81879771528
}
\references{
  Feistel, R., 2008: A Gibbs function for seawater thermodynamics for -6 to 80
  dgC and salinity up to 120 g/kg. Deep-Sea Res. I, 55, 1639-1671.

  McDougall, T.J., R. Feistel, F. J. Millero, D. R. Jackett, D. G. Wright,
  B. A. King, G. M. Marion, C.-T. A. Chen and P. Spitzer, 2009: Calculation of
  the Thermophysical Properties of Seawater, Global Ship-based Repeat
  Hydrography Manual, IOCCP Report No. 14, ICPO Publication Series no. 134.
}
\keyword{utilities}

