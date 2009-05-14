\name{convert_StoCl}
\alias{convert_StoCl}
\title{Salinity-Chlorinity Conversion}
\description{Estimates the chlorinity concentration based on salinity
}
\usage{convert_StoCl(S = 35)}
\arguments{
  \item{S }{salinity}
}
\value{
  The chlorinity concentration, g/kg
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
convert_StoCl(20)
}
\references{Cox, R.A., Culkin, F., Riley, J.P., 1967. The electrical
conductivity -- chlorinity relationship in natural seawater.
Deep--Sea Research 14, 203--220.
}
\keyword{utilities}

