\name{convert_p}
\alias{convert_p}
\title{Conversion Between Different Barometric Units}
\description{
  The function converts between different units of pressure.
}
\usage{
convert_p(x, unit = c("Pa", "bar", "at", "atm", "torr"))
}

\arguments{
  \item{x}{vector of given pressure values,}
  \item{unit}{measurement unit of the given value(s).}
}

\value{
  A data frame with converted values.
}
\references{\url{http://en.wikipedia.org/wiki/Bar_(unit)}}

%\seealso{}
\examples{
convert_p(1, "atm")
convert_p(c(2, 3, 4.5), "bar")

convert_p(1, "atm")$Pa

}
\keyword{ misc }

