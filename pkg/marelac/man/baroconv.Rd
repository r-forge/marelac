\name{baroconv}
\alias{baroconv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Conversion Between Different Barometric Units}
\description{
  The function converts between different units of pressure.
}
\usage{
baroconv(x, unit = c("Pa", "bar", "at", "atm", "torr"))
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
baroconv(1, "atm")
baroconv(c(2, 3, 4.5), "bar")

baroconv(1, "atm")$Pa

}
\keyword{ misc }

