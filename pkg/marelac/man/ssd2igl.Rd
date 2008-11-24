\name{ssd2igl}
\alias{ssd2igl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert sunshine duration to global radiation}
\description{
  The function converts values of sunshine duration (in hours) to global radiation
  (in J/cm2/d)
}
\usage{
ssd2igl(S, doy, a = 0.25, b = 0.5, rho = 50.29)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Sunshine duration (hours)}
  \item{doy}{Julian day (for northern hemisphere only)}
  \item{a, b, rho}{site specific conversion parameters, must be fitted to
    measured data}
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(S, doy, a=0.25, b=0.5, rho=50.29) {
  xi <- 0.0172 * doy - 1.39
  S0 <- 12.3 + sin(xi) * (4.3 + (rho - 51))/6
  Rex <- 245 * (9.9 + 7.08 * sin(xi) + 0.18 * (rho - 51) * (sin(xi) - 1))
  Rex * (a + b * S/S0)
  }
}
\keyword{ misc }
