\name{AtmComp}
\alias{AtmComp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
AtmComp(x = c("He", "Ne", "N2", "O2", "Ar", "Kr", "CH4", "CO2", "N2O", "H2", "Xe", "CO", "O3"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ ~~Describe \code{x} here~~ }
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
function(x=c("He", "Ne", "N2", "O2", "Ar", "Kr", "CH4", "CO2","N2O", "H2", "Xe", "CO", "O3")) {
  x <- match.arg(x, several.ok = TRUE)
  ret <- list(He = 5.24e-06, Ne = 1.818e-05, N2 = 0.78084, O2 = 0.20946,
    Ar = 0.00934, Kr = 1.14e-06, CH4 = 1.745e-06, CO2 = 0.000365,
    N2O = 3.14e-07, H2 = 5.5e-07, Xe = 8.7e-08, CO = 5e-08, O3 = 1e-08)
  ret[x]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
