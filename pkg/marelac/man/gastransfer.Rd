\name{gastransfer}
\alias{gastransfer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gastransfer(T = 25, u10 = 1, gas = "O2", method = c("Liss", "Nightingale", "Wanninkhof1", "Wanninkhof2"), Schmidt = schmidt(T = T, gas = gas))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T}{ ~~Describe \code{T} here~~ }
  \item{u10}{ ~~Describe \code{u10} here~~ }
  \item{gas}{ ~~Describe \code{gas} here~~ }
  \item{method}{ ~~Describe \code{method} here~~ }
  \item{Schmidt}{ ~~Describe \code{Schmidt} here~~ }
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
function (T=25,u10=1,gas="O2",method=c("Liss","Nightingale","Wanninkhof1","Wanninkhof2"),
           Schmidt=schmidt(T=T,gas=gas))
{
 method <- match.arg(method)
 S600<-Schmidt/600
 tr <- switch(method, Liss =
 ifelse (u10<3.6,0.17*u10*Schmidt^(-2/3),ifelse (u10<=13, (u10-3.4)*2.8/sqrt(S600),
 (u10-8.4)*5.9/sqrt(S600))), Nightingale = (0.33*u10+0.222*u10*u10)/sqrt(S600),
  Wanninkhof1=0.31*u10^2*sqrt(S600),Wanninkhof2=0.0283*u10^3/sqrt(S600))
 attr(tr, "unit") = "cm/hr"
 tr
  }
}

\keyword{ misc }

