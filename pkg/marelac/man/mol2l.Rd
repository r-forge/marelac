\name{mol2l}
\alias{mol2l}
\title{mol to liter conversion for a gas}
\description{Converts from liter to moles for a gas}
\usage{mol2l(x=1, T=25, P=1, gas=NULL, a=0, b=0)}
\arguments{
  \item{x }{mol of the gas}
  \item{T }{temperature, dg celsius}
  \item{P }{pressure, atmosphere}
  \item{gas}{the gas, if NULL then the coefficients a and b are used, if specified, one of "O2","He","N2" or "H2"}
  \item{a }{Species-specific coefficient, \eqn{dm^6*bar/mol^2}}
  \item{b }{Species-specific coefficient, \eqn{dm^3/mol}}
}
\value{
  volume of the gas, liter
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{#22.4 Liters of Oxygen = 1.00 moles of Oxygen
mol2l(a=1.382,b=0.03186,T=0)

mol2l(T=0,gas="O2")

# the same but with standard coefficients
mol2l(T=0)}

\seealso{\code{\link{l2mol}} the inverse function}
\note{the coefficients a and b are species-specific; values of 0 assume an ideal gas and
in general give good estimates

  For the gasses He, H2, N2 and O2, the following values of a and b are used:

       \tabular{rlll}{
           \tab He \tab 0.034598   \tab 0.023733  \cr
           \tab H2 \tab 0.24646    \tab 0.026665   \cr
           \tab N2 \tab 1.3661     \tab 0.038577  \cr
           \tab oxygen \tab 1.3820     \tab 0.03186
       }}
\keyword{utilities}

