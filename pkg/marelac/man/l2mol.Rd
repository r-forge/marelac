\name{l2mol}
\alias{l2mol}
\title{liter to mol conversion for a gas}
\description{Converts from liter to moles for a gas}
\usage{l2mol(x=1, T=25, P=1, gas=NULL, a=0, b=0)}
\arguments{
  \item{x }{Volume of the gas, liter}
  \item{T }{Temperature, dg celsius}
  \item{P }{Pressure, atmosphere}
  \item{gas}{the gas, if NULL then the coefficients a and b are used, if specified, one of "O2","He","N2" or "H2"}
  \item{a }{Van der Waals constant a, a species-specific coefficient, \eqn{dm^6*bar/mol^2}}
  \item{b }{Van der Waals constant b, a species-specific coefficient, \eqn{dm^3/mol}}
}
\value{
   mol 
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{#22.4 Liters of Oxygen = 1.00 moles of Oxygen
l2mol(x=8,a=1.382,b=0.03186,T=0)*1000
l2mol(x=8,gas="O2",T=0)*1000
plot(x<-seq(1,8,by=0.1),l2mol(x=x,a=1.382,b=0.03186,T=0)*1000,
     xlab="ml/l",type="l",ylab="mmol/m3",main="Oxygen conversion")}
\seealso{\code{\link{mol2l}} the inverse function}
\note{the coefficients a and b are species-specific; values of 0 assume an ideal gas and
in general give good estimates

The values of the van der Waals constants are from
Weast. R. C. (Ed.)1972. Handbook of Chemistry and Physics (53rd Edn.), Cleveland:Chemical Rubber Co.

  For instance, for the gasses He, CO2, N2O and O2, the following values of a and b are used:

       \tabular{rlll}{
           \tab He \tab 0.03457   \tab 0.02370  \cr
           \tab CO2 \tab  3.64000    \tab 0.04267   \cr
           \tab N2O \tab 3.83200    \tab  0.04415 \cr
           \tab oxygen \tab 1.37800     \tab 0.03183
       }}
\keyword{utilities}

