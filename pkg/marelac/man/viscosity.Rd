\name{viscosity}
\alias{viscosity}
\title{shear viscosity of water}
\description{Calculates the shear viscosity of water, in centipoise (0.01 g/cm/sec).
  Valid for 0<T<30 \eqn{^\circ}{degrees }C,  0<S<36, 1<P<1000 bars.
 \cr Based on the code "CANDI" by B.P. Boudreau}
\usage{viscosity(T=25, S=35, P=1.013253)}
\arguments{
  \item{T }{Temperature, \eqn{^\circ}{degrees }C}
  \item{S }{Salinity, -}
  \item{P }{Pressure, bar}
}
\value{
  shear visocisity in centipoise
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{plot(0:30,viscosity(S=35,T=0:30,P=1),xlab="temperature",ylab="centipoise", 
       main="shear viscosity of water",type="l")
  lines(0:30,viscosity(S=0,T=0:30,P=1),col="red")
  lines(0:30,viscosity(S=35,T=0:30,P=100),col="blue")
  legend("topright",col=c("black","red","blue"),lty=1,
        legend=c("S=35,P=1","S=0,P=1","S=35,P=100"))}
\references{Based on the FORTRAN implementation of the diagenetic model "CANDI" of B.P. Boudreau:
\cr Boudreau BP,
A method-of-lines code for carbon and nutrient diagenesis in aquatic sediments. 
COMPUTERS & GEOSCIENCES 22 (5): 479-496 JUN 1996

Kulkula, DJ, Gebhart, B and Mollendorf, JC 1987. Thermodynamic and transport properties of pure
and saline water. Adv. Heat transfer 18, 325-363.}
\details{The details given in the original code by B. Boudreau are repeated here:
\cr Uses the equation given by Kukulka et al. (1987).}
\keyword{utilities}

