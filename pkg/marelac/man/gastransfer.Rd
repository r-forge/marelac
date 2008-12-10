\name{gastransfer}
\alias{gastransfer}
\title{The gas transfer coefficient, in cm/hr}
\description{The gas transfer coefficient, in cm/hr, for certain gases in seawater (S=35)}
\usage{gastransfer(T=25,u10=1,gas="O2",
      method=c("Liss","Nightingale","Wanninkhof1","Wanninkhof2"),
      Schmidt=schmidt(T=T,gas=gas))
}
\arguments{
  \item{T }{Temperature, dgC}
  \item{u10 }{wind speed, in m/sec at a nominal height of 10 m above sea level}
  \item{gas }{The gas, one of "He","Ne","N2","O2","Ar","Kr","Rn","CH4","CO2","N2O","CCl2F2","CCL3F","SF6","CCl4"}
  \item{method}{one of "Liss", for Liss and Merlivat, 1986; "Nightingale", for Nightingale et al., 2000;
  "Wanninkhof1", for Wanninkhof 1992, or "Wanninkhof2" for Wanninkhof and McGills 1999.}
  \item{Schmidt}{the Schmidt number, when given this overrules the arguments \code{gas} and \code{T}}
}
\value{
The gas transfer velocity, for seawater, in cm/hr }

\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
useq <- 0:15
plot(useq,gastransfer(u10=useq),type="l",lwd=2,xlab="u10,m/s",
     ylab="cm/hr", main="O2 gas transfer velocity",ylim=c(0,75))
lines(useq,gastransfer(u10=useq,method="Nightingale"),lwd=2,lty=2)
lines(useq,gastransfer(u10=useq,method="Wanninkhof1"),lwd=2,lty=3)
lines(useq,gastransfer(u10=useq,method="Wanninkhof2"),lwd=2,lty=4)

legend("topleft",lty=1:4,lwd=2,legend=c("Liss and Merlivat 1986","Nightingale et al. 2000",
 "Wanninkhof 1992","Wanninkhof and McGills 1999"))
 }
\references{

Sarmiento J.L. and N. Gruber, 2006. Ocean Biogeochemical Dynamics. Princeton University Press, Princeton. p 85.
Liss P.S. and L. Merlivat, 1986. Air-sea gas exchange rates: introduction and synthesis.
In: the role of air-sea exchange in Geochemical cycling, edited by P. Buat-Menard, pp 113-127.
D. Reidel, Dordrecht, the Netherlands.

Nightingale et al. 2000. In situ evaluation of air-sea gas exchnage prameterizations using novel
conservatie and volatile tracers. Global biogeochemical cycles 14, 373-387.

Wanninkhof R. 1992. Relationship between wind speed and gas exchange over the ocean. Journ. Geophys. Res. 97, 7373-7383.

Wanninkhof R. and W. McGillis 1999. A cubic relationshp between ari-sea CO2 exchange and wind speed.
Geophys. Res. Lett. 26, 1889-1892. }

\seealso{\code{\link{schmidt}}}

\keyword{utilities}

