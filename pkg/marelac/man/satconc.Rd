\name{satconc}
\alias{satconc}
\title{Saturated concentrations of gasses in water}
\description{Calculates the saturated concentration of several gasses in water for s given temperature, salinity and pressure}
\usage{satconc(S = 35, T = 25, P = 1, gas = "O2", atm = atmComp(gas))}
\arguments{
  \item{S }{Salinity, -}
  \item{T }{Temperature, dgC}
  %\item{K }{Temperature, dgK}
  \item{P }{Pressure, atm}
  \item{gas }{The gas, one of "He","Ne","N2","O2","Ar","Kr","CH4","CO2","N2O"}
  \item{atm }{the number of moles of the gas per unit mole of air, the "mixing ratio".
  When present, this overrules the \code{gas} argument. When
  unspecified, the value from \code{AtmComp} for the gas is taken.}
}
\value{
   the saturated concentration of the gass in \eqn{mmol m^{-3}}{mmol/m3}.
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{satconc(gas="O2")
Temp<-seq(from=0,to=30,by=0.1)
Sal <- seq(from=0,to=35,by=0.1)

mf  <-par(mfrow=c(1,2))
plot(Temp,satconc(T=Temp,gas="N2"),xlab="temperature",log="y",
     ylab="mmol/m3",main="Saturated conc (S=35)",type="l",lwd=2,ylim=c(1,700))
lines(Temp,satconc(T=Temp,gas="CO2"),col="red",lwd=2)
lines(Temp,satconc(T=Temp,gas="O2"),col="blue",lwd=2)
lines(Temp,satconc(T=Temp,gas="CH4"),col="green",lwd=2)
lines(Temp,satconc(T=Temp,gas="N2O"),col="yellow",lwd=2)

plot(Sal,satconc(S=Sal,gas="N2"),xlab="salinity",log="y",
     ylab="mmol/m3",main="Saturated conc (T=20)",type="l",lwd=2,ylim=c(1e-3,700))
lines(Sal,satconc(S=Sal,gas="CO2"),col="red",lwd=2)
lines(Sal,satconc(S=Sal,gas="O2"),col="blue",lwd=2)
lines(Sal,satconc(S=Sal,gas="CH4"),col="green",lwd=2)
lines(Sal,satconc(S=Sal,gas="N2O"),col="yellow",lwd=2)
legend("right",col=c("black","red","blue","green","yellow"),lwd=2,
       legend=c("N2","CO2","O2","NH4","N2O"))

par(mfrow=mf)

# generate table 3.2.4 from Sarmiento and Gruber
Temp <- seq (0,30,by=5)
# saturated concentrations in mmol/m3
A <- data.frame(rbind( T=Temp,
            N2=satconc(T=Temp,gas="N2"),
            O2=satconc(T=Temp,gas="O2"),
            CO2=satconc(T=Temp,gas="CO2"),
            Ar=satconc(T=Temp,gas="Ar")))
format(A,digits=4)
# table values at 20 dg C:
# 425.7   230.5  11.61  11.29
# note the deviations for CO2 !

# saturated concentrations in micromol/m3
AA <- data.frame(rbind( T=Temp,
            N2O=satconc(T=Temp,gas="N2O")*1000,
            Ne=satconc(T=Temp,gas="Ne")*1000,
            Kr=satconc(T=Temp,gas="Kr")*1000,
            CH4=satconc(T=Temp,gas="CH4")*1000,
            He=satconc(T=Temp,gas="He")*1000))
format(AA,digits=4)
# table values at 20 dgC:
# 7.16  6.94  2.50  2.12  1.70
# Note: different for N2O

}
\note{Compared to the table in Sarmiento and Gruber, there is a slight deviation for N2O ,and He;
      CO2 is OK for Temperature 0 only}
\references{Sarmiento J.L. and N. Gruber, 2006. Ocean Biogeochemical Dynamics. Princeton University Press, Princeton. p 74

who cite:

for he and Ne: Weiss, R. 1971. Solubility of helium and neon in water and seawater. Journ. Chem. Eng. Data 16, 235-241.

N2, O2 and Ar: Weiss, R. 1970. The solubility of nitrogen, oxygen, and argon in water and seawater. Deep-Sea Res. 17: 721-35.

Kr: Weiss, R.  and T.K. Kyser, 1978. Silubility of Krypton in water and seawater. Journ. Chem. Eng. Data 23, 69-72.

Rn: Hackbusch 1979. Eine methode zur Bestimmung der Diffusions, Loslichkeits un Permeabilitats Konstanten von Radon-222 in Wasser und Meereswasser. PhD thesis, Univeristy of Heidelberg, Germany.

CH4: Wiesenburg, D.A. and J.N.L. Guinasso 1979. Equilibrium solubilities of methane, carbon monoxide and hydrogen in water and sea water. Journ. Chem. Eng. Data 24, 256-360.

CO2 and N2O: Weiss, R. and B.A. Price, 1980. Nitrous oxide solubility in wate and sewater. Mar. Chem. 8, 347-359.

CFC-11 and CFC-12: Warner, M.J. and R. Weiss, 1985. Solubilities of chlorofluorocarbons 11 and 12 in water and sewater. Deep-Sea Res. 32, 1485-1497.

SF6: Bullister et al 2002. The solubility of sulfur hexafluroide in water and sewater. Deep-Sea Res. I, 49: 175-188.

CCl4: Bullister J.L. and D.P. Wisegarver 1998.  The solubility of carbon tetrachloride in water and seawater. Deep-Sea Res. I, 1285-1302.
}
\keyword{utilities}

