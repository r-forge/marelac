\name{solubility}
\alias{solubility}
\title{Solubility parameters}
\description{Solubility parameters SA, mmol/m3/atm, calculated from the Bunsen solubility coefficients
and the volumetric solubility coefficients.
}
\usage{solubility(S = 35, T = NULL, K = NULL, P = 1, gas = "O2")}
\arguments{
  \item{S }{Salinity, -}
  \item{T }{Temperature, dgC}
  \item{K }{Temperature, dgK}
  \item{P }{Pressure, atm}
  \item{gas }{The gas, one of "He","Ne","N2","O2","Ar","Kr","Rn","CH4","CO2","N2O","CCl2F2","CCL3F","SF6","CCl4"}
}
\value{
  the solubility, mmol/m3/atm            }
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{solubility(T=1:20,S=35,P=1,gas="CO2")
solubility(T=0:5,S=35,P=1,gas="O2")

Temp<-seq(from=0,to=30,by=0.1)
mf<-par(mfrow=c(1,2))

plot(Temp,solubility(T=Temp,gas="CCl4"),xlab="temperature",
     ylab="mmol/m3/atm",main="solubility (S=35)",type="l",lwd=2,ylim=c(0,100000))
lines(Temp,solubility(T=Temp,gas="CO2"),col="red",lwd=2)
lines(Temp,solubility(T=Temp,gas="N2O"),col="blue",lwd=2)
lines(Temp,solubility(T=Temp,gas="Rn"),col="green",lwd=2)
lines(Temp,solubility(T=Temp,gas="CCl2F2"),col="yellow",lwd=2)

legend("topright",col=c("black","red","blue","green","yellow"),lwd=2,
       legend=c("CCl4","CO2","N2O","Rn","CCl2F2"))

plot(Temp,solubility(T=Temp,gas="Kr"),xlab="temperature",
     ylab="mmol/m3/atm",main="solubility (S=35)",type="l",lwd=2,ylim=c(0,4000))
lines(Temp,solubility(T=Temp,gas="CH4"),col="red",lwd=2)
lines(Temp,solubility(T=Temp,gas="Ar"),col="blue",lwd=2)
lines(Temp,solubility(T=Temp,gas="O2"),col="green",lwd=2)
lines(Temp,solubility(T=Temp,gas="N2"),col="yellow",lwd=2)
lines(Temp,solubility(T=Temp,gas="Ne"),col="grey",lwd=2)

legend("topright",col=c("black","red","blue","green","yellow","grey"),lwd=2,
       legend=c("Kr","CH4","Ar","O2","N2","Ne"))

par(mfrow=mf)

# generate table 3.2.3 from Sarmiento and Gruber
Temp <- seq (0,30,by=5)
# saturated concentrations in mmol/m3
A <- data.frame(rbind( T=Temp,
            He=solubility(T=Temp,gas="He"),
            Ne=solubility(T=Temp,gas="Ne"),
            N2=solubility(T=Temp,gas="N2"),
            O2=solubility(T=Temp,gas="O2"),
            Ar=solubility(T=Temp,gas="Ar"),
            Kr=solubility(T=Temp,gas="Kr"),
            Rn=solubility(T=Temp,gas="Rn"))  )
format(A,digits=4)
# table values at 20 dg C:
# 332.9   390.7   557.9  1126  1236  2241  14130
# note the slight deviations for Rn

# saturated concentrations in micromol/m3
AA <- data.frame(rbind( T=Temp,
            CH4=solubility(T=Temp,gas="CH4"),
            CO2=solubility(T=Temp,gas="CO2"),
            N2O=solubility(T=Temp,gas="N2O"),
            CCL2F2=solubility(T=Temp,gas="CCl2F2"),
            CCL3F=solubility(T=Temp,gas="CCl3F"),
            SF6=solubility(T=Temp,gas="SF6"),
            CCl4=solubility(T=Temp,gas="CCl4")))
format(AA,digits=4)

# Table values at 20 dgC:
# 1241  33110 23870  2566  9242   195.8, 30307
# Note: there are slight deviations for CO2, N2O and He!

}
\note{The molar volume used for the Bensen coefficient conversion
is the ideal gas value of 22.4136 l/mol.

These coefficients are to be used with pAmoist, the partial
pressure of the gas in moist air.

To convert them for use with partial pressure in dry air, divide by (1-vapor(S,T))}
\references{Sarmiento J.L. and N. Gruber, 2006. Ocean Biogeochemical Dynamics. Princeton University Press, Princeton. p 74

who cite:

for He and Ne: Weiss, R. 1971. Solubility of helium and neon in water and seawater. Journ. Chem. Eng. Data 16, 235-241.

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

