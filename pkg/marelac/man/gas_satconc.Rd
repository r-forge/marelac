\encoding{latin1}
\name{gas_satconc}
\alias{gas_satconc}
\title{Saturated Concentrations of Gases in Water}
\description{
  Calculates the saturated concentration of several gases in water for a
  given temperature, salinity and pressure.
}
\usage{
gas_satconc(S = 35, t = 25, P = 1.013253,
        x =c("He","Ne","N2","O2","Ar","Kr","CH4","CO2","N2O",
        "H2","Xe", "CO", "O3"), atm = atmComp(x))}
\arguments{
  \item{S }{Salinity (dimensionless),}
  \item{t }{Temperature, \eqn{^\circ}{degrees }C,}
  %\item{K }{Temperature, dgK}
  \item{P }{True pressure, bar}
  \item{x }{The gas }
  \item{atm }{the number of moles of the gas per unit mole of air in the
    atmosphere, the "mixing ratio".  When present, this overrules the
    \code{x} argument. When unspecified, the value from \code{atmComp}
    for the gas is taken.
  }
}
\value{
   The saturated concentration of the gass in \eqn{mmol m^{-3}}{mmol/m3}.
}
\author{Karline Soetaert < k.soetaert@nioo.knaw.nl > }
\examples{
gas_satconc(x = "O2")
Temp <- seq(from = 0,to = 30,by = 0.1)
Sal  <- seq(from = 0,to = 35,by = 0.1)

mf  <- par(mfrow = c(1,2))
plot(Temp, gas_satconc(t = Temp, x = "N2"), xlab = "temperature",
  log = "y", ylab = "mmol/m3", main = "Saturated conc (S=35)",
  type = "l", lwd = 2, ylim = c(1, 700))
lines(Temp,gas_satconc(t = Temp, x = "CO2"), col = "red",    lwd = 2)
lines(Temp,gas_satconc(t = Temp, x = "O2"),  col = "blue",   lwd = 2)
lines(Temp,gas_satconc(t = Temp, x = "CH4"), col = "green",  lwd = 2)
lines(Temp,gas_satconc(t = Temp, x = "N2O"), col = "yellow", lwd = 2)

plot(Sal, gas_satconc(S = Sal, x = "N2"), xlab = "salinity",
  log = "y", ylab = "mmol/m3", main = "Saturated conc (t=20)",
  type = "l", lwd = 2, ylim = c(1e-3, 700))

lines(Sal,gas_satconc(S = Sal, x = "CO2"), col = "red",    lwd = 2)
lines(Sal,gas_satconc(S = Sal, x = "O2"),  col = "blue",   lwd = 2)
lines(Sal,gas_satconc(S = Sal, x = "CH4"), col = "green",  lwd = 2)
lines(Sal,gas_satconc(S = Sal, x = "N2O"), col = "yellow", lwd = 2)
legend("right",col = c("black", "red", "blue", "green", "yellow"), lwd = 2,
       legend = c("N2", "CO2", "O2", "NH4", "N2O"))

par(mfrow = mf)

## generate table 3.2.4 from Sarmiento and Gruber
Temp <- seq (0, 30, by = 5)
## saturated concentrations in mmol/m3
A <- data.frame(rbind( t = Temp,
            N2 = gas_satconc(t = Temp, P=1, x = "N2"),
            O2 = gas_satconc(t = Temp, P=1, x = "O2"),
            CO2 = gas_satconc(t = Temp, P=1, x = "CO2"),
            Ar = gas_satconc(t = Temp, P=1, x = "Ar")))
format(A,digits = 4)
## table values at 20 dg C:
## 425.7   230.5  11.61  11.29
## note the deviations for CO2 !

## saturated concentrations in micromol/m3
AA <- data.frame(rbind( t = Temp,
            N2O = gas_satconc(t = Temp, P=1, x = "N2O")*1000,
            Ne = gas_satconc(t = Temp, P=1, x = "Ne")*1000,
            Kr = gas_satconc(t = Temp, P=1, x = "Kr")*1000,
            CH4 = gas_satconc(t = Temp, P=1, x = "CH4")*1000,
            He = gas_satconc(t = Temp, P=1, x = "He")*1000))
format(AA, digits = 4)
## table values at 20 dgC:
## 7.16  6.94  2.50  2.12  1.70
## Note: different for N2O
}
\note{
  Compared to the table in Sarmiento and Gruber, there is a slight
  deviation for N2O, and He CO2 is OK for temperature 0 only.
}
\references{
  Sarmiento J.L. and N. Gruber, 2006. Ocean Biogeochemical
  Dynamics. Princeton University Press, Princeton. p 74.

  who cite:

  for he and Ne: Weiss, R. 1971. Solubility of helium and neon in water
  and seawater. Journ. Chem. Eng. Data 16, 235-241.

  N2, O2 and Ar: Weiss, R. 1970. The solubility of nitrogen, oxygen, and
  argon in water and seawater. Deep-Sea Res. 17: 721-35.

  Kr: Weiss, R.  and T.K. Kyser, 1978. Silubility of Krypton in water
  and seawater. Journ. Chem. Eng. Data 23, 69-72.

  Rn: Hackbusch 1979. Eine Methode zur Bestimmung der Diffusions-,
  L\enc{ö}{oe}slichkeits un Permeabilitats Konstanten von Radon-222 in
  Wasser und Meereswasser. Dissertation, University of Heidelberg,
  Germany.

  CH4: Wiesenburg, D.A. and J.N.L. Guinasso 1979. Equilibrium
  solubilities of methane, carbon monoxide and hydrogen in water and sea
  water. Journ. Chem. Eng. Data 24, 256-360.

  CO2 and N2O: Weiss, R. and B.A. Price, 1980. Nitrous oxide solubility
  in wate and sewater. Mar. Chem. 8, 347-359.

  CFC-11 and CFC-12: Warner, M.J. and R. Weiss, 1985. Solubilities of
  chlorofluorocarbons 11 and 12 in water and sewater. Deep-Sea Res. 32,
  1485-1497.

  SF6: Bullister et al 2002. The solubility of sulfur hexafluroide in
  water and sewater. Deep-Sea Res. I, 49: 175-188.

  CCl4: Bullister J.L. and D.P. Wisegarver 1998.  The solubility of
  carbon tetrachloride in water and seawater. Deep-Sea Res. I,
  1285-1302.
}
\keyword{utilities}

