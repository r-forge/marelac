\name{O2sat}
\alias{O2sat}

\title{Saturation concentration of oxygen in water}
\description{
  Empirical formulae that can be used to compute saturation concentration of
  oxygen in water in \eqn{mg/L}{mg L^{-1}}
}
\usage{
O2sat(Tc = T - 273.15, T = Tc + 273.15, NN = 0, S = 0, method = c("APHA", "Weiss", "Paul"))
}

\arguments{
  \item{Tc}{Temperature in degrees centigrade,}
  \item{T}{Temperature in degrees Kelvin,}
  \item{NN}{height above sea level (for method "Paul" only),}
  \item{S}{salinity in PSU (for method "Weiss" only),}
  \item{method}{formula to be used, see references for correct application.}
}
\details{
  Method APHA is the standard formula in Limnology, method Weiss may be the 
  standard formula in marine sciences. The method after Paul is a 
  simple formula fitted on available tables.
}
\value{
  Vector with oxygen saturation concentration in mg/L.
}
\references{

American Public Health Association, Inc. (1985): Standard Methods for the 
Examination of Water and Wastewater. 18th edition, 1992.

Benson B. B. and Krause, D. (1980): The concentration and isotopic fractionation 
of gases dissolved in freshwater in equilibrium with the atmosphere. I. Oxygen. 
Limnology and Oceanography 25, 662-671.

Brown, L.C. and Barnwell, T. O. Jr. (1987): The Enhanced Stream Water Quality 
Models QUAL2E and QUAL2E-UNCAS: Documentation and User Manual. U.S. Environmental 
Protection Agency, Athens, Georgia. EPA/600/3-87/007, p. 41. http://www.epa.gov)

DIN 38408-23, Ausgabe:1987-11: Deutsche Einheitsverfahren zur Wasser-, Abwasser- 
und Schlammuntersuchung; Gasf{\"o}rmige Bestandteile (Gruppe G); Bestimmung des
Sauerstoffs{\"a}ttigungsindex (G 23). \url{http://www.din.de}

Paul, L. (1985). Das thermische Regime der Talsperre Saidenbach und einige 
Beziehungen zwischen abiotischen und biotischen Komponenten. Dissertation, 
TU Dresden, Fakult{\"a}t Bau-, Wasser- und Forstwesen. 84 pp.

Weiss, R. (1970). The solubility of nitrogen, oxygen, and argon in water and
seawater. Deep-Sea Res. 17: 721-35.

Wagner, R. (1979): Die Praxis der Bestimmung des biochemischen Sauerstoffbedarfs - 
Ergebnisse einer Umfrage (Berichtigung und Erg{\"a}nzung zur Ver{\"o}ffentlichung).
Vom Wasser 53, S. 283-285. 
}


\examples{
  Tc <- seq(0, 30, 0.1)
  plot(Tc, O2sat(Tc), type="l", ylim=c(0,15))
  lines(Tc, O2sat(Tc, S=30, method="Weiss"))
}
\keyword{ misc }

