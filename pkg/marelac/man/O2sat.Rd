\encoding{latin1}
\name{O2sat}
\alias{O2sat}

\title{Saturation concentration of oxygen in water}
\description{
  Empirical formulae that can be used to compute saturation concentration of
  oxygen in water in \eqn{mg/L}{mg L^{-1}}
}
% Note the default value of S. Should we change it to 35?
\usage{
O2sat(S = 0, T = 25, NN = 0, method = c("Weiss", "APHA", "Paul"))
}

\arguments{
  \item{S}{salinity (dimensionless, for method "Weiss" only),}
  \item{T}{Temperature in \eqn{^\circ}{degrees }C,}
  %\item{K}{absolute Temperature in degrees Kelvin,}
  \item{NN}{height above sea level (in m, for method "Paul" only),}
  \item{method}{formula to be used, see references for correct application.}
}
% Karline, correct this statement if wrong
\details{
  Method APHA is the standard formula in Limnology, method Weiss the
  standard formula in marine sciences. The method after Paul is a 
  simple formula fitted on available tables.
  To avoid confusion between the arguments (S, T, NN) it is advisable to use
  named arguments in general, e.g. \code{O2sat(T=4)}.
}
\value{
  Vector with oxygen saturation concentration in \eqn{mg L^{-1}}{mg/L}.
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
und Schlammuntersuchung; Gasf\enc{ö}{oe}rmige Bestandteile (Gruppe G); Bestimmung des
Sauerstoffs\enc{ä}{ae}ttigungsindex (G 23).

Paul, L. (1985). Das thermische Regime der Talsperre Saidenbach und einige 
Beziehungen zwischen abiotischen und biotischen Komponenten. Dissertation, 
TU Dresden, Fakult\enc{ä}{ae}t Bau-, Wasser- und Forstwesen. 84 pp.

Weiss, R. (1970). The solubility of nitrogen, oxygen, and argon in water and
seawater. Deep-Sea Res. 17: 721-35.

Wagner, R. (1979): Die Praxis der Bestimmung des biochemischen Sauerstoffbedarfs - 
Ergebnisse einer Umfrage (Berichtigung und Erg\enc{ä}{ae}nzung zur
Ver\enc{ö}{oe}ffentlichung). Vom Wasser 53, S. 283-285.
}

\seealso{\code{\link{satconc}}}

\examples{
  O2sat(T = 20)
  T <- seq(0, 30, 0.1)
  plot(T, O2sat(T = T, method="APHA"), type = "l", ylim = c(0, 15), lwd=2)
  lines(T, O2sat(S = 0, T = T, method = "Weiss"), col = "yellow", lwd = 2, lty = "dashed")
  lines(T, O2sat(S = 35, T = T, method = "Weiss"), col = "red", lwd = 2)
}
\keyword{ misc }

