\name{marelac-package}
\alias{marelac-package}
\alias{marelac}
\docType{package}
\title{
  Datasets, Constants, Conversion Factors, and Utilities for the
  MArine, Riverine, Estuarine, LAcustrine and Coastal Sciences
}
\description{
  
  R-package marelac has been designed as a tool for use by scientists
  working in the MArine, Riverine, Estuarine, LAcustrine and Coastal
  sciences.

  It contains:
  \itemize{
    \item chemical and physical constants, e.g. atomic weights, gas constants.
    \item conversion factors, e.g. from salinity to chlorinity, from mol to gram,
       etc..
    \item utility functions, e.g. to estimate concentrations of conservative
      substances as a function of salinity, ...
  }
  
  
 About the symbols used.
 
 
 Here we adopt the symbolism as in McDougall et al.,, 2009:
 
 \itemize{
 \item S for absolute salinity, (g/kg)
 \item P for absolute (total) pressure (bar)
 \item p for sea pressure (also called gauge or applied pressure (bar), the
 pressure relative to P0, one standard atmosphere (=1.01325 bar)
 \item t for temperature in \eqn{^\circ}{degrees }C
 \item T for absolute temperature, in \eqn{^\circ}{degrees }K; T = t+273.15

 }
  Many of the functions are from the UNESCO 1983 paper, or from Feistel, 2008.
  Note that in these papers, pressure is expressed in dbar.
}
\details{
  \tabular{ll}{
    Package: \tab marelac\cr
    Type: \tab Package\cr
    Version: \tab 1.4\cr
    Date: \tab 2009-01-12\cr
    License: \tab GNU Public License 2 or above\cr
  }
}
\author{
  Karline Soetaert (Maintainer)

  Filip Meysman
  
  Thomas Petzoldt
}
\references{
  
  For solubilities, atmospheric composition, gas exchange coefficients:

  Sarmiento J.L. and N. Gruber, 2006. Ocean Biogeochemical
  Dynamics. Princeton University Press, Princeton. p 74

  For diffusion coefficients, viscosity

  Boudreau BP, 1996.  A method-of-lines code for carbon and nutrient
  diagenesis in aquatic sediments.  COMPUTERS & GEOSCIENCES 22 (5):
  479-496.

  For many other fundamental properties of seawater: either the UNESCO report (1983)

  Fofonoff, P. and
  Millard, R.C. Jr. Unesco 1983. Algorithms for computation of
  fundamental properties of seawater.  Unesco Tech. Pap. in Mar. Sci.,
  No. 44, 53 pp.
  
  or the more recent report and papers:

  McDougall, T.J., R. Feistel, F. J. Millero, D. R. Jackett, D. G. Wright, B. A. King, G. M. Marion, C.-T. A. Chen
and P. Spitzer, 2009: Calculation of the Thermophysical Properties of Seawater, Global Ship-based Repeat
Hydrography Manual, IOCCP Report No. 14, ICPO Publication Series no. 134.

Feistel, R., 2008: A Gibbs function for seawater thermodynamics for -6 to 80 dgC and salinity up
to 120 g/kg. Deep-Sea Res. I, 55, 1639-1671.

  Millero, F.J., Feistel, R., Wright, D.G., and McDougall, T.J., 2008: The composition of Standard
Seawater and the definition of the Reference-Composition Salinity Scale, Deep-Sea Res. I,
55, 50-72.


}

\examples{
\dontrun{
## show examples (see respective help pages for details)
example(Constants)
example(mol.vol)

## open the directory with documents
browseURL(paste(system.file(package="marelac"), "/doc", sep=""))

## open the UNEXCO 1983 document
browseURL(paste(system.file(package="marelac"), "/doc/UnescoMarine.pdf", sep=""))
}
}

\keyword{ package }
