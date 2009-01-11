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
  
  Many of the functions are from the UNESCO 1983 paper. Note that in
  this paper, pressure is hydrostatic or gauge pressure, i.e. pressure
  relative to atmospheric pressure, and expressed in dbar. In the
  functions of \code{marelac}, pressure is either true pressure, or
  hydrostatic pressure, both expressed in bar.
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

  For many other fundamental properties of seawater: Fofonoff, P. and
  Millard, R.C. Jr. Unesco 1983. Algorithms for computation of
  fundamental properties of seawater.  Unesco Tech. Pap. in Mar. Sci.,
  No. 44, 53 pp.
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
