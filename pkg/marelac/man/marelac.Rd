\name{marelac-package}
\alias{marelac-package}
\alias{marelac}
\docType{package}
\title{
Datasets, constants, conversion factors, and utilities for the
MArine, Riverine, Estuarine, LAcustrine and Coastal sciences
}
\description{
R-package marelac has been designed as a tool for use by scientists working in
the MArine, Riverine, Estuarine, LAcustrine and Coastal sciences.

It contains:

\item chemical and physical constants, e.g. atomic weights, gas constants.

\item conversion factors, e.g. from salinity to chlorinity, from mol to gram, etc..

\item functions, e.g. to estimate concentrations of conservative substances as a function of salinity, ...

}
\details{
  \tabular{ll}{
    Package: \tab marelac\cr
    Type: \tab Package\cr
    Version: \tab 1.4\cr
    Date: \tab 2008-11-031\cr
    License: \tab GNU Public License 2 or above\cr
  }

}
\author{
  Karline Soetaert (Maintainer)

  Thomas Petzoldt

  Filip Meysman
}
%\references{
%% to be filled
%}

\examples{
\dontrun{
## show examples (see respective help pages for details)
example(Constants)
example(mol.vol)

## open the directory with documents
browseURL(paste(system.file(package="marelac"), "/doc", sep=""))


}
}

\keyword{ package }
