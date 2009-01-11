\name{atmComp}
\alias{atmComp}
\title{Atmospheric Gas Composition}
\description{
  The atmospheric composition, given in units of the moles of each gas
  to the total number of moles of gas in dry air (the mixing ratio).
}
\usage{
atmComp(x = c("He", "Ne", "N2", "O2", "Ar", "Kr", "CH4", "CO2",
  "N2O", "H2", "Xe", "CO", "O3"))
}
\arguments{
  \item{x}{character vector with gas(ses) whose composition should be
    estimated.
  }
}
\value{
  A vector with the composition of the gases.
}
\references{
  Sarmiento J.L. and N. Gruber, 2006. Ocean Biogeochemical
  Dynamics. Princeton University Press, Princeton. p 74.  who cite Weast
  and Astle (1982) for all gasses except CO2, CH4 and N2O, which are
  1998 concentrations taken from Ramaswamy et al., 2001.
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}

\examples{
atmComp("O2")
atmComp(c("O2", "CH4", "CO2", "N2O"))
atmComp()
}
\keyword{datasets}
