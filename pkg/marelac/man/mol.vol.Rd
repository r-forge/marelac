\name{mol.vol}
\alias{mol.vol}
\title{Mol to Liter Conversion for a Gas}
\description{Converts from liter to moles for a gas.}
\usage{mol.vol(t = 25, P = 1.013253,
   x = c("ideal","Ar","CO2","CS2","CO","CCl4","Cl2","C2H6S","C2H5OH",
         "C6H5F","CH3F","CH4","CH3OH","C5H12","C3H8","H2O",
         "He","H2","HBr","HCl","H2S","Hg","Kr","NH3","Ne",
         "NO","N2","NO2","N2O","O2","PH3","SiH4","SiF4","SO2","Xe"),
    q = 1, a = 0, b = 0)}
\arguments{
  \item{t }{temperature, \eqn{^\circ}{degrees }C}
  \item{P }{True pressure, bar}
  \item{x}{the gas, if NULL then the coefficients a and b are used}
  \item{q }{quantity, mol, of the gas}

  \item{a }{Van der Waals constant a, a species-specific coefficient,
    \eqn{dm^6*bar/mol^2}
  }
  \item{b }{Van der Waals constant b, a species-specific coefficient,
    \eqn{dm^3/mol}
  }
}
\value{
  volume of the gas, liter
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl> }
\examples{
#molecular volume of an ideal gas.

mol.vol(x=NULL) # 24.46559
plot(0:30, mol.vol(t = 0:30, x=NULL),
    xlab = "Temperature, dgC", ylab = "Molar volume")

#
mol.vol(a = 1.382, b = 0.03186, x=NULL, t = 0)

mol.vol(t = 0, x = "O2")

# the same but for all gasses
mol.vol(t = 0)

# table for different pressures
mol.vol(P = 1:5, x = "O2")

# the inverse function
1/mol.vol(x="O2")
}

\seealso{\code{\link{mol.weight}} the molecular weight}
\note{
  The coefficients a and b are species-specific; values of 0 assume an
  ideal gas and in general give good estimates.

  The values of the van der Waals constants are from
  Weast. R. C. (Ed.)1972. Handbook of Chemistry and Physics (53rd Edn.),
  Cleveland:Chemical Rubber Co.
  
  as found in \url{http://en.wikipedia.org/wiki/Van_der_Waals_constants_(data_page)}
  
  Use 1/mol.vol to convert from liter to moles.
  
  the default calculates the molar volume of an ideal gas
}
\keyword{utilities}

