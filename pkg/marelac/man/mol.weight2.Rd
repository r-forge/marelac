\name{mol.weight2}
\alias{mol.weight2}
\title{Molecular weight of a species}
\description{Calculates the molecular weight of chemical species}
\usage{mol.weight2(species)}
\arguments{
  \item{species }{chemical species whose molecular weight should be estimated}
}
\value{
  the molecular weight in g/mol
}
\details{
  Molecular weights of chemical elements may vary due to different isotope
  compositions, depending on geology, industrial processes or biological activity.
  Please consult the IUPAC Technical report about the details.
  The function returns NA for elements (and their compounds) which have no
  stable isotopes (except U, Th, Pa).
}
\references{Wieser, M.E. (2006). Atomic weights of the elements 2005 (IUPAC Technical Report).
 Pure Appl. Chem., Vol. 78, No. 11, pp. 2051--2066, 2006. doi:10.1351/pac200678112051
}
\author{Thomas Petzoldt}
\note{
  This function is not restricted to a pre-defined set of biologically active species.
  As it uses text parsing of chemical formula, it is strictly case sensitive (in contrast
  to the original function in this package)!
}
\examples{
mol.weight("CO3")
mol.weight2("CO3")

mol.weight("HCO3")
mol.weight2("HCO3")

mol.weight("H")
mol.weight2("H")

mol.weight("H3PO4")
mol.weight2("H3PO4")

mol.weight2("H2O")

mol.weight2("NaHCO3")
mol.weight2("C6H12O6")
mol.weight2("Ca(HCO3)2")

mol.weight2("Pb(NO3)2")
mol.weight2("(NH4)2SO4")
mol.weight2("C2H5OH")
mol.weight2("CH3CH2CHCHCH2CHCHCH2CHCHCH2CHCHCH2CHCH(CH2)3COOH")
mol.weight2("C20H30O2")

mol.weight("SiOH4")
mol.weight2("SiOH4")

# the old function is not case sensitive
mol.weight("Co2")
mol.weight("CO2")

## but the new function is (and must be!)
mol.weight2("Co") # cobalt
mol.weight2("CO") # carbon monoxide

}

\keyword{utilities}

