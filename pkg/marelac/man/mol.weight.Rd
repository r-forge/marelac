\name{mol.weight}
\alias{mol.weight}
\title{Molecular Weight of a Chemical Species}
\description{Calculates the molecular weight of chemical species.}
\usage{
mol.weight(species)
}
\arguments{
  \item{species }{character vector with chemical species whose molecular
    weight is requested.
  }
}
\value{
  Vector with the molecular weights in g/mol.
}
\details{
  Molecular weights of chemical elements may vary due to different
  isotope compositions, depending on geology, industrial processes or
  biological activity.  Please consult the IUPAC Technical report about
  the details.  The function returns NA for elements (and their
  compounds) which have no stable isotopes (except U, Th, Pa).
}
\references{
  Wieser, M.E. (2006). Atomic weights of the elements 2005 (IUPAC
  Technical Report).  Pure Appl. Chem., Vol. 78, No. 11, pp. 2051--2066,
  2006. doi:10.1351/pac200678112051
}
\author{Thomas Petzoldt}
\note{
  This function uses text parsing of chemical formula, it is
  strictly case sensitive.
}
\seealso{
  \code{\link{AtomicWeight}}: Atomic weight data of stable elements.
}
\examples{
mol.weight("CO2")
mol.weight("HCO3")
mol.weight("H")
mol.weight("H3PO4")

mol.weight("CH3CH2CHCHCH2CHCHCH2CHCHCH2CHCHCH2CHCH(CH2)3COOH")
mol.weight("C20H30O2")

## works also with vectors
mol.weight(c("C2H5OH", "CO2", "H2O"))
mol.weight(c("SiOH4", "NaHCO3", "C6H12O6", "Ca(HCO3)2", "Pb(NO3)2", "(NH4)2SO4"))

## Note that chemical formulae are case-sensitive
mol.weight("Co") # cobalt
mol.weight("CO") # carbon monoxide


## from gram to mol
1/mol.weight("CO3")

}
\keyword{utilities}

