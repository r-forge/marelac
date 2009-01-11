\name{g2mol}
\alias{g2mol}
\title{Conversion from Gram to }
\description{Converts Gram to Mol of a chemical species.}
\usage{
g2mol(species)
}
\arguments{
  \item{species }{character vector with chemical species whose
    weight should be be converted to mol.}
}
\value{
  mol/g
}
\seealso{
  \code{\link{mol.weight}} the inverse function and
  \code{\link{AtomicWeight}} for atomic weights of the elements.
}
\examples{
g2mol("CO3")
g2mol("HCO3")
g2mol("H")
g2mol("H3PO4")

## works also with vectors
g2mol(c("C2H5OH", "CO2", "H2O"))

## Note that chemical formulae are case-sensitive
g2mol("Co") # cobalt
g2mol("CO") # carbon monoxide
}
\note{See \code{\link{mol.weight}} for more information}
\keyword{utilities}

