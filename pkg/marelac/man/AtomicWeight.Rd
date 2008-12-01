\name{AtomicWeight}
\docType{data}
\alias{AtomicWeight}
\title{The atomic weights of chemical elements}
\description{Atomic weights of chemical elements useful for marine and lacustrine sciences}
\usage{AtomicWeight}
\format{
  a list
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl> and Thomas Petzoldt}
\details{
  Molecular weights of chemical elements may vary due to different isotope
  compositions, depending on geology, industrial processes or biological activity.
  Please consult the IUPAC Technical report about the details.
  The data set contains NA for elements which have no stable isotopes
  (except U, Th, Pa).
}
\seealso{\code{\link{mol.weight}}, \code{\link{g2mol}}}
\examples{AtomicWeight
AtomicWeight$H
(W_H2O<- with (AtomicWeight, 2 * H + O))
}
\references{Wieser, M.E. (2006). Atomic weights of the elements 2005 (IUPAC Technical Report).
 Pure Appl. Chem., Vol. 78, No. 11, pp. 2051--2066, 2006. doi:10.1351/pac200678112051
}
\keyword{datasets}


