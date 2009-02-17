\name{AtomicWeight}
\docType{data}
\alias{AtomicWeight}
\title{The Atomic Weights of Chemical Elements}
\description{
  Atomic weights of chemical elements useful for marine and lacustrine sciences.
}
\usage{AtomicWeight}
\format{
  Atomic weights of elements are given as list. The list index is the atomic number.
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl> and Thomas Petzoldt}
\details{
  Molecular weights of chemical elements may vary due to different isotope
  compositions, depending on geology, industrial processes or biological activity.
  Please consult the IUPAC Technical report about the details.
  The data set contains \code{NA} for elements which have no stable isotopes
  (except U, Th, Pa).
}
\seealso{
  \code{\link{mol.weight}}
}
\examples{
AtomicWeight
AtomicWeight$H
(W_H2O<- with (AtomicWeight, 2 * H + O))
}
\references{
  Wieser, M.E. (2006). Atomic weights of the elements 2005 (IUPAC
  Technical Report). Pure Appl. Chem., Vol. 78, No. 11, pp. 2051--2066,
  2006. \url{http://dx.doi.org/10.1351/pac200678112051}
}
\keyword{datasets}


