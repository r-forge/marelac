\name{AtomicWeight}
\docType{data}
\alias{AtomicWeight}
\alias{atomicweight}
\title{The Atomic Weights of Chemical Elements}
\description{
  Atomic weights of chemical elements useful for marine and lacustrine sciences.
}
\usage{
AtomicWeight
atomicweight
}
\format{
  Data frame with the IUPAC table of atomic weights.
  The data frame contains the following colums: AtomicNumber, Name,
  Symbol, AtomicWeight, Footnotes. See IUPAC report for details.
  The lower case version is in list-only format that allows symbolic
  computations with elements.
}
\author{Karline Soetaert and Thomas Petzoldt}
\details{
  Molecular weights of chemical elements may vary due to different isotope
  compositions, depending on geology, industrial processes or biological activity.
  Please consult the IUPAC Technical report about the details.
  The data set contains \code{NA} for elements which have no stable isotopes
  (except U, Th, Pa).
}
\seealso{
  \code{\link{molweight}} for molecular weight
}
\examples{
}
\references{
  Wieser, M.E. (2006). Atomic weights of the elements 2005 (IUPAC
  Technical Report). Pure Appl. Chem., Vol. 78, No. 11, pp. 2051--2066,
  2006. \url{http://dx.doi.org/10.1351/pac200678112051}
}
\keyword{datasets}


