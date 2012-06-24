\name{ImageOcean}
\docType{data}
\alias{ImageOcean}
\title{Image of the Ocean bathymetry}
\description{Realtively crude Image of the Ocean bathymetry, at 1 dg intervals.
Uses dataset \link{Bathymetry} from the R-package \code{marelac}.
}
\usage{
ImageOcean(...)

}
\arguments{
  \item{\dots}{ arguments passed to function \link{Image}. }
}  
  \author{Karline Soetaert <karline.soetaert@nioz.nl>}

\examples{

ImageOcean()
ImageOcean( asp = TRUE, xlab = "dg E", ylab = "dg N", add.contour = TRUE)
ImageOcean( xlim = c(-90, 90), NAcol = grey(0.4))

}
\keyword{hplot}

