\name{marelacTeaching-package}
\alias{marelacTeaching-package}
\alias{marelacTeaching}
\docType{package}
\title{
 Datasets and tutorials for use in the
 MArine, Riverine, Estuarine, LAcustrine and Coastal sciences
}
\description{
  R-package marelacTeaching has been designed as a tool to help in teaching
  students from the MArine, Riverine, Estuarine, LAcustrine and Coastal sciences.

  It contains:
  \itemize{
    \item a tutorial for novice R-users, focussing on the use of R in scientific
      computing.
    \item data sets used in science teaching.
  }
}
\details{
  \tabular{ll}{
    Package: \tab marelacTeaching\cr
    Type: \tab Package\cr
    Version: \tab 1.1\cr
    Date: \tab 2009-11-11\cr
    License: \tab GNU Public License 2 or above\cr
  }


To open the lecture notes, copy-paste this sentence to R:

\code{
browseURL(paste(system.file(package="marelacTeaching"), "/doc/UsingR.pdf", sep=""))}

which will open a PDF generated with WORD (Microsoft)

or type \code{vignette("UsingR")} which will open a LaTeX-generated file.
}
\author{
  Karline Soetaert (Maintainer)

  Filip Meysman
}
%\references{
%% to be filled
%}

\examples{
\dontrun{
## show examples (see respective help pages for details)
example(Bathymetry)
example(SCOC)
example(Zoogrowth)

## open the directory with documents
browseURL(paste(system.file(package="marelacTeaching"), "/doc", sep=""))


## show document with the introduction to R - one file, a pdf generated with WORD:
browseURL(paste(system.file(package="marelacTeaching"), "/doc/UsingR.pdf", sep=""))

## the introduction to R - a vignette
vignette("UsingR")

## and the answers
vignette("Answers")

edit(vignette("UsingR"))
edit(vignette("Answers"))
}
}

\keyword{ package }