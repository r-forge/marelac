\name{redfield}
\alias{redfield}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Redfield Ratio Calculator}
\description{
  Estimate elemental composition of biomass (or media) according to the Redfield ratio.
}
\usage{
redfield(element, x, method = c("mol", "mass"))
}
\arguments{
  \item{element}{The element that is given ("C", "H", "O", "N", "P"),}
  \item{x}{amount of substance of that element (in mol or mass units),}
  \item{method}{measurement unit ("mol" or "mass").}
}
\details{
  The Redfield ratio is C106 H80 O45 N16 P1.
}
\value{
  Estimated ratio of the main elements.
}

% Add **original** reference here, and not only N:P
\references{Redfield A.C. (1934). On the proportions of organic derivations
  in sea water and their relation to the composition of plankton.
  In James Johnstone Memorial Volume. (ed. R.J. Daniel).
  University Press of Liverpool, 177-192.
}

\seealso{\code{\link{mol.weight}}}
\examples{
# Redfield ratio
redfield("P", 1)
# returns the molar Redfield ratio, rescaled to nitrogen
redfield("N", 1)
# how many mass units are related to 2 mass units (e.g. mg) P
redfield("P", 2, "mass")

# mass percentage of elements
x <- redfield("P", 1, "mass")
x / sum(x)

# approximate fresh mass (in mg) formed by 1 microgram P
redfield("P", 1, "mass")["C"] * 2 * 10 / 1000

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc}

