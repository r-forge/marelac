\name{redfield}
\alias{redfield}

\title{Redfield Ratio Calculator}
\description{
  Estimate elemental composition of biomass (or media) according to the
  Redfield ratio.
}
\usage{
redfield(x, element, method = c("mol", "mass"))
}
\arguments{
  \item{x}{amount of substance of that element (in mol or mass units),}
  \item{element}{The element that is given ("C", "H", "O", "N", "P"),}
  \item{method}{measurement unit ("mol" or "mass").}
}
\details{
  The Redfield ratio is \eqn{C_{106}H_{80}O_{45}N_{16}P_{1}}{C106 H80
  O45 N16 P1}.
}
\value{
  A data frame with the estimated ratio of the main elements.
}

% Add **original** reference here, and not only N:P
\references{
  Redfield A.C. (1934). On the proportions of organic derivations in sea
  water and their relation to the composition of plankton.  In: James
  Johnstone Memorial Volume. (ed. R.J. Daniel).  University Press of
  Liverpool, 177-192.
}
\seealso{\code{\link{molweight}}, the molecular weight}
\examples{
## Redfield ratio
redfield(1, "P")
## returns the molar Redfield ratio, rescaled to nitrogen
redfield(1, "N")
## how many mass units are related to 2 mass units (e.g. mg) P
redfield(2, "P", "mass")
redfield(c(1, 2, 3), "N", "mass")

## mass percentage of elements
x <- redfield(1, "P", "mass")
x / sum(x)

## rule of thumb for fresh mass (in mg) formed by 1 microgram P
redfield(1, "P", "mass")$C * 2 * 10 / 1000
}
\keyword{ misc}

