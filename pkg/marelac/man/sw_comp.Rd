\name{sw_comp}
\alias{sw_comp}
\title{Reference sea salt Composition}
\description{
  The sea salt composition definition for reference salinity of the standard
  ocean at 25 dgC and 101325 Pa (atmospheric pressure), given in mass fractions).
}
\usage{
sw_comp(x=c("Na","Mg","Ca","K","Sr","Cl","SO4","HCO3",
            "Br","CO3","BOH4","F","OH","BOH3","CO2"))
}
\arguments{
  \item{x}{character vector with components whose composition should be
    estimated.
  }
}
\value{
  A vector with the mass fractions.
}
\references{
Millero, F.J., Waters, J., Woosley, R., Huang, F. and Chanson, M., 2008: The effect of
composition of the density of Indian Ocean waters, Deep-Sea Res. I, 55, 960-470. }
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}

\examples{
sw_comp("CO2")
sw_comp()
sum(sw_comp())
}
\keyword{utilities}
