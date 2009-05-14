\name{convert_salinity}
\alias{convert_PStoAS}
\alias{convert_AStoPS}
\title{Practical - Absolute Salinity Conversions}
\description{Conversion from practical to absolute salinity and vice versa.
}
\usage{
convert_PStoAS(S = 35, p = max(0,P-1.013253), P=1.013253,
  lat=NULL, lon=NULL)

convert_AStoPS(S=35, p=max(0,P-1.013253), P=1.013253,
  lat=NULL, lon=NULL)

}
\arguments{
  \item{S }{Salinity, either practical salinity (\code{convert\_PStoAS}),
    dimensionless or absolute salinity (\code{convert\_AStoPS}, g/kg)
  }
  \item{p }{gauge or applied pressure, pressure relative to the
            local atmospheric pressure, bar
            }
  \item{P }{true pressure, bar
  }
  \item{lat }{latitude (-90 to +90)
  }
  \item{lon }{longitude (0 to 360)
  }

}

\value{
  The absolute salinity (\code{convert\_PStoAS}) or practical salinity
  (\code{convert\_AStoPS}), or a matrix with conversion factors for the
  respective gauge pressure (\code{convert\_Sfactor}).
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
convert_PStoAS(S=35)
convert_AStoPS(S=35)
convert_PStoAS(S=35, lat=10, lon=10, p=0)
}
\references{
  Millero, F.J., Feistel, R., Wright, D.G., and McDougall, T.J., 2008:
  The composition of Standard Seawater and the definition of the
  Reference-Composition Salinity Scale, Deep-Sea Res. I, 55, 50-72.
  
  Uses the Fortran code written by David Jackett.
  \url{http://www.marine.csiro.au/~jackett/TEOS-10/}
}
\keyword{utilities}

