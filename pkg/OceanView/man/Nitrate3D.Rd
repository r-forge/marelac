\name{Nitrate3D}
\docType{data}
\alias{Nitrate3D}
\title{
  Yearly averaged Nitrate concentration from the NODC World Ocean
  Atlas 2005 (WOA05).
}

\description{
  Nitrate concentrations from the NODC World Ocean Atlas 2005. mmol m-3
  
  Originally made available by CSIRO.
  
  The surface values are gridded in 1dg * 1dg longitude latitude sets, 
  and there are 14 depth intervals. 
  The original data were averaged over the 4 seasons.

}

\usage{
data(Nitrate3D)
}

\format{
 list with 
  \itemize{
   \item lon, the longitude (dgE), at 1dg resolution, 360 values 
   \item lat, the latitude (dg N), at 1dg resolution, 180 values
   \item depth, the water depth, m.
   \item val, the value.   
    \code{val} is an array of dimension (360, 180, 4), (lon, lat, season).
   \item name, the long name of the variable
   \item units, the units of measurement
  }
}
\details{
  The "objectively analyzed climatology" has been used to extract these data.
}
\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}

\examples{
# load the data
  data(Nitrate3D)

# plot surface data
  Image(z = Nitrate3D$val[,,1], x = Nitrate3D$lon, y = Nitrate3D$lat,
        main = "surface nitrate data 2005")

# plot depth profiles
  lat <- Nitrate3D$lat
  Image.array(z = Nitrate3D$val, margin = c(2, 3), x = Nitrate3D$lat, 
            y = Nitrate3D$depth, subset = (lat > 19 & lat < 22),
            ylim = c(500,0), NAcol = "black")

}

\references{
  \url{http://www.nodc.noaa.gov}

  \url{ftp://ftp.nodc.noaa.gov/pub/data.nodc/woa/WOA05nc/seasonal/}
  
  Mark A. Collier and Paul J. Durack, 2005. CSIRO netCDF version of the
   NODC World Ocean Atlas 2005. CSIRO Marine and Atmospheric Research Paper 015.
   December 2006
}

\seealso{
  \link{Image} for plotting.
}

\keyword{datasets}

