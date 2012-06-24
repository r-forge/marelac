\name{Oxsat3D}
\docType{data}
\alias{Oxsat3D}
\title{
  Yearly averaged oxygen saturation from the NODC World Ocean
  Atlas 2005 (WOA05).
}

\description{
  Percentage Oxygen Saturation from the NODC World Ocean Atlas 2005. 
  
  Originally made available by CSIRO.
  
  The surface values are gridded in 2dg * 2dg longitude - latitude sets, 
  and there are 33 depth intervals. 
  The original data were averaged over the 4 seasons, and converted to half 
  the resolution for latitude and longitude.
  The longitude was converted to the European view 
  (original from (0,360) converted to (-180, 180))
  
}

\usage{
data(Oxsat3D)
}

\format{
 list with 
  \itemize{
   \item lon, the longitude (dgE), at 2dg resolution, 180 values 
   \item lat, the latitude (dg N), at 2dg resolution, 90 values
   \item depth, the water depth, m.
   \item val, the value.   
    \code{val} is an array of dimension (180, 90, 33), (lon, lat, depth).
   \item name, the long name of the variable
   \item units, the units of measurement
  }
}
\details{
  The "objectively analyzed climatology" has been used to extract these data.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{
# load the data
  data(Oxsat3D)

# plot surface data
  Image(z = Oxsat3D$val[,,1], x = Oxsat3D$lon, y = Oxsat3D$lat,
        main = "surface oxygen saturation \% 2005")

# plot depth profiles; input is an array
  lon <- Oxsat3D$lon
  Image (z = Oxsat3D$val, margin = c(2, 3), x = Oxsat3D$lat, 
         y = Oxsat3D$depth, subset = (lon > 18 & lon < 23),
         ylim = c(5500, 0), NAcol = "black", zlim = c(0, 110))

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

