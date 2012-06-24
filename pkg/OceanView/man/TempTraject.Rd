\name{TrackProf}
\docType{data}
\alias{TrackProf}
\title{
  Temperature profiles made along ship track.
}

\description{
  Profiles of temperature made along ship track.
  
  Originally made available by US NOAA NODC.
  
  The data were merged from 29 input files named \code{gtspp_103799_xb_111.nc}
  till \code{gtspp_103827_xb_111.nc}.

  These data were acquired from the US NOAA National Oceanographic Data Center 
  (NODC) on 9/06/2012 from \url{http://www.nodc.noaa.gov/gtspp/}. 
}

\usage{
data(TrackProf)
}

\format{
 list with 
  \itemize{
   \item meta, a data.frame with the metadata, containing for each of the 
    29 profiles the following:
   \itemize{
    \item station, the number of the station (part of the filename)
    \item filename, the original name of the NetCDF file.
    \item date, the date of sampling
    \item the time of sampling, a number relative to 1-1-1900 0 hours
    \item longitude, dg E 
    \item latitutde, dg N
   }
   \item temp, the seawater temperature, at the \code{depth} of the 
     measurement in dg C. A matrix of dimension \code{(29, 93)} for the 
     29 profiles and (at most) 93 depth values; \code{NA} means no measurement.
   \item depth, the depth of the measurement in \code{temp}, in metres, 
     positive downward. A matrix of dimension \code{(29, 93)} for the 
     29 profiles and (at most) 93 depth values; \code{NA} means no measurement.
  }
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{

  pm <- par(mfrow = c(2, 2))
   
## =============================================================================
## show the metadata
## =============================================================================
  print(TrackProf$meta)

## =============================================================================
## display the cruisetrack on the Ocean Bathymetry data  
## =============================================================================  
  # 1. plots the ocean's bathymetry and add samplings
  ImageOcean(xlim = c(-50, 50), ylim = c(-50, 50), 
             main = "cruise track")  
  points(TrackProf$meta$longitude, TrackProf$meta$latitude, pch = "+")
  # mark starting point
  points(TrackProf$meta$longitude[1], TrackProf$meta$latitude[1], 
         pch = 18, cex = 2, col = "green")

## =============================================================================    
## image plots of raw data
## =============================================================================    

  Image(z = TrackProf$depth, main = "raw depth values",  
        xlab = "station nr", ylab = "sample nr", zlab = "depth")
  Image(z = TrackProf$temp, main = "raw temperature values",  
        xlab = "station nr", ylab = "sample nr", zlab = "dgC")

## =============================================================================    
## image plots of temperatures at correct depth
## =============================================================================    
  
  # water depths to which data set is interpolated
  depth <- 0 : 809    

  # map from "sigma" to "depth" coordinates
  Temp_Depth <- mapsigma (TrackProf$temp, TrackProf$depth, depth = depth)$z

  # image with depth increasing downward  
  Image(z = Temp_Depth, main = "Temperature-depth",
        ylim = c(809, 0), y = depth, NAcol ="black", resfac = 2,
        xlab = "station nr", ylab = "depth, m", zlab = "dgC")
  
  par(mfrow = pm)
   
## =============================================================================    
## scatterplot on Ocean Bathymetry
## =============================================================================    

  SurfTemp <- cbind(TrackProf$meta[ ,c("longitude", "latitude")], 
                    Temp = TrackProf$temp[,1])

  mr <- par("mar")
  par(mar = mr + c(0, 0, 0, 2))
  
  # No colors, but add contours
  ImageOcean(xlim = c(-30, 30), ylim = c(-40, 40), 
           main = "cruise track", col = NA, add.contour = TRUE) 
  
  scatterplot(SurfTemp, zvar = "Temp", add = TRUE)
  
  par(mar = pm)

}

\references{
  \url{http://www.nodc.noaa.gov/gtspp/}

  U.S. National Oceanographic Data Center

  Global Temperature and Salinity Profile Programme (GTSPP) Data.
  project = Joint IODE/JCOMM Global Temperature-Salinity Profile Programme

}

\seealso{
  \link{Image} for plotting.
}

\keyword{datasets}

