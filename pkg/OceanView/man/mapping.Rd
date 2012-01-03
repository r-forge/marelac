\name{mapxy}
\alias{mapxy}
\alias{mapsigmatoz}
\title{
 \code{mapxy} maps 2-D data to higher or lower resolution
 
 \code{mapsigmatoz} maps 2-D data in sigma coordinates to z-coordinates
}
\description{
  \code{mapxy} maps a matrix 'z' from (x, y) to (xto, yto) by 2-D 
  linear interpolation.

  \code{mapsigtoz} maps a matrix containing (x, depth) values from 
  sigma-depth values to z-depth coordinates.
  The result is a matrix which will contain \code{NA}s where the z-coordinates
  extend beyond the sigma values.
}
\usage{
  mapxy (M, x, y, xto, yto)
  mapsigmatoz (M, sigma, zcoord = NULL)
}

\arguments{
  \item{M }{Matrix or array to be mapped to lower or higher resolution (\code{mapxy}) 
    or to be mapped to the z-depth coordinates (\code{mapsigmatoz}).
    }
  \item{x }{Vector with original x-coordinates of the matrix \code{M} to be mapped.
    Length should be = number of rows of \code{z}.  
    }
  \item{y }{Vector with original y-coordinates of the matrix \code{M} to be mapped. 
    Length should be = number of columns of \code{z}.  
    }
  \item{xto }{Vector with x-coordinates to which the matrix \code{M} should be mapped. 
    The elements in \code{xto} should be embraced by the elements in \code{x}
    (it is not allowed to extrapolate outside of the region).
    }
  \item{yto }{Vector with y-coordinates to which the matrix \code{M} should be mapped. 
    The elements in \code{yto} should be embraced by the elements in \code{y}
    (it is not allowed to extrapolate outside of the region).
    }
  \item{sigma }{The sigma coordinates, a matrix with the same dimension as \code{M}, 
    a matrix with dimension (nrow(M), ncol(M)). 
    }
  \item{zcoord }{The z coordinates to which matrix \code{M} has to be mapped. 
    If \code{NULL} then \code{seq(min(sigma), max(sigma), length.out = ncol(M))}. 
    }
}
\details{
  \code{mapxy} can be used to increase or decrease the resolution of a matrix,
  or to zoom in on a certain area.
  If \code{M} is an array, then \code{xto} and \code{yto} should be one value.
  In this case, the value returned will be a vector.  
  
  \code{mapsigmatoz} should be used to make images from data that are in sigma 
  coordinates
}

\value{
  \code{mapxy}: The higher or lower resolution matrix, 
   with dimension = (length(x.to), length(y.to)).
  
  \code{mapsigmatoz}: 
    \itemize{
     \item{M} {A matrix with columns in z-coordinates.}
     \item{zcoord} {The z-coordinates.}
    }
}
\seealso{
  \link{Image} for a plotting function
}

\examples{

## =======================================================================
##  mapxy: changing the resolution and zooming
## =======================================================================
 pm <- par(mfrow = c(2,2))
 x  <- 1 : nrow(volcano)
 y  <- 1 : ncol(volcano)
 Image(x = x, y = y, volcano, main = "original")

# increasing the resolution
 x2 <- seq(1, 87, 0.25)
 y2 <- seq(1, 61, 0.25)

 VOLC1 <- mapxy(volcano, x = x, y = y, xto = x2, yto = y2)
 Image(x = x2, y = y2, z = VOLC1, main = "high resolution")

# low resolution
 xb <- seq(1, 87, 2)
 yb <- seq(1, 61, 3)
 VOLC2 <- mapxy(volcano, x, y, xb, yb)
 Image(VOLC2, main = "low resolution")

# zooming in high resolution
 xc <- seq(10, 40, 0.1)
 yc <- seq(10, 40, 0.1)

 VOLC3 <- mapxy(volcano,x, y, xc, yc) 
 Image(VOLC3, main = "zoom")

# Get one value or a grid of values
 mapxy(volcano, x, y, xto = 2.5,     yto = 5) 
 mapxy(volcano, x, y, xto = c(2, 5), yto = c(5, 10)) 

## =======================================================================
##  mapsigmatoz: changing from sigma coordinates into z-coordinates
## =======================================================================
 val <- t(matrix (nrow = 10, ncol = 50, data = 1:10))
 Image(val, main = "values in sigma coordinates")

# The depth at each 'column' 
 Depth <- approx(x = 1:5, y = c(10, 4, 5, 6, 4), 
                 xout = seq(1,5, len = 50))$y

# Sigma coordinates
 sigma <- t(matrix(nrow = 10, ncol = 50, data = Depth, byrow = TRUE) * 
                seq(from = 0, to = 1, length = 10))
 matplot(sigma, type = "l", main = "sigma coordinates", 
         ylim = c(10, 0))

# Mapping to the default z coordinates
 valz <- mapsigmatoz(M = val, sigma)
 Image(valz$M, y = valz$zcoord, NAcol = "black", ylim = c(10, 0), 
       main = "z-coord, low resolution")

# Mapping to z-coordinates at higher resolution of zcoord
 valz <- mapsigmatoz(val, sigma, zcoord = seq(0, 10, length = 100))
 Image(valz$M, y = valz$zcoord, NAcol = "black", ylim = c(10, 0), 
       main = "z-coord, high resolution")

# reset plotting parameters
 par(mfrow = pm)
}
\keyword{ hplot }

