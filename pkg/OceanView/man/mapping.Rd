\name{mapxy}
\alias{mapxy}
\alias{transectxy}
\alias{mapsigma}

\title{
   mapping and extracting from 2-D or 3D data
}
\description{
  \code{mapxy} maps a variable from a matrix (or array) from (x, y) coordinates 
  to a matrix or array with first and second dimension given by xto, and yto respectively.
  Extrapolation to all combinations of xto and yto. Uses 2-D linear interpolation. Returns a matrix or array.

  \code{transectxy} maps a variable from a matrix or array from (x, y) to the xy coordinate pair xyto
   by linear interpolation. Returns a vector or matrix.

  \code{mapsigma} maps a matrix \code{z} containing values defined at (x, sigma) coordinates 
  to (x, depth) coordinates.
  The depths corresponding to the sigma values in \code{z} are in input matrix sigma.
  The result is a matrix which will contain \code{NA}s where the depth-coordinates
  extend beyond the sigma values.

}
\usage{
  mapxy (z, x, y, xto = NULL, yto = NULL)
  transectxy (z, x, y, xyto)
  mapsigma (z, sigma, depth = NULL)
}

\arguments{
  \item{z }{Matrix or array with z-values to be mapped to lower or higher resolution (\code{mapxy}) 
    or to be mapped to the depth coordinates (\code{mapsigma}).
    }
  \item{x }{Vector with original x-coordinates of the matrix \code{z} to be mapped.
    Length should be = number of rows of \code{z}.  
    For \code{mapxy} \code{x} can also be a matrix, of same dimension as \code{z}.
    }
  \item{y }{Vector with original y-coordinates of the matrix \code{z} to be mapped. 
    Length should be = number of columns of \code{z}.  
    For \code{mapxy} y can also be a matrix, of same dimension as \code{z}.
    }
  \item{xto }{Vector with x-coordinates to which the matrix \code{z} should be mapped.
    The elements in \code{xto} should be embraced by the elements in \code{x}
    (it is not allowed to extrapolate outside of the region).
    If \code{NULL} then the range of \code{x} is covered, with the same number of points.
    }
  \item{yto }{Vector with y-coordinates to which the matrix \code{z} should be mapped. 
    The elements in \code{yto} should be embraced by the elements in \code{y}
    (it is not allowed to extrapolate outside of the region).
    If \code{NULL} then the range of \code{y} is covered, with the same number of points.
    }
  \item{xyto }{Two-columned matrix, with first and second column specifying the 
    x- respectively y-coordinates to which the matrix \code{z} should be mapped. 
    The elements should be embraced by the elements in \code{x} (first column)
    and \code{y} (second column)
    (it is not allowed to extrapolate outside of the region).
    }
  \item{sigma }{The sigma coordinates, a matrix with the same dimension as \code{z}, 
    i.e. (nrow(z), ncol(z)). 
    }
  \item{depth }{The depth (or z) coordinates to which matrix \code{z} has to be mapped. 
    If \code{NULL} then \code{seq(min(sigma), max(sigma), length.out = ncol(z))}. 
    }
}
\details{
  \code{mapxy} can be used to increase or decrease the resolution of a matrix,
  or to zoom in on a certain area. If \code{x} or \code{y} are a matrix, then the 
  function \link{interp} from the \code{akima} package is used.

  If \code{z} is an array, then the dimensions of the return element of 
  \code{mapxy} will be equal
  to c(length(\code{xto}), length(\code{yto}), dim(z)[3]). 
  Thus, if xto and yto are one value, the value returned will be a vector.  
  if xto = one value and yto a vector,then a matrix willl be returned, else it will 
  be an array.

  Function \code{transectxy} will return an object with one dimension less 
  than \code{mapxy}
    
  \code{mapsigma} should be used to make images from data that are in sigma 
  coordinates
}

\value{
  \code{mapxy}: 
    \itemize{
     \item{z} {The higher or lower resolution matrix or array, 
          with dimension = c(length(xto), length(yto), dim(z)[3])}
     \item{x} {The x coordinates, corresponding to first dimension of \code{z}}
     \item{y} {The y coordinates, corresponding to second dimension of \code{z}}
    }
  
  
  
  \code{transectxy}:
    \itemize{
     \item{z} {The higher or lower resolution object, 
   with dimension = c(nrow(xyto), dim(z)[3]).}
     \item{xy} {The pairs of (x,y) coordinates}
    }
  
  \code{mapsigma}: 
    \itemize{
     \item{z} {A matrix with columns in depth-coordinates.}
     \item{depth} {The depth-coordinates (also known as 'z'-coordinates.}
    }
}
\seealso{
  \link{Image} for a plotting function
}

\examples{

## =======================================================================
##  mapxy: changing the resolution and zooming of a matrix
## =======================================================================
 pm <- par(mfrow = c(2, 2))
 nr <- nrow(volcano)
 nc <- ncol(volcano)
 
 x  <- 1 : nr
 y  <- 1 : nc
 Image(x = x, y = y, volcano, main = "original")

# increasing the resolution
 x2 <- seq(from = 1, to = nr, by = 0.25)
 y2 <- seq(from = 1, to = nc, by = 0.25)

 VOLC1 <- mapxy(volcano, x = x, y = y, xto = x2, yto = y2)$z
 Image(x = x2, y = y2, z = VOLC1, main = "high resolution")

# low resolution
 xb <- seq(from = 1, to = nr, by = 2)
 yb <- seq(from = 1, to = nc, by = 3)
 VOLC2 <- mapxy(volcano, x, y, xb, yb)$z
 Image(VOLC2, main = "low resolution")

# zooming in high resolution
 xc <- seq(10, 40, 0.1)
 yc <- seq(10, 40, 0.1)

 VOLC3 <- mapxy(volcano,x, y, xc, yc)$z 
 Image(VOLC3, main = "zoom")

# Get one value or a grid of values
 mapxy(volcano, x, y, xto = 2.5,     yto = 5) 
 mapxy(volcano, x, y, xto = c(2, 5), yto = c(5, 10)) 

# Crosssection or transect
 transectxy(volcano, x, y, xyto = cbind(c(2, 5), c(5, 10)))
 
## =======================================================================
## take a cross section or transect of volcano
## =======================================================================

Image(volcano, x = 1:nr, y = 1:nc, mfrow = NULL)
xyto <- cbind(seq(from = 1,  to = nr, length.out = 20), 
              seq(from = 20, to = nc, length.out = 20))
points(xyto[,1], xyto[,2], pch = 16)              
              
(Crossection <- transectxy (volcano, x = 1:nr, y = 1:nc, 
            xyto = xyto))

plot(Crossection$z, type = "b")                         

## =======================================================================
##  mapxy: changing the resolution and zooming of an *array*
## =======================================================================
# volcano now as an array
u <- array(dim = c(dim(volcano), 2), 
           data = c(volcano, sqrt(volcano)))
nr <- nrow(volcano)
nc <- ncol(volcano)

mapxy (u, 1:nr, 1:nc, xto = 1.5, yto = 1.5:4.5)
Image(mapxy (u, x = 1:nr, y = 1:nc, 
             xto = seq(1, nr, by = 2), 
             yto = seq(1, nc, by = 2))$z,
             main = c("original", "sqrt"))

Crossection <- transectxy (u, x = 1:nr, y = 1:nc, xyto = xyto)
head(cbind(Crossection$xy, Crossection$z))

## =======================================================================
##  Mapping with x and y a matrix
## =======================================================================

# x- and y-coordinates of 'volcano'
volcx <- matrix(nrow = 87, ncol = 61, data = 1:87)
volcx <- volcx + matrix(nrow = 87, ncol = 61, 
         byrow = TRUE, data = seq(0., 15, length.out = 61))

volcy <- matrix(ncol = 87, nrow = 61, data = 1:61)
volcy <- t(volcy + matrix(ncol = 87, nrow = 61, 
         byrow = TRUE, data = seq(0., 25, length.out = 87)))


par(mfrow = c(1, 1))

# remap tilted volcano on square grid
tilted <- mapxy(volcano, volcx, volcy)
Image(tilted$z, x = tilted$x, y = tilted$y, main = "tilted volcano")

# add grid points
points(volcx, volcy, pch = ".", cex = 1.5)  

## =======================================================================
##  mapsigma: changing from sigma coordinates into depth-coordinates
## =======================================================================

 par(mfrow = c(2, 2))
 val <- t(matrix (nrow = 10, ncol = 50, data = 1:10))
 Image(val, main = "values in sigma coordinates")

# The depth at each 'column' 
 Depth <- approx(x = 1:5, y = c(10, 4, 5, 6, 4), 
                 xout = seq(1,5, len = 50))$y

# Sigma coordinates
 sigma <- t(matrix(nrow = 10, ncol = 50, data = Depth, byrow = TRUE) * 
                seq(from = 0, to = 1, length = 10))
 matplot(sigma, type = "l", main = "sigma coordinates", 
         ylab = "depth", ylim = c(10, 0))

# Mapping to the default depth coordinates
 valz <- mapsigma(z = val, sigma)
 Image(valz$z, y = valz$depth, NAcol = "black", ylim = c(10, 0), 
       main = "depth-coord, low resolution")

# Mapping to z-coordinates at higher resolution of depth coordinates
 valz <- mapsigma(val, sigma, depth = seq(0, 10, length.out = 100))
 Image(valz$z, y = valz$depth, NAcol = "black", ylim = c(10, 0), 
       main = "depth-coord, high resolution")

# reset plotting parameters
 par(mfrow = pm)


}
\keyword{ hplot }

