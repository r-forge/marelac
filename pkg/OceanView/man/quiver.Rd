\name{quiver}
\alias{quiver}
\alias{quiver.matrix}
\alias{quiver.array}
\alias{flowpath}
\title{
  Plots velocitys as arrows or as trajectory plots.
}
\description{
  Functions \code{quiver} displays velocity vectors as arrows.
  
  Function \code{flowpath} displays the flow paths of particles, based on 
  velocity vectors.
  
}
\usage{
  quiver(u, ...)

  \method{quiver}{matrix}(u, v, x = NULL, y = NULL, scale = 1, 
         arr.max = 0.2, arr.min = 0, by = NULL, add = FALSE, ...) 

  \method{quiver}{array} (u, v, margin = c(1, 2), subset, ask = NULL, ...)
                    
  flowpath(u, v, x = NULL, y = NULL, startx = NULL, starty = NULL,
           scale = 1, plot = TRUE, add = FALSE, 
           numarr = 0, arr.length = 0.2, maxstep = 1000, ...)
}

\arguments{
  \item{u }{A matrix (\code{quiver}) or array (\code{quiver.array})
    with velocities in x-direction. 
    For \code{quiver} 
    the number of rows should be = Nx or Nx+1 (Nx = length(x), if x given).  
    Number of columns should be = Ny or Ny+1 (Ny = length(y), if y given). 
    }
  \item{v }{A matrix (\code{quiver}) or array (\code{quiver.array})
    with velocities in y-direction. 
    For \code{quiver} Number of rows should be = Nx or Nx+1;  
    Number of columns should be = Ny or Ny+1.  
    }
  \item{x }{Vector with x-coordinates of the velocities. If \code{NULL}, it
    is a sequence between (0, 1). Length = Nx. 
    }
  \item{y }{Vector with y-coordinates of the velocities. If \code{NULL}, it
    is a sequence between (0, 1). Length = Ny. 
    }
  \item{startx }{Vector with start position in x-direction of flow paths. 
    Length > =1. If not specified, then all combinations of x and y at 
    the outer margins will be used as starting point.
    }
  \item{starty }{Vector with start position in y-direction of flow paths. 
    Length = length of \code{startx}. 
    }
  \item{scale }{Scaling factor for the arrows. 
    }
  \item{arr.max }{Maximal size of arrowhead, in cm (approximately). 
    The arrows are scaled according to the velocity (sqrt(u^2 + v^2)). 
    \code{arr.max} is associated with the maximal velocity.
    }
  \item{arr.min }{Minimal size of arrowhead, in cm (approximately). 
    Set \code{arr.min} = \code{arr.max} for constant size.
    }
  \item{by }{Number increment for plotting vectors. 
    Useful if the vector density is too high. 
    }
  \item{margin }{a vector giving the subscripts which the plotting 
    function will be applied over. For instance, \code{c(1, 2)},
     indicates rows(x) and columns(y);  \code{c(2, 1)} the same but transposed.
    \code{margin} should be a vector with two numbers inbetween \code{1}, and \code{3}.
    The plotting function will loop over the index that is not in \code{margin}.
    }
  \item{ask }{logical; if \code{TRUE}, the user is asked before each plot, 
     if \code{NULL} the user is only asked if more than one page of plots is 
     necessary and the current graphics device is set interactive, 
     see \link{par}(ask) and \link{dev.interactive}.
    }
  \item{add }{If \code{TRUE}, will add to current plot. 
    Else will start a new plot. 
    }
  \item{plot }{If \code{FALSE}, will not plot the flow paths, 
    but will return the matrix with path values instead.
    }
  \item{numarr }{The number of arrows added on the flow paths. 
    }
  \item{arr.length }{Constant size of arrowhead, in cm (approximately). 
    }
  \item{maxstep }{Maximum number of steps for calculating the flow paths.  
    } 
  \item{\dots}{additional arguments passed to the methods. 
    }
  \item{subset }{logical expression indicating over which elements to loop
     missing values are taken as \code{FALSE}
    }
}
\details{
  \code{quiver} plots vectors specified by \code{u,v} at the coordinates
  \code{x,y}.
  
  \code{flow path} uses the velocities \code{u, v} at the coordinates
  \code{x,y} to create trajectories, starting at points 
  \code{startx, starty}. It can also be used to return the flow path 
  points by setting \code{plot} equal to \code{FALSE}
}

\value{
  \code{flow path} also returns (as \code{invisible}) a 2-column 
  matrix with the x-y coordinates of the flow paths. Separate flow paths
  are separated with \code{NA}.
}
\seealso{
  \link{Image} for an image function
  
  \link{vectorplot} for plotting velocity vectors as spikes

  \link{Arrows} for the arrow function from package \code{shape} on which quiver is based.
}

\examples{
## =======================================================================
##  EXAMPLE 1: 
## =======================================================================
 pm <- par("mfrow")
 par(mfrow = c(2, 2))

# generate arrows
 x <- seq(-1, 1, by = 0.2)
 y <- seq(-1, 1, by = 0.2)
 dx <- outer(x, y , function(x, y) -y)
 dy <- outer(x, y , function(x, y) x)

# velocity plots
 quiver(u = dx, v = dy, x = x, y = y)

# different color for left/right pointing arrows
 Col <- dx
 Col[dx > 0] <- "blue"
 Col[dx < 0] <- "red"
 quiver(u = dx, v = dy, x = x, y = y, 
        arr.max = 0.4, arr.min = 0.1, col = Col)

# different scale
 quiver(u = dx, v = dy, x = x, y = y, by = 2, scale = 2)

# three flow paths
 flowpath(u = dx, v = dy, x = x, y = y, startx = 0.1, starty = 0.1)
 flowpath(u = dx, v = dy, x = x, y = y, 
           startx = c(0.9, -0.9), starty = c(0.0, 0.0), col = "red",
           numarr = 2, add = TRUE)

## =======================================================================
##  EXAMPLE 2: 
## =======================================================================
 par(mfrow = c(1, 1))
 x <- seq(-2, 2, by = 0.2)
 y <- seq(-1, 1, by = 0.2)
 z <- outer (x, y, function(x, y) x^3 - 3*x -2*y^2)
 contour(x, y, z, col = femmecol(10))

# gradients in x- and y-direction (analytical)
 dX <- outer(x, y, function(x,y) 3*x^2 - 3)
 dY <- outer(x, y, function(x,y) -4*y)

 quiver(u = dX, v = dY, x = x, y = y, scale = 1, add = TRUE, by = 1)
 flowpath(u = dX, v = dY, x = x, y = y, startx = c(-2, 1.1), 
            starty = c(-1, -1), add = TRUE, arr.len = 0.5,
            col = "darkgreen", lwd = 3, numarr = 1)

## =======================================================================
##  EXAMPLE 3: 
## =======================================================================
 
 x <- y <- 1:20
 u <- outer (x, y, function (x, y) cos(2*pi*y/10))
 v <- outer (x, y, function (x, y) cos(2*pi*x/10))

 quiver(x = x, y = y, u = u, v = v, col = "grey")

# flowpaths using all combinations of x and y at edges
 flowpath(x = x, y = y, u = u, v = v, add = TRUE, 
           lwd = 2, col = "orange")

## =======================================================================
##  EXAMPLE 4: quiver of an array.. 
## =======================================================================
 
 x <- y <- 1:20
 u2 <- outer (x, y, function (x, y) sin(2*pi*y/10))
 v2 <- outer (x, y, function (x, y) sin(2*pi*x/10))

# merge u, u2 and v, v2 to create an "array"
 U <- array(dim = c(dim(u2), 2), data = c(u, u2), 
            dimnames = c("x", "y", "times"))
 V <- array(dim = c(dim(v2), 2), data = c(v, v2), 
            dimnames = c("x", "y", "times"))

 quiver(u = U, v = V, x = x, y = y, main = c("time 1", "time 2"))

# quiver over x and time, for a subset of y-values:
 quiver(u = U, v = V, x = x, y = 1:2, 
        margin = c(1, 3), main = paste("y ", y), 
        subset = y < 10)

\dontrun{
 quiver(u = U, v = V, x = x, y = y, ask = TRUE, 
        mfrow = c(1,1))

 quiver(u = U, v = V, x = x, y = 1:2, ask = TRUE, 
        margin = c(1, 3), main = paste("y ", y),
        mfrow = c(1,1))

}

## =======================================================================
##  EXAMPLE 5: 
## =======================================================================
 par(mfrow = c(1, 1))

 Image(x = 1:nrow(volcano), y = 1:ncol(volcano), 
       volcano, add.contour = TRUE)

# Assume these are streamfunctions, we calculate the velocity field as:
 dx <- dy <- 1
 v <-   (volcano[-1, ] - volcano[-nrow(volcano), ] )/dx
 u <- - (volcano[, -1] - volcano[ ,-ncol(volcano)] )/dy

 quiver(x = 1:nrow(u), y = 1:ncol(v), 
        u = u, v = v, add = TRUE, by = 3)
 
 flowpath(x = 1:nrow(u), y = 1:ncol(v), numarr = 10,
          u = u, v = v, add = TRUE, lwd = 2, col = "grey", 
          startx = 20, starty = 30)

 par(mfrow = pm)
}
\keyword{ hplot }

