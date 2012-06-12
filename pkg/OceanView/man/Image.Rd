\name{Image}
\alias{Image}
\alias{Image.matrix}
\alias{Image.array}
\alias{Image.list}
\title{
 Improved image plot for 2-D data and for 3-D data
}
\description{
  Image is similar to R's \link{image} function, except that it is possible to 
  add a legend and contours, it has a different default coloring scheme, 
  it can deal with decreasing x- and y-values, 
  and it is possible to increase or decrease the resolution by 2-D 
  linear interpolation. 
  
}
\usage{
  Image(z, ...)

  \method{Image}{matrix}(z, col = femmecol(100), add.contour = FALSE, 
                   legend = TRUE, resfac = 1, NAcol = NULL, 
                   zlab = NULL, cex.zlab = par("cex.axis"),...) 
  \method{Image}{array}(z, margin = c(1, 2), subset, ask = NULL, ...) 
  \method{Image}{list}(z, ...) 
}

\arguments{
  \item{z }{Matrix (2-D) or array (3-D) or a list with matrices or arrays,
    to be plotted as an \link{image}. 
    }
  \item{col }{A list of colors to be used for the \link{image} function
    }
  \item{NAcol }{Color to be used for \code{NA} values; default is ``white''. 
    }
  \item{add.contour }{If \code{TRUE} then contours will be added to the image plot.
    }
  \item{legend }{If \code{TRUE} then a color legend will be added in the right margin.
    } 
  \item{zlab }{Only if \code{legend = TRUE}, the label to be written on top of the
    colorlegend. The label will be written at the same level as the main title.
    to lower it, \code{zlab} can be made a vector, with the first values empty 
    strings - see last example. 
    To use different zlabs for different plots, use a list 
    (for \code{Image.array} and \code{Image.list}).
    } 
  \item{cex.zlab }{Only if \code{legend = TRUE}, the size of the label 
    written on top of the colorlegend; default = same as other axes. 
    } 
  \item{resfac }{Resolution factor, one value or a vector of two numbers, for 
     the x and y- values respectively. If >= 1, which will increase the 
    resolution. For instance, if \code{resfac} equals \code{(3} then for each 
    adjacent pair of x- and y-values, z will be interpolated to two intermediary points. 
    This uses simple linear interpolation. If \code{resfac} is one number then
    the resolution will be increased similarly in x and y-direction.  
    }
  \item{margin }{a vector giving the subscripts which the \code{image} 
    function will be applied over. For instance, \code{c(1, 2)},
     indicates rows (x) and columns (y);  \code{c(2, 1)} the same but transposed.
    \code{margin} should be a vector with two numbers inbetween \code{1}, and \code{3}.
    The image function will loop over the index that is not in \code{margin}.
    }
  \item{ask }{logical; if \code{TRUE}, the user is asked before each plot, 
     if \code{NULL} the user is only asked if more than one page of plots is 
     necessary and the current graphics device is set interactive, 
     see \link{par}(ask) and \link{dev.interactive}.
    }
  \item{subset }{logical expression indicating over which elements to loop
     missing values are taken as \code{FALSE}
    }
  \item{\dots}{additional arguments passed to the methods. 
    }
}
\details{
  \code{Image} is an extension to the default \link{image} plot that has
  the possibility to add a legend and contourlines, and to increase the 
  resolution, in order to make smoother images. It also uses a different color
  scheme and can deal with decreasing x- and
  y-values. \code{Image.array} and \code{Image.list} are versions that
  accept a 3 dimensional array resp. a list as their first argument to
  produce multiple plots.
}

\value{
  Upon returning from the \code{image} function, the figure coordinates are set
  to the main figure (excluding the legend). Thus, one can safely add other 
  plotting elements.
}
\seealso{
  \link{mapxy} for changing the resolution of an image.

  \link{image} for the original R function.
  
  \link{quiver} for a quiver function.
}

\examples{

## =======================================================================
##  EXAMPLE 1: 
## =======================================================================
  Image(volcano, x = 1:nrow(volcano), y = 1:ncol(volcano), 
        main = "volcano", zlab = "height, m")
  points(50, 30, pch = 3, cex = 5, lwd = 3)
  abline(v = seq(10, 80, by = 10))
  abline(h = seq(10, 60, by = 10))

## =======================================================================
##  EXAMPLE 2: large data set, input is an array
## =======================================================================
  data(Oxsat3D)
  Image(z = Oxsat3D$val[,,1], x = Oxsat3D$lon, y = Oxsat3D$lat,
        main = "surface oxygen saturation data 2005", NAcol = "black", 
        zlab = c("","","\%"))

## -----------------------------------------------------------------------
  # images at first 9 depths
  
  pm <- par("mfrow")
  Image(z = Oxsat3D$val[ , ,1:9], 
        x = Oxsat3D$lon, y = Oxsat3D$lat,
        margin = c(1, 2), NAcol = "black", 
        xlab = "longitude", ylab = "latitude", 
        zlim = c(0, 115),
        main = paste("depth ", Oxsat3D$depth[1:9], " m"),
        mfrow = c(3, 3))
  par(mfrow = pm)

## -----------------------------------------------------------------------
  # images at latitude - depth section

  pm <- par("mfrow")
  Image(z = Oxsat3D$val[,  Oxsat3D$lat > - 5 & Oxsat3D$lat < 5, ], 
        x = Oxsat3D$lon, y = Oxsat3D$depth,
        margin = c(1, 3), NAcol = "black", 
        resfac = 3, ylim = c(5000, 0))
  # show position of transects 
  Image(z = Oxsat3D$val[ , ,1], 
        x = Oxsat3D$lon, y = Oxsat3D$lat,
        NAcol = "black")
  abline(h = Oxsat3D$lat[Oxsat3D$lat > - 5 & Oxsat3D$lat < 5])      
  par(mfrow = pm)

## =======================================================================
##  EXAMPLE 3: list of matrices
## =======================================================================

  volcano2 <- log(volcano)
  Image(list(volcano, volcano2), 
        x = 1:nrow(volcano), y = 1:ncol(volcano), 
        main = c("volcano", "log(volcano)"), 
        zlab = list("height, m", "log(m)"),
        zlim = list(c(80, 200), c(4.4, 5.5)))

\dontrun{
## =======================================================================
##  EXAMPLE 4: list of arrays
## =======================================================================
# crude conversion from oxsat to oxygen 
  Oxygen <- Oxsat3D$val/100 * 360
  
  pm <- par("mfrow")
  Image(z = list(Oxsat3D$val, Oxygen), 
        x = Oxsat3D$lon, y = Oxsat3D$lat,
        margin = c(1, 2), NAcol = "black", 
        main = c("Oxygen saturation ", " Oxygen concentration"),
        mtext = paste("depth ", Oxsat3D$depth, " m")
        )
  par(mfrow = pm)

}
}
\keyword{ hplot }

