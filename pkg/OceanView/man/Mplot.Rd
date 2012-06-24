\name{Mplot}
\alias{Mplot}
\alias{Msplit}
\alias{Mcommon}
\alias{scatterplot}
\title{
 Functions for plotting matrices
}
\description{
  \code{Mplot} plot data in (a list of) matrices. 
  
  \code{scatterplot} plots a matrix with z-colors.
  
  \code{Msplit} splits a matrix in a list according to factors.
  
  \code{Mcommon} creates a list of matrices that have only common variables.
}
\usage{
  Mplot (M, ..., x = 1, select = NULL, which = select, 
         subset = NULL, ask = NULL, legendpar = list(x = "top")) 

  scatterplot (M, x = 1, select = NULL, which = select, 
               zvar = x, subset = NULL, col = NULL, 
               ask = NULL, discrete = FALSE, 
               legendpar = list(x = "top"), add = FALSE, ... ) 

  Msplit (M, column = 1, subset = NULL) 

  Mcommon (M, ..., verbose = FALSE)      

}

\arguments{
  \item{M }{Matrix or data.frame to be plotted. For \code{Mplot}, \code{M} can be
    a list with matrices or data.frames. 
    }
  \item{x }{Name of number of the column to be used as the x-values. 
    }
  \item{select }{which variable/columns to be selected.  This is added for
    consistency with the R-function \code{subset}.
  }
  \item{which }{the name(s) or the index to the variables that should be
    plotted or selected.  Default = all variables, except \code{time}. }

  \item{subset }{logical expression indicating elements or rows to keep in 
    \code{select}: missing values are taken as \code{FALSE}
  }
  \item{ask }{logical; if \code{TRUE}, the user is \emph{ask}ed before
    each plot, if \code{NULL} the user is only asked if more than one
    page of plots is necessary and the current graphics device is set
    interactive, see \code{\link{par}(ask)} and
    \code{\link{dev.interactive}}.}
  \item{col }{A list of colors to be used for coloring the z variable 
    in \code{scatterplot}.    
    }
  \item{legendpar }{A \code{list} with parameters for the color legend 
    to be added; ignored when \code{discrete} is \code{TRUE}.
    } 
    
  \item{zvar }{Name of the column that holds the variable that defines the
    color in \code{scatterplot}.
    } 
  \item{discrete }{If \code{TRUE}, the z-variable is represented by discrete values;
    if \code{FALSE} it is continuous. 
    } 
  \item{add }{If \code{TRUE}, the dots will be added to an existing plot. 
    } 
  \item{column }{The name or number of the column with the factor 
    according to which the matrix will be split. 
    } 
  \item{verbose }{If \code{TRUE} will write output to the screen. 
    } 
  \item{\dots}{additional arguments passed to the methods. For \code{Mplot}:
    can also be extra matrices to plot. 
    }
}

\value{
  Function \code{Msplit} returns a list with the matrices, split according to 
  the factors; the names of the elements is set by the factor's name.

  Function \code{Mcommon} returns a list with the matrices, which only have
  the common variables.
}
\seealso{
  \link{Image} for plotting an image.

  \link{TrackProf} for a data set to be plotted.
  
}

\examples{

## =============================================================================
## Create three dummy matrices
## =============================================================================

M1 <- matrix(nrow = 10, ncol = 5, data = 1:50)
colnames(M1) <- LETTERS[1:5]

M2 <- M1[, c(1, 3, 4, 5, 2)] 
M2[ ,-1] <- M2[,-1] /2
colnames(M2)[3] <- "CC"  # Different name

M3 <- matrix(nrow = 5, ncol = 4, data = runif(20)*10) 
M3[,1] <- sort(M3[,1])
colnames(M3) <- colnames(M1)[-3]

head(M1); head(M2); head(M3)

# plot results
Mplot(M1, M2, M3, lwd = 2, mtext = "All variables",
      legendpar = list(x = "top", legend = c("M1", "M2", "M3")))

# Only common elements
Mplot(Mcommon(M1, M2, M3), lwd = 2, mtext = "common variables", 
      legendpar = list(x = "top", legend = c("M1", "M2", "M3")))

Mplot(M1, M2, M3, x = "B", select = c("A", "E"), pch = c(NA, 16, 1), 
      type = c("l", "p", "b"), col = "black",
      legendpar = list(x = "right", legend = c("M1", "M2", "M3")))

## =============================================================================
## The iris and Orange data set
## =============================================================================

Irislist <- Msplit(iris, column = "Species")
names(Irislist)

Mplot(Irislist, type = "p", 
      legendpar = list(x = "bottomright", title = "species"))

scatterplot(iris, x = 1, zvar = 5, discrete = TRUE, 
      legendpar = list(x = "right"))
scatterplot(iris, x = 2, which = c(1, 3), zvar = 4)

Mplot(Msplit(Orange,1), lwd = 2,
      legendpar = list(x = "topleft", title = "tree nr"))

## =============================================================================
## Scatterplot on an image
## =============================================================================

# Image of oxygen saturation
zlim <- range(Oxsat3D$val[,,1], na.rm  = TRUE) 
Image(Oxsat3D$val[,,1], x = Oxsat3D$lon, y = Oxsat3D$lat,
      xlab = "longitude", ylab = "latitude", 
      main = "Oxygen saturation", zlim = zlim, zlab = "\%")

# (imaginary) measurements at 5 sites
O2sat <- cbind(lon = c( 11.2,   6.0, 0.9,  -4, -8.8), 
               lat = c(-19.7,-14.45,-9.1,-3.8, -1.5),
               sat = c(   90,    95,  92,  85,  100))

# add to image; use same zrange; avoid adding legend
scatterplot(O2sat, zvar = "sat", zlim = zlim, 
            add = TRUE, cex = 2, legendpar = list())
}
\keyword{ hplot }

