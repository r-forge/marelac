\name{vectorplot}
\alias{vectorplot}
\title{
 Vector velocity plot.
}
\description{
  Displays vectors as segments.
}
\usage{
  vectorplot(x, y, times = NULL,  by = 1, arr = FALSE, tfac = NULL, ...) 
}

\arguments{
  \item{x }{A vector with quantities (velocities) in x-direction. 
    }
  \item{y }{A vector with quantities (velocities) in y-direction. 
    Should have the same length as \code{x}
    }
  \item{times }{A vector with x-values. If \code{0}, everything will 
    be radiating from the origin. 
    }
  \item{by }{Number increment for plotting vectors. 
    Set this to an integer > \code{1} if the vector density is too high. 
    }
  \item{tfac }{Only for \code{times} not \code{NULL}, the proportionality
    factor with which the vectors on the x-axis must be drawn. A value of \code{1}
    means that the distance of one will be drawn as one time unit on the x-axis.
    For a value of \code{2} a distance of 1 will appear as 2 time units on the x-axis. 
    if \code{NULL}, the range on the y-axis is used. In that case, it may be
    necessary to manually set the \code{xlim} of the figure.
    }
  \item{arr }{If \code{TRUE}, then \link{Arrows} will be drawn; 
    if \code{FALSE}, \link{segments} will be drawn. 
    }
  \item{\dots}{additional arguments passed to the plotting methods. 
    }
}

\value{
 none
}
\seealso{
 \link{quiver}, \link{flowpath} for other plotting functions.
}

\examples{
## =======================================================================
##  EXAMPLE: 
## =======================================================================
 u <- cos(seq(0, 2*pi, 0.1))
 v <- sin(seq(0, 2*pi, 0.1)+ 1)

 par(mfrow = c(2, 2))
 vectorplot(x = u, y = v)
 vectorplot(x = u, y = v, times = 1:length(u))
 vectorplot(x = u, y = v, times = 1:length(u), tfac = 1)
 vectorplot(x = u, y = v, times = 1:length(u), tfac = 10)
}
\keyword{ hplot }

