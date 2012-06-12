## =============================================================================
## =============================================================================
## Maps a matrix 'M' from (x,sigma) to (x, z) 
## =============================================================================
## =============================================================================

Approx <- function(x, y, ...) {
  if(diff(range(x, na.rm = TRUE)) == 0)
    warning("Warning in approx: all 'x' values are the same")
  if (any(is.na(c(x, y)))) {
    ii <- unique(c(which(is.na(x)), which(is.na(y))))
    x <- x[-ii]
    y <- y[-ii]
  } 
  approx(x, y, ...)
}

mapsigmatoz <- function (M, 
                         sigma, 
                         zcoord = NULL) {
  Nr <- nrow(M)
  Nc <- ncol(M)
  if (is.null(zcoord))
    zcoord <- seq(min(sigma), max(sigma), length.out = ncol(sigma))
  if(nrow(sigma) != Nr | ncol(sigma) != Nc)
    stop ("'sigma' should be of same dimension as 'M'")
 
  Mnew <- matrix(nrow = Nr, ncol = length(zcoord), data =0)
  
  for (i in 1:Nr) 
    Mnew[i,] <- Approx(x = sigma[i,], y = M[i,], xout = zcoord)$y 
     
  list(M = Mnew, zcoord = zcoord)
}

## Accounts for occurrence of decreasing values...
  FindInterval <- function(x, vec, ...) {

      if (all(diff(vec) < 0)) {# swap
        vec <- rev(vec)
        c(length(vec):1) [findInterval(x, vec, ...)]
      } else findInterval(x, vec, ...)
    
  }


## =============================================================================
## =============================================================================
## Maps a matrix or array 'M' from (x, y, ..) to (xto, yto, ..) by 
## linear 2-D interpolation
## =============================================================================
## =============================================================================

mapxy <- function(M, x, y, xto, yto) {

  Nx <- length(x)
  Ny <- length(y)

  if (min(xto) < min(x) | max(xto) > max(x)) 
    stop("'x' should embrace 'xto'")
  if (min(yto) < min(y) | max(yto) > max(y)) 
    stop("'y' should embrace 'yto'")

  dx  <- c(diff(x), 1)  # 1= for last value
  dy  <- c(diff(y), 1)

  Du <- dim(M)
  if (length(Du) > 3)
    stop ("'M' should be either a matrix or an array of dim 3")
  if (Du[1] != Nx)
    stop("'M' and 'x' not compatible: 1st dimension not equal to length (x)")
  if (Du[2] != Ny)
    stop("'M' and 'y' not compatible: 2nd dimension not equal to length (y)")


  Transf <- function (xto, yto, u, type = 1) {

  # find embracing values : first interval
    ix <- FindInterval(xto, x )
    iy <- FindInterval(yto, y)

  # next interval
    ixp1 <- pmin(ix+1, Nx)
    iyp1 <- pmin(iy+1, Ny)

  # interpolation factor
    xfac <- (xto-x[ix])/dx[ix]
    yfac <- (yto-y[iy])/dy[iy]

  # interpolate
    if (type == 1)
     (1-yfac)*((1-xfac)*u[cbind(ix,iy)  ]+xfac*u[cbind(ixp1,iy)]) +
         yfac*((1-xfac)*u[cbind(ix,iyp1)]+xfac*u[cbind(ixp1,iyp1)])
    else if (type == 2)
     (1-yfac)*((1-xfac)*u[ix,iy,  ]+xfac*u[ixp1,iy,]) +
         yfac*((1-xfac)*u[ix,iyp1,]+xfac*u[ixp1,iyp1,])
  } # end Transf
  if(length(Du) == 2)
    outer(xto, yto, FUN = Transf, u = M, type = 1)

  else if (length(Du) == 3 & length(xto) == 1 & length(yto) == 1) 
    Transf(xto, yto, M, type = 2) 

  else if (length(Du) == 3)  {
     MM <- NULL
     i3 <- dim(M)[3]
     for (i in 1:i3) { 
       MM <- c(MM, outer(xto, yto, FUN = Transf, u = M[ , ,i], type = 1))
     }
     MM <- array(dim = c(length(xto), length(yto), i3), data = MM)
     MM             
    }

  else stop('cannot run mapxy; xto and yto not of correct type')
}

## =============================================================================
## =============================================================================
## Takes a transect across a matrix or array 'M' from (x, y, ..) to 
## cbind(xto, yto) by linear 2-D interpolation
## =============================================================================
## =============================================================================

transectxy <- function(M, x, y, xyto) {

  Nx <- length(x)
  Ny <- length(y)

  if (ncol(xyto) != 2)
    stop("'xyto' should be a two-columned matrix")
  
  xto <- xyto[ ,1]
  yto <- xyto[ ,2]
    
  if (min(xto) < min(x) | max(xto) > max(x)) 
    stop("'x' should embrace elements in first column of 'xyto'")
  if (min(yto) < min(y) | max(yto) > max(y)) 
    stop("'y' should embrace elements in second column of 'xyto'")

  dx  <- c(diff(x), 1)  # 1= for last value
  dy  <- c(diff(y), 1)

  Du <- dim(M)
  if (length(Du) > 3)
    stop ("'M' should be either a matrix or an array of dim 3")
  if (Du[1] != Nx)
    stop("'M' and 'x' not compatible: 1st dimension not equal to length (x)")
  if (Du[2] != Ny)
    stop("'M' and 'y' not compatible: 2nd dimension not equal to length (y)")

 # find embracing values : first interval
    ix <- FindInterval(xto, x )
    iy <- FindInterval(yto, y)

  # next interval
    ixp1 <- pmin(ix+1, Nx)
    iyp1 <- pmin(iy+1, Ny)

  # interpolation factor
    xfac <- (xto-x[ix])/dx[ix]
    yfac <- (yto-y[iy])/dy[iy]

  # interpolate
    if (is.matrix(M))
     (1-yfac)*((1-xfac)*M[cbind(ix,iy)  ]+xfac*M[cbind(ixp1,iy)]) +
         yfac*((1-xfac)*M[cbind(ix,iyp1)]+xfac*M[cbind(ixp1,iyp1)])
    else {
     MM <- NULL
     i3 <- dim(M)[3]
     for (i in 1:i3) { 
       MM <- c(MM, 
         (1-yfac)*((1-xfac)*M[cbind(ix,iy, i)  ]+xfac*M[cbind(ixp1,iy, i)]) +
         yfac    *((1-xfac)*M[cbind(ix,iyp1, i)]+xfac*M[cbind(ixp1,iyp1, i)])
       
       )
     }
     MM <- matrix(nrow = nrow(xyto), ncol = i3, data = MM)
     MM             
    }
}

