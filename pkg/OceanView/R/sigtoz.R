## =============================================================================
## =============================================================================
## Maps a matrix 'M' from (x,sigma) to (x, z) 
## =============================================================================
## =============================================================================

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
    Mnew[i,] <- approx(x = sigma[i,], y = M[i,], xout = zcoord)$y 
     
  list(M = Mnew, zcoord = zcoord)
}


## =============================================================================
## =============================================================================
## Maps a matrix or array 'M' from (x,y) to (xto, yto) by 
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
    stop("'M' and 'x' not compatible")
  if (Du[2] != Ny)
    stop("'M' and 'y' not compatible")

  if (length(Du) == 3 & (length(xto) > 1 | length(yto) > 1))
    stop ("if 'M' is an array, then xto and yto should be one value")

  Transf <- function (xto, yto, u) {

  # find embracing values : first interval
    ix <- findInterval(xto, x )
    iy <- findInterval(yto, y)

  # next interval
    ixp1 <- pmin(ix+1,Nx)
    iyp1 <- pmin(iy+1,Ny)

  # interpolation factor
    xfac <- (xto-x[ix])/dx[ix]
    yfac <- (yto-y[iy])/dy[iy]

  # interpolate
    if (is.matrix(u))
     (1-yfac)*((1-xfac)*u[cbind(ix,iy)  ]+xfac*u[cbind(ixp1,iy)]) +
         yfac*((1-xfac)*u[cbind(ix,iyp1)]+xfac*u[cbind(ixp1,iyp1)])
    else
     (1-yfac)*((1-xfac)*u[ix,iy,  ]+xfac*u[ixp1,iy,]) +
         yfac*((1-xfac)*u[ix,iyp1,]+xfac*u[ixp1,iyp1,])
  } # end Transf

  if (length(Du) == 3) Transf(xto, yto, M) 
  else outer(xto, yto, FUN = Transf, u = M)
}
