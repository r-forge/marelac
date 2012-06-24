
## =============================================================================
## =============================================================================
## QUIVER FUNCTIONS
## =============================================================================
## =============================================================================

checkinput <- function(u, v, x = NULL, y = NULL, scale = 1, by = 1) {
  
  if (is.null(x)) 
     x <- seq(0, 1, length.out = nrow(u))
  if (is.null(y)) 
     y <- seq(0, 1, length.out = ncol(v))

  if (is.array(x)) {
    if (length(dim(x)) ==1) 
     x <- matrix(ncol = ncol(u), nrow = length(x), data = x)
    else if (length(dim(x)) > 2) 
     stop ("'x' cannot be an array with more than 2 dimension")
  }

  if (is.array(y)) {
    if (length(dim(y)) ==1) 
     y <- matrix(nrow = nrow(v), ncol = length(y), data = y, byrow = TRUE)
    else if (length(dim(y)) > 2) 
     stop ("'y' cannot be an array with more than 2 dimension")
  }

  if (is.vector(x)) {
     x <- matrix(ncol = ncol(u), nrow = length(x), data = x)
  }
  if (is.vector(y)) {
     y <- matrix(nrow = nrow(v), ncol = length(y), data = y, byrow = TRUE)
  }

  
  
  
  div <- dim(u) - dim(v)
  
# u and v may have a dimension 1 different, if they are defined on the edges.       
  if (div[1] == -1 & div[2] == 1) {
     u <- 0.5*(u[,-1] + u[,-ncol(u)])
     v <- 0.5*(v[-1,] + v[-nrow(v),])
  } else if (div[1] == 1 & div[2] == -1)  {
     u <- 0.5*(u[-1,] + u[-nrow(u),])
     v <- 0.5*(v[,-1] + v[,-ncol(v)])
  } else if (div[1] != 0 | div[2] != 0)
    stop("dimensions of 'u' and 'v' not compatible") 

# Check x- and y
    Dim <- dim(u)
    
     if (is.matrix(x)) {
       if(ncol(x) == Dim[2] +1)
         x <- 0.5*(x[,-1] + x[,-ncol(x)])
       if(nrow(x) == Dim[1] +1)
         x <- 0.5*(x[-1,] + x[-nrow(x),])
       if(ncol(x) != Dim[2] | nrow(x) != Dim[1])
         stop ("'x' not compatible with u or v")
     }
     if (is.matrix(y)) {
       if(ncol(y) == Dim[2] +1)
         y <- 0.5*(y[,-1] + y[,-ncol(y)])
       if(nrow(y) == Dim[1] +1)
         y <- 0.5*(y[-1,] + y[-nrow(y),])
       if(ncol(y) != Dim[2] | nrow(y) != Dim[1])
         stop ("'x' not compatible with u or v")
     }

  if (is.null(x)) 
     x <- (row(u)-0.5)/nrow(u)
  if (is.null(y)) 
     y <- (col(v)-0.5)/ncol(v)

# ------------------------------------------------------------------------------
# select the elements of x, y, u and v
# ------------------------------------------------------------------------------    
  if (! is.null(by)) {                     
     by <- rep(by, length = 2)
     ix <- seq(1, nrow(u), by = by[1])
     iy <- seq(1, ncol(u), by = by[2])
     u <- u[ix,iy]
     x <- x[ix,iy]
     v <- v[ix,iy]
     y <- y[ix,iy]
  } else {
    ix <- 1:nrow(u)
    iy <- 1:ncol(u)
  }
  isna <- which (is.na(u))
  u[is.na(u)] <- 0
  v[is.na(v)] <- 0
    
# ------------------------------------------------------------------------------
# size of the arrows
# ------------------------------------------------------------------------------    
  speed <- sqrt(u^2 + v^2)
  maxspeed <- max(speed)
  xr <- diff (range(x)) / max(dim(u))
  yr <- diff (range(y)) / max(dim(u))


  if (!is.null(scale)) {
    u <- u * scale / maxspeed * xr 
    v <- v * scale / maxspeed * yr 
  }
  list(x = x, y = y, u = u, v = v, 
       speed = speed, maxspeed = maxspeed, 
       isna = isna, ix = ix, iy = iy)
}

## =============================================================================
## quiver function - uses some ideas from: 
## from http://tolstoy.newcastle.edu.au/R/help/01c/2711.html
## Robin Hankin Tue 20 Nov 2001 - 13:10:28 EST
## =============================================================================
quiver <- function(u, ...) UseMethod ("quiver")
quiver.default <- function (u, ...) quiver.matrix(u, ...)
quiver.matrix  <- function(u, v, x = NULL, y = NULL,           #   pol = NULL,
                    scale = 1, arr.max = 0.2, arr.min = 0, 
                    by = NULL, add = FALSE, ...)  {
# ------------------------------------------------------------------------------
# check input
# ------------------------------------------------------------------------------
  MM <- checkinput(u, v, x, y, scale, by = by) 

# matrices have rows for y and columns for x - should be transposed
  x <- t(MM$x)
  y <- t(MM$y)
  u <- t(MM$u)
  v <- t(MM$v)

  xto <- x + u
  yto <- y + v

  dotlist  <- splitpardots( list(...) )
  dp <- dotlist$points
  NAcol <- dp$NAcol
  dp$NAcol <- NULL
  # 

  # transpose dp elements that are matrices
  dp <- lapply(dp, FUN = function(x) if (is.matrix(x)) t(x[MM$ix, MM$iy]) else x)
  dm <- dotlist$main
#
  if (is.null(dm$xlab)) dm$xlab <- "x"
  if (is.null(dm$ylab)) dm$ylab <- "y"

  if (!add) do.call("matplot", c(alist(rbind(x, xto), rbind(y, yto)), 
                     type = "n",dm))

  if (is.null(dp$arr.type))  dp$arr.type <- "triangle"  
  if (is.null(dp$arr.lwd))   dp$arr.lwd <- 1
  dp$arr.length <- t(MM$speed) / MM$maxspeed * (arr.max - arr.min) + arr.min
  
  do.call("Arrows", c(alist(x, y, xto, yto),  dp))
}

## =============================================================================

quiver.array <- function (u, v, margin = c(1, 2), subset, ask = NULL, ...) {

  DD <- dim(u)
  if (length(DD) != 3)
    stop ("Can only make quiver of 3-D array, 'u' has dimension ", length(DD))
  if (length(dim(v)) != 3)
    stop ("Can only make quiver of 3-D array, 'v' has dimension ", length(dim(v)))

  if (length(margin) != 2)
    stop ("'margin' should contain two numbers, the x, y subscripts of which to make quiver plots")
   
  if ( max(margin) > 3 | min (margin) < 1)
    stop ("indexes in 'margin' should be inbetween 1 and 3")

  index <- (1:3) [- margin]

  if (index > 3 || index <1)
    stop ("'index' to loop over should be inbetween 1 and 3")
  
  x <- 1:DD[index]
  
  if (!missing(subset)){
      e <- substitute(subset)
      r <- eval(e, as.data.frame(x), parent.frame())
      if (!is.logical(r))
          stop("'subset' must evaluate to logical")
      isub <- r & !is.na(r)
      isub <- which(isub)
  } else isub <- x

  np     <- length(isub)
  ldots  <- list(...)

  ## Set par mfrow and ask
  ask <- setplotpar(ldots, np, ask)
  if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
  }
  if (is.null(ldots$main)) title <- isub
  else   title <- rep(ldots$main, length.out = length(isub))
  i1 <- 1
  for (i in isub) {
    LL <- ldots
    if (index == 1) { 
      LL$u <- u[i, , ] 
      LL$v <- v[i, , ] 
      }
    else if (index == 2) {
      LL$u <- u[ ,i , ] 
      LL$v <- v[ ,i , ] 
    }  
    else  {
      LL$u <- u[ ,, i ]
      LL$v <- v[ ,, i ]
    }  
    if (margin[2] < margin[1]){
      LL$u <- t(LL$u)
      LL$v <- t(LL$v)
    }  
    LL$main <- title[i1]
    i1 <- i1 +1
    do.call(quiver, LL)
    #quiver(u = uu, v = vv, ...) 
  }
  
}

## =============================================================================
## =============================================================================
## flowpaths based on flow fields
## =============================================================================
## =============================================================================

flowpath <- function (u, v, x = NULL, y = NULL, startx = NULL, starty = NULL,
                       scale = 1, plot = TRUE, add = FALSE, 
                       numarr = 0, arr.length = 0.2, maxstep = 1000, ...) {

  if (is.null(x)) 
     x <- seq(0, 1, length.out = nrow(u))
  if (is.null(y)) 
     y <- seq(0, 1, length.out = ncol(v))

  MM <- checkinput(u, v, x, y, scale, by = 1) 
  
  x <- MM$x
  y <- MM$y
  u <- MM$u
  v <- MM$v

  Nx <- nrow(x)
  Ny <- ncol(y)

  xx <- x[,1]
  yy <- y[1,]
  dx <- c(diff(xx), xx[Nx]-xx[Nx-1])
  dy <- c(diff(yy), yy[Ny]-yy[Ny-1])

  if(is.null(startx) ) {
    startx <- c(rep(xx[1], Ny),rep(xx[Nx], Ny), xx, xx) 
    starty <- c(yy, yy, rep(yy[1], Nx), rep(yy[Ny], Nx))
  } 
  
  lx <- length(startx)
  ly <- length(starty)
  if (lx != ly) {
    startx <- rep(startx, len = ly)
    starty <- rep(starty, len = length(startx))
  }
  
  xto <- startx
  yto <- starty
  
  if (min(xto) < min(x)) stop("'x' should embrace 'startx'")
  if (min(yto) < min(y)) stop("'y' should embrace 'starty'")
    
  interp <- function (xto, yto, u) {

  # find embracing values : first interval
    ix <- findInterval(xto, xx)
    iy <- findInterval(yto, yy)

  # next interval
    ixp1 <- pmin(ix+1,Nx)
    iyp1 <- pmin(iy+1,Ny)

  # interpolation factor
    xfac <- (xto-xx[ix])/dx[ix]
    yfac <- (yto-yy[iy])/dy[iy]

  # interpolate
    (1-yfac)*((1-xfac)*u[cbind(ix,iy)]+xfac*u[cbind(ixp1,iy)]) +
    yfac*((1-xfac)*u[cbind(ix,iyp1)]+xfac*u[cbind(ixp1,iyp1)])

  } # end Transf

  if (plot) {
    dotlist  <- splitpardots( list(...) )
    dp <- dotlist$points
    dm <-  dotlist$main

    if (is.null(dm$xlab)) dm$xlab <- "x"
    if (is.null(dm$ylab)) dm$ylab <- "y"
 
    if (!add) do.call("matplot", c(alist(x, y), 
                      type = "n",dm))
    dpa <- dp                   
    dp[grep("arr",names(dp))] <- NULL
    if ( numarr >0) {
      if (is.null(dpa$arr.type))  dpa$arr.type <- "triangle"  
      dpa$arr.length <- arr.length
    }
  }
#  XY <- cbind 
  xr <- range(xx)
  yr <- range(yy)
  mm <- min(diff(xx), diff(yy))*1e-12
  ij <- length(startx)

  xy <- NULL
  
  for (j in 1:ij) {
    xto <- startx[j]
    yto <- starty[j]       
    xy <- rbind(xy, c(xto, yto))
    for (i in 2: maxstep) {
      u2 <- interp(xto, yto, u) 
      v2 <- interp(xto, yto, v) 
      xto <- xto + u2
      yto <- yto + v2
      xy <- rbind(xy, c(xto, yto))

      if (xto < xr[1] | xto > xr[2] | yto < yr[1] | yto > yr[2] | (abs(u2) < mm & abs(v2) < mm)) break
    }
    if (plot) {
      do.call ("lines",c(alist(xy[1:i, ]),dp))                     

      if( numarr > 0 & i > 2) {
        i1 <- 1/( numarr +1)
        ii <- pmax(2,as.integer(seq(i*i1, i *(1-i1), len =  numarr )))
        do.call("Arrows", 
                c(alist(xy[ii-1,1], xy[ii-1,2], xy[ii,1], xy[ii,2]), dpa, cex = 0.5))
      }
      xy <- NULL

    } else {
      xy <- rbind(xy, c(NA, NA))
    }
  }
  invisible(xy)
}

