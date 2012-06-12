## =============================================================================
## =============================================================================
## Image S3 functions
## =============================================================================
## =============================================================================

Image <- function(z, ...) UseMethod ("Image")
Image.default <- function (z, ...) Image.matrix(z, ...)
Image.matrix <- function (z, col = femmecol(100), add.contour = FALSE, 
                   legend = TRUE, resfac = 1, NAcol = NULL, 
                   zlab = NULL, cex.zlab = par("cex.axis"), ...) {

  if (legend) {
    parplt <- par("plt") - c(0, 0.08, 0, 0)
    parleg <- c(parplt[2] + 0.02, parplt[2] + 0.05, parplt[3], parplt[4])
    plt.or <- par(plt = parplt)
    on.exit(par(plt = plt.or))
  }  
    dots <- list(...)
    
  # plotting parameters : split in plot parameters and point parameters
   imnames <- c("xlab", "ylab", "xlim", "ylim", "main", "sub", "log", "asp",
                 "ann", "axes", "frame.plot", "panel.first", "panel.last",
                 "cex.lab", "cex.axis", "cex.main")

  # plot.default parameters
   cntnames <- c("nlevels", "levels", "labels", "labcex", "drawlabels",
                    "method", "vfont", "lty", "lwd")	
   x <- dots[["x"]]
   y <- dots[["y"]]
   if (is.null (x)) x <- seq(0, 1, length.out = nrow(z))
   if (is.null (y)) y <- seq(0, 1, length.out = ncol(z))
   
   resfac<- rep(resfac, length.out = 2)  
   if (any(resfac > 1)) {
       diffx <- c(diff(x), 0)
       diffy <- c(diff(y), 0)
       XX <- x
       YY <- y 
       RX <- 1/resfac[1]
       RY <- 1/resfac[2]
       if (resfac[1] > 1)
         for (i in 1: resfac[1])  
           XX <- c(XX, x + diffx * i*RX)
       if (resfac[2] > 1)
         for (i in 1: resfac[2])  
           YY <- c(YY, y + diffy * i*RY)

       XX <- unique(sort(XX))
       YY <- unique(sort(YY))
       z <- mapxy(z, x = x, y = y, xto = XX, yto = YY)
       dots[["x"]] <- XX
       dots[["y"]] <- YY
    }
# Check for decreasing values of x and y    
    if (! is.null(dots[["x"]])) {
      x <- dots[["x"]]
      if (all(diff(x)<0)) {# swap
        if (is.null(dots$xlim)) dots$xlim <- rev(range(x))
        dots[["x"]] <- rev(dots[["x"]])
        z <- z[nrow(z):1, ]
      }
    } else dots[["x"]] <- x
    if (! is.null(dots[["y"]])) {
      y <- dots[["y"]]
      if (all(diff(y)<0)) {# swap
        if (is.null(dots$ylim)) dots$ylim <- rev(range(y))
        dots[["y"]] <- rev(dots[["y"]])
        z <- z[, (ncol(z):1)]
      }
    } else dots[["y"]] <- y
    if (is.null(dots[["xlab"]])) dots[["xlab"]] <- "x"
    if (is.null(dots[["ylab"]])) dots[["ylab"]] <- "y"

    ii <- which(names(dots) %in% cntnames)  
    dotcont <- dots[ii]
    
  # contour parameters
    ip <- c(!names(dots) %in% cntnames) 
    dotimage <- dots[ip]

 
     if (is.null(dots[["zlim"]])) {
       if (length(which(!is.na(z))) == 0)
         dotimage$zlim <- c(0,1)
       else dotimage$zlim <- range(z, na.rm = TRUE)
     } 
     if (!is.null(dots[["x"]])) dotimage$x <- dots[["x"]]
     if (!is.null(dots[["y"]])) dotimage$y <- dots[["y"]]
     
     LegendZlim <- dotimage$zlim
     LegendCol <- col
     if (any (is.na(z)) & ! is.null(NAcol) ) {
        N <- length(col)
        col <- c(NAcol, col)
        rr <- diff(dotimage$zlim)
        if (rr == 0) rr <- dotimage$zlim[1] *1e-3
        z[is.na(z)] <- dotimage$zlim[1] - 1/N*rr
        dotimage$zlim [1] <- dotimage$zlim[1] - 1/N*rr
     }

     do.call("image", c(alist(z = z, col = col), dotimage))
     box()
     if (add.contour) {
       if (!is.null(dots[["x"]])) dotcont$x <- dots[["x"]]
       if (!is.null(dots[["y"]])) dotcont$y <- dots[["y"]]
       do.call("contour", c(list(z = z, add = TRUE), dotcont))
     }  
     if(legend) drawlegend(parleg, LegendCol, LegendZlim, zlab, cex.zlab) 
}

## =============================================================================

Image.array <- function (z, margin = c(1, 2), subset, ask = NULL, ...) {
  
  DD <- dim(z)
  if (length(DD) != 3)
    stop ("Can only make image of 3-D array, 'z' has dimension ", length(DD))

  if (length(margin) != 2)
    stop ("'margin' should contain two numbers, the x, y subscripts of which to make images")
   
  if ( max(margin) > 3 | min (margin) < 1)
    stop ("indexes in 'margin' should be inbetween 1 and 3")

  index <- (1:3) [- margin]

  if (index > 3 || index <1)
    stop ("'index' to loop over should be inbetween 1 and 3")
  
  x <- 1:DD[index]
  
  if (!missing(subset)){
      e <- substitute(subset)
      r <- eval(e, as.data.frame(z), parent.frame())
      if (!is.logical(r))
          stop("'subset' must evaluate to logical")
      isub <- r & !is.na(r)
      isub <- which(isub)
      if (length(isub) == 0)
        stop("cannot continue: nothing selected - check 'subset'")
  } else isub <- x

  np     <- length(isub)
  ldots  <- list(...)
#  Ldots <- setdots(ldots, np)   # removed xlim and ylim first
  Ldots <- ldots 
  ## Set par mfrow and ask
  ask <- setplotpar(ldots, np, ask)
  if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
  }

  if (is.null(ldots$main)) 
    title <- isub
  else 
    title <- rep(ldots$main, length.out = length(isub))

  i1 <- 1
  for (i in isub) {
    if (index == 1) 
      zz <- z[i, , ] 
    else if (index == 2)
      zz <- z[ ,i , ] 
    else
      zz <- z[ ,, i ]
    if (margin[2] < margin[1])
      zz <- t(zz)

    LL <- c(list(z = zz), Ldots)
    LL$main <- title[i1]
    i1 <- i1+1
    do.call(Image, LL)

  }
  
}

## =============================================================================
## Image of a list of matrices or arrays
## =============================================================================

Image.list <- function (z, ...) {
  
# check z: list with similar matrices or arrays of dimension at most 3
  if (! is.list(z)) 
    stop("'z' should be a list containing the matrix or array to plot")

  nz     <- length(z)
  classz <- class(z[[1]])

  if (! classz %in% c("matrix", "array"))
    stop ("'z' should be a list with either matrices or arrays")
    
  DD <- dim(z[[1]])
  if (length(DD) > 3 | length(DD) < 2)
      stop ("Can only make image of 2-D or 3-D array, 'z' has dimension ", length(DD))

  for (i in 2 : nz)
    if (any(dim(z[[i]]) - DD != 0))
      stop("elements of 'z' should have the same dimension, check element", i)
  
# Set the mfrow argument
   if ("matrix" %in% classz)  {
     nc <- min(ceiling(sqrt(nz)), 3)
     nr <- min(ceiling(nz/nc), 3)
   } else { # differs from default in that it is not limited to 3
     nc <- ceiling(sqrt(nz))
     nr <- ceiling(nz/nc)
   }
   mfrow <- c(nr, nc)
   par(mfrow = mfrow)

   mfrow <- c(nr, nc)
   par(mfrow = mfrow)

# Plotting arguments
   Ldots <- list(...) 
   Ldots$mfrow <- mfrow
   if (!is.null(Ldots$main)) {
     main <- rep(Ldots$main, length.out = nz)
     Ldots$main <- NULL
   } else main <- 1:nz

  ask <- Ldots$ask
  if (is.null(ask)) ask <- TRUE
  Ldots$ask <- NULL

  # ylim and xlim can be lists and are at least two values
  yylim  <- expanddotslist(Ldots$ylim, nz)
  xxlim  <- expanddotslist(Ldots$xlim, nz)
  zzlim  <- expanddotslist(Ldots$zlim, nz)
  zzlab  <- expanddotslist(Ldots$zlab, nz)
   
  if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
  }

# Display the images
   if ("matrix" %in% classz)  {
     for (i in 1:nz) {
      Ldots$main <- main[i]
      Ldots$xlim <- xxlim[[i]]
      Ldots$ylim <- yylim[[i]]
      Ldots$zlim <- zzlim[[i]]
      Ldots$zlab <- zzlab[[i]]
      
      LL <- c(list(z = z[[i]]), Ldots)
      do.call(Image, LL)
     }  
   
   } else {  # array
     margin <- Ldots$margin
     Ldots$margin <- NULL
     
     if (is.null(margin)) margin <- 1:2
     if (length(margin) != 2)
       stop ("'margin' should contain two numbers, the x, y subscripts with which to make images")
     if ( max(margin) > 3 | min (margin) < 1)
       stop ("indexes in 'margin' should be inbetween 1 and 3")
     index <- (1:3) [- margin]

     subset <- Ldots$subset
     Ldots$subset <- NULL
     if (!is.null(subset)){
        e <- substitute(subset)
        r <- eval(e, as.data.frame(z), parent.frame())
        if (!is.logical(r))
          stop("'subset' must evaluate to logical")
        isub <- r & !is.na(r)
        isub <- which(isub)
        if (length(isub) == 0)
          stop("cannot continue: nothing selected - check 'subset'")
     } else isub <- 1:DD[index]
     
     nisub     <- length(isub)

     # number of empty plots 
     noplot <- prod(mfrow) - nz
     if (noplot == 0) noplot <- NULL else noplot <- 1:noplot

     # outer margin text
     Mtext <- Ldots$mtext
     Ldots$mtext <- NULL
     
     if(! is.null(Mtext)) 
       Mtext <- rep(Mtext, length.out = nisub)
     else
       Mtext <- isub
     pline <- par("oma")[3]-1  
     # loop first over margin, then over data sets
     for (jj in 1:nisub) {
       j <- isub[jj]      
       for (i in 1:nz) {
        if (index == 1) 
         zz <- z[[i]][j, , ] 
        else if (index == 2)
         zz <- z[[i]][ ,j , ] 
        else
         zz <- z[[i]][ ,, j ]
        if (margin[2] < margin[1])
         zz <- t(zz)
       
        Ldots$main <- main[i]
        Ldots$xlim <- xxlim[[i]]
        Ldots$ylim <- yylim[[i]]
        Ldots$zlim <- zzlim[[i]]
        Ldots$zlab <- zzlab[[i]]
        LL <- c(list(z = zz), Ldots)
        do.call(Image, LL)
       }
       # to make sure all figures are drawn
       for (i in noplot) emptyplot()  
       mtext(text = Mtext[jj], side = 3, outer = TRUE, line = pline)
     }  
   }
}


## =============================================================================
## Checking and expanding arguments in dots (...) with default
## =============================================================================

expanddots <- function (dots, default, n) {
  dots <- if (is.null(dots)) default else dots
  rep(dots, length.out = n)
}

# lists: e.g. xlim and ylim....
expanddotslist <- function (dots, n) {
  if (is.null(dots)) return(dots)
  dd <- if (!is.list(dots )) list(dots) else dots
  rep(dd, length.out = n)
}

## allow multiple titles for color key legend
##  if (is.null(Ldots$key.title)) {
##    key.title <- ""
##  } else {
##    key.title <- Ldots$key.title
##  }
##  if (np %% (length(key.title)) != 0) 
##   warning("length of key.title is not a multiple of elements in 'z'")
  
## make key.title same length as nz
##  key.title <- rep(key.title, length.out = np)

