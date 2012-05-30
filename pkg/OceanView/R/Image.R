

## =============================================================================
## =============================================================================
## Image functions
## =============================================================================
## =============================================================================
#Image <- function(...) UseMethod ("Image")
#Image.matrix <- function(...) Image.matrix(...)

Image <- function (z, col = femmecol(100), add.contour = FALSE, 
                   legend = TRUE, resfac = 1, NAcol = NULL, key.title= NULL, ...) {
  
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
   cntnames <- c("nlevels","levels","labels","labcex","drawlabels",
                    "method","vfont","lty","lwd")	
   x <- dots[["x"]]
   y <- dots[["y"]]
     
   if (resfac > 1) {
       if (is.null (x)) x <- 1: nrow(z)
       if (is.null (y)) y <- 1: ncol(z)
       diffx <- c(diff(x), 0)
       diffy <- c(diff(y), 0)
       XX <- x
       YY <- y 
       RR <- 1/resfac
       for (i in 2: resfac)  {
         XX <- c(XX, x + diffx * i*RR)
         YY <- c(YY, y + diffy * i*RR)
       }
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
    } 
    if (! is.null(dots[["y"]])) {
      y <- dots[["y"]]
      if (all(diff(y)<0)) {# swap
        if (is.null(dots$ylim)) dots$ylim <- rev(range(y))
        dots[["y"]] <- rev(dots[["y"]])
        z <- z[, (ncol(z):1)]
      }
    }     
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
     if(legend) drawlegend(parleg, LegendCol, LegendZlim, key.title) 
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

  if (is.null(ldots$main)) title <- isub
  else title <- rep(ldots$main, length.out = length(isub))
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
#    Image(z = zz,  ...)# extractdots(Ldots, i))) 

  }
  
}

## =============================================================================

Image.list <- function (z, ...) {
  
  if (! is.list(z)) 
    stop("'z' should be a list containing the matrix or array to plot")
  nz <- length(z)
  classz <- class(z[[1]])
  if (! classz %in% c("matrix","array"))
    stop ("'z' should be a list with either matrices or vectors")
  
    
  D1 <- dim(z[[1]])
  for (i in 2:nz)
    if (any(dim(z[[i]]) - D1 != 0))
      stop("elements of 'z' should have the same dimension, check element", i)
  
# Set the mfrow!
   nc <- min(ceiling(sqrt(nz)), 3)
   nr <- min(ceiling(nz/nc), 3)
   mfrow <- c(nr, nc)
   par(mfrow = mfrow)
   Ldots <- list(...) 
   Ldots$mfrow <- NULL
   if (!is.null(Ldots$main)) {
     main <- rep(Ldots$main, length.out = nz)
     Ldots$main <- NULL
   } else main <- 1:nz
   if (classz == "matrix")  {
     for (i in 1:nz) {
      Ldots$main <- main[i]
      LL <- c(list(z = z[[i]]), Ldots)
      do.call(Image, LL)
     }  
   } else stop("Not yet implemented for list of arrays")  
}
