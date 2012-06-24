## =============================================================================
## =============================================================================
## Matrix plotting - most of these function were slightly modified from 
## similar functions in the deSolve and rootSolve packages.
## =============================================================================
## =============================================================================

## ============================================================================
## Finds common variables in a set of matrices
## ============================================================================

Mcommon <- function(M, ..., verbose = FALSE) {     
# creates a list with subsets of matrices that have only common elements

  Colnames <- function(M) 
    if (is.null (cn <- colnames(M))) return (1:ncol(M)) else return(cn)

  LL <- list(...)
  if (length(LL) == 0)
    return(list(M))

  if (is.list(M)) {   
    if (! class(M[[1]]) %in% c("matrix", "data.frame"))
      stop ("elements in list 'M' should be either a 'matrix' or 'data.frame' ")
     LL <- c(M, LL)
  } else {
    if (! class(M) %in% c("matrix", "data.frame"))
      stop ("'M' should be either a 'matrix' or 'data.frame' or a 'list'")
     LL <- c(list(M), LL)
  }
  cn <- Colnames(LL[[1]])
  for (i in 2:length(LL))  {
    if (! class(LL[[i]]) %in% c("matrix", "data.frame"))
      stop ("elements in '...' should be either a 'matrix' or 'data.frame'")
    cn <- cn[cn %in% Colnames(LL[[i]])]
  }
  if (verbose) 
    print(paste("common variable names: ", paste(cn, collapse = ", ")))
  if (length(cn) ==  0) {
    if (verbose) 
      warning("No variable names in common - returning 'NULL'")
    LL <- NULL
  } else {
    if (length(cn) ==  1 & verbose) 
    warning("Only one variable name in common - returning list of vectors")
    for (i in 1:length(LL))  
       LL[[i]] <- LL[[i]][, cn]
  }
  return(LL)
}

## ============================================================================
## Splits a matrix according to values in 'column'; the result is a list
## ============================================================================

Msplit <- function(M, 
                   column = 1, 
                   subset = NULL) {
  LL <- list()
  
  # quick and dirty
  if (!missing(subset)){
      e <- substitute(subset)
      r <- eval(e, as.data.frame(M), parent.frame())
      if (!is.logical(r))
          stop("'subset' must evaluate to logical")
      isub <- r & !is.na(r)
    M <- M[isub, ]
  }  

  ux <- unique(M[ , column])
  if (length(column) == 1) ux <- matrix(ncol = 1, data = ux)
  lux <- nrow(ux)

  isel <- 1 : ncol(M)
  if (is.numeric(column)) 
    isel <- isel[-column]
  else
    isel <- isel[ -which(colnames(M) %in% column)]

  for (i in 1:lux) {
    Sel <- M 
    for (j in 1:length(column))
      Sel <- Sel[Sel[ ,column[j]] == ux[i,j], ]
 
    LL[[i]] <- Sel[isel] 
  }  

  lnames <- ux[,1]
  if (length(column) > 1) for (i in 2:length(column))  
     lnames <- paste(lnames, ux[,i])

  names(LL) <- lnames

  LL
}

## =============================================================================
## Plot a (list of) matrices
## =============================================================================

Mplot <- function (M, ..., 
                   x = 1, 
                   select = NULL, which = select, 
                   subset = NULL, ask = NULL, 
                   legendpar = list(x = "top")) {

  getnames <- function(x) {
    if (is.null (cn <- colnames(x))) return (1:ncol(x)) else return(cn)
  }
                   
  # The ellipsis
  ldots   <- list(...)
  
  mtext <- ldots$mtext
  ldots$mtext <- NULL
  
  Dots    <- splitdots(ldots, M)
  x2      <- Dots$x2
  nother  <- Dots$nother
  nx      <- nother + 1 # total number of objects to be plotted
  varnames <- getnames(x2[[1]])

# x-variable  
  xPos <- vector()
  for (i in 1: length(x2))
    xPos[i] <- selectvar(x, getnames(x2[[i]]))
  xisfactor <- is.factor(x2[[1]][,xPos[1]])
  xname <- varnames[xPos[1]]

 # variables to be plotted
  Which <- which
  if (is.null(Which)) {
    for (i in 1: length(x2))
      Which <- c(Which,getnames(x2[[i]])[- xPos[i]])
    Which <- unique(Which)
  }

  np      <- length(Which)  
  
  # Position of variables to be plotted in "M" and other matrices
  xWhich <- list()

  for (i in 1: length(x2))
    xWhich[[i]] <- selectvar(Which, getnames(x2[[i]]))

  if (! is.character(Which)) 
    Which <- varnames[xWhich[[1]]]

  # number of figures in a row and interactively wait if remaining figures
  ask <- setplotpar(ldots, np, ask)                   
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  Dotmain <- setdots(Dots$main, np)  # expand to np for each plot

  # these are different from the default
  Dotmain$xlab <- expanddots(ldots$xlab, xname      , np)
  Dotmain$ylab <- expanddots(ldots$ylab, ""         , np)
  Dotmain$main <- expanddots(ldots$main, Which      , np)

  # ylim and xlim can be lists and are at least two values
  yylim  <- expanddotslist(ldots$ylim, np)
  xxlim  <- expanddotslist(ldots$xlim, np)

  Dotpoints <- setdots(Dots$points, nx)   # expand all dots to nx values

  # these are different from default
  Dotpoints$type <- expanddots(ldots$type, "l", nx)
  Dotpoints$lty  <- expanddots(ldots$lty, 1:nx, nx)
  Dotpoints$pch  <- expanddots(ldots$pch, 1:nx, nx)
  Dotpoints$col  <- expanddots(ldots$col, 1:nx, nx)
  Dotpoints$bg   <- expanddots(ldots$bg,  1:nx, nx)

  if (!missing(subset)){
    isub <- list()
    for (i in 1:nx) {
      e <- substitute(subset)
      r <- eval(e, as.data.frame(x2[[i]]), parent.frame())
      if (!is.logical(r))
          stop("'subset' must evaluate to logical")
      isub[[i]] <- r & !is.na(r)
    }  
  } else isub <- rep(TRUE, nx)


  # LOOP for each output variable (plot)
   
  for (ip in 1 : np) {

    # plotting parameters for matrix output 1 (opens a plot)
    dotmain   <- extractdots(Dotmain, ip)
    dotpoints <- extractdots(Dotpoints, 1)  # 1st dotpoints

    Xlog <- Ylog <- FALSE
    if (! is.null(dotmain$log)) {
      Ylog  <- length(grep("y",dotmain$log))
      Xlog  <- length(grep("x",dotmain$log))
    }
                                                           
    # first object plotted (new plot created)
    ix <- xWhich[[1]][[ip]]      # position of variable in 'x'

    if (is.null(yylim[[ip]]))
      dotmain$ylim <- SetRange(yylim[[ip]], x2, isub, xWhich, ip, Ylog)
    else
      dotmain$ylim <- yylim[[ip]]

    if (is.null(xxlim[[ip]])) {
      dotmain$xlim <- SetRange(xxlim[[ip]], x2, isub, xPos, 1, Xlog)
      if(xisfactor) dotmain$xlim <- dotmain$xlim + c(-0.5, 0.5)
    } else
      dotmain$xlim <- xxlim[[ip]]
    
    do.call("plot", c(alist(x2[[1]][isub, xPos[1]], x2[[1]][isub[[1]], ix]), dotmain, dotpoints))

    if (nother > 0)        # if other outputs
      for (j in 2:nx) {
        ix <- xWhich[[j]][[ip]]      # position of variable in 'x2'
        if (!is.na(ix))
        do.call("lines", c(alist(x2[[j]][isub[[j]], xPos[j]], x2[[j]][isub[[j]], ix]),
                extractdots(Dotpoints, j)) )
      }
  }
  
  if(! is.null(mtext))
    mtext(outer = TRUE, side = 3, mtext, line = par()$oma[3]-1, 
          cex = par()$cex*1.2)
     
  # Add legend, if legendpar not equal to NULL
  if(length(legendpar) > 0) {
    if(!is.list(legendpar)) 
      stop ("'legendpar' should be a list or NULL")

    if(is.null(legendpar$col))
      legendpar$col <- Dotpoints$col

    if(is.null(legendpar$pt.bg))
      legendpar$pt.bg <- Dotpoints$pt.bg

    if(is.null(legendpar$lwd))
      legendpar$lwd <- Dotpoints$lwd

    if(is.null(legendpar$pch)) {
      legendpar$pch <- Dotpoints$pch
      legendpar$pch[Dotpoints$type == "l"] <- NA
    }
    if(is.null(legendpar$lty)) {
      legendpar$lty <- Dotpoints$lty
      legendpar$lty[Dotpoints$type == "p"] <- NA
      if (all(is.na(legendpar$lty)))
        legendpar$lty <- NULL

    }

    if(is.null(legendpar$legend))
      legendpar$legend <- names(x2)
    if(is.null(legendpar$x))
      legendpar$x <- "top"
    do.call("legend", legendpar)

  }
}

## =============================================================================
## Plot a matrix with z-colors
## =============================================================================

scatterplot <- function (M, x = 1, 
                        select = NULL, which = select, 
                        zvar = x, subset = NULL, 
                        col = NULL, ask = NULL, discrete = FALSE,
                        legendpar = list(x = "top"),
                        add = FALSE, ... ) {

  getnames <- function(x) {
    if (is.null (cn <- colnames(x))) return (1:ncol(x)) else return(cn)
  }

  # The ellipsis
  ldots   <- list(...)
  mtext <- ldots$mtext
  ldots$mtext <- NULL
  
  Dots    <- splitdots(ldots, M)

  varnames <- getnames(M)

# x-variable  
  xPos <- selectvar(x, varnames)
  xname <- varnames[xPos]

# color-variable  
  cPos <- selectvar(zvar, varnames)
  legend <- TRUE
  if (length(legendpar) == 0) legend <- FALSE

  if(discrete) {
    x2 <- unique(M[, cPos])
    if (is.null(col))
      col <- femmecol(length(x2))
  } else {    
    if (is.null(col)) col <- femmecol(100)
    if (legend) {
     if (!add) 
       parplt <- par("plt") - c(0, 0.08, 0, 0)
     else
       parplt <- par("plt") 
     parleg <- c(parplt[2] + 0.02, parplt[2] + 0.05, parplt[3], parplt[4])
     plt.or <- par(plt = parplt)
     on.exit(par(plt = plt.or))
   }
 }
 # variables to be plotted
  Which <- which
  if (is.null(Which)) 
    Which <- unique(getnames(M)[- c(cPos, xPos)])

  np      <- length(Which)  
  
  if(add & np > 1)
    stop ("'add' only makes sense if only one variable should be plotted")
    
  # Position of variables to be plotted  
  xWhich <- selectvar(Which, varnames)
  if (! is.character(Which)) 
    Which <- varnames[xWhich]
 
  # number of figures in a row and interactively wait if remaining figures
  ask <- setplotpar(ldots, np, ask)                   
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  Dotmain <- setdots(Dots$main, np)  # expand to np for each plot

  # these are different from the default
  Dotmain$xlab <- expanddots(ldots$xlab, xname      , np)
  Dotmain$ylab <- expanddots(ldots$ylab, ""         , np)
  Dotmain$main <- expanddots(ldots$main, Which      , np)

  # ylim and xlim can be lists and are at least two values
  yylim  <- expanddotslist(ldots$ylim, np)
  xxlim  <- expanddotslist(ldots$xlim, np)
  zzlim  <- expanddotslist(ldots$zlim, np)
  
  Dotpoints <- Dots$points
  Dotpoints$zlim <- NULL

  # these are different from default
  Dotpoints$type <- "p"
  Dotpoints$lty  <- expanddots(ldots$lty, 1, 1)
  Dotpoints$pch  <- expanddots(ldots$pch, 18, 1)
  Dotpoints$bg   <- expanddots(ldots$bg,  1, 1)

       
  if (!missing(subset)){
      e <- substitute(subset)
      r <- eval(e, as.data.frame(M), parent.frame())
      if (!is.logical(r))
          stop("'subset' must evaluate to logical")
      isub <- r & !is.na(r)
  } else isub <- TRUE


  # LOOP for each output variable (plot)
  if (is.null(Dotpoints$zlab))
    zlab <- varnames[cPos]
  else
    zlab <- Dotpoints$zlab

  cex.zlab <- Dotpoints$cex.zlab
  Dotpoints$zlab <- NULL 
  Dotpoints$cex.zlab <- NULL 
  if (length(legendpar) > 0) {
    legendpos <- legendpar$x
#    if (is.null(legendpos)) legendpos <- "top"
    legendpos <- rep(legendpos , length.out  = np)
  } else legendpos <- NULL  
  for (ip in 1 : np) {

    # plotting parameters for matrix output 1 (opens a plot)
    dotmain   <- extractdots(Dotmain, ip)
    dotpoints <- extractdots(Dotpoints, 1)  # 1st dotpoints

    Xlog <- Ylog <- FALSE
    if (! is.null(dotmain$log)) {
      Ylog  <- length(grep("y",dotmain$log))
      Xlog  <- length(grep("x",dotmain$log))
    }
                                                           
    # first object plotted (new plot created)
    ix <- xWhich[[ip]]      # position of variable in 'x'
    dotmain$ylim <- Range(NULL, M[isub, xWhich[ip]], Ylog)
    dotmain$xlim <- Range(NULL, M[isub, xPos], Xlog)
    # rescale from color variable to 1: length(col)
    Xcol <- as.double(M[,cPos])

    if (! is.null(zzlim[[ip]]))
     zlim <- zzlim[[ip]]
    else  
     zlim <- range(Xcol, na.rm = TRUE)    

    xcol <- 1 + (Xcol - zlim[1])*(length(col)-1)/diff(zlim) 
    
    dcol <- col[xcol]
    dcol [is.na(dcol )] <- "black"
    dotpoints$col <- dcol

    if (add) 
      do.call("points", c(alist(M[isub, xPos], M[isub, ix]), dotpoints))
    
    else
      do.call("plot", c(alist(M[isub, xPos], M[isub, ix]), dotmain, dotpoints))
 
    if (discrete & ! is.null(legendpos[ip])) {
    legendpar <- list()

    legendpar$col <- unique(dcol)
    legendpar$pt.bg <- Dotpoints$pt.bg
    legendpar$lwd <- Dotpoints$lwd
    legendpar$pch <- Dotpoints$pch
    legendpar$lty <- NULL
    legendpar$legend <- as.character(x2)
    legendpar$x <- legendpos[ip]
    
    if(is.null(legendpos))
      legendpar$x <- "top"
    do.call("legend", legendpar)

    } else
      if(legend) drawlegend (parleg, col, zlim, zlab = zlab, cex.zlab = cex.zlab)  
  }

  if(! is.null(mtext))
    mtext(outer = TRUE, side = 3, mtext, line = par()$oma[3]-1, 
          cex = par()$cex*1.2)
  
}


## =============================================================================
## Update range, taking into account neg values for log transformed values
## =============================================================================

Range <- function(Range, x, log) {
   if (is.null(x)) return(Range)
   if (log)
      x[x <= 0] <- min(x[x>0])  # remove zeros
   
   RR <- range(Range, as.double(x), na.rm = TRUE)
   RR[is.infinite(RR)]<- NA

   return( RR )
}


SetRange <- function(lim, x2, isub, xWhich, ip, Log) {

  nx <- length (x2)
  if ( is.null (lim)) {
    yrange <- NULL
      for (j in 1:nx){
        ix <- xWhich[[j]][ip]
        yrange <- Range(yrange, x2[[j]][isub[[j]],ix], Log)
      }  
  } else
     yrange  <- lim

  return(yrange)
}


## =============================================================================
## Set the mfrow parameters and whether to "ask" for opening a new device
## =============================================================================

setplotpar <- function(ldots, nv, ask) {
  nmdots <- names(ldots) 
  # nv = number of variables to plot
  if (!any(match(nmdots, c("mfrow", "mfcol"), nomatch = 0))) {
     nc <- min(ceiling(sqrt(nv)), 3)
     nr <- min(ceiling(nv/nc), 3)
     mfrow <- c(nr, nc)
  } else if ("mfcol" %in% nmdots)
     mfrow <- rev(ldots$mfcol)
  else mfrow <- ldots$mfrow

  if (! is.null(mfrow))  mf <- par(mfrow = mfrow)

  ## interactively wait if there are remaining figures
  if (is.null(ask))
    ask <- prod(par("mfrow")) < nv && dev.interactive()

  return(ask)
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
## =============================================================================

repdots <- function(dots, n)
  if (is.function(dots)) dots else rep(dots, length.out = n)

setdots <- function(dots, n) lapply(dots, repdots, n)

## =============================================================================
## Extracting element 'index' from dots  (...)
## =============================================================================

extractdots <- function(dots, index) {
  ret <- lapply(dots, "[", index)
  ret <- lapply(ret, unlist) # flatten list
  return(ret)
}

## =============================================================================
## find a variable  - and keep the ordering
## =============================================================================

selectvar <- function (Which, var, NAallowed = TRUE) {
  if (!is.numeric(Which)) {
    ln <- length(Which)

    ## the loop is necessary so as to keep ordering...
    Select <- NULL
    for ( i in 1:ln) {
      ss <- which(Which[i] == var)
      if (length(ss) ==0 & ! NAallowed)
        stop("variable ", Which[i], " not in variable names")
      else if (length(ss) == 0)
        Select <- c(Select, NA)
      else
        Select <- c(Select, ss)
    }

  } else {
    Select <- Which  # "Select" now refers to the column number
    if (max(Select) > length(var))
        stop("index in 'which' too large: ", max(Select)-1)
    if (min(Select) < 1)
        stop("index in 'which' should be > 0")
  }
  return(Select)
}

## ============================================================================
## create several lists: x2:   other matrix objects,
##                       dotmain, dotpoints: remaining (plotting) parameters
## ============================================================================

splitdots <- function(ldots, x){
  x2      <- list()
  nother <- 0
  islist <- (! is.data.frame(x) & is.list(x))
  
  if (! islist) {
    x2[[1]] <- x
    names(x2)[1] <-"M"
    }
  else {
   for(i in 1:length(x))
    x2[[i]] <- x[[i]]
   names(x2) <- names(x)
   nother <- length(x) - 1
   
  }

  dots   <- list()
  nd     <- 0
  ndots <- names(ldots)
    
  if (length(ldots) > 0)
    for ( i in 1:length(ldots))
      if ("matrix" %in% class(ldots[[i]]) | "data.frame" %in% class(ldots[[i]])) { # a deSolve object
        nother <- nother + 1        
        x2[[nother + 1]] <- ldots[[i]]
        if (is.null(ndots[i]))
          names(x2)[nother] <- nother 
        else 
          names(x2)[nother] <- ndots[i]
        # a list of matrix objects
      } else if (is.list(ldots[[i]]) & 
        ("matrix" %in% class(ldots[[i]][[1]]) | 
         "data.frame" %in% class(ldots[[i]][[1]]))) {
        for (j in 1:length(ldots[[i]])) {
          nother <- nother + 1        
          x2[[nother+1]] <- ldots[[i]][[j]]
          nn <- names(ldots[[i]])[[j]]
          if (is.null(nn)) 
            nn <- nother
          names(x2)[nother] <- nn
        }
      } else if (! is.null(ldots[[i]])) {  # a graphical parameter
        dots[[nd <- nd+1]] <- ldots[[i]]
        names(dots)[nd] <- ndots[i]
      }

  nmdots <- names(dots)

  # plotting parameters : split in plot parameters and point parameters
  plotnames <- c("xlab", "ylab", "xlim", "ylim", "main", "sub", "log", "asp",
                 "ann", "axes", "frame.plot", "panel.first", "panel.last",
                 "cex.lab", "cex.axis", "cex.main")

  # plot.default parameters
  ii <- names(dots) %in% plotnames
  dotmain <- dots[ii]

  # point parameters
  ip <- !names(dots) %in% plotnames
  dotpoints <- dots[ip]
  list(points = dotpoints, main = dotmain, nother = nother, x2 = x2)
}

