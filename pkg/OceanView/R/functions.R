## =============================================================================
## =============================================================================
## General utility functions
## =============================================================================
## =============================================================================
ImageOcean <- function(...) {
  dots <- list(...)
  if (is.null(dots$xlab)) 
    dots$xlab  <- "latitude"
  if (is.null(dots$ylab)) 
    dots$ylab  <- "longitude"
  if (is.null(dots$zlab)) 
    dots$zlab  <- "depth, m"
  if (is.null(dots$NAcol)) 
    dots$NAcol  <- "black"
  zz       <- Bathymetry$z
  zz[zz>0] <- NA

  
  do.call("Image", c(alist(zz, x = Bathymetry$x, y = Bathymetry$y), dots))
}


## =============================================================================
## Expanding arguments in dots  (...)
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
## set mfrow and ask
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
## Split plotting parameters in general (main) and point parameters
## =============================================================================

splitpardots <- function(dots) {

  nmdots <- names(dots)

  # plotting parameters : split in plot parameters and point parameters
  plotnames <- c("xlab", "ylab", "xlim", "ylim", "main", "sub", "log", "asp",
                 "ann", "axes", "frame.plot", "panel.first", "panel.last",
                 "cex.lab", "cex.axis", "cex.main","col.lab", "col.axis", "col.main")

  # plot.default parameters
  ii <- names(dots) %in% plotnames
  dotmain <- dots[ii]

  # point parameters
  ip <- !names(dots) %in% plotnames
  dotpoints <- dots[ip]
  list(points = dotpoints, main = dotmain)
}

## =============================================================================
## Legend on right hand side of plot
## =============================================================================

drawlegend <- function (parleg, col, zlim, zlab = NULL, cex.zlab = NULL) {
  Plt <- par(plt = parleg)
  PP <- par()
  par(new = TRUE)
  usr <- par("usr")
  ix <- 1
  minz <- zlim[1]
  maxz <- zlim[2]
  binwidth <- (maxz - minz)/64
  iy <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
  iz <- matrix(iy, nrow = 1, ncol = length(iy))
  if (! is.numeric(cex.zlab)) cex.zlab <- 1.
  image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "",
        ylab = "", col = col, main = zlab, cex.main = cex.zlab)
#  cex <- par("cex")
#  mtext(key.title, 3, line = 1, cex = cex)
  
  do.call("axis", list(side = 4, mgp = c(3, 1, 0), las = 2, cex.axis = cex.zlab))
  box() # thpe: to avoid different thicknes of axis and surrounding box

  par(plt = Plt)
  par(usr = usr)
  par(xlog = PP$xlog)
  par(ylog = PP$ylog)
  par(new = FALSE)  
}
