## =============================================================================
## =============================================================================
## General utility functions
## =============================================================================
## =============================================================================


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

drawlegend <- function (parleg, col, zlim) {
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

  image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "",
        ylab = "", col = col)

  do.call("axis", list(side = 4, mgp = c(3, 1, 0), las = 2))

  par(plt = Plt)
  par(usr = usr)
  par(xlog = PP$xlog)
  par(ylog = PP$ylog)
  par(new = FALSE)  
}
