vectorplot <- function(x, y, times = NULL, by = 1, arr = FALSE, tfac = NULL, ...) {
  
  dots <- splitpardots(list(...))
  dm <- dots$main
  dp <- dots$points
  
  if (is.null(dm$xlab)) dm$xlab <- ""
  if (is.null(dm$ylab)) dm$ylab <- ""
  
  ii <- seq(1, length(x), by = by)
  ll <- length(ii)

# the plot
  
  if (is.null(times)) {
    if (is.null(dm$ylim)) dm$ylim <- range(y[ii])
    if (is.null(dm$xlim)) dm$xlim <- range(x[ii])
    x0 <- rep(0, ll)
    y0 <- rep(0, ll)
    xe <- x[ii]
    ye <- y[ii]
    LL <- c(alist(0, type = "n"), dm)
    do.call("plot",LL)
  } else {
    ii <- seq(1, length(times), by = by)
    y0 <- rep(0, length(times[ii]))
    ye <- y[ii]
    if (is.null(tfac))
      tfac <- diff(range(times))/diff(range(c(y0, ye)))
    x0 <- times[ii]
    xe <- times[ii] + x[ii]*tfac
    if (is.null(dm$ylim)) dm$ylim <- range(c(y0, ye))
    if (is.null(dm$xlim)) dm$xlim <- range(c(x0, xe))
    LL <- c(alist(0, type = "n"), dm)
    do.call("plot",LL)
    pusr <- par("usr")
    if (is.null(tfac))
     tfac <- diff(pusr[1:2])/diff(pusr[3:4])
    xe <- times[ii] + x[ii]*tfac
  }

# the segments/arrows
  if (arr) {
    if (is.null(dp$arr.length)) dp$arr.length <- 0.1
    if (is.null(dp$arr.type)) dp$arr.type <- "triangle"
  }
  Ls <- c(alist(x0, y0, xe, ye), dp) 
  if (arr)  do.call("Arrows", Ls)
  else do.call("segments", Ls)
}


