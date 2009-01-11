vertmean <- function(depth, vari, level, top, bot, vol) {
  idepth  <- which(depth >= top & depth <= bot)
  ndepth  <- c(depth[idepth])
  nvari   <- c(vari[idepth])
  mdepth  <- c(top, (ndepth[-1] + ndepth[-length(ndepth)])/2, bot)
  nvol    <- vol(mdepth, level)
  dvol    <- abs(diff(nvol))
  sum(nvari * dvol) / sum(dvol)
}

