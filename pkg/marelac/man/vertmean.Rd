\name{vertmean}
\alias{vertmean}

\title{Vertical Volumetric Weighted Mean of Matter Concentrations in Water Bodies}
\description{
  Calculate vertical mean values which respect depths of different layers and
  lake morphometry.
}
\usage{
vertmean(depth, vari, level, top, bot, vol)
}

\arguments{
  \item{depth}{sorted vector of sampled depths}
  \item{vari}{measurements corresponding to \code{depth} (concentration, temperature, \dots)}
  \item{level}{surface water level (above ground or m NN, depending on bathymetric function used}
  \item{top}{top water level of requested layer}
  \item{bot}{bottom water level of requested layer}
  \item{vol}{hypsographic function to be used (e.g. \code{vol.depth})}
}

\value{
  Volumetric average value of the requested layer between depths \code{top} and \code{bottom}.
}

\seealso{\code{\link{vol.tsb}}, \code{\link{vol.depth}}}
\examples{
## define a bathymetric formula for a given lake or basin
vol.tsb <- function(z, level) {
  ## Ground level 154m NN
  zz  <- level - 154 - z
  ## interpolate at very small
  vol <- ifelse(zz < 0.5,
    approx(c(0, 0.5), c(0, 18020.74), zz)$y,
    ## Polynomial (Thesis U. Miersch)
    vol <- 44331 - 202176 * zz + 363281 * zz^2 - 150342 * zz^3 + 47275.8 * zz^4 -
          6791.125 * zz^5 + 503.9106 * zz^6 - 18.95684 * zz^7 + 0.2871075 * zz^8
  )
  vol
}

## dummy formula for depth-weighted averaging
## (water column, instead of bathymetric curve)
vol.depth <- function(z, level) {z}


## Plot of bathymetric curve
x<-seq(0, 155, 0.01)
y<-vol.tsb(0, x)
plot(x, y, xlim = c(154, max(x)))

## Test Data
level <- 167.5

depth <- c(0, 1, 6, 7, 8, 9, 10, 11)
vari  <- c(18.8, 17, 12.1, 11.9, 11.8, 11.7, 11.4, 11.4)

vertmean(depth, vari, level, top=0, bot=4, vol.depth)
vertmean(depth, vari, level, top=4, bot=13.5, vol.depth)

vertmean(depth, vari, level, top=0, bot=4, vol.tsb)
vertmean(depth, vari, level, top=4, bot=13.5, vol.tsb)

vertmean(depth, vari, level, top=0, bot=3.7, vol.depth)
vertmean(depth, vari, level, top=3.7, bot=13.5, vol.depth)

vertmean(depth, vari, level, top=0, bot=3.7, vol.tsb)
vertmean(depth, vari, level, top=3.7, bot=13.5, vol.tsb)

}
\keyword{misc}

