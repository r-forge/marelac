## -----------------------------------------------------------------------------
## the Redfield ratio
## -----------------------------------------------------------------------------
redfield <- function(q, species, method=c("mol", "mass")) {
  if (!is.numeric(q)) stop("q must be numeric")
  method   <- match.arg(method)
  ratio    <- c(C=106, H=180, O=45, N=16, P=1)
  speciess <- names(ratio)
  species  <- match.arg(species, speciess)

  if (method == "mass") {
    ratio <- with(atomicweight,{
      ratio * c(C, H, O, N, P)
    })
  }
  p <- match(species, speciess)
  as.data.frame(q %o% ratio / ratio[p])
}

