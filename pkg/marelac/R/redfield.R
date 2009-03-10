## -----------------------------------------------------------------------------
## the Redfield ratio
## -----------------------------------------------------------------------------
redfield <- function(x, element, method=c("mol", "mass")) {
  if (!is.numeric(x)) stop("x must be numeric")
  method   <- match.arg(method)
  ratio    <- c(C=106, H=180, O=45, N=16, P=1)
  elements <- names(ratio)
  element  <- match.arg(element, elements)

  if (method == "mass") {
    ratio <- with(atomicweight,{
      ratio * c(C, H, O, N, P)
    })
  }
  p <- match(element, elements)
  as.data.frame(x %o% ratio / ratio[p])
}

