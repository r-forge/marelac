redfield <- function(element, x, method=c("mol", "mass")) {
  method   <- match.arg(method)
  ratio    <- c(C=106, H=180, O=45, N=16, P=1)
  elements <- names(ratio)
  element  <- match.arg(element, elements)

  if (method == "mass") {
    ratio <- with(AtomicWeight,{
      ratio * c(C, H, O, N, P)
    })
  }
  p <- match(element, elements)
  x * ratio / ratio[p]
}

