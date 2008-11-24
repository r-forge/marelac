`rhoH2O` <-
function(S=35, T=25, P=0, method=c("Millero", "Chen", "marine", "limnological")) {
   method <- match.arg(method)
   rho <- switch(method,
     Millero      = rho(S, T, P),               # from package seacarb
     marine       = rho(S, T, P),               # from package seacarb
     Chen         = 1000 *rhoH2O_Chen(S, T, P), # limnological range
     limnological = 1000 *rhoH2O_Chen(S, T, P), # limnological range
   )
   attr(rho, "unit") = "(kg/m3)"
   return(rho)

}

