## -----------------------------------------------------------------------------
## Density of Freshwater
## -----------------------------------------------------------------------------

rhoH2O_Chen <- function(S=0, t=25, p) {
  # Chen, Ch.-T. and F.J. Millero (1986) - Precise thermodynamic
  # properties of natural waters covering only the limnological range.
  # Limnol. Oceanogr. 31 No. 3, 657 - 662
  # rho in g/cm^3
  # P in bar - surface = 0 bar.
  # t in °C
  P <- p
  rho0 <- 0.9998395 + 0.000067914 * t -
          0.0000090894 * t * t + 0.00000010171 * t * t * t -
          0.0000000012846 * t ^ 4 + 0.000000000011592 * t ^ 5 -
          5.0125E-14 * t ^ 6 + (0.0008181 - 0.00000385 * t + 0.0000000496 * t * t) * S
  K    <- 19652.17 + 148.113 * t * -2.293 * t * t + 0.01256 * t ^ 3 - 0.0000418 * t ^ 4 +
       (3.2726 - 0.0002147 * t + 0.0001128 * t ^ 2) * P + (53.238 - 0.313 * t + 0.005728 * P) * S
  rho0 * (1 - P / K) ^ -1
}


rhoH2O <- function(S=35, t=25, p=max(0,P-1.013253), P=1.013253,
          method=c("Millero", "Chen", "marine", "limnological")) {
   method <- match.arg(method)
   P   <- p      # hydrostatic pressure; make surface = 0 bar
   rho <- switch(method,
     Millero      = rho(S=S, T=t, P=p),               # from package seacarb
     marine       = rho(S=S, T=t, P=p),               # from package seacarb
     Chen         = 1000 *rhoH2O_Chen(S=S, t=t, p=p), # limnological range
     limnological = 1000 *rhoH2O_Chen(S=S, t=t, p=p), # limnological range
   )
   attr(rho, "unit") = "(kg/m3)"
   return(rho)
}
