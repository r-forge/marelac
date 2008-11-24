rhoH2O_Chen <- function(S=0, T=25, P=0) {
  # Chen, Ch.-T. and F.J. Millero (1986) - Precise thermodynamic
  # properties of natural waters covering only the limnological range.
  # Limnol. Oceanogr. 31 No. 3, 657 - 662
  # rho in g/cm^3
  # P in bar
  # T in °C
  rho0 <- 0.9998395 + 0.000067914 * T -
          0.0000090894 * T * T + 0.00000010171 * T * T * T -
          0.0000000012846 * T ^ 4 + 0.000000000011592 * T ^ 5 -
          5.0125E-14 * T ^ 6 + (0.0008181 - 0.00000385 * T + 0.0000000496 * T * T) * S
  K    <- 19652.17 + 148.113 * T * -2.293 * T * T + 0.01256 * T ^ 3 - 0.0000418 * T ^ 4 +
       (3.2726 - 0.0002147 * T + 0.0001128 * T ^ 2) * P + (53.238 - 0.313 * T + 0.005728 * P) * S
  rho0 * (1 - P / K) ^ -1
}


rhoH2O <- function(S=35, T=25, P=0, method=c("Millero", "Chen", "marine", "limnological")) {
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

