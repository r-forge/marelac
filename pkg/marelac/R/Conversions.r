##########################################################################
# molar volume of ideal gas
##########################################################################
mol.vol <- function(T = 25,      # temperature, dg celsius
                    P = 1.013253 # pressure, bar
                   )
{
  #  T in degrees Kelvin
  TK <- 273.15 + T
  R  <- 0.082058  # liter-atm / mole-K
  return(R * TK / P*1.013253)  # molar volume of an ideal gas
}

##########################################################################
# mol to liter and liter to mol conversion for a gas
##########################################################################
mol2l <- function(x=1,      # mol of the gas
                  T=25,     # temperature, dg celsius
                  P=1.013253,  # pressure, bar
                  gas=NULL,
                  a=0,      # dm^6*bar/mol^2
                  b=0)      # dm3/mol
# the van der Waals coefficients...
{
  if (! is.null(gas))
  {
    ## simple indexed vectors are 5 times faster than row name indexing
    ## in data frame and even slightly faster than a matrix stored in /data
    Waals.a <- c(1.363, 3.640, 11.77, 1.505, 19.7483, 6.579, 13.04, 12.18,
      20.19, 4.692, 2.283, 9.649, 19.26, 8.779, 5.536, 0.03457,
      0.2476, 4.510, 3.716, 4.490, 8.200, 2.349, 4.225, 0.2135,
      1.358, 1.408, 5.354, 3.832, 1.378, 4.692, 4.377, 4.251, 6.803, 4.250)

    Waals.b <- c( 0.03219, 0.04267, 0.07685, 0.03985, 0.1281,  0.05622, 0.09213,
      0.08407, 0.1286,  0.05264, 0.04278, 0.06702, 0.146,   0.08445, 0.03049,
      0.0237, 0.02661, 0.04431, 0.04081, 0.04287, 0.01696, 0.03978, 0.03707,
      0.01709, 0.02789, 0.03913, 0.04424, 0.04415, 0.03183, 0.05156, 0.05786,
      0.05571, 0.05636, 0.05105)

    Waals.names <- c("Ar", "CO2", "CS2", "CO", "CCl4", "Cl2", "C2H6S", "C2H5OH",
      "C6H5F", "CH3F", "CH4", "CH3OH", "C5H12", "C3H8", "H2O", "He",
      "H2", "HBr", "HCl", "H2S", "Hg", "Kr", "NH3", "Ne",
      "NO", "N2", "NO2", "N2O", "O2", "PH3", "SiH4", "SiF4", "SO2", "Xe")

    a <- Waals.a[gas == Waals.names]
    b <- Waals.b[gas == Waals.names]
    if (is.na(a)) stop(paste("do not have a and b values for the gas", gas))
  }

  #  T in degrees Kelvin
  TK <- 273.15 + T
  R  <- 0.082058

  #V=xRT/P ;    R=0.0821 liter-atm / mole-K
  if (a==0 & b==0) {
    V =  x * R * TK / P * 1.013253
  } else {
    V <- NULL
    for (TT in TK) {
      for (PP in P) {
        for (xx in x) {
          # V = c(V, uniroot(fun<- function (V)((PP+xx*xx*a/(V^2))*(V/xx-b)-R*TT),c(-10,1e6))$root)
          V <- c(V,
            uniroot(function (V) ((PP* 1.013253 + xx * xx * a/(V^2)) *
              (V/xx - b) - R * TT), c(-10, 1e6))$root)
        }
      }
    }
  }
  return(V)  # volume, liter
}

# The reverse: liter to mol
l2mol <- function(x=1, # litre of the gas
                  T=25,# temperature, dg celsius
                  P= 1.013253, # pressure, bar
                  gas=NULL,
                  a=0, # dm^6*bar/mol^2
                  b=0) # dm3/mol
  x/mol2l(1, T=T, P=P, a=a, b=b, gas=gas)     # mole of gas



