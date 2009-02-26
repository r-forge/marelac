
##########################################################################
# mol to liter and liter to mol conversion for a gas
##########################################################################
molvol <- function(t=25,     # temperature, dg celsius
                  P=1.013253,  # pressure, bar
                  x=c("ideal","Ar", "CO2", "CS2", "CO", "CCl4", "Cl2", "C2H6S", "C2H5OH",
      "C6H5F", "CH3F", "CH4", "CH3OH", "C5H12", "C3H8", "H2O", "He",
      "H2", "HBr", "HCl", "H2S", "Hg", "Kr", "NH3", "Ne",
      "NO", "N2", "NO2", "N2O", "O2", "PH3", "SiH4", "SiF4", "SO2", "Xe"),
                  q=1,      # quantity (mol) of the gas
                  a=0,      # dm^6*bar/mol^2
                  b=0)      # dm3/mol
# the van der Waals coefficients...
{

  Names <- eval(formals(sys.function(sys.parent()))$x)

  if (! is.null(x))
  {
   x <- match.arg(x, several.ok = TRUE) # check if valid input...
   ii <- pmatch(x,Names) # position of x

    ## simple indexed vectors are 5 times faster than row name indexing
    ## in data frame and even slightly faster than a matrix stored in /data
    # a in L2bar/mol2
    Waals.a <- c(0,1.363, 3.640, 11.77, 1.505, 19.7483, 6.579, 13.04, 12.18,
      20.19, 4.692, 2.283, 9.649, 19.26, 8.779, 5.536, 0.03457,
      0.2476, 4.510, 3.716, 4.490, 8.200, 2.349, 4.225, 0.2135,
      1.358, 1.408, 5.354, 3.832, 1.378, 4.692, 4.377, 4.251, 6.803, 4.250)

    Waals.b <- c(0, 0.03219, 0.04267, 0.07685, 0.03985, 0.1281,  0.05622, 0.09213,
      0.08407, 0.1286,  0.05264, 0.04278, 0.06702, 0.146,   0.08445, 0.03049,
      0.0237, 0.02661, 0.04431, 0.04081, 0.04287, 0.01696, 0.03978, 0.03707,
      0.01709, 0.02789, 0.03913, 0.04424, 0.04415, 0.03183, 0.05156, 0.05786,
      0.05571, 0.05636, 0.05105)

    aa <- Waals.a[ii]
    bb <- Waals.b[ii]
  } else
  {aa<-a
  bb<-b}

  #  t in degrees Kelvin
  R  <- 0.082058
  TK <- 273.15 + t

  VV <- NULL
  what <-  NULL
  for (i in 1:length(aa))
  {
  a <- aa[i]
  b <- bb[i]
  
  V <- NULL
  for (TT in TK) {
    for (PP in P) {
      for (xx in q) {
       if (a==0 & b==0) {
          V    <- c(V, q * R * TT / PP * 1.013253)            # CHECK FACTOR 1.013253
          } else {
           V <- c(V,
            uniroot(function (V) ((PP* 1.013253 + xx * xx * a/(V^2)) *
              (V/xx - b) - R * TT), c(-10, 1e6))$root)
          }
        }
      }
    }
  VV <- cbind(VV,V)
  }
  colnames(VV) <-x
  if(nrow(VV) ==1) {VV<- as.vector(VV);names(VV)<-x}
  return(VV)  # volume, liter
}

