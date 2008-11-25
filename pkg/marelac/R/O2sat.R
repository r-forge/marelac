O2sat <- function(S=0, T = K - 273.15, K = T+273.15, NN = 0, method=c("Weiss", "APHA", "Paul")) {
   log10 <- function(x) log(x)/log(10)
   method <- match.arg(method)
   if ((S !=0) & (method != "Weiss")) warning("Salinity value ignored by this method!")
   if ((NN !=0) & (method != "Paul")) warning("Sea level height ignored by this method!")
   ret <- switch(method,
     # American Public Health Association
     APHA  = exp(-139.34411 + (157570.1/K) - (66423080/K^2) + (12438000000/K^3)-
             (862194900000/K^4)),
     # Weiss, R. (1970). "The solubility of nitrogen, oxygen, and argon 
     # in water and seawater". Deep-Sea Res. 17: 721-35.
     Weiss = 1.426 * exp(-173.4292 + 249.6339 * 100 / K +
             143.3483 * log(K / 100) - 21.8492 * K / 100 +
             S * (-0.033096 + 0.014259 * K / 100 + - 0.001700 * (K / 100)^2)),
     # Paul, L. approximation that respects height above see level
     Paul = (1012-0.12 * NN)/1013 * (14.674 - 13.644 * log10(1 + K/12.8))
   )
   attr(ret, "unit") = "(g/m3)"
   ret
}

