
sw_comp <- function(x=c("Na","Mg","Ca","K","Sr","Cl","SO4","HCO3","Br","CO3",
                        "BOH4","F","OH","BOH3","CO2")) {
  x <- match.arg(x, several.ok = TRUE)
  ret <- c(Na= 0.3065958,Mg= 0.0365055,Ca=0.0117186,K=0.0113495,
Sr=0.0002260,Cl=0.5503396,SO4=0.0771319,HCO3=0.0029805,
Br=0.0019134,CO3=0.0004078,BOH4=0.0002259,F=0.0000369,
OH=0.0000038,BOH3=0.0005527,CO2=0.0000121 )
  ret[x]
}



