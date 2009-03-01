## -----------------------------------------------------------------------------
## Saturated Concentration of Gasses in Seawater
## -----------------------------------------------------------------------------

gas_satconc <- function(S = 35, t = 25, P = 1.013253,
          x = c("He", "Ne", "N2", "O2", "Ar", "Kr", "CH4", "CO2","N2O"),
          atm = atmComp(x)) {

  if (is.null(atm)) atm <- atmComp(x)
  Vapor <- vapor(t = t, S = S)
  if (length(Vapor) ==1)
    gas_solubility(S = S, t = t, x = x) * P * atm * (1 - Vapor)
  else {
    gs <- gas_solubility(S = S, t = t, x = x)
    gs * P * matrix(nc=ncol(gs),nr=nrow(gs),data=atm,byrow=TRUE) *
        matrix(nc=ncol(gs),nr=nrow(gs),data=(1 - Vapor))
  }
}
