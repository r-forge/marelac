gas_satconc <- function(S = 35, t = 25, P = 1.013253,
x = c("He", "Ne", "N2", "O2", "Ar", "Kr", "CH4", "CO2","N2O", "H2", "Xe", "CO", "O3"), atm = atmComp(x))
{
  if (is.null(atm)) atm <- atmComp(x)
  gas_solubility(S = S, t = t, P = P, x = x) * P * atm * (1 - vapor(t = t, S = S))
}
