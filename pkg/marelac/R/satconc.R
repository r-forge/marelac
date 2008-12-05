satconc <- function(S=35,T=25,K=273.15+T,P=1,gas="O2",atm=atmComp(gas))
{
if (is.null(atm)) atm <- atmComp(gas)
solubility(S=S,K=K,P=P,gas=gas)*P* atm * (1-vapor(T=T,S=S))
}
