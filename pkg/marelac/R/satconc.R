satconc <- function(S=35,T=25,K=273.15+T,P=1,gas="O2",atm=AtmComp[[gas]])
{
if (is.null(atm)) atm <- AtmComp[[gas]]
solubility(S=S,K=K,P=P,gas=gas)*P* atm
}
