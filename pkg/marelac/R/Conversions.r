##########################################################################
# molar volume of ideal gas
##########################################################################

mol.vol <- function(T=25,     # temperature, dg celsius
                    P=1)      # pressure, atmosphere

{
#  T in degrees Kelvin
T = 273.15 + T
R = 0.082058  # liter-atm / mole-K

return(R* T / P)  # molar volume of an ideal gas
}

##########################################################################
# mol to liter and liter to mol conversion for a gas
##########################################################################

mol2l <- function(x=1,      # mol of the gas               
                  T=25,     # temperature, dg celsius
                  P=1,      # pressure, atmosphere
                  gas=NULL,
                  a=0,      # dm^6*bar/mol^2
                  b=0)      # dm3/mol
#     a         b
# He 0.034598   0.023733
# H2 0.24646    0.026665
# N2 1.3661     0.038577
# O2 1.3820     0.03186

{
if (! is.null(gas))
{
 if (gas == "He"){a=0.034598; b=0.023733} else
 if (gas == "H2"){a=0.24646 ; b=0.026665} else
 if (gas == "N2"){a=1.3661  ; b=0.038577} else
 if (gas == "O2"){a=1.3820  ; b=0.03186} else stop(paste("do not have a and b values for the gas",gas))
}

#  T in degrees Kelvin
T = 273.15 + T
R = 0.082058

#V=xRT/P ;    R=0.0821 liter-atm / mole-K
if (a==0 & b==0) V=  x*R* T / P else
{
V <- NULL
for (TT in T) {for (PP in P) {for (xx in x)
 V = c(V,uniroot(fun<- function (V)((PP+xx*xx*a/(V^2))*(V/xx-b)-R*TT),c(-10,1e6))$root)
 }}
}
return(V)  # volume, liter
}

# The reverse: liter to mol
l2mol <- function(x=1, # litre of the gas             
                  T=25,# temperature, dg celsius    
                  P=1, # pressure, atmosphere       
                  gas=NULL,
                  a=0, # dm^6*bar/mol^2
                  b=0) # dm3/mol                    
x/mol2l(1,T=T,P=P,a=a,b=b,gas=gas)     # mole of gas



