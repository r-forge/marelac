# --------------------------------------------------
# Salt, chlorinity and conductivity functions
# -------------------------------------------------------

# Compute salinity from conductivity, temperature, and pressure

convert_RtoS <- function (R=1, t=25, P=1.013253, p=max(0,P-1.013253))
{
#=========================================================================
# using UNESCO 1983 polynomial.
#
# units in K/bar (!)
# REFERENCES:
#    Fofonoff, P. and Millard, R.C. Jr
#    Unesco 1983. Algorithms for computation of fundamental properties of
#    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
#=========================================================================

    P   <- p*10    # uses decibar in calculations, hydrostatic pressure
    C_P <- (2.070e-5 + (-6.370e-10 + 3.989e-15*P)*P)*P

    DT  <- t - 15.0
    R_T <- 0.6766097 + (2.00564e-2 + (1.104259e-4 +
                        (-6.9698e-7 + 1.0031e-9*t)*t)*t)*t
    A_T <- 4.215e-1 + -3.107e-3*t
    B_T <- 1.0 + (3.426e-2 + 4.464e-4*t)*t

    RT  <- R/(R_T*(1.0 + C_P/(B_T + A_T*R)))
    RT  <- sqrt(abs(RT))

    DS = (DT / (1+0.0162*DT) ) *
         (0.0005 + (-0.0056 + (-0.0066 + (-0.0375 + (0.0636 + -0.0144*RT)*RT)*RT)*RT)*RT)

    return (0.0080 + (-0.1692 + (25.3851 + (14.0941 +
           (-7.0261 + 2.7081*RT)*RT)*RT)*RT)*RT + DS)
}

# -------------------------------------------------
#Compute conductivity ratio from salinity, temperature, and pressure

convert_StoR <- function(S=35, t=25, P=1.013253, p=max(0,P-1.013253))
#=========================================================================
# using UNESCO 1983 polynomial.
#
# REFERENCES:
#    Fofonoff, P. and Millard, R.C. Jr
#    Unesco 1983. Algorithms for computation of fundamental properties of
#    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
#=========================================================================
{
   fun <- function(x)
     convert_RtoS(x, t=t, p=p) - S
   cond <- uniroot(fun, c(0,5),tol=1e-10)$root
   return(cond)
}


# -------------------------------------------------
# salinity-chlorinity conversion
# -------------------------------------------------

convert_StoCl <- function(S=35) # salinity
          S/1.80655      # chlorinity, in g/kg
# in g/kg
# from Wooster et al.,1969

# -------------------------------------------------
# concentrations of conservative species in seawater
# -------------------------------------------------

# µmol/kg solution
sw_conserv<- function(S=35)
{
 Borate   = 4.16e2*S/35                # Millero 95
 Calcite  = 0.01028e6 *S/35            # Millero 95
 Sulphate = convert_StoCl(S)*0.14e6*1/96.062  # Morris & Riley 1966. Deep-See Res. 13:699-705
 Fluoride = convert_StoCl(S)*67/18.9984       # Riley 1965.Deep-Sea Res. 12:219-220
 return(data.frame(Borate=Borate,
                   Calcite=Calcite,
                   Sulphate=Sulphate,
                   Fluoride=Fluoride))
}
