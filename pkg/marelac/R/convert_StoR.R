# -------------------------------------------------
#Compute conductivity ratio from salinity, temperature, and pressure

convert_StoR <- function(S=35, t=25, p=max(0,P-1.013253), P=1.013253)
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

