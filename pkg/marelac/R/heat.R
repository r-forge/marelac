
##########################################################################
#  Heat Capacity of sea water
##########################################################################

cp <- function (T=25, # dgC
                S=35 ,# -
                P=1.013253,
                hydroP=max(0,P-1.013253))  # bar

#=========================================================================
# using UNESCO 1983 polynomial.
#
# units in J kg-1 dgC^-1
# REFERENCES:
#    Fofonoff, P. and Millard, R.C. Jr
#    Unesco 1983. Algorithms for computation of fundamental properties of
#    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
#=========================================================================
# Check: CPSW = 3849.500 J/(kg dg.C) for S = 40, T=40, P=1000
{
P     <- hydroP  # hydrostatic pressure is called "P" here...
S3_2  <- S*sqrt(S)

# eqn 26 p.32
# specific heat for P=0

Cpst0 =  4217.4    - 3.720283 *T + 0.1412855*T^2 -2.654387e-3*T^3 + 2.093236e-5*T^4 +
        (-7.64357  + 0.1072763*T -1.38385e-3*T^2)*S +
      	(0.1770383 -4.07718e-3*T +  5.148e-5*T^2)*S*sqrt(S)

# eqn 28 p.33
# pressure and temperature terms for S=0

del_Cp0t0 =  (-4.9592e-1 + 1.45747e-2*T - 3.13885e-4*T^2 + 2.0357e-6*T^3 + 1.7168e-8*T^4)*P    +
	           ( 2.4931e-4 - 1.08645e-5*T + 2.87533e-7*T^2 - 4.0027e-9*T^3 + 2.2956e-11*T^4)*P^2 +
             (- 5.422e-8 +  2.6380e-9*T - 6.5637e-11*T^2 + 6.136e-13*T^3)*P^3

# eqn 29 p.34
# pressure and temperature terms for S>0


del_Cpstp <-((4.9247e-3  -1.28315e-4*T + 9.802e-7*T^2 + 2.5941e-8*T^3 -2.9179e-10*T^4)*S +
             (-1.2331e-4 -  1.517e-6*T + 3.122e-8*T^2)*S3_2)*P                           +
	          ((-2.9558e-6 +1.17054e-7*T -2.3905e-9*T^2 + 1.8448e-11*T^3)*S                +
	            9.971e-8*S3_2)*P^2                                                         +
            ((5.540e-10 - 1.7682e-11*T + 3.513e-13*T^2)*S                                +
	                       -1.4300e-12*T*S3_2)*P^3

# specific heat
cp <- Cpst0 + del_Cp0t0 + del_Cpstp

return(cp)
}

##########################################################################
# Compute adiabatic temperature gradient
##########################################################################

 adtgrad <- function (T=25, S=35, P=1.013253, hydroP=max(0,P-1.013253))
#=========================================================================
# using UNESCO 1983 polynomial.
#
# units in K/bar (!)
# REFERENCES:
#    Fofonoff, P. and Millard, R.C. Jr
#    Unesco 1983. Algorithms for computation of fundamental properties of
#    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
#=========================================================================
  {
    P <- hydroP*10 #hydrostatic pressure, in dbar
    val<-  3.5803e-5 + (8.5258e-6 + (-6.836e-8 + 6.6228e-10*T)*T)*T  +
          (1.8932e-6   -4.2393e-8*T)*(S-35)                          +
	    ( (1.8741e-8 + (-6.7795e-10 + (8.733e-12 - 5.4481e-14*T)*T)*T) +
          (-1.1351e-10 +2.7759e-12*T)*(S-35))*P                      +
         (-4.6206e-13 + (1.8676e-14 -2.1687e-16*T)*T )*P*P
    return(val*10)        # K/dbar->K/bar
    }
    
# ---------------------------------------------------------------
##########################################################################
# Potential temperature
##########################################################################

 temppot<- function (T=25, S=35, hydroP, hydroPref=0)
#=========================================================================
# using UNESCO 1983 polynomial.
#
# units in K
# REFERENCES:
#    Fofonoff, P. and Millard, R.C. Jr
#    Unesco 1983. Algorithms for computation of fundamental properties of
#    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
#=========================================================================
 {
    P    <- hydroP
    Pref <- max(0,hydroPref)
    H  <- 10*(Pref-P)
    XK <- H*adtgrad(T=T,S=S,hydroP=P)/10

    T  <- T + 0.5*XK
    Q  <- XK
    P  <- P + 0.05*H
    XK <- H*adtgrad(T=T,S=S,hydroP=P)/10

    T  <- T + 0.29289322*(XK-Q)
    Q  <- 0.58578644*XK + 0.121320344*Q
    XK <- H*adtgrad(T=T,S=S,hydroP=P)/10

    T  <- T + 1.707106781*(XK-Q)
    Q  <- 3.414213562*XK - 4.121320344*Q
    P  <- P + 0.05*H
    XK <- H*adtgrad(T=T,S=S,hydroP=P)/10

    return (T + (XK-2.0*Q)/6.0)
 }
 
##########################################################################
# Compute freezing temperature of sea water
##########################################################################

tempfreeze <- function(S=35, P=1.013253, hydroP=max(0,P-1.013253))
#=========================================================================
# using UNESCO 1983 polynomial.
#
# units in dg
# REFERENCES:
#    Fofonoff, P. and Millard, R.C. Jr
#    Unesco 1983. Algorithms for computation of fundamental properties of
#    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
#=========================================================================

{
    return(-0.0575*S + 1.710523e-3*S^1.5 -2.154996e-4*S^2 -7.53e-3*hydroP)
}
