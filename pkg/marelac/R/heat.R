
##########################################################################
#  Heat Capacity of sea water
##########################################################################

sw_cp <- function (S=35 ,# -
                t=25, # dgC
                p=P-1.013253,
                P=1.013253,
                UNESCO=FALSE)  # bar

#=========================================================================
# using UNESCO 1983 polynomial.
#
# units in J kg-1 dgC^-1
# REFERENCES:
#    Fofonoff, P. and Millard, R.C. Jr
#    Unesco 1983. Algorithms for computation of fundamental properties of
#    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
#=========================================================================
# Check: CPSW = 3849.500 J/(kg dg.C) for S = 40, t=40, P=1000
{
if (UNESCO)
{
P     <- p  # hydrostatic pressure is called "P" here...
S3_2  <- S*sqrt(S)

# eqn 26 p.32
# specific heat for P=0

Cpst0 =  4217.4    - 3.720283 *t + 0.1412855*t^2 -2.654387e-3*t^3 + 2.093236e-5*t^4 +
        (-7.64357  + 0.1072763*t -1.38385e-3*t^2)*S +
      	(0.1770383 -4.07718e-3*t +  5.148e-5*t^2)*S*sqrt(S)

# eqn 28 p.33
# pressure and temperature terms for S=0

del_Cp0t0 =  (-4.9592e-1 + 1.45747e-2*t - 3.13885e-4*t^2 + 2.0357e-6*t^3 + 1.7168e-8*t^4)*P    +
	           ( 2.4931e-4 - 1.08645e-5*t + 2.87533e-7*t^2 - 4.0027e-9*t^3 + 2.2956e-11*t^4)*P^2 +
             (- 5.422e-8 +  2.6380e-9*t - 6.5637e-11*t^2 + 6.136e-13*t^3)*P^3

# eqn 29 p.34
# pressure and temperature terms for S>0


del_Cpstp <-((4.9247e-3  -1.28315e-4*t + 9.802e-7*t^2 + 2.5941e-8*t^3 -2.9179e-10*t^4)*S +
             (-1.2331e-4 -  1.517e-6*t + 3.122e-8*t^2)*S3_2)*P                           +
	          ((-2.9558e-6 +1.17054e-7*t -2.3905e-9*t^2 + 1.8448e-11*t^3)*S                +
	            9.971e-8*S3_2)*P^2                                                         +
            ((5.540e-10 - 1.7682e-11*t + 3.513e-13*t^2)*S                                +
	                       -1.4300e-12*t*S3_2)*P^3

# specific heat
cp <- Cpst0 + del_Cp0t0 + del_Cpstp
} else {

cp = -(t+273.15)*sw_gibbs(S,t,p,dS=0,dt=2,dp=0)
}
return(cp)
}

##########################################################################
# Compute adiabatic temperature gradient
##########################################################################

 sw_adtgrad <- function (S=35, t=25, p=P-1.013253, P=1.013253)
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
    P <- p*10 #hydrostatic pressure, in dbar
    val<-  3.5803e-5 + (8.5258e-6 + (-6.836e-8 + 6.6228e-10*t)*t)*t  +
          (1.8932e-6   -4.2393e-8*t)*(S-35)                          +
	    ( (1.8741e-8 + (-6.7795e-10 + (8.733e-12 - 5.4481e-14*t)*t)*t) +
          (-1.1351e-10 +2.7759e-12*t)*(S-35))*P                      +
         (-4.6206e-13 + (1.8676e-14 -2.1687e-16*t)*t )*P*P
    return(val*10)        # K/dbar->K/bar
    }
    
# ---------------------------------------------------------------
##########################################################################
# Potential temperature
##########################################################################

sw_tpot<- function (S=35, t=25, p, pref=0)
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
    P    <- p
    Pref <- max(0,pref)
    H  <- 10*(Pref-P)
    XK <- H*sw_adtgrad(t=t,S=S,p=P)/10

    t  <- t + 0.5*XK
    Q  <- XK
    P  <- P + 0.05*H
    XK <- H*sw_adtgrad(t=t,S=S,p=P)/10

    t  <- t + 0.29289322*(XK-Q)
    Q  <- 0.58578644*XK + 0.121320344*Q
    XK <- H*sw_adtgrad(t=t,S=S,p=P)/10

    t  <- t + 1.707106781*(XK-Q)
    Q  <- 3.414213562*XK - 4.121320344*Q
    P  <- P + 0.05*H
    XK <- H*sw_adtgrad(t=t,S=S,p=P)/10

    return (t + (XK-2.0*Q)/6.0)
 }
 
##########################################################################
# Compute freezing temperature of sea water
##########################################################################

sw_tfreeze <- function(S=35, p=P-1.013253, P=1.013253)
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
    return(-0.0575*S + 1.710523e-3*S^1.5 -2.154996e-4*S^2 -7.53e-3*p)
}
