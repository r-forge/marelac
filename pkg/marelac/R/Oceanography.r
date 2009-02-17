
##########################################################################
# Coriolis factor as a function of latitude
##########################################################################

coriolis <- function (lat)  # latitude in degrees north (-90:+90)
  2* 7.2921e-5 *sin(pi/180*lat)

#=========================================================================
# Coriolis force: 2*omega*sin(lat), omega=7.292e-5 radians/sed, lat
# coriolis factor is in per second
# reference
#   S. Pond & G.Pickard  2nd Edition 1986
#   Introductory Dynamical Oceanogrpahy
#   Pergamon Press Sydney.  ISBN 0-08-028728-X
#
# omega from:
#
#   A.E. Gill 1982. p.597
#   "Atmosphere-Ocean Dynamics"
#   Academic Press: New York.  ISBN: 0-12-283522-0
#=========================================================================

##########################################################################
# velocity of sound
##########################################################################

sw_svel <- function (S=35, t=25, p=P-1.013253, P=1.013253, UNESCO=FALSE)
#=========================================================================
# using UNESCO 1983 polynomial.
#
# units in m/sec
# REFERENCES:
#    Fofonoff, P. and Millard, R.C. Jr
#    Unesco 1983. Algorithms for computation of fundamental properties of
#    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
#=========================================================================
{

if (UNESCO) {
    P = p  # P is used in code as synonym for hydroP
    T = t
    P2 <- P*P
    P3 <- P2*P
    Cw <- 1402.388 + (5.03711   + (-5.80852e-2 + (3.3420e-4 + (-1.47800e-6 + 3.1464e-9*T)*T)*T)*T)*T   +
       + (0.153563 + (6.8982e-4 + (-8.1788e-6  + (1.3621e-7 -6.1185e-10*T)*T)*T)*T)*P          +
       + (3.1260e-5 +(-1.7107e-6+ (2.5974e-8 + (-2.5335e-10 + 1.0405e-12*T)*T)*T)*T)*P2        +
       + (-9.7729e-9+(3.8504e-10 -2.3643e-12*T)*T)*P3

    A <- 1.389      + (-1.262e-2  + (7.164e-5    + (2.006e-6  -3.21e-8*T)*T)*T)*T        +
      + (9.4742e-5  + (-1.2580e-5 + (-6.4885e-8  + (1.0507e-8 -2.0122e-10*T)*T)*T)*T)*P  +
      + (-3.9064e-7 + (9.1041e-9  + (-1.6002e-10 + 7.988e-12*T)*T)*T)*P2                 +
      + (1.100e-10  + (6.649e-12  -3.389e-13*T)*T)*P3

    B <- -1.922e-2 -4.42e-5*T + (7.3637e-5 + 1.7945e-7*T)*P

    D <- 1.727e-3 + -7.9836e-6*P
    svel <- Cw + A*S + B*S**1.5 + D*S**2
} else {

n0 = 0; n1 = 1; n2 = 2;

g_tt <- sw_gibbs(S=S,t=t,p=p,dS=0,dt=2,dp=0)
g_2  <- sw_gibbs(S=S,t=t,p=p,dS=0,dt=0,dp=1)
g_tp <- sw_gibbs(S=S,t=t,p=p,dS=0,dt=1,dp=1)
svel <- sw_gibbs(S=S,t=t,p=p,dS=0,dt=0,dp=1) *
             sqrt(g_tt/(g_tp*g_tp - g_tt*sw_gibbs(S=S,t=t,p=p,dS=0,dt=0,dp=2)))

}
    return (svel)
}

# ----------------------------------------------------------------
gravity <- function (lat=0, UNESCO = FALSE)
{
# Compute gravity from latitude

    X <- sin(lat*pi/180.)
    X <- X*X
    X2 <- (sin(2*lat*pi/180.))^2
   if (UNESCO)
      grav <-  9.780318 * (1.0 + (5.2788e-3 + 2.36e-5*X)*X)
   else grav <-  9.780327*(1.+0.0053024*X -0.0000058*X2)

   return (grav)
}
# ----------------------------------------------------------------
sw_depth <- function (p=P-1.013253, P=1.013253, lat=0)
{
# Compute depth from hydrostatic pressure and latitude

    P <- p*10    # P=hydrostatic pressure, in dbar
    denom = gravity(lat)+ 1.092e-6*P
    nom = (9.72659 + (-2.2512e-5 + (2.279e-10 -1.82e-15*P)*P)*P)*P

    return (nom / denom)
}




# thermal expansion coefficient
sw_alpha <- function ( S=35, t=25, p=P-1.013253, P=1.013253)
    sw_gibbs(S,t,p,dS=0,dt=1,dp=1)/sw_gibbs(S,t,p,dS=0,dt=0,dp=1)

# haline contraction coefficient
sw_beta <- function ( S=35, t=25, p=P-1.013253, P=1.013253)
   -sw_gibbs(S,t,p,dS=1,dt=0,dp=1)/sw_gibbs(S,t,p,dS=0,dt=0,dp=1)

# specific enthalpy
sw_enthalpy <- function ( S=35, t=25, p=P-1.013253, P=1.013253)
  sw_gibbs(S,t,p,dS=0,dt=0,dp=0)-(t+273.15)*sw_gibbs(S,t,p,dS=0,dt=1,dp=0)

# specific enthropy
sw_entropy <- function ( S=35, t=25, p=P-1.013253, P=1.013253)
  -sw_gibbs(S,t,p,dS=0,dt=1,dp=0)
  
# isentropic compressibility
sw_kappa <- function ( S=35, t=25, p=P-1.013253, P=1.013253)
{
n0 = 0; n1 = 1; n2 = 2;

g_tt = sw_gibbs(S,t,p,dS=0,dt=2,dp=0)
g_tp = sw_gibbs(S,t,p,dS=0,dt=1,dp=1)

1e5 * (g_tp * g_tp - g_tt *sw_gibbs(S,t,p,dS=0,dt=0,dp=2))/
   (sw_gibbs(S,t,p,dS=0,dt=0,dp=1)*g_tt)
}


# isothermal compressibility
sw_kappa_t <- function ( S=35, t=25, p=P-1.013253, P=1.013253)
  -1e4 * sw_gibbs(S,t,p,dS=0,dt=0,dp=2)/sw_gibbs(S,t,p,dS=0,dt=0,dp=1)
