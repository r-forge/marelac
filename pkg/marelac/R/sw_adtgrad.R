
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
    
