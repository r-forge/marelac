
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
