##########################################################################
# Coriolis factor as a function of latitude
##########################################################################

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

coriolis <- function (lat)  # latitude in degrees north (-90:+90)
  2* 7.2921e-5 *sin(pi/180*lat)

