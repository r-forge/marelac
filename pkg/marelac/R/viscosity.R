######################################################################
# Calculates water viscosity
# From the FORTRAN implementation of B. Boudreau
# Boudreau BP
# A method-of-lines code for carbon and nutrient diagenesis in aquatic sediments 
# COMPUTERS & GEOSCIENCES 22 (5): 479-496 JUN 1996 
######################################################################


viscosity <- function (S=35,     # Salinity, ppt
                       t=25,     # Temperture, degrees C
                       P=1.013253)      # Pressure, bar
#---------------------------------------------------------------------
#  VISCO      Calculates the shear viscosity of water using the equation
#             given by Kukulka et al. (1987).
#             Calculated viscosity is in centipoise.
#
#             Valid for 0<t<30 and 0<S<36.
#---------------------------------------------------------------------

      1.7910 - t*(6.144e-02 - t*(1.4510e-03 - t*1.6826e-05))            +
      - 1.5290e-04*P + 8.3885e-08*P*P + 2.4727e-03*S                    +
      + (6.0574e-06*P - 2.6760e-09*P*P)*t + (t*(4.8429e-05              +
      - t*(4.7172e-06 - t*7.5986e-08)))*S

    # end viscosity

