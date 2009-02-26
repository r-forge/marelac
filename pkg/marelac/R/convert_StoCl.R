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
