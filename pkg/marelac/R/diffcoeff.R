## -----------------------------------------------------------------------------
## Molecular diffusion coefficients
## -----------------------------------------------------------------------------

diffcoeff <- function(S=35,      # Salinity, ppt
                      t=25,      # Temperture, degrees C
                      P=1.013253,       # Pressure, atm
                      species=c("O2","CO2","NH3","H2S","CH4","HCO3","CO3","NH4",
      "HS","NO3","H2PO4","HPO4","PO4","H","OH","Ca","Mg","Fe","Mn",
      "SO4","H3PO4","BOH3","B0H4","H4SiO4")
                      ) {

  TK  <- t + 273.15   # Temperature, Kelvin
  species   <- match.arg(species, several.ok = TRUE)

  #  The viscosity for the true sample conditions.
  H2OViscosity <- viscosity(S,t,P)

  #  Diffusion coefficients in pure water at sample temperature.
  V0 <- viscosity(S=0,t=t,P=1)

  #  To correct for salinity, the Stokes-Einstein relationship is used.
  #  This is not quite accurate, but is at least consistent.
  #  Also corrected from cm2/sec to cm2/hour
  fac <- V0/H2OViscosity*3600.

  #  Diffusion coefficient of Water
  A <- 12.5e-09*exp(-5.22e-04*P)
  B <- 925.0*exp(-2.6e-04*P)
  T0 <- 95.0 + 2.61e-02*P
  D_H2O <- A*sqrt(TK)*exp(-B/(TK-T0))*1.0e+04 *fac

  #  Dissolved gases : from Wilke and Chang (1955)
  #  D = (A/VM^0.6)*7.4e-8*fac
  #    note: 1) mw = molecular weight of water
  #          2) VM = molar volumes (cm**3/mol) (Sherwood et al., 1975)
  #  The factor phi is reduced to 2.26 as suggested by Hayduk and Laudie (1974).

  phi <- 2.26
  mw  <- 18.0
  A   <- sqrt(phi*mw)*TK/V0
  fac2 <- 7.4e-08*fac

  D_O2  <- (A/(25.6^0.6)) *fac2        #  Oxygen
  D_CO2 <- (A/(34.0^0.6)) *fac2        #  CO2
  D_NH3 <- (A/(25.8^0.6)) *fac2        #  NH3
  D_H2S <- (A/(32.9^0.6)) *fac2        #  H2S
  D_CH4 <- (A/(37.7^0.6)) *fac2        #  CH4

  #  The coefficients in pure water for the following species are
  #  calculated by linear functions of temperature (deg C)
  #  coefficients as in Boudreau (1996).

  #  i.e. NO3-,HS-,H2PO4-,CO3=,SO4=,Ca++,Mg++,Mn++,Fe++,NH4+,H+ & OH-,
  #         HCO3-, HPO4=, PO4(3-)

  fac2    <- 1.0e-06*fac
  D_HCO3  <- (5.06 + 0.275*t)*fac2
  D_CO3   <- (4.33 + 0.199*t)*fac2
  D_NH4   <- (9.5  + 0.413*t)*fac2
  D_HS    <- (10.4 + 0.273*t)*fac2
  D_NO3   <- (9.50 + 0.388*t)*fac2
  D_H2PO4 <- (4.02 + 0.223*t)*fac2
  D_HPO4  <- (3.26 + 0.177*t)*fac2
  D_PO4   <- (2.62 + 0.143*t)*fac2
  D_H     <- (54.4 + 1.555*t)*fac2
  D_OH    <- (25.9 + 1.094*t)*fac2
  D_Ca    <- (3.60 + 0.179*t)*fac2
  D_Mg    <- (3.43 + 0.144*t)*fac2
  D_Fe    <- (3.31 + 0.150*t)*fac2
  D_Mn    <- (3.18 + 0.155*t)*fac2
  D_SO4   <- (4.88 + 0.232*t)*fac2

  # H3PO4 : Least (1984) determined D(H3PO4) at 25 deg C and 0 S.
  #         Assume that this value can be scaled by the Stokes-Einstein
  #         relationship to any other temperature.

  D_H3PO4 <- 0.87e-05
  TS      <- 25.0
  SS      <- 0.0
  V25     <- viscosity(SS,TS,1)
  VTK     <- viscosity(SS,t,1)
  D_H3PO4 <- D_H3PO4*V25/298.15*TK/VTK*fac

  #  B(OH)3 : Mackin (1986) determined D(B(OH)3) at 25 deg C and
  #           about 29.2 S.
  #           Assume that this value can be scaled by the Stokes-Einstein
  #           relationship to any other temperature.

  D_BOH3 <- 1.12e-05
  TS     <- 25.0
  SS     <- 29.2
  V25   <- viscosity(SS,TS,1)
  VTK   <- viscosity(0,t,1)

  D_BOH3 = D_BOH3*V25/298.15*TK/VTK*fac

  #  B(OH)4 : No information on this species ! Boudreau and
  #           Canfield (1988) assume it is 12.5% smaller than B(OH)3.
  #
  D_B0H4 <- 0.875*D_BOH3

  #  H4SiO4 : Wollast and Garrels (1971) found D(H4SiO4) at 25 deg C
  #           and 36.1 ppt S.
  #           Assume that this value can be scaled by the Stokes-Einstein
  #           relationship to any other temperature.
  D_H4SiO4 <- 1.0E-05
  TS <- 25.0
  SS <- 36.1
  V25 <- viscosity(SS,TS,1)
  VTK <- viscosity(0,t,1)
  D_H4SiO4 <- D_H4SiO4*V25/298.15*TK/VTK*fac

  fac <- 1/(3600*1e4)       # from cm2/hr to m2/s
  diffc <-data.frame(       #
     O2  = D_O2*fac,         # Oxygen
     CO2 = D_CO2*fac,
     NH3 = D_NH3*fac,
     H2S = D_H2S*fac,
     CH4 = D_CH4*fac,
     HCO3  = D_HCO3*fac,
     CO3   = D_CO3*fac,
     NH4   = D_NH4*fac,
     HS    = D_HS*fac,
     NO3   = D_NO3*fac,
     H2PO4 = D_H2PO4*fac,
     HPO4  = D_HPO4*fac,
     PO4   = D_PO4*fac,
     H     = D_H*fac,
     OH    = D_OH*fac,
     Ca    = D_Ca*fac,
     Mg    = D_Mg*fac,
     Fe    = D_Fe*fac,
     Mn    = D_Mn*fac,
     SO4   = D_SO4*fac,
     H3PO4 = D_H3PO4*fac,
     BOH3 = D_BOH3*fac,
     B0H4 = D_B0H4*fac,
     H4SiO4 = D_H4SiO4*fac)

  return(diffc[species])
}
