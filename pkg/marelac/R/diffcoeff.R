## -----------------------------------------------------------------------------
## Molecular diffusion coefficients
## -----------------------------------------------------------------------------

diffcoeff <- function(S=35,        # Salinity, ppt
                      t=25,        # Temperature, degrees C
                      P=1.013253,  # Pressure, bar
                      species=c("H2O","O2","CO2","H2","CH4","DMS",
      "He","Ne","Ar","Kr","Xe","Rn",
      "N2","H2S","NH3","NO","N2O","CO","SO2",
      "OH","F","Cl","Br","I",
      "HCO3","CO3","H2PO4","HPO4","PO4",
      "HS","HSO3","SO3","HSO4","SO4","IO3","NO2","NO3",
      "H","Li","Na","K","Cs","Ag","NH4",
      "Ca","Mg","Fe","Mn",
      "Ba","Be","Cd","Co","Cu","Hg","Ni","Sr","Pb","Ra","Zn","Al","Ce","La","Pu",
      "H3PO4","BOH3","BOH4","H4SiO4")
      ) {
  
   
  species   <- match.arg(species, several.ok = TRUE)

  TK  <- t + 273.15   # Temperature, Kelvin
  Patm <- 1.013253

  #  The viscosity in pure water at atmospheric pressure and sample temperature.
  mu_0 <- viscosity(S=0,t=t,P=1.013253)

  #  Diffusion coefficient of Water at S = 0, P = P=1.013253, t = t
  #  Cohen MH and Turnbull D (1959). 
  #  Krynicki K, Green CD and Sawyer DW (1978). 
  
  A <- 12.5e-09*exp(-5.22e-04*P)
  B <- 925.0*exp(-2.6e-04*P)
  T0 <- 95.0 + 2.61e-02*P
  D_H2O <- A*sqrt(TK)*exp(-B/(TK-T0))
  D_H2O <- D_H2O*viscosity(S=0,t,P)/mu_0

  
  #  Diffusion coefficient of O2 and CO2 
  #  Boudreau (1997)
  A <- 0.2604
  B <- 0.006383
  D_O2 <- (A+B*(TK/mu_0))*1E-09

  A <- 0.1954
  B <- 0.005089
  D_CO2 <- (A+B*(TK/mu_0))*1E-09
  
  #  Other dissolved gases
  #  Boudreau (1997)
  #  TK <- 298.15

  Arrhenius <- function (A,Ea) A*exp(-(Ea*1000)/(8.314472*TK))*1.0E-09

  D_H2  <- Arrhenius(A= 3338, Ea =16.06)  # Jahne et al (1987)
  D_CH4 <- Arrhenius(A= 3047, Ea =18.36)  # Jahne et al (1987)
  D_DMS <- Arrhenius(A= 2000, Ea =18.10)  # Saltzman et al. (1993)

  D_He  <- Arrhenius(A=  818, Ea =11.70)  # Jahne et al (1987)
  D_Ne  <- Arrhenius(A= 1608, Ea =14.84)  # Jahne et al (1987)
  D_Ar  <- Arrhenius(A= 7238, Ea =19.81)  # Ohsumi and Horibe (1984)
  D_Kr  <- Arrhenius(A= 6393, Ea =20.20)  # Jahne et al (1987)
  D_Xe  <- Arrhenius(A= 9007, Ea =21.61)  # Jahne et al (1987)
  D_Rn  <- Arrhenius(A=15877, Ea =23.26)  # Jahne et al (1987)
 

  #  Other dissolved gases
  #  from Wilke and Chang (1955) as modified by Hayduk and Laudie (1974)

  WilkeChang <- function (Vb) 4.72E-07*TK/(mu_0*Vb^0.6)*1E-04

  D_N2 <- WilkeChang(Vb=34.7)
  D_H2S <- WilkeChang(Vb=35.2)
  D_NH3 <- WilkeChang(Vb=24.5)
  D_NO <- WilkeChang(Vb=23.6)
  D_N2O <- WilkeChang(Vb=36.0)
  D_CO <- WilkeChang(Vb=34.5)
  D_SO2 <- WilkeChang(Vb=43.8)
  
   
  #  The coefficients in pure water for the following species are
  #  calculated by linear functions of temperature (deg C)
  #  coefficients as in Boudreau (1997).

  Boudreau <- function (m0,m1) (m0+m1*t)*1.0e-10

  D_OH    <- Boudreau(25.9,1.094)
  D_Br    <- Boudreau(10.0,0.441)
  D_Cl    <- Boudreau(9.60,0.438)
  D_F     <- Boudreau(6.29,0.343)
  D_I     <- Boudreau(9.81,0.432)

  D_HCO3  <- Boudreau(5.06,0.275)
  D_CO3   <- Boudreau(4.33,0.199)

  D_H2PO4 <- Boudreau(4.02,0.223)
  D_HPO4  <- Boudreau(3.26,0.177)
  D_PO4   <- Boudreau(2.62,0.143)

  D_HS    <- Boudreau(10.4,0.273)
  D_HSO3  <- Boudreau(6.35,0.280)
  D_SO3   <- Boudreau(4.82,0.266)
  D_HSO4  <- Boudreau(5.99,0.307)
  D_SO4   <- Boudreau(4.88,0.232)

  D_IO3   <- Boudreau(4.66,0.252)
  D_NO2   <- Boudreau(10.3,0.331)
  D_NO3   <- Boudreau(9.50,0.388)

  D_H     <- Boudreau(54.4,1.555)
  D_Li    <- Boudreau(4.43,0.241)
  D_Na    <- Boudreau(6.06,0.297)
  D_K     <- Boudreau(6.06,0.297)
  D_Cs    <- Boudreau(10.3,0.416)
  D_Ag    <- Boudreau(7.82,0.359)
  D_NH4   <- Boudreau(9.50,0.413)

  D_Ca    <- Boudreau(3.60,0.179)
  D_Mg    <- Boudreau(3.43,0.144)
  D_Fe    <- Boudreau(3.31,0.150)
  D_Mn    <- Boudreau(3.18,0.155)

  D_Ba    <- Boudreau(4.06,0.176)
  D_Be    <- Boudreau(2.57,0.140)
  D_Cd    <- Boudreau(3.31,0.152)
  D_Co    <- Boudreau(3.31,0.152)
  D_Cu    <- Boudreau(3.39,0.158)
  D_Hg    <- Boudreau(3.63,0.208)
  D_Ni    <- Boudreau(3.36,0.130)
  D_Sr    <- Boudreau(3.69,0.169)
  D_Pb    <- Boudreau(4.46,0.198)
  D_Ra    <- Boudreau(3.91,0.199)
  D_Zn    <- Boudreau(3.31,0.151)

  D_Al    <- Boudreau(2.79,0.172)
  D_Ce    <- Boudreau(2.95,0.131)
  D_La    <- Boudreau(2.78,0.136)
  D_Pu    <- Boudreau(2.71,0.120)

  # H3PO4 : Least (1984) determined D(H3PO4) at 25 deg C and 0 S.
  #         Assume that this value can be scaled by the Stokes-Einstein
  #         relationship to any other temperature.

  D_H3PO4 <- 0.87e-09
  tS      <- 25.0
  SS      <- 0.0
  mu_S   <- viscosity(SS,tS,Patm)
  D_H3PO4 <- D_H3PO4*(mu_S/mu_0)*(TK/(tS + 273.15))

  #  B(OH)3 : Mackin (1986) determined D(B(OH)3) at 25 deg C and
  #           about 29.2 S.
  #           Assume that this value can be scaled by the Stokes-Einstein
  #           relationship to any other temperature.

  D_BOH3 <- 1.12e-09
  tS     <- 25.0
  SS     <- 29.2
  mu_S   <- viscosity(SS,tS,Patm)
  D_BOH3 = D_BOH3*(mu_S/mu_0)*(TK/(tS + 273.15))

  #  B(OH)4 : No information on this species ! Boudreau and
  #           Canfield (1988) assume it is 12.5% smaller than B(OH)3.
  #
  D_BOH4 <- 0.875*D_BOH3

  #  H4SiO4 : Wollast and Garrels (1971) found D(H4SiO4) at 25 deg C
  #           and 36.1 ppt S.
  #           Assume that this value can be scaled by the Stokes-Einstein
  #           relationship to any other temperature.
  D_H4SiO4 <- 1.0E-09
  tS <- 25.0
  SS <- 36.1
  mu_S   <- viscosity(SS,tS,Patm)
  D_H4SiO4 <- D_H4SiO4*(mu_S/mu_0)*(TK/(tS + 273.15))

  # SALINITY AND PRESSURE CORRECTION
  #  To correct for pressure and salinity, the Stokes-Einstein relationship 
  #  is used. This is not quite accurate, but it is at least consistent.
  
  mu <- viscosity(S,t,P)   #  viscosity at sample conditions
  fac <- (mu_0/mu)

  diffc <- data.frame(      
#    "H2O","O2","CO2","H2","CH4","DMS",
     H2O  = D_H2O*fac,        
     O2  = D_O2*fac,        
     CO2 = D_CO2*fac,
     H2  = D_H2*fac,        
     CH4 = D_CH4*fac,
     DMS = D_DMS*fac,        
#    "He","Ne","Ar","Kr","Xe","Rn",
     He  = D_He*fac,        
     Ne  = D_Ne*fac,        
     Ar  = D_Ar*fac,        
     Kr  = D_Kr*fac,        
     Xe  = D_Xe*fac,        
     Rn  = D_Rn*fac,        
#    "N2","H2S","NH3","NO","N2O","CO","SO2",
     N2  = D_N2*fac,        
     H2S = D_H2S*fac,        
     NH3 = D_NH3*fac,        
     NO  = D_NO*fac,        
     N2O = D_N2O*fac,        
     CO  = D_CO*fac,        
     SO2 = D_SO2*fac,        
#    "OH","F","Cl","Br","I",
     OH  = D_OH*fac,        
     F  = D_F*fac,        
     Cl  = D_Cl*fac,        
     Br  = D_Br*fac,        
     I  = D_I*fac,        
#    "HCO3","CO3","H2PO4","HPO4","PO4",
     HCO3  = D_HCO3*fac,
     CO3   = D_CO3*fac,
     H2PO4 = D_H2PO4*fac,
     HPO4  = D_HPO4*fac,
     PO4   = D_PO4*fac,
#    "HS","HSO3","SO3","HSO4","SO4","IO3","NO2","NO3",
     HS    = D_HS*fac,
     HSO3  = D_HSO3*fac,
     SO3   = D_SO3*fac,
     HSO4  = D_HSO4*fac,
     SO4   = D_SO4*fac,
     IO3   = D_IO3*fac,
     NO2   = D_NO2*fac,
     NO3   = D_NO3*fac,
#    "H","Li","Na","K","Cs","Ag","NH4",
     H     = D_H*fac,
     Li    = D_Li*fac,
     Na    = D_Na*fac,
     K     = D_K*fac,
     Cs    = D_Cs*fac,
     Ag    = D_Ag*fac,
     NH4   = D_NH4*fac,
#    "Ca","Mg","Fe","Mn",
     Ca    = D_Ca*fac,
     Mg    = D_Mg*fac,
     Fe    = D_Fe*fac,
     Mn    = D_Mn*fac,
#    "Ba","Be","Cd","Co","Cu","Hg","Ni","Sr","Pb","Ra","Zn","Al","Ce","La","Pu",
     Ba    = D_Ba*fac,
     Be    = D_Be*fac,
     Cd    = D_Cd*fac,
     Co    = D_Co*fac,
     Cu    = D_Cu*fac,
     Hg    = D_Hg*fac,
     Ni    = D_Ni*fac,
     Sr    = D_Sr*fac,
     Pb    = D_Pb*fac,
     Ra    = D_Ra*fac,
     Zn    = D_Zn*fac,
     Al    = D_Al*fac,
     Ce    = D_Ce*fac,
     La    = D_La*fac,
     Pu    = D_Pu*fac,
#    "H3PO4","BOH3","B0H4","H4SiO4")
     H3PO4 = D_H3PO4*fac,
     BOH3 = D_BOH3*fac,
     BOH4 = D_BOH4*fac,
     H4SiO4 = D_H4SiO4*fac)

  return(diffc[species])
}
       