## -----------------------------------------------------------------------------
## Molecular Diffusion Coefficients  - more species
## -----------------------------------------------------------------------------

Diffusion <- function(S=35, t=25, P=0, x=c("H2","He","NO","N2O","N2","NH3",
     "O2","CO","CO2","SO2","H2S","Ar","Kr","Ne","CH4","CH3Cl","C2H6","C2H4",
     "C3H8","C3H6","C4H10","Urea","Methanol","Ethanol","Glycine","Alanine","Serine","Valine","Leucine",
     "Proline","Dextrose","Sucrose","Formic acid","Acetic acid",
     "Proprionic acid","Butyric acid","H2O","H3PO4","BOH3","BOH4","H4SiO4",
     "H+","D+","Li+","Na+","K+","Rb+","Cs+","Ag+","NH4+",
     "Ba++","Be++","Ca++","Cd++","Co++","Cu++","Fe++","Hg++","Mg++",
     "Mn++","Ni++","Sr++","Pb++","Ra++","Zn++","Al++","Ce+++","La+++","Pu+++",
     "OH-","OD-","AlOH4-","Br-","Cl-","F-","HCO3-","H2PO4-","HS-",
     "HSO3-","HSO4-","I-","IO3-","NO2-","NO3-","Acetate-","Lactate-",
     "CO3--","HPO4--","SO3--","SO4--","S2O3--","S2O4--","S2O6--","S2O8--","Malate--",
     "PO4---","Citrate---")) {

  Names <- eval(formals(sys.function(sys.parent()))$x)
  if (is.null(x)) {
    ii <-1:length(Names);x<-Names
  } else ii <- pmatch(x,Names)  # position of element

  x <- match.arg(x, several.ok = TRUE) # check if valid input...
    
  TK <- t + 273.15
  DifGas <- DifIon <- NULL

  iGas <- ii[ii <=36]

  #  The viscosity for the true sample conditions, in centipoise
  H2OViscosity <- viscosity(t=t,S=S,P=P)

  #  Diffusion coefficients in pure water at sample temperature.
  V0 <- viscosity(S=0,t=t,P=1)

  #  To correct for salinity, the Stokes-Einstein relationship is used.
  #  This is not quite accurate, but is at least consistent.
  #  Also corrected from cm2/sec to cm2/hour
  fac <- V0/H2OViscosity*3600.

  if (length(iGas) > 0) {
    Vb <- c(28.5,31.9,23.6,36,34.7,24.5,27.9,34.5,37.3,43.8,35.2,29.2,34.9,
    16.8,37.7,50.6,53.5,49.4,74.5,69,96.6,58.0,42.5,62.6,78,100,108,145,
    167,127,166,325,41.1,63.8,86.0,108)

  #  Dissolved gases : from Wilke and Chang (1955)
  #  D = (A/VM^0.6)*7.4e-8*fac
  #    note: 1) mw = molecular weight of water
  #          2) VM = molar volumes (cm**3/mol) at boiling temp (Sherwood et al., 1975)
  #  The factor phi is reduced to 2.26 as suggested by Hayduk and Laudie (1974).

    phi <- 2.26
    mw  <- 18.0
    A   <- sqrt(phi*mw)*TK/V0
    fac2 <- 7.4e-08*fac/(3600*1e4)
    DifGas <- A/(Vb[iGas]^0.6)*fac2
    names(DifGas) <- Names[iGas]       #m2/sec
  }
  if ("H2O" %in% x) {
  #  Diffusion coefficient of Water : from Cohen and Turnbull (1959) and Krynicki et al. (1978)
    A <- 12.5e-09*exp(-5.22e-04*P)
    B <- 925.0*exp(-2.6e-04*P)
    T0 <- 95.0 + 2.61e-02*P
    D_H2O <- A*sqrt(TK)*exp(-B/(TK-T0))*1.0e+04 *fac
    DifGas <- c(DifGas, H2O= D_H2O)
  }
  if ("H3PO4" %in% x) {
  # H3PO4 : Least (1984) determined D(H3PO4) at 25 deg C and 0 S.
  #         Assume that this value can be scaled by the Stokes-Einstein
  #         relationship to any other temperature.

    D_H3PO4 <- 0.87e-05
    TS      <- 25.0
    SS      <- 0.0
    V25     <- viscosity(SS,TS,1)
    VTK     <- viscosity(SS,t,1)
    D_H3PO4 <- D_H3PO4*V25/298.15*TK/VTK*fac
    DifGas <- c(DifGas, H3PO4= D_H3PO4)
  }
  if (any(c("BOH3","BOH4") %in% x)) {
  #  B(OH)3 : Mackin (1986) determined D(B(OH)3) at 25 deg C and
  #           about 29.2 S.
  #           Assume that this value can be scaled by the Stokes-Einstein
  #           relationship to any other temperature.
    D_BOH3 <- 1.12e-05
    TS     <- 25.0
    SS     <- 29.2
    V25   <- viscosity(SS,TS,1)
    VTK   <- viscosity(0,t,1)

    D_BOH3 <- D_BOH3*V25/298.15*TK/VTK*fac
    if ("BOH3" %in% x) DifGas <- c(DifGas, BOH3= D_BOH3)
  #  B(OH)4 : No information on this species ! Boudreau and
  #           Canfield (1988) assume it is 12.5% smaller than B(OH)3.
      if ("BOH4" %in% x) DifGas <- c(DifGas,BOH3=0.875*D_BOH3)
  }
  
  if ("H4SiO4" %in% x) {
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
    DifGas <- c(DifGas,H4SiO4=D_H4SiO4)
  }

  fac2 <- 1.0e-06*fac   *1/(3600*1e4)  # second part:cm2/hour -> m2/sec

  iiIon <- ii[ii >41]
  if (length(iiIon) > 0) {
    iIon <- iiIon - 41

    #Gives units of 10-5cm2/sec. All values from Hayduk and Laudie (1974), as in
    #Boudreau 1998.
    m0<-c(54.4,31.9,4.43,6.06,9.55,10.2,10.3,7.82,9.5,4.06,2.57,3.6,3.31,
    3.31,3.39,3.31,3.63,3.43,3.18,3.36,3.69,4.46,3.91,3.31,2.79,2.95,2.78,2.71,
    25.9,15.3,4.46,10,9.6,6.29,5.06,4.02,10.4,6.35,5.99,9.81,4.66,10.3,
    9.5,4.8,4.41,4.33,3.26,4.53,4.88,4.82,3.95,5.3,4.88,3.36,2.62,2.67)

    m1<-c(1.555,1.402,0.241,0.297,0.409,0.424,0.416,0.359,0.413,0.176,0.140,
    0.179,0.152,0.152,0.158,0.15,0.208,0.144,0.155,0.13,0.169,0.198,0.199,0.151,
    0.172,0.131,0.136,0.12,1.094,0.667,0.243,0.441,0.438,0.343,0.275,0.223,
    0.273,0.280,0.307,0.432,0.252,0.331,0.388,0.245,0.241, 0.199,0.177,0.249,
    0.232,0.266,0.187,0.291,0.267,0.183,0.143,0.146)
    DifIon <-(m0[iIon]+m1[iIon]*t)*fac2
    names(DifIon) <- Names[iiIon]
  }
  return(c(DifGas,DifIon))
}
