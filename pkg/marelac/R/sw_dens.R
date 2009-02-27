## -----------------------------------------------------------------------------
## Seawater Density
## -----------------------------------------------------------------------------

sw_dens <- function(S=35, t=25, p=max(0,P-1.013253), P=1.013253, UNESCO=FALSE) {
  if (UNESCO) {
    dens <- rho(S=S, T=t, P=p) # use seacarb function
    attributes(dens) <-NULL    # remove "unit" attribute
  }
  else      # S = absolute salinity
    dens <- 1/ sw_gibbs(S,t,p,dS=0,dt=0, dp=1)

  return(dens)
}
