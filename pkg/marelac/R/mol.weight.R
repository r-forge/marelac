##########################################################################
# molar weight
# Original version of K. Soetaert, limited to a small set of species
##########################################################################

#mol.weight <- function(species)
## converts from mol to g
#{
#  SP <- toupper(species)
#  switch(SP,
#         C     = AtomicWeight$C,
#         O     = AtomicWeight$O,
#         H     = AtomicWeight$H,
#         N     = AtomicWeight$N,
#         S     = AtomicWeight$S,
#         P     = AtomicWeight$P,
#         SI    = AtomicWeight$Si,
#         B     = AtomicWeight$B,
#         CO2   = AtomicWeight$C+2*AtomicWeight$O,
#         CO3   = AtomicWeight$C+3*AtomicWeight$O,
#         HCO3  = AtomicWeight$H+AtomicWeight$C+3*AtomicWeight$O,
#         NO3   = AtomicWeight$N+3*AtomicWeight$O,
#         NH3   = AtomicWeight$N+3*AtomicWeight$H,
#         NH4   = AtomicWeight$N+4*AtomicWeight$H,
#         HS    = AtomicWeight$S+AtomicWeight$H,
#         SO4   = AtomicWeight$S+4*AtomicWeight$O,
#         HSO4  = AtomicWeight$H+AtomicWeight$S+4*AtomicWeight$O,
#         PO4   = AtomicWeight$P+4*AtomicWeight$O,
#         HPO4  = AtomicWeight$H+AtomicWeight$P+4*AtomicWeight$O,
#         H2PO4 = 2*AtomicWeight$H+AtomicWeight$P+4*AtomicWeight$O,
#         H3PO4 = 3*AtomicWeight$H+AtomicWeight$P+4*AtomicWeight$O,
#         SIOH4 = AtomicWeight$Si+AtomicWeight$O+4*AtomicWeight$H,
#         SIO2  = AtomicWeight$Si+2*AtomicWeight$O,
#         BOH3  = AtomicWeight$B+AtomicWeight$O+3*AtomicWeight$H,
#         NA     )
#}

##########################################################################
# converts from mol to g, works with vectors and simple lists
# general parser for arbitrary species
# chemical species: whose molecular weight should be estimated
# (slower than the original limited version)
# ToDo: ??? add speedup for certain common species
##########################################################################
mol.weight <- function(species) {
  if (!is.vector(species)) stop("species must be a vector")
  if (!is.character(unlist(species))) stop("species must be character")
  molweight <- function(species) {
    with(AtomicWeight, {
      ## insert * before number (with one or more digits)
      s1 <- gsub("([0-9]+)", "*\\1+", species)
      ## insert + after capital letters
      s1 <- gsub("([a-z,A-Z])", "\\1+", s1, perl=TRUE)
      ## remove + before lower case letters
      s1 <- gsub("\\+([a-z])", "\\1", s1, perl=TRUE)
      ## replace +* with only *
      s1 <- gsub("+*", "*", s1, fixed = TRUE)
      ## remove trailing +
      s1 <- gsub("\\+$", "", s1)
      ## remove trailing +)
      s1 <- gsub("\\+\\)", ")", s1)
      ## calculate molar mass
      eval(parse(text=s1))
    })
  }
  sapply(species, molweight)
}

