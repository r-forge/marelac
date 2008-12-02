
.onLoad <- function(lib, pkg) {
  cat("marelacTeaching package loaded\n")
  if (.Platform$GUI == "Rgui") addMarelacMenu() #else print("Marelac menue works only in RGui")
}

## this is called by:  detach(package:marelacTeaching)
.Last.lib <- function(libpath) remMarelacMenu()