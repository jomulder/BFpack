.onAttach <- function(libname, pkgname) {
  version <- read.dcf(
    file = system.file("DESCRIPTION", package = pkgname),
    fields = "Version"
  )
  packageStartupMessage(
    "\n",
    "This is ", paste(pkgname, version),".", "\n",
    "BFpack is free software. Please report any bugs."
  )
}


# .onAttach <- function(libname, pkgname) {
#   packageStartupMessage("test")
# }
