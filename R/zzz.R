.onAttach <- function(libname, pkgname) {
  cli::cli_alert_info("BFpack ({utils::packageVersion('BFpack')}) loaded")
  cli::cli_text("Citation: https://doi.org/10.18637/jss.v100.i18")
}

# .onAttach <- function(libname, pkgname) {
#   version <- read.dcf(
#     file = system.file("DESCRIPTION", package = pkgname),
#     fields = "Version"
#   )
#   packageStartupMessage(
#     "\n",
#     "This is ", paste(pkgname, version),".", "\n",
#     "BFpack is free software. Please report any bugs."
#   )
# }

# .onAttach <- function(libname, pkgname) {
#   packageStartupMessage("test")
# }
