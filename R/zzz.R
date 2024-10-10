.onAttach <- function(libname, pkgname) {
  version <- read.dcf(
    file = system.file("DESCRIPTION", package = pkgname),
    fields = "Version"
  )
  packageStartupMessage(
    "\n",
    "This is ", paste(pkgname, version),".", "\n",
    "Updates on defaults:","\n",
    "- For standard (exploratory) tests, the default prior probability for a zero, negative,", "\n",
    "and positive effect are 0.5, 0.25, and 0.25, respectively (change using argument 'prior.hyp.explo').", "\n",
    "- For linear regression, ANOVA, t-tests, the default Bayes factor is the fractional","\n",
    "Bayes factor (change using argument 'BF.type')."
  )
}


# .onAttach <- function(libname, pkgname) {
#   packageStartupMessage("test")
# }
