.onAttach <- function(libname, pkgname) {
  version <- read.dcf(
    file = system.file("DESCRIPTION", package = pkgname),
    fields = "Version"
  )
  packageStartupMessage(
    "\n",
    "This is ", paste(pkgname, version),".", "\n",
    "Updates on default settings:","\n",
    "- For standard (exploratory) tests, the default prior probability for a zero, negative,", "\n",
    "and positive effect are 0.5, 0.25, and 0.25, respectively. The previous default was 1/3 for each","\n",
    "hypothesis. Changing these prior probabilities can be done using the argument 'prior.hyp.explo'.", "\n",
    "- For linear regression, ANOVA, t-tests, the fractional Bayes factor ('FBF') is the new default.","\n",
    "To change this to the adjusted fractional Bayes factor (the previous default), users can set","\n",
    "the argument: BF.type='AFBF'."
  )
}


# .onAttach <- function(libname, pkgname) {
#   packageStartupMessage("test")
# }
