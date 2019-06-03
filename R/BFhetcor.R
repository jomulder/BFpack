


#' @importFrom MCMCpack rinvgamma
#' @method BF hetcor
#' @export
BF.hetcor <- function(x,
                       hypothesis = NULL,
                       prior = NULL,
                       ...){
  P <- nrow(x$std.errors)
  estimates <- x$correlations
  stderr <- x$std.errors[lower.tri(diag(P))]
  errcov <- diag(stderr**2)

  #for prior part use



}










