#BF method for survreg classes


#' @method BF survreg
#' @export
BF.survreg <- function(x,
                   hypothesis = NULL,
                   prior = NULL,
                   ...){

  #Extract summary statistics
  n <- length(x$y)
  sigma <- x$var
  estimate <- x$coefficients

  out <- BF_Gaussian(estimate, sigma, n, hypothesis, prior)
  out$model <- x
  out$call <- match.call()
  out

}


