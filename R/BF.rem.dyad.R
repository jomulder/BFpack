#BF method for rem.dyad class objects


#' @method BF rem.dyad
#' @export
BF.rem.dyad <- function(x,
                        hypothesis = NULL,
                        prior = NULL,
                        ...){

  #Extract summary statistics
  n <- x$m
  sigma <- x$cov
  estimate <- x$coef

  out <- BF_Gaussian(estimate, sigma, n, hypothesis, prior)
  out$model <- x
  out$call <- match.call()
  out

}


