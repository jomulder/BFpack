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

  out <- BF(x$coef, hypothesis, prior, sigma=sigma, n=n)
  out$model <- x
  out$call <- match.call()
  out

}


