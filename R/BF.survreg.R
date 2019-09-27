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

  out <- BF(x$coefficients, hypothesis, prior, sigma=sigma, n=n)
  out$model <- x
  out$call <- match.call()
  out

}


