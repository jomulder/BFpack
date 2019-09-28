#BF method for rem class objects


#' @method BF rem
#' @export
BF.rem <- function(x,
                   hypothesis = NULL,
                   prior = NULL,
                   ...){

  #Extract summary statistics
  n <- x$df.null
  sigma <- x$cov

  out <- BF(x$coef, hypothesis, prior, sigma=sigma, n=n)
  out$model <- x
  out$call <- match.call()
  out

}
