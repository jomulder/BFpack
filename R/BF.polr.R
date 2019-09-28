#BF method for polr classes


#' @method BF polr
#' @export
BF.polr <- function(x,
                       hypothesis = NULL,
                       prior = NULL,
                       ...){

  #Extract summary statistics
  n <- nrow(x$fitted.values)
  sigma <- vcov(x)

  out <- BF(c(coef(x),x$zeta), hypothesis, prior, sigma=sigma, n=n)
  out$model <- x
  out$call <- match.call()
  out

}


