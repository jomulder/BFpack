#BF method for zeroinfl classes


#' @method BF zeroinfl
#' @export
BF.zeroinfl <- function(x,
                    hypothesis = NULL,
                    prior = NULL,
                    ...){

  #Extract summary statistics
  n <- length(x$residuals)
  sigma <- vcov(x)

  out <- BF(c(coef(x),x$zeta), hypothesis, prior, sigma=sigma, n=n)
  out$model <- x
  out$call <- match.call()
  out

}


