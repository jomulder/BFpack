#BF method for zeroinfl classes


#' @method BF zeroinfl
#' @export
BF.zeroinfl <- function(x,
                    hypothesis = NULL,
                    prior = NULL,
                    ...){

  #Extract summary statistics
  n <- nrow(x$fitted.values)
  sigma <- vcov(x)
  estimate <- c(coef(x),x$zeta)

  out <- BF_Gaussian(estimate, sigma, n, hypothesis, prior)
  out$model <- x
  out$call <- match.call()
  out

}


