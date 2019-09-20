#BF method for glm classes


#' @method BF glm
#' @export
BF.glm <- function(x,
                   hypothesis = NULL,
                   prior = NULL,
                   ...){

  #Extract summary statistics
  n <- nobs(x)
  sigma <- vcov(x)
  estimate <- coef(x)

  out <- BF_Gaussian(estimate, sigma, n, hypothesis, prior)
  out$model <- x
  out$call <- match.call()
  out

}


