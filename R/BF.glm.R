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

  out <- BF(coef(x), hypothesis, prior, sigma=sigma, n=n)
  out$model <- x
  out$call <- match.call()
  out

}


