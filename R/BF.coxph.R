#BF method for coxph class objects


#' @method BF coxph
#' @export
BF.coxph <- function(x,
                      hypothesis = NULL,
                      prior = NULL,
                      ...){

  #Extract summary statistics
  n <- x$nevent
  sigma <- vcov(x)

  out <- BF(coef(x), hypothesis, prior, sigma=sigma, n=n)
  out$model <- x
  out$call <- match.call()
  out
}




