#BF method for rem.dyad class objects


#' @method BF rem.dyad
#' @export
BF.rem.dyad <- function(x,
                        hypothesis = NULL,
                        prior = NULL,
                        ...){

  #Extract summary statistics
  n <- x$m
  covmN <- x$cov
  meanN <- x$coef

  out <- Gaussian_estimator(meanN, covmN, n, hypothesis, prior)
  out$model <- x
  out

}


