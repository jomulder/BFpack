#BF method for rem class objects


#' @method BF rem
#' @export
BF.rem <- function(x,
                   hypothesis = NULL,
                   prior = NULL,
                   ...){

  #Extract summary statistics
  n <- x$df.null
  covmN <- x$cov
  meanN <- x$coef

  out <- Gaussian_estimator(meanN, covmN, n, hypothesis, prior)
  out$model <- x
  out

}
