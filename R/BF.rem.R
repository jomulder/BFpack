#BF method for rem class objects


#' @method BF rem
#' @export
BF.rem <- function(x,
                   hypothesis = NULL,
                   prior = NULL,
                   ...){

  #Extract summary statistics
  n <- x$df.null
  covmN <- x$cov #test
  meanN <- x$coef

  out <- BF_Gaussian(meanN, covmN, n, hypothesis, prior)
  out$model <- x
  out$call <- match.call()
  out

}
