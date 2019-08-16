#BF method for glm classes


#' @method BF glm
#' @export
BF.glm <- function(x,
                   hypothesis = NULL,
                   prior = NULL,
                   ...){

  #Extract summary statistics
  n <- nobs(x)
  covmN <- vcov(x) #test
  meanN <- coef(x)

  out <- BF_Gaussian(meanN, covmN, n, hypothesis, prior)
  out$model <- x
  out$call <- match.call()
  out

}


