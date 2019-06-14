#BF method for glm classes


#' @method BF glm
#' @export
BF.glm <- function(x,
                   hypothesis = NULL,
                   prior = NULL,
                   ...){

  #Extract summary statistics
  n <- nobs(x)
  covmN <- vcov(x)
  meanN <- coef(x)

  out <- Gaussian_estimator(meanN, covmN, n, hypothesis, prior)
  out$model <- x
  out

}
