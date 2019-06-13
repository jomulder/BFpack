#BF method for glmerMod class objects
#Currently specifying hypotheses doesn't work because of problems with the method 'isSymmetric'


#' @method BF glmerMod
#' @export
BF.glmerMod <- function(x,
                      hypothesis = NULL,
                      prior = NULL,
                      ...){

  warning("NB! The functionality for glmerMod class objects is in development!")
  #Extract summary statistics
  n <- nobs(x)
  covmN <- vcov(x)
  meanN <- x@beta
  names(meanN) <- rownames(covmN)

  out <- Gaussian_estimator(meanN, covmN, n, hypothesis, prior)
  out$model <- x
  out

}
