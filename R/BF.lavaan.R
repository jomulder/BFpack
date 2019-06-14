#BF method for lavaan class objects

#Note, this method currently has very limited support for specifying hypotheses
#because the parse_hypothesis() function from bain cannot parse parameters with
#the structure "y4~~y4" or "ind60=~x2"


#' @method BF lavaan
#' @export
BF.lavaan <- function(x,
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
