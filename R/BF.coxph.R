#BF method for coxph class objects


#' @method BF coxph
#' @export
BF.coxph <- function(x,
                      hypothesis = NULL,
                      prior = NULL,
                      ...){

  #Extract summary statistics
  n <- x$nevent
  covmN <- vcov(x)
  meanN <- coef(x)

  out <- BF_Gaussian(meanN, covmN, n, hypothesis, prior)
  out$model <- x
  out$call <- match.call()
  out
}




